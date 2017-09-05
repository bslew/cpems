// #include <time.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_sf_legendre.h>
// #include <gsl/gsl_sf_result.h>
// #include <gsl/gsl_errno.h>
// //#include <fftw3.h>
// //#include "cpeds-consts.h"
// //#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"

mscsMap& mscsMap::SH_synthesis(const mscsAlms& a, long lmax, const mscsWindowFunction& beamtf, const mscsWindowFunction&  pixtf, double smoothGauss, long method) {
	// int mscsMap::SHT_alm2map(int how, mscsAlms* a, long lmax, string pixtfFile, string smoothing, double thFWHM, int save_partial, string partial_files_prefix, string plms_file, string show_process) {
	// int show_on_fly; if (show_process=="show") { show_on_fly = 1; } else { show_on_fly = 0; }
	// int Plms_from_file;
	// bool smooth_map;
	//  mscsAlm b, aplm;
	//long int iteration_num,iteration_num_th,iteration_num_phi,iter;
	// double tmpd2,skippedms=0; // to be automated somehow ????
	/*   filenamestr speed_log_file = "speed_log_inverse"; */
	// FILE * tmpf = NULL;
	time_t eend,SStart;
	double total_time=0;
	string str;
	// Plms_tab * Plms;
	// bool was_masked = false;
	//---------------------------------------------------------------------------------


	msgs->say("Shperical Harmonic synthesis.  Method"+msgs->toStr(method)+". lmax = "+msgs->toStr(lmax),High);

	//
	// prepare stuff
	//

	// check ranges and resolutions
	if (lmax>a.lmax()) lmax=a.lmax();
	if (nside()==0) {
		msgs->warning("Nside is zero. Will use default value of 512 ",High);
		set_nside(512);
	}
	if (ordering()==unknownMapOrdering) set_ordering(nested);
	if (!mapLoaded()) { makekill_space_manager("make","T",1); clear_map();     } // if no map is loaded yet then allocate the map memory and zero it
	if (!coordLoaded()) { set_map_coord(0,0);  }

	// definition of mask for calculation: if there is non defined we do for full sky, fully transparent mask is created and deleted when it's done
	// if (!maskLoaded()) { makekill_space_manager("make","m",-1); clean_mask(); was_masked = false; } else { was_masked = true;}


	// prepare gaussian kernel
	mscsWindowFunction smGbl("gaussian smoothing",lmax,getVerbosityLevel());
	if (smoothGauss==0) smGbl.make_unit_kernel();
	else smGbl.make_gaussian_kernel(smoothGauss);

	mscsFunction w("weight",NULL,0);
	w.setVerbosityLevel(getVerbosityLevel());






	// // print info for this SHT
	// msgs->say("approximated pixel size [deg]: "+msgs->toStr(PI180inv * cpeds_pix_size_healpix(nside)),Medium);
	// msgs->say("pixel size [ster]: "+msgs->toStr(fourPI/(double)pix_num),Medium);
	// msgs->say("number of iterations per map per m = "+msgs->toStr(pix_num),Medium);
	// msgs->say("gaussian map smoothing status: "+msgs->toStr(smooth_map),Medium);
	// if (smooth_map) { msgs->say("gaussian smoothing kernel (FWHM) [deg]: "+msgs->toStr(PI180inv*thFWHM),Medium); } //changed during transition into version-1.0




	/*   //-------------------------------------------------------------------------------------------------------------------------- */

	if (method == 1) {

		long8 i,ring_num, *ThBreaks;
		long8 l,m;
#ifdef GO128BIT
		flt8 *ThVals, *phi0, *dphi;
#else
		flt8 *ThVals, *phi0, *dphi;
#endif
		//    mask_map_merge();
		// pixel *mapRING = new pixel; // ring version of the map
		// *mapRING = clone_map_space(map,pix_num);
		conv_nest2ring();
		ring_num = ringNum(); // 2(2*nside-1)+1
		msgs->say("number of rings in map: "+msgs->toStr(ring_num),Low);

		ThVals = new flt8[ring_num]; // tabulated sin and cos  functions
		ThBreaks = new long8[ring_num]; // this holds the ring number of each pixel in the map
		phi0 = new flt8[ring_num];
		dphi = new flt8[ring_num];

		// /* converting to ring ordering from  optimalization reasons */
		// copy_map(map,mapRING,pix_num,NULL);
		// conv_nest2ring(mapRING);
		// precalculate_fourier_stuff(ring_num, mapRING,ThVals,ThBreaks, phi0,dphi);

		precalculate_fourier_stuff(ring_num, *this,ThVals,ThBreaks, phi0,dphi);
		// lm_num = alm2num(0,0); lm_num_start = lm_num;
		SStart=time(NULL);

		SHT_FFTW_COMPLEX *alm_fftw = new SHT_FFTW_COMPLEX[a.almsNum()]; //for (i=0; i<alms_num;i++) { alm_fftw[i].re = alm[i].R; alm_fftw[i].im = alm[i].I; }
		coordStruct coords;
		coords.nPix = pixNum();
		coords.nThetaVals = ring_num;
		coords.thetaVals=ThVals;
#ifdef GO128BIT
		coords.thetaBreaks=(long8*)ThBreaks;
#elif GO64BIT
		coords.thetaBreaks=(long8*)ThBreaks;
#else
		coords.thetaBreaks=(long8*)ThBreaks;
#endif
		coords.phi0=phi0;
		coords.deltaPhi=dphi;
		coords.nGaps=0;
		coords.gaps=NULL;
		coords.pixSize=4*PI/(double)pixNum();
		SHT_FFTW_COMPLEX * SHTout = new SHT_FFTW_COMPLEX[pixNum()];


		// smoothing in fourier space
		msgs->say("smoothing in SH space:",Low);
		i=0;

		w=beamtf * smGbl * pixtf;
		w.setPointsNum(lmax+1);
		double w0;
		for (m=-lmax;m<=lmax;m++) { 	  //start = clock();
			for (l=abs(m); l<=lmax; l++) {
				w0=w.getY(l);
				alm_fftw[i][0] = a.get(i).real() * w0;
				alm_fftw[i][1] = a.get(i).imag() * w0;
//							if (cpeds_isnan(alm_fftw[i][0])) printf("alm(%li) re l=%li,m=%li, is nan\n",i,l,m);
//							if (cpeds_isnan(alm_fftw[i][1])) printf("alm(%li) im l=%li,m=%li, is nan\n",i,l,m);
				i++;
			}
		}

		msgs->say("doing SHT",Low);
		backwardSHT(alm_fftw,coords,lmax,SHTout,1);
		long pix_num=pixNum();
		for (i=0;i<pix_num;i++) { T(i) = SHTout[i][0];  
//			if (cpeds_isnan(T(i))) printf("T(%li) is nan\n",i);

		} // for real SHT we only use real part of the output
		// for (i=0;i<pix_num;i++) { T(i) = SHTout[i].im; }
		destroyCoordsCPP(&coords);
		delete [] alm_fftw;
		delete [] SHTout;

		//converting to nested ordering
		// copy_map(mapRING,map,pix_num,NULL); // this conversion could be done just in here !!! but in case the conv_tabs are not loaded this is safer
		// kill_map_space(mapRING);   delete mapRING;
		conv_ring2nest();

		eend = time(NULL); total_time = ((double)eend - (double)SStart) / 60.0;
		msgs->say("SHT done, total time: "+msgs->toStr(total_time)+" [min]",High);   //   fprintf(tmpf,"%i %lf\n",l,total_time);

	}

	/*   -------------------------------------------------------------------------------------------------------------------------- */

	if (method == 2) {
#ifndef INCLUDE_SPHEREPACK
		msgs->say("spherepack method for SHT was chosen and INCLUDE_SPHEREPACK compilation flag was not used at compilation time. Will stop here.", High);
		exit(0);
#else

		conv_nest2ring();
		long nlat= ringNum(); // 2(2*nside-1)+1
		msgs->say("number of rings in map: "+msgs->toStr(nlat),Low);

		SStart=time(NULL);

		fftw_complex *alm_fftw = new fftw_complex[a.almsNum()]; //for (i=0; i<alms_num;i++) { alm_fftw[i].re = alm[i].R; alm_fftw[i].im = alm[i].I; }
		coordStruct coords;
		coords.nPix = pixNum();
		coords.nThetaVals = ring_num;
		coords.thetaVals=ThVals;
#ifdef GO128BIT
		coords.thetaBreaks=(long long*)ThBreaks;
#elif GO64BIT
		coords.thetaBreaks=(long*)ThBreaks;
#else
		coords.thetaBreaks=(int*)ThBreaks;
#endif
		coords.phi0=phi0;
		coords.deltaPhi=dphi;
		coords.nGaps=0;
		coords.gaps=NULL;
		coords.pixSize=4*PI/(double)pixNum();
		fftw_complex * SHTout = new fftw_complex[pixNum()];


		// smoothing in fourier space
		msgs->say("smoothing in SH space:",Low);
		i=0;
		w=beamtf * smGbl * pixtf;
		w.setPointsNum(lmax+1);
		double w0;
		for (m=-lmax;m<=lmax;m++) { 	  //start = clock();
			for (l=abs(m); l<=lmax; l++) {
				w0=w.getY(l);
				alm_fftw[i].re = a.get(i).real() * w0;
				alm_fftw[i].im = a.get(i).imag() * w0;
				i++;
			}
		}

		msgs->say("doing SHT",Low);
		backwardSHT(alm_fftw,coords,lmax,SHTout,1);
		long pix_num=pixNum();
		for (i=0;i<pix_num;i++) { T(i) = SHTout[i].re;  } // for real SHT we only use real part of the output
		// for (i=0;i<pix_num;i++) { T(i) = SHTout[i].im; }
		destroyCoords(&coords);
		delete [] alm_fftw;
		delete [] SHTout;

		//converting to nested ordering
		conv_ring2nest();

		eend = time(NULL); total_time = ((double)eend - (double)SStart) / 60.0;
		msgs->say("SHT done, total time: "+msgs->toStr(total_time)+" [min]",High);   //   fprintf(tmpf,"%i %lf\n",l,total_time);



#endif
	}

	/*   //-------------------------------------------------------------------------------------------------------------------------- */

	if (method == 7) { // my own calculation with GSL (but not for long) -- calculate all m's  and calculate only for the pixels on sky -- this does not average the transform at points that belongs to the same pixel since the transform is calculated only for the directions of pixels in the map; but also precalculates the Plm for every theta in the map at cost of larger memory usage; corrected for the d\Omega stuff = 1/pix_num in integration; acclelerated with precalculated Plm's for the rings in the map.
		// as a matter of fact the below form of the fourier transform forces the antisymmetrized alms beacuse the negaitve part of them
		//is not used at all. the antisymmetrization in the generation of the alms is only for the exactness purposes.
		// this method uses precomputed Plms if required

		//     // setting graphical visualization of the progress
		//     if (show_on_fly == 1) show_progress("map_start","");

		//     pixel *mapl = new pixel; pixel *ptr = map; // map of a single multipole
		//     pixel *mapRING = new pixel; // ring version of the map
		//     *mapl = clone_map_space(map,pix_num);
		//     *mapRING = clone_map_space(map,pix_num);

		// /*     /\*     seed the coordinates in the map -- if not yet loaded   *\/ */
		// /*     if (coord_loaded == 0) {  printf("settig coordinates for the map\n"); set_map_coord(0,0);    } */

		//     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*nside-1)+1
		//     printf("  -- number of rings in map= %li \n",ring_num);
		//     Ptab = new double[ring_num]; // table keeps the associated Legendre polynominals values calculated at all rings
		//     ringtab = new long int[pix_num]; // this holds the ring number of each pixel in the map
		//     ringpixtab = new long int[pix_num]; // this keeps the number of the pixel (-1) that starts the next ring
		//     COS = new double[ring_num]; // tabulated sin and cos  functions
		//     SIN = new double[ring_num];

		//     /* converting to ring ordering from  optimalization reasons */
		//     //for (i=0;i<pix_num;i++) { mapRING[i] = map[i];} // this is done so because the different names of te types maptmp and mapRING belong to. probably this could be done better by reorganizing the structures of the object with relations with other libraries like CPEDS etc.
		//     copy_map(map,mapRING,pix_num,NULL);
		//     conv_nest2ring(mapRING);    // the furtther calculations are performed on mapl and summed to mapRING in RING ordering

		//     precalculate_fourier_stuff(COS,SIN,mapRING,ring_num,ringtab,ringpixtab);
		//     rmax = 2*nside-1; // maximal ring number in northern hemisphere
		//     //calculate_Plm_harmonics(COS,ring_num-1,Ptab,1,1);
		//     /*     inverse Fourier transformation  itself   */
		//     Start=clock();
		//     for (l=0;l<=lmax;l++) { 	  start = clock();
		//       clear_map(mapl);
		//       for (m=0; m<=l; m++) { // this is because it's enough the calculate half of them since it's the output is to be real and the al-m = (-1)^m * alm^*. notice the 2 factor in the transformation// also Y_l-m = (-1)^m*Y_lm^*
		//       //for (m=-l; m<=l; m++) { // full version
		// 	//setting Plm_harmonic function values accordingly to the pixelization system
		// 	if (strcmp(plms_file,"") != 0) {  calculate_Plm_harmonics(ring_num,Ptab,Plms);} else { calculate_Plm_harmonics(COS,ring_num,rmax,Ptab,l,m); }

		// 	lm_num = alm2num(l,m);
		// 	philm = alm[lm_num].P; Dlm = alm[lm_num].D; if (m == 0) { Dlm = 0.5*Dlm; } // this is because in the transform there is a factor 2 which can be applied only for sumation from m=1...l a_l0, should be summed only once,! hence 0.5 factor here
		// // ***** EXPERIMENTAL *****
		// /* 	Rlm = sqrt(C_l[l][1]/2)*cpeds_random_gauss_number(0,1,1000,1);	Ilm = sqrt(C_l[l][1]/2)*cpeds_random_gauss_number(0,1,1000,1); */
		// /* 	Dlm = sqrt(pow(Rlm,2)+pow(Ilm,2)); philm = cpeds_cart2sph(1,Rlm,Ilm,0.0); */
		// /* 	printf("Rlm=%lE Ilm=%lE, Dlm=%lE, philm=%.3lf |||| sqrt(Rlm^2+Ilm^2)=%lE  ||| C_l[][0] = %lE, C_l[][1] = %lE\n",Rlm,Ilm,Dlm,philm,sqrt(Rlm*Rlm+Ilm*Ilm),C_l[l][0],C_l[l][1]); */
		// /* 	//philm = cpeds_random_uniform_number(0,twoPI); Dlm = C_l[l][1]/2*sqrt(pow(cpeds_random_gauss_number(0,1,1000,1),2)+pow(cpeds_random_gauss_number(0,1,1000,1),2)); */
		// /* 	if (m==0) { Dlm *=2; philm=0; } */
		// // ***** EXPERIMENTAL *****

		// 	for (i=0; i<pix_num;i++) {
		// 	  Plm = Ptab[ringtab[i]]; //printf("l=%i, m=%i, i=%li, phi=%lE pi  cos(th)=%lE ringpixtab[%li]=%li, plm = %lE, ringtab[i]=%li\n", l,m,i,phi/PI,COS[ringtab[i]],i,ringpixtab[i],Plm, ringtab[i]);
		// 	  if (Plm == 0.0) { //printf("l=%i, m=%i, i=%li  ",l,m,i);
		// 	    i = ringpixtab[i]; //printf("i=%li\n",i);
		// 	  }
		// 	  else {
		// 	    phi = mapRING->n[i].l;
		// 	    Ylm_re = Plm * cos( (double)m * phi + philm ); // we take here the real part of a_lm * Y_lm	    //mapl->T[i] = mapl->T[i] + Dlm * Ylm_re; // T(n) = \Sum_lm D_lm * P_lm(cos(\th) * cos(m \phi + \phi_lm) = \Sum_lm Re( a_lm * Y_lm(n) )
		// 	    //mapl->T[i] = mapl->T[i] + 2 * Dlm * Ylm_re  * mapRING[i].m * gauss_sm[l]; // T(n) = \Sum_lm D_lm * P_lm(cos(\th) * cos(m \phi + \phi_lm) = \Sum_lm Re( a_lm * Y_lm(n) )
		// 	    mapl->T[i] += 2 * Dlm * Ylm_re  * mapRING->m[i] * bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l);
		// 	    //mapl->T[i] = mapl->T[i] + Dlm * Ylm_re * map[i].m; // full version
		// 	    //printf("pix = %li pix_num = %li\n",pix, pix_num);
		// 	  }
		// 	}
		//       }

		//       for (i=0; i< pix_num;i++) { mapRING->T[i] += mapl->T[i]; } // add the whole multipole to the summed map

		//       // and from here on some useless stuff: print out info, save partial files, save speed log, and visualize the progress
		//       tmpd2 = (double)pix_num*(2*(double)l+1);
		//       printf("l=%li, skipped points for this l: %.0lf/%.0lf (%lf %%) ",l,skippedms,tmpd2,skippedms/tmpd2*100);

		//       end = clock();         total_time_p=total_time;       total_time = ((double)end - (double)start) / (double)CLOCKS_PER_SEC;
		//       cpu_time_used = ((double)end - (double)Start) / (double)CLOCKS_PER_SEC;  tot = 0.5*(total_time-total_time_p)*pow((double)lmax,2);       //tot = 0.5*()*pow((double)lmax,2);
		//       printf("time per l: %.2lf [s] TOT(l=%li)=%.2lf [h] ETA = %0.2lf [s]\n",total_time,lmax,tot/3600.0,cpu_time_used);      //fprintf(tmpf,"%i %lf\n",l,total_time);
		//       if (show_on_fly == 1) {map = mapRING; plot_map(1,0,0,0,50,"color","",14); map=ptr; }// visualize map

		//       if (save_partial >= 1) {  // this saves the individal temperature maps for each multipole and for sum of multipoles upto l
		// 	map=mapRING; sprintf(flatmap_file_name,"%s_sum_lmax%li",partial_files_prefix,l);  savebinT(flatmap_file_name,0); map=ptr; // saves sum of multipoles
		// 	if (save_partial == 2) {
		// 	  map = mapl; sprintf(flatmap_file_name,"%s_l%li",partial_files_prefix,l); savebinT(flatmap_file_name,0); map = ptr; // saves one multipole
		// 	}
		//       }
		//       //vis_xtab[l] = (float)l; vis_ytab[l] = (float)total_time; //vis_x = (float)l; vis_y = (float)total_time;
		//       //vis_iter++; //show_progress("xytabFourier2"); //vis_x_prev=vis_x;       vis_y_prev=vis_y;
		//     }

		//     //converting to nested ordering
		//     //for (i=0;i<pix_num;i++) { map[i] = mapRING[i]; } // if you want that the transform adds only the output values to the existing map then this line should be changed accordingly
		//     copy_map(mapRING,map,pix_num,NULL); // this conversion could be done just in here !!!
		//     conv_ring2nest(map);

		//     kill_map_space(mapRING);     kill_map_space(mapl);    delete mapRING; delete mapl;  //exit(0);
		//     delete SIN; delete COS; delete ringtab; delete ringpixtab; delete Ptab; //delete gauss_sm;
		//     if (flat_coord_loaded == 1) { killflatcoord(); }
		//     //fclose(tmpf);
		//     if (show_on_fly == 1) show_progress("end","");
		// /*     if (show_on_fly == 1) cpgclos(); */
	}

	/*   //-------------------------------------------------------------------------------------------------------------------------- */

	if (method == 8) { // this is my own but improved calculation by the partial separation of variables. Instead of N^4 complexity this should have 2*N^3 of so.
		// the accuracy level should be the same as for method 7. no approximations are made.
		// this method uses precomputed Pmthls if required

		//     pixel *mapRING = new pixel; // ring version of the map
		//     *mapRING = clone_map_space(map,pix_num);

		//     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*nside-1)+1
		//     printf("  -- number of rings in map= %li \n",ring_num);
		//     Ptab = new double[lmax+1]; // table keeps the associated Legendre polynominals values calculated at all rings
		//     ringtab = new long int[pix_num]; // this holds the ring number of each pixel in the map
		// /*     ringpixtab = new long int[pix_num]; // this keeps the number of the pixel (-1) that starts the next ring  */
		//     COS = new double[ring_num]; // tabulated sin and cos  functions
		// /*     SIN = new double[ring_num]; */
		//     count = ring_num*(lmax+1);
		//     b = new a_lm[count];    for (i=0;i<count;i++) { b[i].R = 0; b[i].I = 0; } count = 0; // here \Sum_l a_lm*P_lm*B_l is kept for all th knots and m values from 0 to lmax
		//     aplm = new a_lm[alms_num];
		//     /* converting to ring ordering from  optimalization reasons */
		//     printf("  -- making  ring ordering\n");
		//     copy_map(map,mapRING,pix_num,NULL);
		//     conv_nest2ring(mapRING);    // the furtther calculations are performed on mapl and summed to mapRING in RING ordering
		//     precalculate_fourier_stuff(COS,mapRING,ring_num,ringtab);
		//     //rmax = 2*nside-1; // maximal ring number in northern hemisphere
		//     lm_num = alm2num(0,0); lm_num_start = lm_num;
		//     //printf("alm2num(lmax,lmax) = %li, alm2num(0,0) = %li",alm2num(lmax,lmax),alm2num(0,0));
		//     /*     inverse Fourier transformation  itself   */
		//     SStart=time(NULL);


		//     // smoothing in fourier space
		//     //if (smooth_map == 1) {
		//     printf("  -- smoothing in fourier space: \n");
		//     i=lm_num;
		//     for (m=0;m<=lmax;m++) { 	  //start = clock();
		//       //for (l=m; l<=lmax; l++) {	aplm[i].R = alm[i].R*gauss_sm[l]; aplm[i].I = alm[i].I*gauss_sm[l]; i++;      }
		//       for (l=m; l<=lmax; l++) {
		// 	aplm[i].R = alm[i].R * bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l);
		// 	aplm[i].I = alm[i].I * bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l);
		// 	i++;
		//       }
		//     }
		//     //}

		//     // this calculates the b_l(th) coefficients

		//     printf("  -- calculating b_m(th) coefficients: \n");
		//     for (m=0;m<=lmax;m++) { 	  //start = clock();
		//       lm_num_start = lm_num;
		//       for (i=0; i<ring_num;i++) {
		// 	lm_num = lm_num_start;
		// 	if (Plms_from_file == 1) {  skip = (lmaxP-lmax); /* skip=0; */ calculate_Pmth_harmonics(m,skip,Ptab,Plms);  } else { calculate_Pmth_harmonics(Ptab,m,COS[i]); }
		// 	for (l=m; l<=lmax; l++) {
		// 	  Plm = Ptab[l-m];
		// 	  if (Plm != 0.0) { b[count].R += aplm[lm_num].R*Plm; b[count].I += aplm[lm_num].I*Plm; }
		// 	  //printf("i=%li, aplm[i].R=%lE, aplm[i].I=%lE, l=%i, m=%i, cos(th) = %lE Plm = %lE, count=%li, b[count].R=%lE, b[count].I=%lE, lm_num = %li\n",i,aplm[i].R,aplm[i].I,l,m,COS[i],Plm,count,b[count].R, b[count].I,lm_num);
		// 	  lm_num++;
		// 	}
		// 	count++;
		//       }
		//       /*       eend = time(NULL); total_time = (double)eend - (double)SStart; tot = total_time*(double)lmax/(double)m; */
		//       /*       printf("  -- mode number: %i of %i, TOT = %.2lf [min] ETA = %.2lf [min]\n",m,lmax,tot/60.0,(tot-total_time)/60.0); */
		//     }

		//     eend = time(NULL);  total_time = ((double)eend - (double)SStart) / 60.0;
		//     printf("  -- time for b_l(th): %.2lf [min]\n",total_time); //     fprintf(tmpf,"%i %lf\n",l,total_time);

		//     // filling the rest

		//     count = 0;
		//     for (m=0; m<=lmax; m++) {
		//       for (i=0; i<ring_num; i++) {
		// 	b[count].D = sqrt(pow(b[count].R,2)+pow(b[count].I,2)); b[count].P = cpeds_cart2sph(1,b[count].R,b[count].I,0);
		// 	if (m == 0) { b[count].D /= 2; }
		// 	count++;
		//       }
		//     }

		//     // now we do the inverse fourier transform like FFT in phi direction---> \Sum_m=-lmax^lmax b_m(th_i) * e^(im\phi_i) * M(th_i,phi_i)
		//     // this could be done by GSL function of FFTW library in paralel implementation since there are many rings to transform.

		//     printf("  -- doing inverse FFT in phi direction  --  summing over the pixels in the map (almost there...:) \n");
		//     for (i=0; i<pix_num; i++) {
		//       count = ringtab[i];
		//       phi = mapRING->n[i].l;
		//       for (m=0; m<=lmax; m++) {
		// 	Dlm = b[count].D; philm = b[count].P;
		// 	mapRING->T[i] += 2 * Dlm * cos( (double)m * phi + philm )  * mapRING->m[i]; // T(n) = \Sum_lm D_lm * P_lm(cos(\th) * cos(m \phi + \phi_lm) = \Sum_lm Re( a_lm * Y_lm(n) )
		// 	count += ring_num;
		//       }
		//     }

		//     EEnd = time(NULL);   total_time = ((double)EEnd - (double)SStart) / 60.0;
		//     printf("total time: %.2lf [min] \n",total_time);   //   fprintf(tmpf,"%i %lf\n",l,total_time);
		//     //converting to nested ordering
		//     copy_map(mapRING,map,pix_num,NULL); // this conversion could be done just in here !!! but in case the conv_tabs are not loaded this is safer
		//     kill_map_space(mapRING);   delete mapRING;
		//     conv_ring2nest(map);

		//     // free space

		//     delete COS; delete ringtab; delete Ptab; //delete gauss_sm;
		//     delete b; delete aplm;
		//     if (flat_coord_loaded == 1) { killflatcoord(); }
		//     //fclose(tmpf);
		// /*     if (show_on_fly == 1) cpgclos(); */
	}


	/*   //-------------------------------------------------------------------------------------------------------------------------- */

	if (method == 9) { // like in 8 but uses fftw for integration over phi
		//     lmax_loc = lmax;    lmax = 2*nside; // cheat it

		//     // this method uses precomputed Pmthls if required

		//     pixel *mapRING = new pixel[pix_num]; // ring version of the map
		//     *mapRING = clone_map_space(map,pix_num);

		//     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*nside-1)+1
		//     printf("  -- number of rings in map= %li \n",ring_num);
		//     Ptab = new double[lmax+1]; // table keeps the associated Legendre polynominals values calculated at all rings
		//     ringtab = new long int[pix_num]; // this holds the ring number of each pixel in the map
		// /*     ringpixtab = new long int[pix_num]; // this keeps the number of the pixel (-1) that starts the next ring  */
		//     COS = new double[ring_num]; // tabulated sin and cos  functions
		// /*     SIN = new double[ring_num]; */
		//     count = ring_num*(lmax+1);
		//     b = new a_lm[count];    for (i=0;i<count;i++) { b[i].R = 0; b[i].I = 0; } count = 0; // here \Sum_l a_lm*P_lm*B_l is kept for all th knots and m values from 0 to lmax
		//     aplm = new a_lm[alms_num];
		//     /* converting to ring ordering from  optimalization reasons */
		//     printf("  -- making  ring ordering\n");
		//     copy_map(map,mapRING,pix_num,NULL);
		//     conv_nest2ring(mapRING);    // the furtther calculations are performed on mapl and summed to mapRING in RING ordering
		//     precalculate_fourier_stuff(COS,mapRING,ring_num,ringtab);
		//     //rmax = 2*nside-1; // maximal ring number in northern hemisphere
		//     lm_num = alm2num(0,0); lm_num_start = lm_num;
		//     //printf("alm2num(lmax,lmax) = %li, alm2num(0,0) = %li",alm2num(lmax,lmax),alm2num(0,0));
		//     /*     inverse Fourier transformation  itself   */
		//     SStart=time(NULL);

		//     // smoothing in fourier space
		//     //if (smooth_map == 1) {
		//       printf("  -- smoothing in fourier space: \n");
		//       i=lm_num;
		//       for (m=0;m<=lmax;m++) { 	  //start = clock();
		// 	//for (l=m; l<=lmax; l++) {	aplm[i].R = alm[i].R*gauss_sm[l]; aplm[i].I = alm[i].I*gauss_sm[l]; i++;      }
		// 	for (l=m; l<=lmax; l++) {
		// 	  aplm[i].R = alm[i].R * bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l);
		// 	  aplm[i].I = alm[i].I * bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l);
		// 	  i++;
		// 	}
		//       }
		//       //}


		//     // this calculates the b_m(th) coefficients
		//     printf("  -- calculating b_m(th) coefficients: \n");
		//     for (m=0;m<=lmax;m++) { 	  //start = clock();
		//       lm_num_start = lm_num;
		//       for (i=0; i<ring_num;i++) {
		// 	lm_num = lm_num_start;
		// 	if (Plms_from_file == 1) {  skip = (2*nside-(long)lmax); calculate_Pmth_harmonics(m,skip,Ptab,Plms);  } else { calculate_Pmth_harmonics(Ptab,m,COS[i]); }
		// 	for (l=m; l<=lmax; l++) {
		// 	  Plm = Ptab[l-m];
		// 	  if (Plm != 0.0) { b[count].R += aplm[lm_num].R*Plm; b[count].I += aplm[lm_num].I*Plm; }
		// 	  //printf("i=%li, aplm[i].R=%lE, aplm[i].I=%lE, l=%i, m=%i, cos(th) = %lE Plm = %lE, count=%li, b[count].R=%lE, b[count].I=%lE, lm_num = %li\n",i,aplm[i].R,aplm[i].I,l,m,COS[i],Plm,count,b[count].R, b[count].I,lm_num);
		// 	  lm_num++;
		// 	}
		// 	count++;
		//       }
		// /*       eend = time(NULL); total_time = (double)eend - (double)SStart; tot = total_time*(double)lmax/(double)m; */
		// /*       printf("  -- mode number: %i of %i, TOT = %.2lf [min] ETA = %.2lf [min]\n",m,lmax,tot/60.0,(tot-total_time)/60.0); */
		//     }

		//     eend = time(NULL);  total_time = ((double)eend - (double)SStart) / 60.0;
		//     printf("  -- time for b_l(th): %.2lf [min]\n",total_time); //     fprintf(tmpf,"%i %lf\n",l,total_time);

		//     // filling the rest

		//     count = 0;
		//     for (m=0; m<=lmax; m++) {
		//       for (i=0; i<ring_num; i++) {
		// 	b[count].D = sqrt(pow(b[count].R,2)+pow(b[count].I,2)); b[count].P = cpeds_cart2sph(1,b[count].R,b[count].I,0);
		// 	if (m == 0) { b[count].D /= 2; }
		// 	count++;
		//       }
		//     }

		//     printf("  -- doing inverse FFT in phi direction  --  FFTW \n");


		// /*  FFTW phi direction transform  */

		// /*     fftw_complex *inf, *outf, *inb, *outb; */
		// /*     fftw_plan pf,pb; */
		// /*     long xpix = 2*nside; // this is such in case the lmax is not set to some the maximal possible value, in which case FFTW would not work OK. */
		// /*     //inf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); */
		// /*     inb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); */
		// /*     //outf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); */
		// /*     outb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); */
		// /*     //pf = fftw_plan_dft_1d(lmax, inf, outf, FFTW_FORWARD, FFTW_MEASURE); */
		// /*     pb = fftw_plan_dft_1d(lmax, inb, outb, FFTW_BACKWARD, FFTW_MEASURE); */
		// /*     count2 = 0; */
		// /*     for (i=0; i<ring_num;i++) { */
		// /*       count = ring_num+i; // we start copy from m=1 not m=0 for FFT */
		// /*       b0th = b[i].R; // this is the transform for m=0 for all theta */
		// /*       for (m=0; m<xpix; m++) { // loop for all m = 2n_s values */
		// /* 	inb[m][0] = b[count].R; inb[m][1] = b[count].I; */
		// /* 	//inf[m][0] = b[count].R; inf[m][1] = -b[count].I; */
		// /* 	count +=ring_num; */
		// /*       } // load array of b_m(th) for a given theta */
		// /*       fftw_execute(pb); // do the transform */
		// /*       pixinring = cpeds_get_pixnum_in_ring_healpix(nside,i); */
		// /*       for (m=0; m<pixinring; m++) { */
		// /* 	phi = mapRING->n[count2].l; x = (long)round(phi*(xpix-1)/twoPI)+1; // the region m (1...lmax) is mapped onto phi <0,twoPI> */
		// /* 	x=m; if (x>=xpix) x= xpix-1; */
		// /* 	mapRING->T[count2] = outb[x][0] + outb[x][0] + b0th;  */
		// /* 	if (m==pixinring-1 || x==xpix)  mapRING->T[count2] = outb[0][0] + outb[0][1] + b0th; // cyclic condition ?!?! */
		// /* 	count2++; */
		// /*       } */
		// /*     } */
		// /*     printf(" count2=%li \n",count2); */
		//     //fftw_destroy_plan(pf);
		// /*     fftw_destroy_plan(pb); */
		// /*     fftw_free(inb); fftw_free(outb); */
		//     //fftw_free(inf); fftw_free(outf);

		// /*  ----------------------------  */


		//     EEnd = time(NULL);   total_time = ((double)EEnd - (double)SStart) / 60.0;
		//     printf("total time: %.2lf [min] \n",total_time);   //   fprintf(tmpf,"%i %lf\n",l,total_time);


		//     copy_map(mapRING,map,pix_num,NULL); // this conversion could be done just in here !!! but in case the conv_tabs are not loaded this is safer
		//     kill_map_space(mapRING);   delete mapRING;
		//     conv_ring2nest(map);

		//     // free space

		//     delete COS; delete ringtab; delete Ptab; //delete gauss_sm;
		//     delete b; delete aplm;
		//     if (flat_coord_loaded == 1) { killflatcoord(); }
		//     //fclose(tmpf);
		// /*     if (show_on_fly == 1) cpgclos(); */


	}

	// if (was_masked == false) { makekill_space_manager("kill","m",-1); }
	// kill_window_functions();
	//loaded.m = true;
	// if (how != 1) if (strcmp(plms_file,"") != 0) {  delete Plms; }
	loaded.T=true;
	return *this;
}
/* **************************************************************************************************** */
mscsMap& mscsMap::SH_synthesis(const mscsAlms& a, long lmax, const mscsWindowFunction& beamtf, string pixtfType, double smoothGauss, long method) {
	//
	// prepare pixel transfer function
	//
	mscsWindowFunction pixtf("pixtf",lmax,getVerbosityLevel());
	cpedsStatusCodes res;
	string str;
	if (pixtfType=="healpix") {
		str=MSCS_DATA_DIR+MSCS_GLOBAL__HEALPIX_PIXTF_PREF+msgs->toStr(nside())+MSCS_GLOBAL__HEALPIX_PIXTF_SUFF;
		pixtf=readPixTf(str,&res);
	}
	else {
		pixtf=readPixTf("",&res);
		if (res!=cpedsSuccess) {
			pixtf.setPointsNum(lmax+1);
			pixtf.make_unit_kernel();
		}
	}


	//
	// do SHT
	//
	return SH_synthesis(a,lmax,beamtf,pixtf,smoothGauss,method);
}




























const mscsAlms mscsMap::SH_analysis(long lmax, const mscsWindowFunction& beamtf, string pixtfType, double desmoothGauss, long method) {

	//
	// prepare pixel transfer function
	//
	mscsWindowFunction pixtf("pixtf",lmax,getVerbosityLevel());
	cpedsStatusCodes res;
	string str;
	if (pixtfType=="healpix") {
		str=MSCS_DATA_DIR+MSCS_GLOBAL__HEALPIX_PIXTF_PREF+msgs->toStr(nside())+MSCS_GLOBAL__HEALPIX_PIXTF_SUFF;
		pixtf=readPixTf(str,&res);
	}
	else {
		pixtf=readPixTf("",&res);
		if (res!=cpedsSuccess) {
			pixtf.setPointsNum(lmax+1);
			pixtf.make_unit_kernel();
		}
	}


	//
	// do SHT
	//
	return SH_analysis(lmax,beamtf,pixtf,desmoothGauss,method);
}






//************************************************************************
//                      CALCULATES THE ASYMETRICAL SINE FORWARD FOURIER TRANSFORMATION ( Real data --> Complex data )
// how  - the method of calculations
// how = 1 - ccSHTx
// how = 2 - GSL skipping m<0 for l>30
// save_partial = 1 - save the partial files with alms according to names given by partial_files_prefix
// save_partial = 0 - do not save partial files
// partial_files_prefix - prefix for saving the partial files for maps with increasing number of multipoles included and individual multipoles maps
// resolution - not used
// !!!! 2006/11/24 -- redefinition of the unsued so far parameter resolution !!!!
// if resolution parameter is given other than 0, then the pixel transfer function is used for that resolution. this parameter specifies the ns resolution of the Healpix  pix system.

// pix_system - 1 - HEALPIX pixelization system - the transform is performed using given pixelization system
// pix_system - 2 - TETRAHEDRON pixelization system (not implemented yet)

const mscsAlms mscsMap::SH_analysis(long lmax, const mscsWindowFunction& beamtf, const mscsWindowFunction& pixtf, double desmoothGauss, long method) {

	//int mscsMap::calculate_transformF(int how, double resolution, strarg smoothing, double thFWHM, int save_partial, strarg partial_files_prefix, strarg plms_file) {  // makes alms from a map
	//long int pix,thi, iter,iteration_num,iteration_num_th,iteration_num_phi;
	//double calc_level; // to be automated somehow ????
	//double tmpd,tmpd2,skippedms=0;

	long8 i,ring_num;
	long8 l,m, *ThBreaks;
#ifdef GO128BIT
	flt8 *dphi,  *phi0;
	flt8  *ThVals;
#else
	flt8 *dphi,  *phi0;
	flt8  *ThVals;
#endif	
	double total_time=0;
	time_t eend,SStart;
	//filenamestr speed_log_file = "speed_log_forward";
	bool was_masked = false;


	msgs->say("Shperical Harmonic analysis.  Method"+msgs->toStr(method)+". lmax = "+msgs->toStr(lmax),High);


	//
	// prepare stuff
	//

	// check ranges and resolutions
	mscsAlms a("alms",lmax);
	if (nside()==0) {
		msgs->error("Nside is zero. Will not continue",High);
		return a;
	}
	if (!coordLoaded()) { set_map_coord(0,0);  }
	//definition of mask for calculation: if there is non defined we do for full sky, fully transparent mask is created and deleted when it's done
	if (!maskLoaded()) { makekill_space_manager("make","m",-1); clean_mask(); was_masked = false; } else { was_masked = true;}


	// prepare gaussian kernel
	mscsWindowFunction smGbl("gaussian smoothing",lmax,getVerbosityLevel());
	if (desmoothGauss==0) smGbl.make_unit_kernel();
	else smGbl.make_gaussian_kernel(desmoothGauss);

	mscsFunction w;



	//   // print info for this SHT
	// printf("  -- approximated pixel size [deg]: %lf\n",180/PI * cpeds_pix_size_healpix(nside));
	// printf("  -- number of iterations per map per m = %li \n",pix_num);
	// printf("  -- map (de)smoothing in fourier sapce is: "); if (smooth_map == 1) { printf("ON, gaussian beam (FWHM) [deg]: %lE\n",180/PI*thFWHM); } else printf("OFF\n");

	//************************************************************************
	if (method == 1) { // this uses the ccSHTx software for calculation

		mask_map_merge();

		conv_nest2ring();
		ring_num = ringNum(); // 2(2*nside-1)+1
		msgs->say("number of rings in map= "+msgs->toStr(ring_num),Medium);

		ThVals = new flt8[ring_num]; // tabulated sin and cos  functions
		ThBreaks = new long8[ring_num]; // this holds the ring number of each pixel in the map
		phi0 = new flt8[ring_num];
		dphi = new flt8[ring_num];


		precalculate_fourier_stuff(ring_num, *this, ThVals,ThBreaks, phi0,dphi);

		SStart=time(NULL);

		SHT_FFTW_COMPLEX *alm_fftw = new SHT_FFTW_COMPLEX[a.almsNum()];
		coordStruct coords;
		coords.nPix = pixNum();
		coords.nThetaVals = ring_num;
		coords.thetaVals=ThVals;
#ifdef GO128BIT
		coords.thetaBreaks=(long8*)ThBreaks;
#elif GO64BIT
		coords.thetaBreaks=(long8*)ThBreaks;
#else
		coords.thetaBreaks=(int*)ThBreaks;
#endif
		coords.phi0=phi0;
		coords.deltaPhi=dphi;
		coords.nGaps=(long)0;
		coords.gaps=NULL;
		coords.pixSize=(double)(4*PI/(double)pixNum());
		long pix_num=pixNum();
		flt8 * SHTin = new flt8[pixNum()]; for (i=0; i<pix_num;i++) { SHTin[i] = get_T(i);  }

		//
		// do SHT
		//
		msgs->say("doint  SHT:",Top);
		forwardSHT(SHTin,0,coords,lmax,alm_fftw);

		//
		// (de)smoothing in fourier space
		//
		msgs->say("de-smoothing in SH space:",Low);

		i=0;
		w=beamtf * smGbl * pixtf;
		for (m=-lmax;m<=lmax;m++) { 	  //start = clock();
			for (l=abs(m); l<=lmax; l++) {
				a[i]=mscsAlm( alm_fftw[i][0] / w(l), alm_fftw[i][1] / w(l) );
				i++;
			}
		}
//		a.printtxtAlms();
		destroyCoordsCPP(&coords);
		delete [] alm_fftw;
		delete [] SHTin;

		//converting to nested ordering
		/*     copy_map(mapRING,map,pix_num,NULL); // this conversion could be done just in here !!! but in case the conv_tabs are not loaded this is safer */
		// kill_map_space(mapRING);   delete mapRING;
		conv_ring2nest();
		eend = time(NULL); total_time = ((double)eend - (double)SStart) / 60.0;
		msgs->say("done, total time: "+msgs->toStr(total_time)+" [min]",High);   //   fprintf(tmpf,"%i %lf\n",l,total_time);

	}

	//************************************************************************

	if (method == 7) { // my own calculation only for the pixels in the map
		// this method uses precomputed Plms if required

		//     // setting graphical visualization of the progress
		//     //show_progress("start2");

		//     dOm = 4*PI/(double)pix_num;//*(1.0-f_sky);

		//     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*(nside-1)+1)
		//     printf("  -- number of rings in map= %li \n",ring_num);

		//     Ptab = new double[ring_num];
		//     ringtab = new long int[pix_num];
		//     ringpixtab = new long int[pix_num];
		//     COS = new double[ring_num];
		//     SIN = new double[ring_num];
		//     //gauss_sm = new double[lmax+1]; gauss_sm[0] = (double)smooth_map; gauss_sm[1] = sigmaB; // the fourier coefficients of gauss function and just temporairly keep the information about smoothing details in the array have it in the subroutine precalculate_fourier_stuff  !!! WARNING - this is dangerous if lmax = 0;

		//     if (map_ordering == 1) { conv_nest2ring(map);    } // if nested - convert to ring

		//     precalculate_fourier_stuff(COS,SIN,map,ring_num,ringtab,ringpixtab);
		//     rmax = 2*nside-1; // maximal ring number in northern hemisphere

		//     /*     forward Fourier transformation  itself   */
		//     Start=clock();
		//     for (l=0;l<=lmax;l++) { start = clock();
		//       for (m = 0; m<=l; m++) { // from 0 because we have real data. alms yield the reality condition
		// 	if (strcmp(plms_file,"") != 0) {  calculate_Plm_harmonics(ring_num,Ptab,Plms);} // -1 because we want the number from 0 not the amount
		// 	else calculate_Plm_harmonics(COS,ring_num,rmax,Ptab,l,m);

		// 	for (i=0; i<pix_num;i++) {
		// 	  Plm = Ptab[ringtab[i]];
		// 	  if ( Plm == 0 ) {  i = ringpixtab[i];} //printf("i=%li\n",i); // go to the next ring
		//        	  else {
		// 	    phi = map->n[i].l;
		// 	    Ylm_re = Plm * cos((double)m*phi); // we take here the real part of a_lm * Y_lm
		// 	    almR += Ylm_re * map->T[i] * map->m[i] / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l));
		// 	    Ylm_im = -Plm * sin((double)m*phi); // we take here the real part of a_lm * Y_lm
		// 	    almI +=  Ylm_im * map->T[i] * map->m[i] / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l));
		// 	  }
		// 	}

		// 	almR *= dOm; 	almI *= dOm;
		// 	i = alm2num(l,m);
		// /* 	printf("-----------> alms2num: %li\n",i);  alm[i].R = 12.12; exit(0);i=210;  */
		// 	alm[i].R = almR;  alm[i].I = almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0);
		// 	i = alm2num(l,-m); alm[i].R = pow(-1,(double)m)*almR; alm[i].I = pow(-1,(double)(m+1))*almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0);
		// 	//if (m==0) { printf("almR = %lE, almI = %lE\n",almR,almI); }
		// 	almR = almI = 0;

		//       }

		//       /*    normalization  for each l  */
		// /*       if (l == 2) {norm = sqrt(C_2WMAP/calculate_single_C_l(2,1)); printf("generated alms normalization factor: %lf\n",norm);} */
		// /*       for (m=-l; m<=l;m++) {  */
		// /* 	i = alm2num(l,m);  */
		// /* 	alm[i].R = norm * alm[i].R;  */
		// /* 	alm[i].I = norm * alm[i].I;  */
		// /* 	alm[i].D = sqrt(pow(alm[i].R,2) + pow(alm[i].I,2)); */
		// /* 	alm[i].P = cpeds_cart2sph(1,alm[i].R,alm[i].I,0); */
		// /*       } */

		//       if (save_partial == 1) {  // this saves the individal multipole set of alms and a sum of them upto l
		// 	tmpalmnum  = alm2num(l,l);
		// 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,4,l); savetxtF(flatmap_file_name,0); // save one multipole
		// 	tmpl = lmax; lmax = l;
		// 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,0);  // save sum of multipoles upto l'th
		// 	lmax = tmpl;
		//       }

		//       end = clock();         total_time_p=total_time;       total_time = ((double)end - (double)start) / (double)CLOCKS_PER_SEC;
		//       cpu_time_used = ((double)end - (double)Start) / (double)CLOCKS_PER_SEC;  tot = 0.5*(total_time-total_time_p)*pow((double)lmax,2);
		//       //tot = 0.5*()*pow((double)lmax,2);
		//       printf("l= %li time per l: %.2lf [s] TOT(l=%li)=%.2lf [h] ETA = %0.2lf [s]\n",l,total_time,lmax,tot/3600,cpu_time_used);
		//       //fprintf(tmpf,"%li %lf\n",l,total_time); vis_xtab[l] = (float)l; vis_ytab[l] = (float)total_time;
		//       //vis_x = (float)l; vis_y = (float)total_time; vis_iter++;      show_progress("xytabFourier2"); vis_x_prev=vis_x;       vis_y_prev=vis_y;
		//     }
		//     delete SIN; delete COS;   delete ringtab; delete ringpixtab; delete Ptab; //delete gauss_sm;

		//     //converting to nested ordering
		//     conv_ring2nest(map);
		//     //fclose(tmpf); //exit(0);
		//     //cpgclos();
	}

	/*   //-------------------------------------------------------------------------------------------------------------------------- */

	if (method == 8) { // this is my own but improved calculation by the partial separation of variables. Instead of N^4 complexity this should have 2*N^3 of so.
		// the accuracy level should be the same as for method 6. no approximations are made.
		// this method uses precomputed Plms if required


		// /*     tmpf = fopen(speed_log_file,"w"); */

		//     th0 = 0; phi0 = 0;
		// /*     dth = resolution*cpeds_pix_size_healpix(nside); */
		//     printf("  -- approximated pixel size [deg]: %lf\n",180/PI * cpeds_pix_size_healpix(nside));
		//     printf("  -- **** the resolution for this run is always the size of the pixel in the map ****\n");
		//     printf("  -- number of iterations per map per m = %li \n",pix_num);
		//     printf("  -- map (de)smoothing in fourier sapce is: "); if (smooth_map == 1) { printf("ON, gaussian beam (FWHM) [deg]: %lE\n",180/PI*thFWHM); } else printf("OFF\n");
		//     dOm = 4*PI/pix_num;//*(1.0-f_sky); // assumption that the each pixel in the map coners that save area; this is that case for the Healpix but  generally not true.

		//     /* prepare the directon map in a requested pixelization system for this accuracy */

		//     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*(nside-1)+1)
		//     printf("  -- number of rings in map= %li \n",ring_num);

		//     mapmasked = new double[pix_num];
		//     Ptab = new double[ring_num];
		//     ringtab = new long int[pix_num];
		//     ringpixtab = new long int[pix_num];
		//     COS = new double[ring_num];
		//     SIN = new double[ring_num];
		//     //gauss_sm = new double[lmax+1]; gauss_sm[0] = (double)smooth_map; gauss_sm[1] = sigmaB; // the fourier coefficients of gauss function and just temporairly keep the information about smoothing details in the array have it in the subroutine precalculate_fourier_stuff  !!! WARNING - this is dangerous if lmax = 0;
		//     count = ring_num*((long)lmax+1);
		//     b = new a_lm[count];    for (i=0;i<count;i++) { b[i].R = 0; b[i].I = 0; } count = 0; // here \Sum_l a_lm*P_lm*B_l is kept for all th knots and m values from 0 to lmax

		//     if (map_ordering == 1) { conv_nest2ring(map);    } // if nested - convert to ring
		//     for (i=0; i<pix_num; i++) { mapmasked[i] = map->T[i] * map->m[i]; } // make masked version of the map from optimalization reasons

		//     precalculate_fourier_stuff(COS,SIN,map,ring_num,ringtab,ringpixtab);
		//     rmax = 2*nside-1; // maximal ring number in northern hemisphere

		//     /*     forward Fourier transformation  itself   */
		//     SStart=time(NULL);

		//     printf("  -- doing FFT in phi direction (calculating b_m(theta) coefficients)\n");
		//     count = 0;
		//     for (m=0; m<=lmax; m++) {
		//       for (i=0; i<pix_num; i++) { // for each ring actually - look "here"
		// 	ring = ringpixtab[i]+1;
		// 	for (j=i; j< ring; j++) { // integration over phi
		// 	  phi = map->n[j].l;
		// 	  b[count].R += mapmasked[j] * cos( (double)m * phi);       b[count].I += -mapmasked[j] * sin( (double)m * phi);
		// 	}
		// 	count++;
		// 	i = j-1; // "here" - jump with i to the next ring.
		//       }
		//     }

		//     eend = time(NULL);  total_time = ((double)eend - (double)SStart) / 60.0;
		//     printf("  -- time for b_l(th): %.2lf [min]\n",total_time); //     fprintf(tmpf,"%i %lf\n",l,total_time);

		//     printf("  -- doing integration over theta\n");
		//     for (l=0;l<=lmax;l++) { //start = clock();
		//       count = 0;
		//       for (m = 0; m<=l; m++) { // from 0 because we have real data. alms yield the reality condition
		// 	if (Plms_from_file == 1) {  calculate_Plm_harmonics(ring_num,Ptab,Plms); } else { calculate_Plm_harmonics(COS,ring_num,rmax,Ptab,l,m); }
		// 	for (i=0; i<pix_num;i++) { // integration over theta
		// 	  Plm = Ptab[ringtab[i]];
		// 	  if ( Plm == 0 ) {  i = ringpixtab[i]; }//printf("ring=%li l=%li m=%li\n",i,l,m);} // go to the next ring
		//        	  else {
		// 	    almR += b[count].R * Plm;    almI += b[count].I * Plm;
		// 	  }
		// 	  i = ringpixtab[i];
		// 	  count++;
		// 	}

		// 	almR = almR * dOm / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l)); 	almI = almI * dOm / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l));
		// 	i = alm2num(l,m); alm[i].R = almR; alm[i].I = almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0);
		// 	if (m != 0) { i = alm2num(l,-m); alm[i].R = pow(-1.0,(double)m)*almR; alm[i].I = pow(-1.0,(double)(m+1))*almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0); } // this is stupid - why calculate twice almost the same thing ?
		// 	almR = almI = 0;
		//       }

		//       if (save_partial == 1) {  // this saves the individal multipole set of alms and a sum of them upto l
		// 	tmpalmnum  = alm2num(l,l);
		// 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,4,l); savetxtF(flatmap_file_name,0); // save one multipole
		// 	tmpl = lmax; lmax = l;
		// 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,0);  // save sum of multipoles upto l'th
		// 	lmax = tmpl;
		//       }
		//     }

		//     delete SIN; delete COS;   delete ringtab; delete ringpixtab; delete Ptab; //delete gauss_sm;
		//     delete b; delete [] mapmasked;
		//     eend = time(NULL); total_time = ((double)eend - (double)SStart) / 60.0;
		//     printf("total time: %.2lf [min] \n",total_time);   //   fprintf(tmpf,"%i %lf\n",l,total_time);
		//     //converting to nested ordering
		//     conv_ring2nest(map);
		// /*     fclose(tmpf); //exit(0); */

	}


	/*   //-------------------------------------------------------------------------------------------------------------------------- */





	if (method == 9) { // this is like in 8 but uses FFTW for integration over phi
		/*     // this method uses precomputed Plms if required */
		/*     printf("\nCalculating Fourier transform: method %i\n",how); */

		/*     if (strcmp(plms_file,"") != 0) {  Plms = new Plms_tab(nside,lmax,plms_file); Plms_from_file = 1; } else { Plms_from_file = 0; } */
		/*     //------ space allocation in case it's not allocated ---------- */
		/*     makekill_space_manager("make","F",-1); */

		/*     tmpf = fopen(speed_log_file,"w"); */
		/*     // setting graphical visualization of the progress */
		/*     //show_progress("start2"); */

		/*     th0 = 0; phi0 = 0;  */
		/*     dth = resolution*cpeds_pix_size_healpix(nside); */
		/*     printf("  -- approximated pixel size [deg]: %lf, calculation resoltion = %lf [deg]\n",180/PI * cpeds_pix_size_healpix(nside),180/PI * dth); */
		/*     printf("  -- **** the resolution for this run is always the size of the pixel in the map ****\n"); */
		/*     printf("  -- number of iterations per map per m = %li \n",pix_num); */
		/*     printf("  -- map (de)smoothing in fourier sapce is: "); if (smooth_map == 1) { printf("ON, gaussian beam (FWHM) [deg]: %lE\n",180/PI*thFWHM); } else printf("OFF\n"); */
		/*     dphi=dth; */
		/*     //dOm = dth * dphi; */
		/*     dOm = 4*PI/pix_num;//\*(1.0-f_sky); // assumption that the each pixel in the map coners that save area; this is that case for the Healpix but  generally not true. */

		/*     /\*     seed the coordinates in the map -- if not yet loaded   *\/ */
		/*     if (coord_loaded == 0) { set_map_coord(0,0);    } */
		/*     /\* prepare the directon map in a requested pixelization system for this accuracy *\/ */

		/*     ring_num = cpeds_get_ring_num_healpix(nside); // 2(2*(nside-1)+1) */
		/*     printf("  -- number of rings in map= %li \n",ring_num); */

		/*     mapmasked = new pixel[pix_num]; */
		/*     Ptab = new double[ring_num];  */
		/*     ringtab = new long int[pix_num]; */
		/*     ringpixtab = new long int[pix_num]; */
		/*     COS = new double[ring_num]; */
		/*     SIN = new double[ring_num]; */
		/*     //gauss_sm = new double[lmax+1]; gauss_sm[0] = (double)smooth_map; gauss_sm[1] = sigmaB; // the fourier coefficients of gauss function and just temporairly keep the information about smoothing details in the array have it in the subroutine precalculate_fourier_stuff  !!! WARNING - this is dangerous if lmax = 0; */
		/*     count = ring_num*((long)lmax+1);  */
		/*     b = new a_lm[count];    for (i=0;i<count;i++) { b[i].R = 0; b[i].I = 0; } count = 0; // here \Sum_l a_lm*P_lm*B_l is kept for all th knots and m values from 0 to lmax */

		/*     if (map_ordering == 1) { conv_nest2ring(map);    } // if nested - convert to ring */
		/*     for (i=0; i<pix_num; i++) { mapmasked[i].T = map[i].T * map[i].m; } // apply the mask to the map */

		/*     precalculate_fourier_stuff(COS,SIN,map,ring_num,ringtab,ringpixtab); */
		/*     rmax = 2*nside-1; // maximal ring number in northern hemisphere */

		/*     /\*     forward Fourier transformation  itself   *\/ */
		/*     SStart=time(NULL); */

		/*     printf("  -- doing FFT in phi direction (calculating b_m(theta) coefficients) using FFTW \n"); */
		/*     count = 0; */

		/* /\* REPLACE THIS BIT *\/ */


		/* /\*  FFTW phi direction transform  *\/ */

		/* /\*     fftw_complex *inf, *outf, *inb, *outb; *\/ */
		/* /\*     fftw_plan pf,pb; *\/ */
		/* /\*     long xpix = 2*nside+1; // this is to be changed ^^^ above if summing not to l = 2*nside then still the b[] must be of size = ring_num*(4*nside) and the modes with l>lmax should be zeroed *\/ */
		/* /\*     inf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); *\/ */
		/* /\*     inb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); *\/ */
		/* /\*     outf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); *\/ */
		/* /\*     outb = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (xpix)); *\/ */
		/* /\*     pf = fftw_plan_dft_1d(lmax, inf, outf, FFTW_FORWARD, FFTW_MEASURE); *\/ */
		/* /\*     pb = fftw_plan_dft_1d(lmax, inb, outb, FFTW_BACKWARD, FFTW_MEASURE); *\/ */
		/* /\*     count2 = 0;  *\/ */
		/* /\*     for (i=0; i<ring_num;i++) { *\/ */
		/* /\*       count = i; *\/ */
		/* /\*       b0th = b[count].R; *\/ */
		/* /\*       for (m=0; m<xpix; m++) {	 *\/ */
		/* /\* 	inb[m][0] = b[count].R; inb[m][1] = b[count].I;  *\/ */
		/* /\* 	inf[m][0] = b[count].R; inf[m][1] = -b[count].I; *\/ */
		/* /\* 	count +=ring_num;   *\/ */
		/* /\*       } // load array of b_m(th) for a given theta *\/ */
		/* /\*       fftw_execute(pf); // do the transform *\/ */
		/* /\*       fftw_execute(pb); // do the transform *\/ */
		/* /\*       pixinring = cpeds_get_pixnum_in_ring_healpix(nside,i); *\/ */
		/* /\*       for (m=0; m<pixinring; m++) {	 *\/ */
		/* /\* 	phi = mapRING[count2].n.l; x = (long)(phi*(xpix)/twoPI); *\/ */
		/* /\* 	mapRING[count2].T = outb[x][0] + outf[x][0] - b0th; count2++;  *\/ */
		/* /\*       } *\/ */
		/* /\*     } *\/ */
		/* /\*     fftw_destroy_plan(pf); *\/ */
		/* /\*     fftw_destroy_plan(pb); *\/ */
		/* /\*     fftw_free(inb); fftw_free(outb); *\/ */
		/* /\*     fftw_free(inf); fftw_free(outf); *\/ */

		/* /\*  ----------------------------  *\/ */


		/*     for (m=0; m<=lmax; m++) {  */
		/*       for (i=0; i<pix_num; i++) { // for each ring actually - look "here" */
		/* 	ring = ringpixtab[i]+1; */
		/* 	for (j=i; j< ring; j++) { // integration over phi */
		/* 	  phi = map[j].n.l; */
		/* 	  b[count].R += mapmasked[j].T * cos( (double)m * phi);       b[count].I += -mapmasked[j].T * sin( (double)m * phi); */
		/* 	} */
		/* 	count++; */
		/* 	i = j-1; // "here" - jump with i to the next ring. */
		/*       } */
		/*     } */

		/* /\* REPLACE THIS BIT *\/ */

		/*     eend = time(NULL);  total_time = ((double)eend - (double)SStart) / 60.0;   */
		/*     printf("  -- time for b_l(th): %.2lf [min]\n",total_time); //     fprintf(tmpf,"%i %lf\n",l,total_time);  */

		/*     printf("  -- doing integration over theta\n"); */
		/*     for (l=0;l<=lmax;l++) { //start = clock();  */
		/*       count = 0; */
		/*       for (m = 0; m<=l; m++) { // from 0 because we have real data. alms yield the reality condition */
		/* 	if (Plms_from_file == 1) {  calculate_Plm_harmonics(ring_num,Ptab,Plms); } else { calculate_Plm_harmonics(COS,ring_num,rmax,Ptab,l,m); } */
		/* 	for (i=0; i<pix_num;i++) { // integration over theta */
		/* 	  Plm = Ptab[ringtab[i]]; */
		/* 	  if ( Plm == 0 ) {  i = ringpixtab[i];} //printf("i=%li\n",i); // go to the next ring */
		/*        	  else {   */
		/* 	    almR += b[count].R * Plm;    almI += b[count].I * Plm; */
		/* 	  } */
		/* 	  i = ringpixtab[i]; */
		/* 	  count++; */
		/* 	} */

		/* 	almR = almR * dOm / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l)); 	almI = almI * dOm / (bl_tab[8]->get_wl(l) * bl_tab[9]->get_wl(l)); */
		/* 	i = alm2num(l,m); alm[i].R = almR; alm[i].I = almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0); */
		/* 	if (m != 0) { i = alm2num(l,-m); alm[i].R = pow(-1.0,(double)m)*almR; alm[i].I = pow(-1.0,(double)(m+1))*almI; alm[i].D = sqrt(pow(almR,2) + pow(almI,2)); alm[i].P = cpeds_cart2sph(1,almR,almI,0); } */
		/* 	almR = almI = 0; */
		/*       } */

		/*       if (save_partial == 1) {  // this saves the individal multipole set of alms and a sum of them upto l */
		/* 	tmpalmnum  = alm2num(l,l); */
		/* 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,4,l); savetxtF(flatmap_file_name,0); // save one multipole */
		/* 	tmpl = lmax; lmax = l; */
		/* 	sprintf(flatmap_file_name,"%s",partial_files_prefix);  savebinF(flatmap_file_name,0);  // save sum of multipoles upto l'th */
		/* 	lmax = tmpl; */
		/*       } */
		/*     } */

		/*     delete SIN; delete COS;   delete ringtab; delete ringpixtab; delete Ptab; //delete gauss_sm;  */
		/*     delete b; delete mapmasked; */
		/*     eend = time(NULL); total_time = ((double)eend - (double)SStart) / 60.0;   */
		/*     printf("total time: %.2lf [min] \n",total_time);   //   fprintf(tmpf,"%i %lf\n",l,total_time);  */

		/*     //converting to nested ordering */
		/*     conv_ring2nest(map);  */
		/*     fclose(tmpf); //exit(0); */

	}






	if (was_masked == false) { makekill_space_manager("kill","m",-1); }
	// kill_window_functions();

	// if (how !=1) if (strcmp(plms_file,"") != 0) {  delete Plms; }
	// alms_loaded = 1;
	return a;
}




// dedicated to the 1st method of SHT
void mscsMap::precalculate_fourier_stuff(long8 ring_num, const mscsMap& mapRING, flt8 *ThVals,long8 *ThBreaks, flt8* phi0, flt8 * dphi) {
	flt8 th,phi;
	long8 i,j,k;

	// if (pix_system == 1) { // healpix
	j=0;
	for (i=0; i<ring_num;i++) {
		th = PIsnd-mapRING.get_C(j).lat(); phi = mapRING.get_C(j).lon(); k=cpeds_get_pixnum_in_ring_healpix(nside(),i);

		ThVals[i] = th;
		ThBreaks[i]=j;
		phi0[i] = phi;
		dphi[i] = twoPI/(flt8)k;
		//printf("ring: %li, thvals: %lE, thbreaks: %li, phi0: %lE dphi: %lE pixnum: %li pixinring: %li\n",i,PI180inv*th,j,PI180inv*phi0[i],PI180inv*dphi[i],j,k);
		j+=k;
	}
	// }
}


























// //************************************************************************
// // derives the normalized Plm function for Fourier transformations also for m < 0
// // this function makes use of the Klm_prev value if possible for faster calculation using recurence relation
// // Klm is bullshit here - it'a a mistake and shouldn't be calculated.
// // P_l-m = (-1)^m (l-m)!/(l+m)! * P_lm but since we calculate the normalized values it's almost the normalized spherical harmonic
// // without the e^i m phi part so the relation is Y_l-m = (-1)^m Y_lm*
// // generally for the real transforms the negative ms need not to be calculated, hence separate function will be used for these calculations (method 7)
// double mscsMap::Plm_harmonic(long l, long m, double x) {
//   double calc_threshold = 1e-240;
//   double Plm=0, klm;
//   gsl_sf_result res;
//   int status;

// /*   if (l < 100) { // this is temporary arbitrary threshold */
// /*     if (m >= 0 ) {  Plm = gsl_sf_legendre_sphPlm(l,m,x); }  */
// /*     else { m = -m; Plm = pow(-1,(double)m) * klm * gsl_sf_legendre_sphPlm(l,m,x);} */
// /*     //else { m = -m; Plm = pow(-1,(double)m) * gsl_sf_legendre_sphPlm(l,m,x); } // this is because th gsl_sf_legendre_sphPlm - is already normalized */
// /*   } */
// /*   else { */
//   if (m >= 0 ) {
//     status = gsl_sf_legendre_sphPlm_e(l,m,x,&res);
//     if (status == 0) { if (fabs(res.val) < calc_threshold) { Plm = 0; } else { Plm = res.val; }}
// /*     printf("----------1 sphplm = %lE ",res.val); */
//   }
//   else {     // for l >= 100
// /*     if (l > 170) { Plm = 0; }  */
// /*     else {	 */
//     m = -m;
//     klm = 1;//Klm(l,m);
//     status = gsl_sf_legendre_sphPlm_e(l,m,x,&res);
// /*     printf("----------2 sphplm = %lE ",res.val); */
//     if (status != 0) { Plm = 0; }
//     else {
//       Plm = pow(-1,(double)m) * klm * res.val;
//       if ( cpeds_isnan(Plm) == 1 ) { Plm = 0; } else {
// 	if ( fabs(Plm) < calc_threshold) {Plm = 0; }
//       }
//     }
//   }

//   return Plm;
// }

// // Plm_positive (to be exact non negative) this calculates the same as above but it's dedicated to m>0 only optimazed for speed
// double mscsMap::Plmp_harmonic(long l, long m, double x) {
//   double calc_threshold = 1e-240;
//   double Plm;//, klm;
//   gsl_sf_result res;
//   int status;

//   status = gsl_sf_legendre_sphPlm_e(l,m,x,&res);
//   if (status == 0) { if (fabs(res.val) < calc_threshold) { Plm = 0.0; } else { Plm = res.val; }} else { Plm = 0; }
//   return Plm;
// }


// //************************************************************************
// // recurence function for deriving the Klm factor: Klm = (l-m)! / (k+m)!
// // rewritten into iterative form since it's much simpler
// double mscsMap::Klm(long l, long m) {
//   double klm=0,tmp=1;
//   long i;

//   if (ok_prev) {
//     if ((l_prev == l) && (m_prev == m)) { klm = Klm_prev; }
//     if (l_prev < l) { // if we go down
//       for (i=l_prev+1;i<=l;i++) { tmp = tmp * ((double)(i-m_prev)/(double)(i+m_prev)); } klm = tmp * Klm_prev;  Klm_prev = klm;  l_prev = l; // first in ls
// /*       printf("down"); */
//     }

//     tmp=1;
//     if (m_prev < m) { // if we go right
//       for (i=m_prev+1;i<=m;i++) { tmp = tmp / (double)( (l_prev+i)*(l_prev-i+1) ); } klm = tmp * Klm_prev; Klm_prev = klm; m_prev = m;
// /*       printf("right"); */
//     }
//     if (m_prev > m) { // if we go left
//       for (i=m_prev-1;i>=m;i--) { tmp = tmp * (double)( (l_prev-i)*(l_prev+i+1) ); } klm = tmp * Klm_prev; Klm_prev = klm; m_prev = m;
// /*       printf("left"); */
//     }

//     tmp=1;
//     if (l_prev > l) { // if we go up
//       for (i=l_prev-1;i>=l;i--) { tmp = tmp * ((double)(i+1+m_prev)/(double)(i+1-m_prev)); } klm = tmp * Klm_prev;  Klm_prev = klm;  // second in ls
// /*       printf("up"); */
//     }

//   }
//   else { klm = cpeds_factorial(l-m)/cpeds_factorial(l+m);  ok_prev = 1; }//printf("**oszustwo2 klm=%lE",klm);}
//   Klm_prev = klm;  l_prev = l; m_prev = m;
//   if ( cpeds_isnan(klm) == 1) { ok_prev = 0; }
// /*   printf("**rekurencja l=%i m=%i, Klm=%lE Klmprev=%lE lprev=%i mprev=%i rekdepth=%i\n",l,m,klm,Klm_prev,l_prev,m_prev,rek_depthl); */
// /*   if (ok_prev) { */
// /*     if (l_prev < l) { //klm = cpeds_factorial(l-m)/cpeds_factorial(l+m); } // if we go down */
// /*       printf("**licze l"); */

// /*       if (l_prev + 1 < l ) {   rek_depthl++; Klm(l-1,m_prev);       rek_depthl--;} */
// /*       if (l_prev - 1 > l ) {   rek_depthl++; Klm(l+1,m_prev);       rek_depthl--;} */

// /*       if (l_prev + 1 == l) { klm = Klm_prev * ( (double)(l-m_prev)/(double)(l+m_prev) ) ;  printf("**jestem lp+1=m Klm=%lE Klmprev = %lE l=%i m=%i lprev = %i mprev=%i\n ",klm,Klm_prev,l,m,l_prev,m_prev);} */
// /*       if (l_prev - 1 == l) { klm = Klm_prev * ((double)(l+1+m_prev)/(double)(l+1-m_prev) ); printf("**jestem lp-1=m Klm=%lE Klmprev = %lE l=%i m=%i\n ",klm,Klm_prev,l,m);} */

// /*       Klm_prev = klm;  */
// /*     }  */

// /*     if ((rek_depthl == 1) && (m_prev != m)) { l_prev=l; // we go horizontally in ms. */

// /*       if (m_prev + 1 < m ) { printf("***rek-"); Klm(l_prev,m-1); } */
// /*       if (m_prev - 1 > m ) { printf("***rek+"); Klm(l_prev,m+1); } */

// /*       if (m_prev + 1 == m) { klm = Klm_prev / (double)( (l_prev+m)*(l_prev-m+1) ); printf("**jestem mp+1=m Klm=%lE Klmprev = %lE l=%i m=%i\n ",klm,Klm_prev,l,m);} */
// /*       if (m_prev - 1 == m) { klm = Klm_prev * (double)( (l_prev-m)*(l_prev+m+1) ); printf("**jestem mp-1=m Klm=%lE Klmprev = %lE l=%i m=%i\n ",klm,Klm_prev,l,m);} */
// /*       Klm_prev = klm;  */
// /*     }  */

// /*     if (l_prev > l) { //klm = cpeds_factorial(l-m)/cpeds_factorial(l+m); } // if we go up */
// /*       printf("**licze l"); */
// /*       if (l_prev + 1 < l ) {   rek_depthl++; Klm(l-1,m_prev);       rek_depthl--;} */
// /*       if (l_prev - 1 > l ) {   rek_depthl++; Klm(l+1,m_prev);       rek_depthl--;} */

// /*       if (l_prev + 1 == l) { klm = Klm_prev * ( (double)(l-m_prev)/(double)(l+m_prev) ) ;  printf("**jestem lp+1=m Klm=%lE Klmprev = %lE l=%i m=%i lprev = %i mprev=%i\n ",klm,Klm_prev,l,m,l_prev,m_prev);} */
// /*       if (l_prev - 1 == l) { klm = Klm_prev * ((double)(l+1+m_prev)/(double)(l+1-m_prev) ); printf("**jestem lp-1=m Klm=%lE Klmprev = %lE l=%i m=%i\n ",klm,Klm_prev,l,m);} */
// /*     } */


// /*     if ((rek_depthl == 1) && (m_prev == m) && (l_prev == l)) { return Klm_prev; } */
// /*   } */
// /*   else { klm = cpeds_factorial(l-m)/cpeds_factorial(l+m);  ok_prev = 1; printf("**oszustwo2 klm=%lE",klm);} */

// /*   Klm_prev = klm;  l_prev = l; m_prev = m; */
// /*   printf("**rekurencja wyjscie l=%i m=%i, Klm=%lE Klmprev=%lE lprev=%i mprev=%i\n",l,m,klm,Klm_prev,l_prev,m_prev); */
// /*   printf("klm=%lE ",klm); */
//   return klm;
// }


// //************************************************************************
// // derives the Nlm coefficient for Fourier transformations
// double mscsMap::Nlm_harmonic(long l, long m) {

//   //int i;
//   double lmm,lpm; // l minus m, l plus m

//   lmm = cpeds_factorial(l-m);
//   lpm = cpeds_factorial(l+m);
//   return sqrt((2*(double)l+1)/(4*PI) * lmm/lpm);
// }
// //************************************************************************
// // this function pixelates the values of given Plm function on the sphere in a given pizelization system but only for the
// void mscsMap::calculate_Plm_harmonics(double* COS, long int ring_num, long int rmax,double * Ptab, long l, long m) {
//   long int i;//,rmax;//, ring_num;
//   //bool *is;

//   //ring_num = cpeds_get_ring_num_healpix(nside)-1; // this is actuall ring_num -1 // this is the number of rings in the map -1
//   // using the PI/2 symmetry
//   if (pix_system == 1) { // healpix
//     //rmax = 2*nside-1; // maximal ring number in northern hemisphere
//     //for (i=0; i<rmax;i++) {
//     for (i=0; i<ring_num;i++) {
//       Ptab[i] = Plmp_harmonic(l,m,COS[i]);
//       //Ptab[i] = Plm_harmonic(l,m,COS[i]);
//     }
//       //Ptab[ring_num-i] = pow(-1,(double)(l+m))*Ptab[i]; // using the symetry PI/2 // this is under condition that the pix. system is symmetrical relative to the equator
//       //printf("*** harmonic calc: ring_num:%li, l=%i, m=%i, COS[%li] = %lE, Ptab[i] = %lE\n",ring_num,l,m,i,COS[i],Ptab[i]);

//     //Ptab[rmax] = Plmp_harmonic(l,m,COS[rmax]); // this is for the equator

//     // not using the symmetry
// /*     for (i=0; i<ring_num;i++) { Ptab[i] = Plm_harmonic(l,m,COS[i]);     } */
//   }
// }

// // this functions reads the next portion of the alms from a file into a memory
//  void mscsMap::calculate_Plm_harmonics(long int ring_num, double * Ptab, Plms_tab * Plms) {
//   long int i;//,rmax;//, ring_num;
//   //bool *is;

//   //ring_num = cpeds_get_ring_num_healpix(nside)-1; // this is actuall ring_num -1 // this is the number of rings in the map -1
//   // using the PI/2 symmetry
//   if (pix_system == 1) { // healpix
//     //rmax = 2*nside-1; // maximal ring number in northern hemisphere
//     //for (i=0; i<rmax;i++) {
//       for (i=0; i<ring_num;i++) {
// 	Ptab[i] = (*Plms).Pnext();
// 	//printf("**********Ptab[%li]=%lE\n",i,Ptab[i]);
//       }
//       //Ptab[ring_num-i] = pow(-1,(double)(l+m))*Ptab[i]; // using the symetry PI/2 // this is under condition that the pix. system is symmetrical relative to the equator
//       //printf("*** harmonic calc: ring_num:%li, l=%i, m=%i, COS[%li] = %lE, Ptab[i] = %lE\n",ring_num,l,m,i,COS[i],Ptab[i]);

//     //Ptab[rmax] = Plmp_harmonic(l,m,COS[rmax]); // this is for the equator

//     // not using the symmetry
// /*     for (i=0; i<ring_num;i++) { Ptab[i] = Plm_harmonic(l,m,COS[i]);     } */
//   }

// }

// // this functions reads the next portion of the alms from a file into a memory
// void mscsMap::calculate_Pmth_harmonics(long m, long skip, double * Ptab, Plms_tab * Plms) {
//   long l;//,rmax;//, ring_num;
//   //bool *is;

//   if (pix_system == 1) { // healpix
//     for (l=m; l<=lmax;l++) {       Ptab[l-m] = (*Plms).Pnext();    }
//     (*Plms).Pskip(skip);
//   }
// }

//  // this functions calculated the Plm(cos(th)) for a given m and theta for all l's
//  // it uses the "old" Ptab array because there is twice less multipoles then the rings in the map. Some part of the array is not used.
// void mscsMap::calculate_Pmth_harmonics(double * Ptab, long m, double costh) {
//   int l;
//   if (pix_system == 1) { // healpix
//     for (l=m; l<=lmax;l++) {
//       Ptab[l-m] = Plmp_harmonic(l,m,costh);
//     }
//   }
// }

// //************************************************************************
// // this function reads the previously precalculated Plm_harmonic value from a table Ptab for a given pixel number which is converted to a given ring number in the map and returns the right Plm value depending on the calculated ring number for a given pixel number
// // WHAT IS THIS ???!?!?!
// double mscsMap::get_Plm_harmonic(double* Ptab, long int pix) {
// /*   long int row,col,ring,region,ring_num,i; */


// /*   if (pix_system == 1) { // healpix */
// /*     region = cpeds_get_healpix_region(nside,pix); */
// /*     row=cpeds_get_healpix_rownum_in_region(nside,pix); */
// /*     col=cpeds_get_healpix_colnum_in_region(nside,pix); */
// /*     ring = (ring_num-1)/2-(row+col)-1; // ring numbering is from 0 to ring_num from north pole to south pole */
// /*     if ((region >=4) && (region < 8)) { ring = ring + nside; } */
// /*     if ((region >=8) && (region < 12)) { ring = ring + (ring_num+1)/2; } */
// /*   } */

// /*   return Ptab[ring]; */
//   return 0;
// }
// //************************************************************************
// // dedicated to the 7th method of SHT
// void mscsMap::precalculate_fourier_stuff(double* COS, double * SIN, pixel * mapRING,long int ring_num, long int * ringtab, long int *ringpixtab) {
//   double th;
//   long int ring,i,j,k;

//   if (pix_system == 1) { // healpix

//     for (i=0; i<ring_num;i++) { ringpixtab[i] = 0; }

//     ringpixtab[0] = 4-1; // -1 is from the fact that the for loop will add this one when passsing by :)

//     for (i=0; i<pix_num;i++) {
//       th = PIsnd-mapRING->n[i].b;
//       //ring = cpeds_get_ring_healpix_nest(nside,i);
//       ring = cpeds_get_ring_healpix_ring(nside,i);
//       COS[ring] = cos(th);
//       SIN[ring] = sin(th);
//       ringtab[i] = ring; // this table contains the number of the ring in which a given pixel lies
//       //thcoord[i] = PIsnd-mapRING->n[i].b;
//       j=0;    for (k=0;k<=ring;k++) { j=j+cpeds_get_pixnum_in_ring_healpix(nside,k); }
//       ringpixtab[i] = j-1; // this contains the nuber of the pixel which starts the next ring starting from ring 0 at i=0 to ring max// -1 is from the fact that the for loop will add this one when passing by :)
//       //printf("**** cos: ring_ring=%li ringpixtab[%li] = %li\n",ring,i,ringpixtab[i]);
//       //printf("ring=%li ringtab=%li, th=%lf cos(th)=%lf COS(ring)=%lf \n", ring,ringtab[i],th,cos(th),COS[ring]);
//     }
//   }
// }

// // dedicated to the 8,9th methods of SHT
// void mscsMap::precalculate_fourier_stuff(double *COS,pixel *mapRING,long int ring_num,long int *ringtab) {
//   double th;
//   long int ring,i,j,k;

//   if (pix_system == 1) { // healpix
//     for (i=0; i<pix_num;i++) {
//       th = PIsnd-mapRING->n[i].b;
//       //ring = cpeds_get_ring_healpix_nest(nside,i);
//       ring = cpeds_get_ring_healpix_ring(nside,i);
//       COS[ring] = cos(th);
//       ringtab[i] = ring;
//       j=0;    for (k=0;k<=ring;k++) { j = j+cpeds_get_pixnum_in_ring_healpix(nside,k); }
//     }
//   }
// }

// // dedicated to the 1st method of SHT
// void mscsMap::precalculate_fourier_stuff(long ring_num, pixel *mapRING, double *ThVals,long *ThBreaks, double* phi0, double * dphi) {
//   double th,phi;
//   long int ring,i,j,k;

//   if (pix_system == 1) { // healpix
//     j=0;
//     for (i=0; i<ring_num;i++) {
//       th = PIsnd-mapRING->n[j].b; phi = mapRING->n[j].l; k=cpeds_get_pixnum_in_ring_healpix(nside,i);

//       ThVals[i] = th;
//       ThBreaks[i]=j;
//       phi0[i] = phi;
//       dphi[i] = twoPI/(double)k;

//       j+=k;
//     }
//   }
// }

// //************************************************************************
// // precalculates the Legendre polynominals for the rings in the map according to the current nside and lmax values
// // the values are stored in a file where
// //in the first line nside is stored
// // in the second line the lmax is stored
// //next are binary written Plm(cos(th)) values from l=0...lmax, m=0,...l,Plm(cos(th(ring0))),...,Plm(cos(th(ring_last)))
// // coords must be loaded for this calculation - i.e. first execute set_map_coord
// void mscsMap::precalculate_Legendre_polynominals(strarg ordering) {
//   FILE *fplms=NULL;
//   filenamestr tmpch;
//   long l=0,m;
//   long int i,ringnum,pix;
//   //pixel * tmpmap = new pixel[pix_num]; // this is only to have coordinates th
//   double th,plm;
//   time_t sstart,eend;
//   double total_time=0;
// /*   int mapcoord; */
// //  FILE* fff;
//   double * PNlmresA,countD;
//   long count;

// /*   mapcoord=coord_loaded; */
//   printf("precalculating the associated Legendre polynominals: nside=%li, lmax=%li\n",nside,lmax);
// /*   map=tmpmap; */
//   set_map_coord(0,0);
//   conv_nest2ring(map);

//   if (strcmp(ordering,"lmth") == 0) {
//     sprintf(tmpch,"Plms-%li-%li",nside,lmax);
//     fplms = fopen(tmpch,"w");
//     fprintf(fplms,"%li\n",nside);   fprintf(fplms,"%li\n",lmax);

//     ringnum = cpeds_get_ring_num_healpix(nside);
//     //long dupa =0;
//     sstart = time(NULL);
//     for (l=0;l<=lmax;l++) {
//       for (m = 0; m<=l; m++) { pix=0;
// 	for (i=0;i<ringnum;i++) { //dupa++;
// 	  th=PIsnd-map->n[pix].b;
// 	  plm = Plmp_harmonic(l,m,cos(th));
// 	  fwrite(&plm,sizeof(plm),1,fplms);
// 	  //printf("dupa=%li, i=%li l=%i m=%i th=%.10lE cos(th) = %lE plm = %.10lE gsl=%lE\n",dupa,i,l,m,180/PI*th, cos(th),plm,gsl_sf_legendre_sphPlm(l,m,cos(th)));
// 	  pix = pix+cpeds_get_pixnum_in_ring_healpix(nside,i);
// 	}
//       }
//       eend = time(NULL);
//       total_time = (double)eend - (double)sstart;
//       printf("l=%li, time elapsed: %lf\r",l,total_time);
//     }
//   }

//   if (strcmp(ordering,"mthl") == 0) {
//     sprintf(tmpch,"Pmthls-%li-%li",nside,lmax);
//     fplms = fopen(tmpch,"w");
//     fprintf(fplms,"%li\n",nside);   fprintf(fplms,"%li\n",lmax);

//     ringnum = cpeds_get_ring_num_healpix(nside);
//     //long dupa =0;
//     sstart = time(NULL);
//     for (m=0;m<=lmax;m++) { pix = 0;
//       for (i=0;i<ringnum;i++) {
// 	th=PIsnd-map->n[pix].b;
// 	for (l=m; l<=lmax; l++) { //dupa++;
// 	  plm = Plmp_harmonic(l,m,cos(th));
// 	  fwrite(&plm,sizeof(plm),1,fplms);
// 	  //printf("dupa=%li, i=%li l=%i m=%i th=%.10lE plm = %.10lE gsl=%lE\n",dupa, i,l,m,180/PI*th,plm,gsl_sf_legendre_sphPlm(l,m,cos(th)));
// 	}
// 	pix = pix+cpeds_get_pixnum_in_ring_healpix(nside,i);
//       }
//       eend = time(NULL);
//       total_time = (double)eend - (double)sstart;
//       printf("l=%li, time elapsed: %lf [min]\r",m,total_time/60.0);
//     }
//   }

//   if (strcmp(ordering,"thml") == 0) {
//     sprintf(tmpch,"Pthmls-%li-%li",nside,lmax);
//     fplms = fopen(tmpch,"w");
//     fprintf(fplms,"%li\n",nside);   fprintf(fplms,"%li\n",lmax);

//     ringnum = cpeds_get_ring_num_healpix(nside);
//     sstart = time(NULL); pix = 0; PNlmresA = new double[lmax+1];
//     for (i=0;i<ringnum;i++) {
//       th=PIsnd-map->n[pix].b;
//       for (m=0;m<=lmax;m++) {
// 	count=0;
// 	for (l=m; l<=lmax; l++) {
// 	  plm = Plmp_harmonic(l,m,cos(th));
// 	  if (plm!=0) { PNlmresA[count]=plm; count++; }
// 	}
// 	countD=(double)count; fwrite(&countD,sizeof(countD),1,fplms);
// 	fwrite(PNlmresA,sizeof(plm),count,fplms);
//       } //printf("i=%li l=%i m=%i th=%.10lE plm = %.10lE gsl=%lE\n", i,l,m,180/PI*th,plm,gsl_sf_legendre_sphPlm(l,m,cos(th)));
//       pix += cpeds_get_pixnum_in_ring_healpix(nside,i);
//       eend = time(NULL);
//       total_time = (double)eend - (double)sstart;
//       printf("ring=%li, time elapsed: %lf [min]\r",i,total_time/60.0);
//     }
//     delete PNlmresA;
//   }

//   printf("l=%li, time elapsed: %lf\n",l,total_time);
//   fclose(fplms);

// }

