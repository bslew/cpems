#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
/* #include <cpgplot.h> */
/* #include <fitsio.h> */
/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */
/* #include <gsl/gsl_vector.h> */
/* #include <gsl/gsl_sf_legendre.h> */
/* #include <gsl/gsl_sf_result.h> */
/* #include <gsl/gsl_cdf.h> */
/* #include <gsl/gsl_errno.h> */
#include "Mscs-map.h"
//#include "cpeds-consts.h"
//#include "cpeds-math.h"


//************************************************************************
// calculates the number of points below some temperature threshold
void map_class::calculate_minkowski_area(double sigmaDA, long thres_num, double th_min, double th_max) {
  double mu,dmu,minTloc,maxTloc;
  int nn;
  minkowski_level_num_area = thres_num;
  if (th_min == th_max) { calculate_minT();   calculate_maxT(); minTloc = minT; maxTloc = maxT; } else { minTloc = th_min; maxTloc = th_max; }
  printf("calculating minkowski functional on the map: area with temeprature level thresholds:%i \n",minkowski_level_num_area);
  printf("Npix = %li, imin=%li, imax=%li, Tmin=%f, Tmax=%f\n\n",pix_num,iminT,imaxT,minT,maxT);
  mu=minTloc; dmu = (maxTloc-minTloc)/(double)minkowski_level_num_area;
  nn=0;

  for (nn=0;nn<thres_num;nn++) { 
/*   while (mu < maxTloc) {  */
    minkowski_area[nn][0] = mu;
    minkowski_area[nn][1] = ((double)pix_num-Ngrmu(mu,0.0,sigmaDA))/(double)pix_num;
/*     printf("mu=%f, Ngrmu=%f\n",mu,minkowski_area[nn][1]/(double)pix_num); */
    mu=mu+dmu; //nn++;
  }

/*   double * tmptab = new double[minkowski_level_num_area]; for (nn=0;nn<minkowski_level_num_area;nn++) {tmptab[nn] = minkowski_area[nn][1];} */
/*   double * params = new double[12]; */
/*   params[0] = 1E-1; params[1] = 1e+1; params[2] = 400; params[3] = 1000; */
/*   params[4] = 1e-7; params[5] = 1e-3; params[6] = 400; params[7] = 1000; */
/*   params[8] = 0.7*minT; params[9] = 0.7*maxT; params[10] = 1000; params[11] = 1000; */
/*   fit_function_3parameters_by_MCMC(2,(long int)minkowski_level_num_area,tmptab,minT,maxT,2,20,params); // looking for hypothesis to check with Pearson X^2 test and chi-squared distribution */
/*   delete tmptab; delete params; */
}

//************************************************************************
// calculates the number of points between some temperature interval defined by (Tmax-Tmin)/thres_num
// calculates the minkowski curcumference statistics on a temperature map with given numbers of observations of each pixel.
// the C(nu) statistics is performed according to the following formula:
// C(nu) = [ N(nu>) - N( (nu+dnu)>) ]/ [ N(0>) - N(dnu>) ] where N - is the number returned by the function Ngrnu - 
// the number of points with temperatures greater then some temperature threshold.
// sigmaDA === sigma0, thres_num  - number of different temperature levels at which calculate the circumference
void map_class::calculate_minkowski_circ(double sigmaDA, long thres_num, double th_min, double th_max) {
  double mu,dmu, circat0,nsigma=1.0,minTloc,maxTloc;
  int nn;
  minkowski_level_num_circ = thres_num;
  if (th_min == th_max) { calculate_minT();   calculate_maxT(); minTloc = minT; maxTloc = maxT; } else { minTloc = th_min; maxTloc = th_max; }
  printf("calculating minkowski functional on the map: circumference with temeprature level thresholds:%i uncertainties at: %.0lf sigma \n",minkowski_level_num_circ,nsigma);
  printf("Npix = %li, imin=%li, imax=%li, Tmin=%f, Tmax=%f\n\n",pix_num,iminT,imaxT,minT,maxT);
  mu=minTloc; dmu = (maxTloc-minTloc)/(double)minkowski_level_num_circ;  nn=0;
  circat0 = Ngrmu(0.0,0.0,sigmaDA) - Ngrmu(dmu,0.0,sigmaDA);
  //circat0 = 1;
  for (nn=0;nn<thres_num;nn++) { 
/*   while (mu <= maxTloc) {  */
    minkowski_circ[nn][0] = mu; // temperature threshold
    minkowski_circ[nn][1] = (Ngrmu(mu,0.0,sigmaDA) - Ngrmu(mu+dmu,0.0,sigmaDA))/circat0; // circumference at that threshold
    if (sigmaDA > 0) {
      minkowski_circ[nn][2] = dmu/2; // lower error in x
      minkowski_circ[nn][3] = dmu/2; // upper error in x
      minkowski_circ[nn][4] = fabs(minkowski_circ[nn][1]-(Ngrmu(mu,-nsigma,sigmaDA) - Ngrmu(mu+dmu,-nsigma,sigmaDA))/circat0); // lower error in y
      minkowski_circ[nn][5] = fabs(minkowski_circ[nn][1]-(Ngrmu(mu,nsigma,sigmaDA) - Ngrmu(mu+dmu,nsigma,sigmaDA))/circat0); // upper error in y
    }

    //printf("mu=%f, Ngrmu=%f\n",mu,Ngrmutab[nn][2]/(double)Npix);
    mu=mu+dmu; //nn++;
  }

}
//************************************************************************
// this is EXPERIMENTALL -- BAD !!!
void map_class::calculate_minkowski_genus(double sigmaDA, long thres_num, double th_min, double th_max) {
/*   double mu,dmu; */
/*   int nn; */
/*   printf("calculating minkowski functional on the map: genus with temeprature level thresholds:%i \n",minkowski_level_num_genus); */
/*   printf("Npix = %li, imin=%li, imax=%li, Tmin=%f, Tmax=%f\n\n",pix_num,iminT,imaxT,minT,maxT); */

/*   mu=minT; dmu = (maxT-minT)/(double)minkowski_level_num_genus; */
/*   nn=0; */

/*   while (mu < maxT) {  */
/*     minkowski_genus[nn][0] = mu; */
/*     minkowski_genus[nn][1] = (Ngrmu(mu) - Nlemu(mu))/pix_num; */
/*     //printf("mu=%f, Ngrmu=%f\n",mu,Ngrmutab[nn][2]/(double)Npix); */
/*     mu=mu+dmu; nn++; */
/*   } */

  //fit_function_3parameters_by_MCMC(1,minkowski_level_num_circ);

}

//************************************************************************
// tests the 

void map_class::test_minkowski_circ() {
  long nn;
  long int k,i,bin_num,mini,maxi;
  double x,Px,Qx,min,max,bin,bin_num_real,one_le_alph,alph, circat0;
  double *data,*err;
  data = new double[minkowski_level_num_circ]; err = new double[minkowski_level_num_circ];
  for (nn=0;nn<minkowski_level_num_circ;nn++) {data[nn] = minkowski_circ[nn][1]; err[nn] = minkowski_circ[nn][4];}
  double * params = new double[12];
/*   params[0] = 1E-1; params[1] = 1e+1; params[2] = 400; params[3] = 1000; */
/*   params[4] = 1e-7; params[5] = 1e-3; params[6] = 400; params[7] = 1000; */
/*   params[8] = 0.7*minT; params[9] = 0.7*maxT; params[10] = 1000; params[11] = 1000; */

  params[0] = 1E-1; params[1] = 1e6; params[2] = 400; params[3] = 1000;
  params[4] = 1e-7; params[5] = 1e-3; params[6] = 400; params[7] = 1000;
  params[8] = 0.7*minT; params[9] = 0.7*maxT; params[10] = 1000; params[11] = 1000;
  
/*   params[0] = 1E0; params[1] = 1e0; params[2] = 400; params[3] = 1000; */
/*   params[4] = 1e-7; params[5] = 1e-3; params[6] = 400; params[7] = 1000; */
/*   params[8] = 0.0*minT; params[9] = 0.0*maxT; params[10] = 1000; params[11] = 1000; */

  fit_function_3parameters_by_MCMC(1,(long int)minkowski_level_num_circ,data,minT,maxT,2,30,params); // looking for hypothesis to check with Pearson X^2 test and chi-squared distribution
  delete params;
  //plot_minkowski_circ();  
  printf("testing sample for gausianity with chi-squared test: null hypotesis %lE*N(%lE,%lE)\n",mink_fit_info.circ.A,mink_fit_info.circ.m,mink_fit_info.circ.s);

  //for (i=0;i<k;i++) { data[i] = map[i].T; } // copy the temperature map - since some sorting will be done.
  //for (i=0;i<k;i++) { t[i] = minkowski_circ[i][1]; }
  k = (long int)minkowski_level_num_circ;
  cpeds_find_minmax_value(data,k,&min,&max,&mini,&maxi);

  x = cpeds_test_chisq(k,data,err,mink_fit_info.circ.A,mink_fit_info.circ.m,mink_fit_info.circ.s);

  //Px = cpeds_chisq_prob(x);
  //printf("sample size: %li, bin_num = %li, X^2 = reduced chisq: %lE, degr. of freedom=%lE\n\n",k,bin_num,x,bin_num_real-1);
  printf("sample size: %li, chisq^2 = reduced chisq: %lE, degr. of freedom=%li\n\n",k,x,k-3);

  cpeds_print_chisq_confidence(0.05,(double)(k-3),x);
  cpeds_print_chisq_confidence(0.1,(double)(k-3),x);
  cpeds_print_chisq_confidence(0.2,(double)(k-3),x);
  cpeds_print_chisq_confidence(0.33,(double)(k-3),x);
  cpeds_print_chisq_confidence(0.5,(double)(k-3),x);

  mink_fit_info.circ.mean = cpeds_mean_value(data,minkowski_level_num_circ);
  mink_fit_info.circ.skew = cpeds_skewness(data,minkowski_level_num_circ);
  mink_fit_info.circ.kurt = cpeds_kurtosis(data,minkowski_level_num_circ);
  delete data; delete err;
  printf("\nminkowski functional circumference:\n mean: %lE\n skewness: %lE\n kurthosis: %lE\n", mink_fit_info.circ.mean, mink_fit_info.circ.skew, mink_fit_info.circ.kurt);

}







//************************************************************************
// num - number of points for calculating and tabulating the gauss function
// f - function type | must be consistent with the tabulate_function function
// ndelta - accuracy parameter - specifies the half side of the grid to be calculated around given point in parameter space to decide in which direction go next
// chainnum2 - number of chains to be started in second pass ( after the order of magnitude of parameters is decided in pass1 )
// data - address for the table where the data for fitting are kept; to these data the function within a parameters freedom will be fitted
// Xmin,Xmax - the range of the independant variable in the data; within this range the function will be tabulated for least squares calculation
// params -  is the address for table of 4*parameters_number - here 3 - with the min and max ranges of the parameters to be fitted
// the convenetion is : [0] - 1st param min value [1] 1st param max value [2] - numbers of levels for 1st pass [3] - numbers of levels for 2nd pass and so on
// CAUTION !!! currently the third parameter is the only one that can take positive and enegative values and it's the only one that is not fitted in log space.
// so constrain it good
// TO BE CHANGED : change the names of the variables like Astuff ---> P1stuff etc...
// CONVENTION: A - P1, s - P1, m - P3
void map_class::fit_function_3parameters_by_MCMC(int f, long int num, double *data, double Xmin, double Xmax,double ndelta, long int chainnum2, double* params) {
  double *fit, *Xsqtmp;
  double * functab = new double[num];
  double A,dA ,Amin, Amax, Astart, Afit,levA;
  double m,dm, mmin, mmax, mstart, mfit, levm;
  double s,ds, smin, smax, sstart, sfit,levs, sfit1stpass,Afit1stpass,mfit1stpass;
  double Xsq_min, Xsq_min_prev;
  double tmpd;
  //double ndelta=3;
  long int cell_num,i1,i2,i3,side_num;
  long int gridsize;  
  double locAstart, locAend, locmstart, locmend, locsstart, locsend, locA, locm, locs;
  long int i,pass;
  long int imin,imax;
  double iter,iternum;
  long int chain_num,number_of_chains,chains_total;
  fit_info *fits;

  side_num=2*(long int)ndelta+1;
  gridsize = side_num*side_num*side_num; 
  double * Xsq = new double[gridsize];
  double direction_tab[gridsize][3]; // this is a table which holds the values of parameters of function (here 3 for now - for gauss funct.) for different grid cells numbered from 0 to gridsize-1; this is kept in order to pick the next step of the chain right away

  if (f == 1) { // fitting gauss for circumference
    printf("fitting gauss function\n");
    printf("  -- local grid size: cell_num = %li\n",gridsize);
    for (pass=1;pass<=2;pass++) {
      printf(" \n------------\n*** pass %li\n-----------\n",pass);
      if (pass == 1) { // find the right order of magnitude of sigma
/* 	smin = 1E-7; smax = 1E-3;  levs=400;   // half order of magnitude step with 40 levels */
/* 	Amin = 1E-1; Amax = 1E+1;  levA=400; */
/* 	mmin = 0.7*minT; mmax = 0.7*maxT; levm=100;  */
	Amin = params[0]; Amax = params[1]; levA = params[2]; // definition of ranges
	smin = params[4]; smax = params[5]; levs = params[6]; // definition of ranges
	mmin = params[8]; mmax = params[9]; levm = params[10]; // definition of ranges
      }
      if (pass == 2) { // find the right value within the found order of magnitude
	Amin = pow(10,log10(Afit1stpass)-dA); Amax = pow(10,log10(Afit1stpass)+dA); Astart = pow(10,Afit1stpass); levA=1000;  levA = params[3];
	smin = pow(10,log10(sfit1stpass)-ds); smax = pow(10,log10(sfit1stpass)+ds); sstart = pow(10,sfit1stpass); levs=1000;  levs = params[7]; 
	mmin = 0.7*minT; mmax = 0.7*maxT; levm=1000; 
	//mmin = pow(10,log10(mfit1stpass)-dm); mmax = pow(10,log10(mfit1stpass)+dm); mstart = pow(10,mfit1stpass); levs=1000; //-- old
	/*       mmin = params[8]; mmax = params[9]; levsm = params[10]; levm = params[11]; // definition of ranges */
      }

      printf("  -- Amin = %lE mmin = %lE, smin = %lE\n", Amin, mmin,smin);
      printf("  -- Amax = %lE mmax = %lE, smax = %lE\n", Amax, mmax,smax);
  
      for (i=0;i<num;i++) { functab[i] = data[i]; }
      if (pass == 1) {       number_of_chains = (long int)(5*log10(smax/smin)*log10(Amax/Amin)*log10(levm)); }
      if (pass == 2) {       number_of_chains = chainnum2; }
      printf("  -- number of chains: %li\n",number_of_chains);
      fits = new fit_info[number_of_chains];
      chain_num = 0;
      chains_total = 0;
      while (chain_num < number_of_chains) {
	//set the monte carlo markov chain start
	if (pass == 1) {  
	  sstart = pow(10,cpeds_random_uniform_number(log10(smin),log10(smax)));
	  Astart = pow(10,cpeds_random_uniform_number(log10(Amin),log10(Amax)));
	  //mstart = pow(10,cpeds_random_uniform_number(log10(mmin),log10(mmax)));
	} // using log uniform distribution
	if (pass == 2) {  
	  sstart = cpeds_random_uniform_number(smin,smax);  
	  Astart = cpeds_random_uniform_number(Amin,Amax);
	  //mstart = cpeds_random_uniform_number(mmin,mmax);	
	} // using uniform distribution
	mstart = cpeds_random_uniform_number(mmin,mmax);
	// set the chain step size
	if (pass == 1) { 
	  ds = (log10(smax)-log10(smin))/levs; s = sstart; 
	  dA = (log10(Amax)-log10(Amin))/levA; A = Astart; 
	  //dm = (log10(mmax)-log10(mmin))/levm; m = mstart; 
	} // this is d(log s)
	if (pass == 2) { 
	  ds = (smax-smin)/levs; s = sstart; 
	  dA = (Amax-Amin)/levA; A = Astart;
	  //dm = (mmax-mmin)/levm; m = mstart;
	}
	dm = (mmax-mmin)/levm; m = mstart;
/* 	printf("  -- ****** dlog(A) = %lE, dlog(m) = %lE, dlog(s) = %lE\n",dA,dm,ds); */
/* 	printf("  -- ****** Astart = %lE, mstart = %lE, sstart = %lE\n",Astart,mstart,sstart); */
	iter = 0; iternum = levA*levm*levs;
	fit = tabulate_function(1,A,s,m,num,Xmin,Xmax); 
	Xsq_min = calculate_chisq(fit,functab,num);	delete fit;
	
	// initiatiate the walk in chain
	do {
	  Xsq_min_prev = Xsq_min;
	  //printf("***********************************************\n");
	  //sniff around in a small grid
	  if (pass==1) { 
	    locsstart = pow(10,log10(s)-ndelta*ds); locsend = pow(10,log10(s)+ndelta*ds); 
	    locAstart = pow(10,log10(A)-ndelta*dA); locAend = pow(10,log10(A)+ndelta*dA); 
	    //locmstart = pow(10,log10(m)-ndelta*dm); locmend = pow(10,log10(m)+ndelta*dm); 
	  }
	  if (pass==2) { 
	    locsstart = s-ndelta*ds; locsend = s+ndelta*ds; 
	    locAstart = A-ndelta*dA; locAend = A+ndelta*dA; 
	    //locmstart = m-ndelta*dm; locmend = m+ndelta*dm; 
	  }
	  locmstart = m-ndelta*dm; locmend = m+ndelta*dm; 
/* 	  printf("*******locstart = %lE locAstart=%lE, locmstart = %lE  \n",locsstart,locAstart,locmstart); */
/* 	  printf("*******locend = %lE locAend=%lE, locmend = %lE\n\n",locsend,locAend,locmend); */

	  locA = locAstart; cell_num=0;
	  for (i1=0;i1<side_num;i1++) {
	    locm = locmstart;
	    for (i2=0;i2<side_num;i2++) {
	      locs = locsstart;
	      //if (locAstart < Amin) { locAstart = Amin; } if (locsstart < smin) { locsstart = smin; } if (locmstart < mmin) { locmstart = mmin; }
	      //if (locAend > Amax)   { locAend = Amax; }   if (locsend > smax)   { locsend = smax; }   if (locmend > mmax)   { locmend = mmax; }
	      for (i3=0;i3<side_num;i3++) {
		direction_tab[cell_num][0] = locA; direction_tab[cell_num][1] = locm; direction_tab[cell_num][2] = locs;
		fit = tabulate_function(1,locA,locs,locm,num,Xmin,Xmax); 
		Xsq[cell_num] = calculate_chisq(fit,functab,num);      delete fit; 
		//printf("*******%li locs = %lE locA=%lE, locm = %lE chisq=%.30lE \n",cell_num,locs,locA,locm,Xsq[cell_num]);
		if (pass==1) { locs=pow(10,log10(locs)+ds); } //exit(0);
		if (pass==2) { locs=locs+ds;}
		cell_num++;
	      }
	      //if (pass==1) { locs=pow(10,log10(locm)+dm); } //exit(0);
	      //if (pass==2) { locm=locm+dm; }
	      locm=locm+dm;
	    }
	    if (pass==1) { locs=pow(10,log10(locA)+dA); } //exit(0);
	    if (pass==2) { locA=locA+dA; }
	  }

	  // find best direction and 
	  //exit(0);
	  //printf("!!!size = %li mini=%li, maxi=  %li, min = %lE, max = %lE\n",gridsize,imin,imax,Xsq_min,tmpd);
	  cpeds_find_minmax_value(Xsq,gridsize,&Xsq_min,&tmpd,&imin,&imax);	  //      printf("!!!!!!!!!!!!! imin: %li gridsize:%li\n",imin, gridsize);
	  //printf("!!!size = %li mini=%li, maxi=  %li, min = %lE, max = %lE\n",gridsize,imin,imax,Xsq_min,tmpd);
	  
	  if (Xsq_min < Xsq_min_prev) { // go there if it's better place 
	    if ((A <= Amax) && (A >= Amin) && (m >= mmin) && (m <= mmax) && (s <= smax) && (s >= smin)) { 
	      A = direction_tab[imin][0]; m = direction_tab[imin][1]; s = direction_tab[imin][2];	  }}

	  //printf("  -- iteration number:%.0lf / %.0lf -- %lE %% (A=%lE s=%lE m=%lE)\r", iter,iternum, iter/iternum*100,A,s,m);
	  iter++; //exit(0);
	} while ((iter < 2*iternum) && (Xsq_min < Xsq_min_prev));
	//if (iter == 1) { chain_num--; } // don't count bad seeded chains
	Afit = A; sfit = s; mfit = m;
	//printf("\n chain:%li **fitted (CDF)gauss function parameters: A = %lE, s = %lE m = %lE, chisq = %lE\n", chain_num,Afit,sfit,mfit,Xsq_min);
	//printf(" chain:%li **fitted (CDF)gauss function parameters: Astart = %lE, sstart = %lE mstart = %lE\n\n",  chain_num,Astart,sstart,mstart);
      
	// remember the current fit and chain information
	fits[chain_num].circ.Astart = Astart;      fits[chain_num].circ.mstart = mstart;       fits[chain_num].circ.sstart = sstart;
	fits[chain_num].circ.dA = dA;      fits[chain_num].circ.dm = dm;       fits[chain_num].circ.ds = ds;
	fits[chain_num].circ.A = Afit;      fits[chain_num].circ.m = mfit;       fits[chain_num].circ.s = sfit;
	fits[chain_num].circ.xmin = Xmin;      fits[chain_num].circ.xmax = Xmax;       fits[chain_num].circ.num = num;
	fits[chain_num].circ.chisq = Xsq_min;      fits[chain_num].circ.ok = 1;
	chains_total++;
	if (Xsq_min < 1E-1) {chain_num++; chains_total=0;}//minimum accuracy requirment
	if (chains_total > 50) chain_num++;
      }
      Xsqtmp = new double[number_of_chains]; // pick up the lowes chisq value to remember
      for (i=0;i<number_of_chains;i++) { Xsqtmp[i] = fits[i].circ.chisq; }//printf("chain %li, chisq: %lE\n",i,Xsqtmp[i]);}
      //printf("!!!!!!!!!!!!! number of chains (skad 62): %li\n",number_of_chains);
      cpeds_find_minmax_value(Xsqtmp,number_of_chains,&Xsq_min,&tmpd,&imin,&imax); delete Xsqtmp;
      //printf("!!!size = %li mini=%li, maxi=  %li, min = %lE, max = %lE\n",gridsize,imin,imax,Xsq_min,tmpd);
      
      // remember only the best fit
      mink_fit_info.circ.ok = 1; mink_fit_info.circ.chisq = fits[imin].circ.chisq;
      mink_fit_info.circ.A = fits[imin].circ.A;   mink_fit_info.circ.m = fits[imin].circ.m;   mink_fit_info.circ.s = fits[imin].circ.s;   
      mink_fit_info.circ.dA = fits[imin].circ.dA;   mink_fit_info.circ.dm = fits[imin].circ.dm;   mink_fit_info.circ.ds = fits[imin].circ.ds;   
      mink_fit_info.circ.Astart = fits[imin].circ.Astart;   mink_fit_info.circ.mstart = fits[imin].circ.mstart;   mink_fit_info.circ.sstart = fits[imin].circ.sstart;   
      mink_fit_info.circ.xmin = fits[imin].circ.xmin;   mink_fit_info.circ.xmax = fits[imin].circ.xmax;    mink_fit_info.circ.num = fits[imin].circ.num; 
    
      sfit1stpass = mink_fit_info.circ.s;       Afit1stpass = mink_fit_info.circ.A;       mfit1stpass = mink_fit_info.circ.m; 
      printf("\n----------------------------------------------------------------------------------\n");
      printf(" chain:%li **fitted gauss function parameters: A = %lE, s = %lE m = %lE, chisq = %lE\n",imin, fits[imin].circ.A,fits[imin].circ.s,fits[imin].circ.m,fits[imin].circ.chisq);
      printf(" chain:%li **fitted gauss function parameters: Astart = %lE, sstart = %lE mstart = %lE\n\n", imin, fits[imin].circ.Astart,fits[imin].circ.sstart,fits[imin].circ.mstart);
      delete fits;
    }
    //delete Xsq;
  }
    //----------------------------------------------------------------------------------------------------------------------------
  if (f == 2) { // fitting gsussCDF for area 
    printf("fitting gaussCDF function\n");
    for (pass=1;pass<=2;pass++) {
      printf(" \n------------\n*** pass %li\n-----------\n",pass);
      if (pass == 1) { // find the right order of magnitude of sigma
	smin = 1E-10; smax = 1E10; sstart = 0; levs=40; }
      if (pass == 2) { // find the right value within the found order of magnitude
	smin = pow(10,floor(log10(sfit1stpass))); smax = pow(10,ceil(log10(sfit1stpass))); sstart = pow(10,sfit1stpass); levs=1000; }
      Amin = 0.8; Amax = 1.2; Astart = 1.1; levA=1000;
      mmin = Xmin/2; mmax = Xmax/2; mstart = 0; levm=1000;
      printf(" Amin = %lE mmin = %lE, smin = %lE\n", Amin, mmin,smin);
      printf(" Amax = %lE mmax = %lE, smax = %lE\n", Amax, mmax,smax);
  
      for (i=0;i<num;i++) { functab[i] = data[2*i+1]; }
      if (pass == 1) {       number_of_chains = 20; }
      if (pass == 2) {       number_of_chains = chainnum2; }

      fits = new fit_info[number_of_chains];
      chain_num = 0;
      while (chain_num < number_of_chains) {
	//set the monte carlo markov chain start
	Astart = cpeds_random_uniform_number(Amin,Amax);
	mstart = cpeds_random_uniform_number(mmin,mmax);
	if (pass == 1) {  sstart = pow(10,cpeds_random_uniform_number(log10(smin),log10(smax)));} // using log uniform distribution
	if (pass == 2) {  sstart = cpeds_random_uniform_number(smin,smax);  } // using uniform distribution
    
	// set the chain step size
	dA = (Amax-Amin)/levA; A = Astart;
	dm = (mmax-mmin)/levm; m = mstart;
	if (pass == 1) { ds = (log10(smax)-log10(smin))/levs; s = sstart; } // this is d(log s)
	if (pass == 2) { ds = (smax-smin)/levs; s = sstart; }
	printf("****** dA = %lE, dm = %lE, dlog(s) = %lE\n",dA,dm,ds);
	printf("****** Astart = %lE, mstart = %lE, sstart = %lE\n",Astart,mstart,sstart);
	iter = 0; iternum = levA*levm*levs;
	fit = tabulate_function(2,A,s,m,num,Xmin,Xmax); 
	Xsq_min = calculate_chisq(fit,functab,num);	delete fit;
	
	// initiatiate the walk in chain
	do {
	  Xsq_min_prev = Xsq_min;

	  //sniff around in a small grid
	  locAstart = A-ndelta*dA; locAend = A+ndelta*dA; 
	  locmstart = m-ndelta*dm; locmend = m+ndelta*dm;
	  if (pass==1) { locsstart = pow(10,log10(s)-ndelta*ds); locsend = pow(10,log10(s)+ndelta*ds);}
	  if (pass==2) { locsstart = s-ndelta*ds; locsend = s+ndelta*ds;}

	  locA = locAstart; cell_num=0;
	  //while (locA<=locAend) {
	  for (i1=0;i1<side_num;i1++) {
	    locm = locmstart;
	    //while (locm<=locmend) {
	    for (i2=0;i2<side_num;i2++) {
	      //if (locAstart < Amin) { locAstart = Amin; } if (locsstart < smin) { locsstart = smin; } if (locmstart < mmin) { locmstart = mmin; }
	      //if (locAend > Amax)   { locAend = Amax; }   if (locsend > smax)   { locsend = smax; }   if (locmend > mmax)   { locmend = mmax; }
	      locs = locsstart;
	      //while (locs<=locsend) {
	      for (i3=0;i3<side_num;i3++) {
		direction_tab[cell_num][0] = locA; direction_tab[cell_num][1] = locm; direction_tab[cell_num][2] = locs;
		fit = tabulate_function(2,locA,locs,locm,num,Xmin,Xmax); 
		Xsq[cell_num] = calculate_chisq(fit,functab,num);      delete fit;
		if (pass==1) { locs=pow(10,log10(locs)+ds); }
		if (pass==2) { locs=locs+ds;}
		cell_num++;
	      }
	      locm=locm+dm;
	    }
	    locA=locA+dA;
	  }

	  // find best direction and 
	  cpeds_find_minmax_value(Xsq,gridsize,&Xsq_min,&tmpd,&imin,&imax);
	  if (Xsq_min < Xsq_min_prev) { // go there if it's better place
	    if ((A <= Amax) && (A >= Amin) && (m >= mmin) && (m <= mmax) && (s <= smax) && (s >= smin)) {
	      A = direction_tab[imin][0]; m = direction_tab[imin][1]; s = direction_tab[imin][2];	  }}

	  printf("iteration number:%.0lf / %.0lf -- %lE %% (A=%lE s=%lE m=%lE)\r", iter,iternum, iter/iternum*100,A,s,m);
	  iter++;
	} while ((iter < 2*iternum) && (Xsq_min < Xsq_min_prev));
	//if (iter == 1) { chain_num--; } // don't count bad seeded chains
	Afit = A; sfit = s; mfit = m;
	printf("\n chain:%li **fitted (CDF)gauss function parameters: A = %lE, s = %lE m = %lE, chisq = %lE\n", chain_num,Afit,sfit,mfit,Xsq_min);
	printf(" chain:%li **fitted (CDF)gauss function parameters: Astart = %lE, sstart = %lE mstart = %lE\n\n",  chain_num,Astart,sstart,mstart);
      
	// remember the current fit and chain information
	fits[chain_num].area.Astart = Astart;      fits[chain_num].area.mstart = mstart;       fits[chain_num].area.sstart = sstart;
	fits[chain_num].area.dA = dA;      fits[chain_num].area.dm = dm;       fits[chain_num].area.ds = ds;
	fits[chain_num].area.A = Afit;      fits[chain_num].area.m = mfit;       fits[chain_num].area.s = sfit;
	fits[chain_num].area.xmin = Xmin;      fits[chain_num].area.xmax = Xmax;       fits[chain_num].area.num = num;
	fits[chain_num].area.chisq = Xsq_min;      fits[chain_num].area.ok = 1;

	chain_num++;
      }
      Xsqtmp = new double[number_of_chains];
      for (i=0;i<number_of_chains;i++) { Xsqtmp[i] = fits[i].area.chisq;}
      cpeds_find_minmax_value(Xsqtmp,number_of_chains,&Xsq_min,&tmpd,&imin,&imax); delete Xsqtmp;

      // remember only the best fit
      mink_fit_info.area.ok = 1; mink_fit_info.area.chisq = fits[imin].area.chisq;
      mink_fit_info.area.A = fits[imin].area.A;   mink_fit_info.area.m = fits[imin].area.m;   mink_fit_info.area.s = fits[imin].area.s;   
      mink_fit_info.area.dA = fits[imin].area.dA;   mink_fit_info.area.dm = fits[imin].area.dm;   mink_fit_info.area.ds = fits[imin].area.ds;   
      mink_fit_info.area.Astart = fits[imin].area.Astart;   mink_fit_info.area.mstart = fits[imin].area.mstart;   mink_fit_info.area.sstart = fits[imin].area.sstart;   
      mink_fit_info.area.xmin = fits[imin].area.xmin;   mink_fit_info.area.xmax = fits[imin].area.xmax;    mink_fit_info.area.num = fits[imin].area.num; 
    
      sfit1stpass = mink_fit_info.area.s;
      printf("\n----------------------------------------------------------------------------------\n");
      printf(" chain:%li **fitted gauss function parameters: A = %lE, s = %lE m = %lE, chisq = %lE\n",imin, fits[imin].area.A,fits[imin].area.s,fits[imin].area.m,fits[imin].area.chisq);
      printf(" chain:%li **fitted gauss function parameters: Astart = %lE, sstart = %lE mstart = %lE\n\n", imin, fits[imin].area.Astart,fits[imin].area.sstart,fits[imin].area.mstart);
      delete fits;
    }
  }
  delete functab; delete Xsq;
}

//************************************************************************




















// MINKOWSKI FUNCTIONALS HELP FUNCTIONS
//************************************************************************
//returns the number of points in the map of temparatures greater then a given threshold level
// mu. The relevant temperatures in the map are T_i +_ nsig*sigma_i
// sigma_i = sigma_0/sqrt(Nobs_i) where the N_obs is the effective number of observations of a given pixel.
double map_class::Ngrmu(double mu, double nsig, double sigma0) {
  long int i;
  long int Ngr=0;

  for (i=0;i<pix_num;i++) {
    //if((*T)[i] > mu) {Ngr++;}
    if (map[i].m != 0) // we don't calculate the minkowski functionals on the mask
/*       if(map[i].T+nsig*sigma0/sqrt(map[i].N) >= mu) {Ngr++;} */
      if(map[i].T >= mu) {Ngr++;}
  }
  return Ngr;
}
//************************************************************************
double map_class::Nlemu(double mu) {
  long int i;
  long int Nle=0;

  for (i=0;i<pix_num;i++) {
    //if((*T)[i] > mu) {Ngr++;}
    if (map[i].m != 0) // we don't calculate the minkowski functionals on the mask
    if (map[i].T < mu) {Nle++;}
  }
  return Nle;
}
//************************************************************************


// some help routines

//************************************************************************
double map_class::calculate_minkowski_area_max() {
  long int i;
  area_max = minkowski_area[0][1];
  for (i=0;i<minkowski_level_num_area;i++) { if (minkowski_area[i][1] > area_max) { area_max = minkowski_area[i][1]; area_maxi = i; }}
  return area_max;
}
//************************************************************************
double map_class::calculate_minkowski_area_min() {
  long int i;
  area_min = minkowski_area[0][1];
  for (i=0;i<minkowski_level_num_area;i++) { if (minkowski_area[i][1] < area_min) { area_min = minkowski_area[i][1]; area_mini = i; }}
  return area_min;
}
//************************************************************************
double map_class::calculate_minkowski_circ_max() {
  long int i;
  circ_max = minkowski_circ[0][1];
  for (i=0;i<minkowski_level_num_circ;i++) { if (minkowski_circ[i][1] > circ_max) { circ_max = minkowski_circ[i][1]; circ_maxi = i; }}
/*   printf("&&&&&&&&&&&&&&&&&&&&&&& circ_max = %lE",circ_max); */
  return circ_max;
}
//************************************************************************
double map_class::calculate_minkowski_circ_min() {
  long int i;
  circ_min = minkowski_circ[0][1];
  for (i=0;i<minkowski_level_num_circ;i++) { if (minkowski_circ[i][1] < circ_min) { circ_min = minkowski_circ[i][1]; circ_mini = i; }}
  return circ_min;
}
//************************************************************************
double map_class::calculate_minkowski_genus_max() {
  long int i;
  genus_max = minkowski_genus[0][1];
  for (i=0;i<minkowski_level_num_genus;i++) { if (minkowski_genus[i][1] > genus_max) { genus_max = minkowski_genus[i][1]; genus_maxi = i; }}
/*   printf("&&&&&&&&&&&&&&&&&&&&&&& genus_max = %lE",genus_max); */
  return genus_max;
}
//************************************************************************
double map_class::calculate_minkowski_genus_min() {
  long int i;
  genus_min = minkowski_genus[0][1];
  for (i=0;i<minkowski_level_num_genus;i++) { if (minkowski_genus[i][1] < genus_min) { genus_min = minkowski_genus[i][1]; genus_mini = i; }}
  return genus_min;
}
