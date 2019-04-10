#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "Mscs-map.h"
#include "cpeds-direction.h"
//#include "cpeds-consts.h"
//#include "cpeds-math.h"


/******************************************************************************************************************************/
// this does partial derivative with respect to phi of the map u and returns the address of the object with u,phi
// the differentiation is performed directly on the map; not in the SH domain
// the returned map is in an input ordering
// the input map ordering should be nested

mscsMap mscsMap::partial_derivative_phi() {
	double dphi;
	long i,j,k,l,m,ring_num,pixinring_num,kr_beg,pixinring_numle1;
	mscsMap::mapOrderings orig_ordering;
	
//	printf(" * calculating partial derivative of the map with respect to phi\n");
	
	//
	//check and prepare the maps
	//
	
	// check and store the ordering of the u map
	orig_ordering = ordering();
	
	// clone u into a new object uphi
	mscsMap uphi=(*this);
	
	// load the proper coordinates into the uphi object
	if ( not uphi.coordLoaded()) { 
		if ( uphi.ordering()== mscsMap::ring) { uphi.conv_ring2nest(); }
		uphi.set_map_coord(0,0); 
	}
	
	// convert to ring ordering for calculation
	if ( uphi.ordering() == mscsMap::nested) { uphi.conv_nest2ring(); }
	if ( ordering() == mscsMap::nested) { conv_nest2ring(); }
	
	
	//
	// calculate the derivative ring by ring (calculation is performed in ring ordering)
	//
	
	ring_num = cpeds_get_ring_num_healpix(nside());
	k=0; // pixel iterator
	for (i=0;i<ring_num;i++) {
		kr_beg = k; // number of the pixel that starts the ring
		pixinring_num = cpeds_get_pixnum_in_ring_healpix(nside(),i);
		dphi = twoPI/(double)pixinring_num; 
		
		pixinring_numle1=pixinring_num-1;
		for (j=0;j<pixinring_num;j++) {
			// define the pixels (l,m) to be used for derivative calculation
			if (k < kr_beg+pixinring_numle1) l=k+1; else l=kr_beg;
			if (k > kr_beg) m=k-1; else m=kr_beg+pixinring_numle1;
			
			if (map.m[k] == 0) { uphi.map.m[l]=0; uphi.map.m[m]=0; } // mark empty pixels to extend the mask  latter
			if (map.m[l] == 0) { uphi.map.m[k]=0; uphi.map.m[m]=0; } // mark empty pixels to extend the mask  latter
			if (map.m[m] == 0) { uphi.map.m[k]=0; uphi.map.m[l]=0; } // mark empty pixels to extend the mask  latter
			
			/*       (*uphi).map->T[k] = 0.5 * ( (*u).map->T[m] + (*u).map->T[l] - 2.0*((*u).map->T[k]) ) / dphi;  // WHY THIS IS NOT THE SAME AS THE LINE BELOW ???? */ 
			uphi.map.T[k] = 0.5 * ( (map.T[k] - map.T[m]) / dphi + ( map.T[l] - map.T[k] ) / dphi );
			/*       (*uphi).map->T[k] = 0.5*((*u).map->T[k] - (*u).map->T[m])/dphi; */
			k++;
		}
	}
	
	
	// restore the original ordering
	if (orig_ordering == mscsMap::nested) { uphi.conv_ring2nest(); conv_ring2nest(); }
	
//	printf(" * derivative ,phi done\n");
	
	/*   (*uphi).savebinT("derivative-uphi",1); exit(0); */
	return uphi;
}

/******************************************************************************************************************************/
// this performs the covariant derivative  on sphere with respect to phi of the map u and returns the address of the object with u;phi
// the differentiation is performed directly on the map; not in the SH domain
mscsMap mscsMap::cov_derivative_phi() {
	double th;
	mscsMap uphi;
	long i;
	
//	printf(" * calculating covariant derivative of the map with respect to phi\n");
	uphi = partial_derivative_phi();
	
	for (i=0;i<pixNum();i++) {
		th = PIsnd - uphi.map.n[i].b();
		uphi.map.T[i] /= sin(th);
	}
	
	return uphi;
}


/******************************************************************************************************************************/
/******************************************************************************************************************************/
/******************************************************************************************************************************/
// this does partial derivative with respect to theta of the map u and returns the address of the new object with u,th
// the differentiation is performed directly on the map; not in the SH domain
// the map passed should be in nest ordering; if it's not the tables are loaded twice and it takes more time
// the returned map is in an input ordering
// if the map u is given in ring ordering it is converted to nested only for calculation 

mscsMap mscsMap::partial_derivative_th() {
	double dthij,dthik, th,phi, th2,phi2, th3,phi3;
	long i,j,k,l,m, ring_num,pixinring_num,ring,ring_numle1,ring_numle2;
	mscsMap::mapOrderings orig_ordering;
	long *pix_num_cuml;
	double rad225=0.39269908169872; // 22.5*PI/180;
	double rad45=0.78539816339745; // 45*PI/180
	double n45,n90;
	double acc=0.1570796326795/(double)nside(); // ~ 2*pi/(4*nside*10) // working accuracy for resolving the 22.5 67.5 ... etc meridians
	long ns2=2*nside();
	long ns3=3*nside();
	long nsle1=nside()-1;
	long ns3le2=ns3-2;
	
//	printf(" * calculating partial derivative of the map with respect to theta\n");
	
	// check and store the ordering of the u map
	orig_ordering = ordering();
	// clone u into a new object uphi
	mscsMap uth=(*this);
	
	// load the proper coordinates into the uth object
	if ( not uth.coordLoaded()) { 
		if ( uth.ordering()== mscsMap::ring) { uth.conv_ring2nest(); }
		uth.set_map_coord(0,0); 
	}
	
	// convert to ring ordering for calculation
	if ( uth.ordering() == mscsMap::nested) { uth.conv_nest2ring(); }
	if ( ordering() == mscsMap::nested) { conv_nest2ring(); }
	
	
	
	//
	// calculate the derivative (calculation is performed in nested ordering) (u,th for the pixels around the south pole are calculated across with each other)
	//
	
	// precalculate the cumulative number of pixels ring by ring and store in pix_num_cuml
	ring_num = cpeds_get_ring_num_healpix(nside());
	pix_num_cuml = new long[ring_num];
	//pix_num_ring = new long[ring_num];
	//pix_num_ring[0] = 
	pix_num_cuml[0] = 4;  
	
	// EXPERIMENTALL
	//FILE *f; f= fopen("testcoord","w");
	//fprintf(f,"%lE\n",(*utmp).map->n[pix_num_cuml[0]].b*PI180inv);
	// EXPERIMENTALL
	
	for (i=1;i<ring_num;i++) { 
		pixinring_num = cpeds_get_pixnum_in_ring_healpix(nside(),i);
		//pix_num_ring[i] = pixinring_num; 
		pix_num_cuml[i] = pix_num_cuml[i-1] + pixinring_num; 
		
		// EXPERIMENTAL
		//fprintf(f,"%lE\n",(*utmp).map->n[pix_num_cuml[i]].b*PI180inv);
		// experimental
	}
	// experimental
	//fclose(f); exit(0);
	// experimental
	uth.load_conv_tabs(nside(),"r2n");
	
	ring_numle1 = ring_num-1;
	ring_numle2 = ring_num-2;
	
	
	
	//
	// derive the derivative for all pixels
	//
	
	ring = 0;
	for (i=0;i<pixNum();i++) {
		th = PIsnd - uth.map.n[i].b();    
		phi = uth.map.n[i].l(); // th and phi coords of the i'th pixel in nested ordering
		
		// set the current ring value
		//ring = cpeds_get_ring_healpix_nest(nside,i); // get the current ring number of i'th pixel in nested ordering
		if (i>=pix_num_cuml[ring]) ring++;
		
		
		// choose the second pixel for the differencial-ratio calculation
		if (ring < nsle1) { // for north polar regions
			n45=(phi-rad225)/rad45;
			if (fabs(n45-round(n45))<acc) { // for meridians 22.5+n*45 , n=0,1,2,3,4,5,6,7 only
				j = pix_num_cuml[ring+1]-1; 
				if (ring > 0) k=pix_num_cuml[ring-1]-1; else k=j;	
			}
			else { // for all other longitudes
				j = pix_num_cuml[ring+2]-1;
				if (ring > 1) k=pix_num_cuml[ring-2]-1; else k=j;
			}
		}
		else { // for equatorial regions
			if (ring < ns2) { // for north equatorial regions
				j = pix_num_cuml[ring+2]-1; 
				
				if (ring==nside()) {
					n90=phi/PIsnd;
					if (fabs(n90-round(n90))<acc) { // for meridians n*90 , n=0,1,2,3 only
						k = pix_num_cuml[ring-1]-1; 
					}
					else k = pix_num_cuml[ring-2]-1; 
				}
				else k = pix_num_cuml[ring-2]-1; 
				
			}
			else { 
				if (ring < ns3) { // for south equatorial regions
					j = pix_num_cuml[ring-2]-1; 
					j=-j; // here we do the derivarive in reverse direction so we need to correct for it latter
					
					if (ring==ns3le2) {
						n90=phi/PIsnd;
						if (fabs(n90-round(n90))<acc) { // for meridians n*90 , n=0,1,2,3 only
							k = pix_num_cuml[ring+1]-1; 
						}
						else k = pix_num_cuml[ring+2]-1; 
					}
					else k = pix_num_cuml[ring+2]-1; 
					
				}
				else { // for south polar regions
					n45=(phi-rad225)/rad45;
					if (fabs(n45-round(n45))<acc) {  // for meridians 22.5+n*45 , n=0,1,2,3,4,5,6,7 only
						j = pix_num_cuml[ring-1]-1;
						j=-j; // here we do the derivarive in reverse direction so we need to correct for it latter
						if (ring < ring_numle1) k=pix_num_cuml[ring+1]-1; else k=j;
					}
					else { // for all other longitudes
						j = pix_num_cuml[ring-2]-1;
						j=-j; // here we do the derivarive in reverse direction so we need to correct for it latter
						if (ring < ring_numle2) k=pix_num_cuml[ring+2]-1; else k=j;
					}
				}
			}
		}
		
		// get coordinates and dth value of the pixel from the two rings below    
		if (j>=0) th2 = PIsnd - uth.map.n[j].b(); else th2 = PIsnd - uth.map.n[-j].b(); 
		if (k>=0) th3 = PIsnd - uth.map.n[k].b(); else th3 = PIsnd - uth.map.n[-k].b(); 
		dthij = th2-th; 
		dthik = th3-th; 
		phi2 = phi;  
		phi3 = phi;
		
		
		// get the pixel number in ring ordering from the map shifted by dth
		cpeds_ang2pix_healpix(nside(),&l,th2,phi2,1); // with th defined by j --> l
		cpeds_ang2pix_healpix(nside(),&m,th3,phi3,1); // with th defined by k --> m
		l=uth.r2n_conv[l];  // convert it to the ring pixel number //   THIS IS TEMPORARY STUFF  !!!
		m=uth.r2n_conv[m];  // convert it to the ring pixel number //   THIS IS TEMPORARY STUFF  !!!
		
		//calculate the derivative -i.e. the differencial ratio
		uth.map.T[i] = 0.5 * ( (map.T[l] - map.T[i])/dthij + (map.T[m] - map.T[i])/dthik );
		
		if (map.m[i] == 0) { uth.map.m[l]=uth.map.m[m]=0; } // mark empty pixels to extend the mask  latter
		if (map.m[l] == 0) { uth.map.m[i]=uth.map.m[m]=0; } // mark empty pixels to extend the mask  latter
		if (map.m[m] == 0) { uth.map.m[i]=uth.map.m[l]=0; } // mark empty pixels to extend the mask  latter
		
	}
	
	//delete utmp; 
	delete pix_num_cuml;
	
	// get back to original ordering
	if (orig_ordering == mscsMap::nested) { uth.conv_ring2nest(); conv_ring2nest();}
	
//	printf(" * derivarive ,th done\n");
	
	/*   (*uth).savebinT("derivative-uth",1); // save the derivative map as nested  */
	/*   exit(0); */
	
	return uth;
}


// new version but not finished - now it's a mess
/*
mscsMap mscsMap::partial_derivative_th() {
	double dthij,dthik, th,phi, th2,phi2, th3,phi3,dphi,dphisnd;
	long i,j,k,l,m, ring_num,pixinring,pixinringpo,ring,ring_numle1,ring_numle2;
	long orig_ordering;
	long *pix_num_cuml; // DEPRECIATED 
	double rad225=0.39269908169872; // 22.5*PI/180; // DEPRECIATED 
	double rad45=0.78539816339745; // 45*PI/180 // DEPRECIATED 
	double n45,n90;  // DEPRECIATED 
	double acc=0.1570796326795/(double)nside; // ~ 2*pi/(4*nside*10) // working accuracy for resolving the 22.5 67.5 ... etc meridians // DEPRECIATED 
	long ns2=2*nside; // DEPRECIATED 
	long ns3=3*nside(); // DEPRECIATED 
	long ns4=4*nside();
	long nsle1=nside()-1;
	long ns3le2=ns3-2;
	
	double *Xint, *Yint, *X, *Y, *TH;
	long Nint, N;
	matrix <double> Dth;
	cpedsDirection n;
	
	printf(" * calculating partial derivative of the map with respect to theta\n");
	
	//
	//check and prepare the map
	//
	
	
	// clone u into a new object uphi
	mscsMap uth=(*this);
	
	// load the proper coordinates into the uth object
	if ( not uth.coordLoaded()) { 
		if ( uth.ordering()== mscsMap::ring) { uth.conv_ring2nest(); }
		uth.set_map_coord(0,0); 
	}
	
	// convert to ring ordering for calculation
	if ( uth.ordering() == mscsMap::nested) { uth.conv_nest2ring(); }
	if ( ordering() == mscsMap::nested) { conv_nest2ring(); }
	
	
	
  // initiate 2D array and interpolation arrays
  ring_num = cpeds_get_ring_num_healpix(nside());
  matrix <double> Dth(ring_num,ns4);
  TH = new double[ring_num];
  Xint = new double[ns4];

  //set the l values at which we will interpolate 
  dphi=twoPI/(double)ns4; dphisnd=dphi/2.0;
  for (j=0;j<ns4;j++) { 
    Xint[j] = (double)j*dphi + dphisnd; // let's interpolate exactly in the directions towards to the centers of the equatorial pixels in the map
  }


  for (i=0;i<ring_num;i++) { 

    // get ring
    pixinring = cpeds_get_pixnum_in_ring_healpix(nside(),i);
    pixinringpo=pixinring+1;
    X = new double[pixinringpo];
    Y = new double[pixinringpo];
    k=0;
    for (j=0;j<pixinring;j++) { 
      n=uth.get_C(k);
      X[j]=n.l();
      Y[j]=uth.get_T(k);
    }
    TH[i]=PIsnd-n.b();

    // enforce periodicity
    X[pixinring]=X[0]+twoPI;
    Y[pixinring]=Y[0];

    // interpolate
    Yint = cpeds_interpolate(X,Y,pixinringpo,Xint,ns4,"cspline_periodic");

    // populate the maxtrix
    for (j=0;j<ns4;j++) { 
      Dth(i,j)=Yint[j];
    }

    // derive derivative
    cpeds_matrix_derivative(&Dth,TH);
    delete [] TH;
    
    // recast back onto the spherical map structure

    // clean up
    


	
	
	
	
	
	
	// restore the original ordering
	if (orig_ordering == mscsMap::nested) { uth.conv_ring2nest(); conv_ring2nest(); }
	
	
	printf(" * derivarive ,th done\n");
	
	return uth;
}
 */

/******************************************************************************************************************************/
// this performs the covariant derivative on sphere with respect to th of the map u and returns the address of the object with u;th
// the differentiation is performed directly on the map; not in the SH domain
mscsMap mscsMap::cov_derivative_th() {
	mscsMap uth;
//	printf("  * calculating covariant derivative of the map with respect to theta\n");
	uth = partial_derivative_th();
	
	return uth;
}
/******************************************************************************************************************************/
