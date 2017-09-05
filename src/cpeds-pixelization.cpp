/* C/C++ Routines for conversions between the ordered pixel numbers and directions on sky in different pixelizations schemes */
/* Written by Bartosz Lew (2005) */



/* ---------------------------------------------------------------------------- */
/* These routines are a part of the CPEDS library v.0.1 by Bartosz Lew */

/* HEALPIX and tetrahedron-based-pixelization-scheme will be available */
/* The "tetrahedron-based-pixelization-scheme" is not available at the moment */
/* At the moment the RING ordering scheme of the HEALPIX is not supported (only the NESTED one)
for the ang2pix pix2ang procedures */

#define _ISOC99_SOURCE
#include <features.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"

void cpeds_pix2colrow_healpix_ring(long nside,long pix, long *row, long *col) { // fills the tab array with the binary representation of the pixel number for the nested ring ordering and returns the col and row of the pixel in the region
	// first define the region
	long rings=cpeds_get_ring_num_healpix(nside);
	long ring=cpeds_get_ring_healpix_ring(nside,pix);
	long pixSt=cpeds_get_pix_num_above_ring_healpix(nside,ring);
	long pixInRing=cpeds_get_pixnum_in_ring_healpix(nside,ring);
	long region;

	if (ring <= nside) { // top regions
		region=(pix-pixSt)/(pixInRing/4);

	} else
	if (ring >rings-nside) { // bottom regions
		region=(pix-pixSt)/(pixInRing/4);

	}
	else { // equatorial regions or other regions
		region=8+(pix-pixSt)/(pixInRing/4);

	}

	printf(" cpeds error: cpeds_pix2colrow_healpix_ring is not implemented yet\n");
	exit(0);


}

void cpeds_pix2colrow_healpix_nest(long nside,long pix, long *row, long *col) { // fills the tab array with the binary representation of the pixel number for the nested ring ordering and returns the col and row of the pixel in the region
  long int byte_num,i,sum,z,X,Y,nsteps,NS;
  double dbyte_num;
  long *tab_num;
  long int pix_rownum=0, pix_colnum=0, region=0;

  // the number of bytes can be considered as set of pairs in which the first byte is the flag for the column number and the second is the flag for the row number
  // in a healpix region. EG. nside=8 pix = 36 ----> 3 pairs of numbers xx xx xx
  //                                                                    10 01 10 = 36
  //      this means that the pixel is in the                           ^  ^  ^   --- > 1 * 2^(3-1)  +  0 * 2^((3-1)-1)  +  1 * 2^((3-1)-0) = 5 'th column and
  //                                                                     ^  ^  ^  --- > 0 * 2^(3-1)  +  1 * 2^((3-1)-1)  +  0 * 2^((3-1)-0) = 2 'nd row counting from 0 in
  // the region

  region = (long int)floor((double)pix/(double)(nside*nside));
  /* reducing the pixel number to the number in a partucular region */
  pix = pix- nside*nside*region;

  dbyte_num = 2*cpeds_log2((double)nside);
  byte_num = (long)round(dbyte_num); // number of bytes on which the pixel number can be stored in binary system
/*     printf("dbyte_num=%lf, byte_num=%li \n",dbyte_num,byte_num); */
// !!!!!!!!!!!!carefull here !!!!!!!!!!!!!! if there is no round here then if you use optimalization flag -O2 then it doesn't work ok. !!!!!!!!!!!!!!!!!!


  /* setting up the binary representation for finding the pixel x,y position */
/*   tab_num = (long*)malloc(byte_num*sizeof(long)); // new int[byte_num]; */
  tab_num = new long[byte_num*(long)sizeof(long)];
  /*     long int *tab_num = new long int[byte_num*sizeof(long)]; // new int[byte_num]; */
  sum = 0;
  for (i=byte_num-1;i>=0;i--) {
    z = cpeds_pow_int(2,i);
    if ((sum+z) > pix) { tab_num[byte_num-1-i] = 0; } else { tab_num[byte_num-1-i] = 1; sum = sum+z; }
    /*       printf("byte_num=%li i = %i z = %i, sum = %i, tab=%i\n",byte_num,i,z, sum, tab_num[i]); */
  }
  /*     for (i=0;i<byte_num;i++) { */
  /*       printf("\n i = %i tab=%i\n",i,tab_num[i]);  */
  /*     } */


  X = Y = 0;  NS = 1;  pix_colnum = pix_rownum = 0;  nsteps = byte_num/2;

  for (i = 0; i< nsteps; i++) {
    NS = cpeds_pow_int(2,nsteps-i-1);
    if (tab_num[2*i] == 1) { pix_colnum +=  NS; }
    if (tab_num[2*i+1] == 1) { pix_rownum += NS; }
/*       printf("NS=%lf pix_colnum=%li,  pix_rownum = %li\n",NS,pix_colnum,pix_rownum); */
  }
  free(tab_num);
/*     delete tab_num; */
  *row = pix_rownum;
  *col = pix_colnum;
}

double cpeds_check_b(double b) {
  b=PIsnd-b;
  double phi=0;
  cpeds_check_thphi(&b,&phi);
  b=PIsnd-b;
  return b;
}

void cpeds_check_bl(double *b, double *l) {
  *b=PIsnd-*b;
  cpeds_check_thphi(b,l);
  *b=PIsnd-(*b);
  /* return *b; */
}

// the coordinates are in radians
void cpeds_check_thphi(double *th, double *phi) {
  double thl = *th, phil = *phi,tmp;
  double loc_acc=1e-12; // local accuracy

  if (thl == -0) { thl=*th=0; } // numerical correction
  if (thl < 0 ) { thl = -thl; tmp = floor(thl / (PI)+1); thl = tmp*PI - thl; phil += PI; }
  if (fabs(thl - PI) < loc_acc ) { thl = PI; } // numerical correction
  if (thl - PI > loc_acc ) { tmp = floor(thl / (PI)); thl = thl - tmp*PI; phil += PI; } // numerical correction

  if (phil == -0) { phil=*phi=0; } // numerical correction
  if (phil < 0 ) { phil = -phil; tmp = floor(phil / (twoPI)+1); phil = tmp*twoPI - phil;}
  if (fabs(phil - twoPI) < loc_acc ) { phil = 0; } // numerical correction
  if (phil - twoPI >= loc_acc ) { tmp = floor(phil / (twoPI)); phil = phil - tmp*twoPI; } // numerical correction

  *th = thl;
  *phi = phil;
}

double cpeds_check_phi(double phi) {  return cpeds_check_phi(&phi); }

double cpeds_check_phi(double *phi) {
  double phil = *phi,tmp;
  double loc_acc=1e-12; // local accuracy

  if (phil == -0) { phil=*phi=0; } // numerical correction
  if (phil < 0 ) { phil = -phil; tmp = floor(phil / (twoPI)+1); phil = tmp*twoPI - phil;}
  if (fabs(phil - twoPI) < loc_acc ) { phil = 0; } // numerical correction  2007-05-09
  if (phil - twoPI >= loc_acc ) { tmp = floor(phil / (twoPI)); phil = phil - tmp*twoPI; } // numerical correction
  *phi = phil;
  return *phi;
}

// returns the region number given nside and pixel number
long int cpeds_get_healpix_region(long int nside, long int pix) {
  long int region;

  region = (long int)floor((double)pix/(double)(nside*nside));
  //printf("*** region = %li pix = %li, **=%lf\n ",region,pix, floor((double)pix/(double)(nside*nside)));
  return region;
}

// returns the number of the pixel num in ordering in the region where it belongs
long int cpeds_get_healpix_pixnum_in_region(long int nside, long int pix) {
  long int region;

  region = cpeds_get_healpix_region(nside,pix);
  return pix-region*nside*nside;
}

// returns the number of the column in ordering in the region where the piexel num  belongs
// the column ordering is from 0 to nside-1
// THIS IS ACTUALLY FOR THE ORDERING WHICH IS NUMBERRED HERE AS 3 - NOT USEFULL FOR NESTED OR RING ORDERING
long int cpeds_get_healpix_colnum_in_region(long int nside, long int pix) {
  long int pix_num_in_reg;
  long int pix_colnum=0;
/*   double pix_colnumd; */

  pix_num_in_reg = cpeds_get_healpix_pixnum_in_region(nside,pix);
  pix_colnum = (long int)floor((double)pix_num_in_reg/(double)nside);
/*   pix_colnumd = (double)pix_num_in_reg/(double)nside; */
/*   if (fabs(pix_colnumd) < 1e-10) { pix_colnum = (long)round(pix_colnumd); } else { pix_colnum = (long)floor(pix_colnumd); } */
  return pix_colnum;
}

// returns the number of the row in ordering in the region where the piexel num belongs
// the column ordering is from 0 to nside-1
// THIS IS ACTUALLY FOR THE ORDERING WHICH IS NUMBERRED HERE AS 3 - NOT USEFULL FOR NESTED OR RING ORDERING
long int cpeds_get_healpix_rownum_in_region(long int nside, long int pix) {
  long int pix_num_in_reg;
  long int pix_rownum, pix_colnum;

  pix_num_in_reg = cpeds_get_healpix_pixnum_in_region(nside,pix);
  pix_colnum = cpeds_get_healpix_colnum_in_region(nside,pix);
  pix_rownum = pix_num_in_reg - pix_colnum*nside;
  return pix_rownum;
}
// this returns the posiotion of the pixels centers in a region, given the pix number and nside in a given ordering
// this works only in nest ordring for now, but for the ring ordering it's not difficult to implement
long cpeds_pix2xy(long int nside, long int pix, double *x, double *y, int ordering) {
  double half_pix_size = 1.0/(double)nside;
  long int pix_rownum=0, pix_colnum=0;

  if (ordering == 3) { // this is my own ordering approach allowing for nside to be any integer value
    pix_colnum = cpeds_get_healpix_colnum_in_region(nside,pix);
    pix_rownum = cpeds_get_healpix_rownum_in_region(nside,pix);
  }

  if (ordering == 1) { // this is standard nested HEALPIX ordering
    cpeds_pix2colrow_healpix_nest(nside,pix,&pix_rownum,&pix_colnum);
  }

  if (ordering == 0) { // this is standard ring HEALPIX ordering
    cpeds_pix2colrow_healpix_ring(nside,pix,&pix_rownum,&pix_colnum);
  }

/*   printf("*** colnum = %li rownum = %li ",pix_colnum,pix_rownum); */

/* this is the simplified form of the comment below */
  *x = half_pix_size*( (double)(pix_rownum - pix_colnum) );
  *y = half_pix_size*( (double)(pix_rownum + pix_colnum + 1) );

/*   *x = 0.5*( (2*half_pix_size*pix_rownum + half_pix_size) - (2*half_pix_size*pix_colnum + half_pix_size) ); */
/*   *y = 0.5*( (2*half_pix_size*pix_rownum + half_pix_size) + (2*half_pix_size*pix_colnum + half_pix_size) ); */
  return 0;
}
/**********************************************************************************/
/* the two functions for converting pixels numbers to directions in sky and vice versa get and return the  */
/* angles in coordinates system consistent with ang2pix and pix2ang routines from HEALPIX package i.e. */
/* 0 <= th <=PI, 0<=phi <= 2PI */
/* these are converted internally only for calculation purposes */

/* ordering -- 1 - when nested ordering */
/* ordering -- 0 - when ring ordering */
long cpeds_pix2ang_healpix(long int nside, long int pix, double *th, double *phi, int ordering) {
  double x,y, quarter;
  long int region;

  cpeds_pix2xy(nside,pix,&x,&y,ordering);
  region = cpeds_get_healpix_region(nside,pix); //(long int)floor((double)pix/(double)(nside*nside));
/*   printf("*** reg. = %li ",region); */
/*   printf("x = %lf, y = %lf, pix = %li ", x , y, pix); */

  if ( (region == 0) || (region == 1) || (region == 2) || (region == 3) ) { // we're in the north polar region
    quarter = (double)region;
    if (y >= 1) {
      *phi = quarter*PIsnd + PI*(2-y+x)/(4*(2-y));
      *th  = asin(1-(2-y)*(2-y)/3);
/*       printf("|||th=%lf y=%lf arg=%lf",*th,y,1-(2-y)*(2-y)/3 ); */
    }
    if (y < 1) { // we're still in the equatorial region
      *phi = quarter*PIsnd + PI*(1+x)/4;
      *th = asin(2*y/3);
    }
  }

  if ((region == 4) || (region == 5) || (region == 6) || (region == 7)) { // we're in the equatorial region
    quarter = (double)(region-4);
    *phi = quarter*PIsnd + PI*x/4; if (*phi < 0) { *phi = twoPI+*phi; }
    *th = asin(2*(y-1)/3);
  }

  if ((region == 8) || (region == 9) || (region == 10) || (region == 11)) { // we're in the south polar region
    quarter = (double)(region-8);
    y = 2-y;
    if (y >= 1) {
      *phi = quarter*PIsnd + PI*(2-y+x)/(4*(2-y));
      *th = -asin((1-(2-y)*(2-y)/3));
    }
    if (y < 1) { // we're still in the equatorial region
      *phi = quarter*PIsnd + PI*(1+x)/4;
      *th = -asin(2*y/3);
    }
  }
  *th = PIsnd-*th;
/*   printf("th_cpeds = %.15lf, phi_cpeds = %lf\n ", 180/PI*(*th) , 180/PI*(*phi)); */
  return 0;
}

/**********************************************************************************/
// !!!! THIS WORKS ONLY IN NESTED ORDERING FOR NOW !!!!
/* ordering -- 1 - when nested ordering */
/* ordering -- 0 - when ring ordering */
/* the invers formulas cannot be used beyond the area of the initial region since they are no longer valid: exactly the formulas for the polar regions y>=1 */
/* usage is ok only in some special cases in polar regions and it's ok in all equatorial regions including y<1 parts of polar regions */

// WARNING !! A BUG WAS FOUND IN THIS ROUTINE on 2008-01-10:
// optimized on 2008-01-10 when tracing a bug that causes phi=0 th=-7.5312388177E-01 to return pixnum=2359133 which corresponds to l,b [deg] = 4.6142579202E+01, -1.0445091833E+00 !!!!
// the incorrectly derived pixel numbers come from l=0 b in (~-1,-8) deg in the 9'th region along bounary with region 6
// similarly  - in the 10'th regions along boundary with reg. 7

long cpeds_ang2pix_healpix(long int nside, long int *pix, double th, double phi, int ordering) {
  double TH_TROPIC = 0.729727656226966; // ~ 41 deg //PI/2-48.189685109680696*PI/180; // boundary of the polar-equator region in corrdinates where pole has theta=PI/2
  double half_pix_size = 1.0/(double)nside;
  double x=0.0,y=0.0,x_eq,y_eq, x_np, y_np, x_sp, y_sp;
  double phi_np,phi_sp,phi_eq,region=0.0,region_np,region_sp,region_eq, pix_rownumd, pix_colnumd,dbnum;//,divcol,divrow;
  long int pix_rownum, pix_colnum,pix_num,i,nsteps,byte_num, col,row,j;
  long test_pix;
  cpeds_check_thphi(&th,&phi);
  // debug info
  double regionNP,regionSP,xNP,yNP,xSP,ySP;

  th = PIsnd-th; // transformation from Healpix coordinates to the coordinates in which Boud's transformations were calculated
  //----------------------------------------------------------------------------------------------------------------------------
  if (th > TH_TROPIC) { // we're in the north polar region
    region = floor(phi*2.0/PI); phi_np = phi - region* PIsnd; // this is to keep x inside the inital region to avoid mistakes.
    x = (-1.0 + 4.0*phi_np/PI) * sqrt(3.0-3.0*sin(th));// + half_pix_size;
/*     x = x + region*2; */
    y = 2.0 - sqrt(3.0-3.0*sin(th)); //- half_pix_size;
/*     printf("(1) x = %lf, y=%lf \n ", x,y); */
    //printf("(1) ------------> x = %lf, y = %lf region = %lf\n ",x,y,region);

    // debug
    xNP=x; yNP=y; regionNP=region;
  }
  //----------------------------------------------------------------------------------------------------------------------------
  if ((th <= TH_TROPIC) && (th >= -TH_TROPIC)) { // we're basically in the equatorial region but still some of the poles region are possible - so we calculate for the two cases
    region_np = floor(phi*2.0/PI); phi_np = phi - region_np* PIsnd; // this is to keep x inside the inital region to avoid mistakes.
    x_np = 4.0*phi_np/PI - 1.0; y_np = 3.0*sin(th)/2.0; // this is for the polar regions in equatorial part
/*     printf("-------------> x_np = %lf, y_np = %lf region_np = %lf\n ",x_np,y_np,region_np); */
    region_eq = floor((phi+PI/4.0)*2.0/PI); phi_eq = phi - region_eq* PIsnd; if (region_eq == 4.0) { region_eq = 0.0;}  // this is to keep x inside the inital region to avoid mistakes.
    x_eq = 4.0*phi_eq/PI; y_eq = 0.5 * (3.0*sin(th) + 2.0); // this is a guess for equatorial region

    th = -th;
    region_sp = floor(phi*2.0/PI); phi_sp = phi - region_sp* PIsnd; // this is to keep x inside the inital region to avoid mistakes.
    x_sp = 4.0*phi_sp/PI - 1.0; y_sp = 2.0-3.0*sin(th)/2.0; // this is for the polar regions in equatorial part

    //printf("(2) x_eq = %lf, y_eq = %lf, x_np = %lf y_np = %lf x_sp = %lf, y_sp = %lf region_np = %lf region_eq = %lf region_sp = %lf\n ",x_eq,y_eq, x_np, y_np,x_sp, y_sp,region_np,region_eq,region_sp);
    th = -th;

    if ( (y_eq < -fabs(x_eq)+2.0) && (y_eq >= fabs(x_eq)) ) { region = region_eq+4.0; x = x_eq; y = y_eq; }
    else {
      if (y_eq >= -fabs(x_eq)+2.0) { region = region_np; x = x_np; y = y_np; } else { region = region_sp+8.0; x = x_sp; y = y_sp; }
    }
/*     printf("-------------> x = %lf, y = %lf region = %lf\n ",x,y,region); */
  }
  //----------------------------------------------------------------------------------------------------------------------------
  if (th < -TH_TROPIC) { // we're in the south polar region
    th = -th;
    region = floor(phi*2.0/PI)+8.0; phi_sp = phi - (region-8.0)* PIsnd; // this is to keep x inside the inital region to avoid mistakes.
/*     x = (-1.0 + 4.0*phi_sp/PI) * sqrt(3.0-3.0*sin(th));// + half_pix_size; */
/*     y = sqrt(3.0-3.0*sin(th));// + half_pix_size; */

// optimized @ 2008-01-10
    y = sqrt(3.0-3.0*sin(th));// + half_pix_size;
    x = (-1.0 + 4.0*phi_sp/PI) * y;// + half_pix_size;
    // this optimization has an impact on the results in fact !!!

/*     x = x + region*2; */

    //  printf("(3) x = %lf, y=%lf region: %lf\n ", x,y,region); // debug stuff
    th = -th;
/*     printf("-------------> x = %lf, y = %lf region = %lf\n ",x,y,region); */
    // debug
    xSP=x; ySP=y; regionSP=region;
  }
  //----------------------------------------------------------------------------------------------------------------------------

/*   pix_rownumd = round(0.5 * (y+x-half_pix_size)/half_pix_size); if (pix_rownumd < 0) printf("WARNING  pixrownumd: %.20lf th=%.20lf phi=%.20lf\n",pix_rownumd,180/PI*th,180/PI*phi); */
/*   pix_colnumd = round(0.5 * (y-x-half_pix_size)/half_pix_size); // obczaic to !!! */

/*   commented out on 2008-01-11 when searching for The bug "ang2pix_healpix-bug1", see BUGTRAQ for more info */
/*   pix_rownumd = floor(0.5 * (y+x)/half_pix_size); if (pix_rownumd >= (double)nside) pix_rownumd = (double)nside-1.0; // this keeps the pixel within the selected region */
/*   pix_colnumd = floor(0.5 * (y-x)/half_pix_size); if (pix_colnumd >= (double)nside) pix_colnumd = (double)nside-1.0; */
/*   pix_rownum = (long int)pix_rownumd; */
/*   pix_colnum = (long int)pix_colnumd; */
/*   commented out on 2008-01-11 when searching for The bug "ang2pix_healpix-bug1", see BUGTRAQ for more info */

// modified on 2008-01-11
  pix_rownumd = (0.5 * (y+x)/half_pix_size);
  pix_colnumd = (0.5 * (y-x)/half_pix_size);
/*   pix_rownum = (long int)round(pix_rownumd); // commented out on 2008-02-08 -- this was buggy as comparred with healpix routines: floor should be used instead*/
/*   pix_colnum = (long int)round(pix_colnumd); // commented out on 2008-02-08 -- this was buggy as comparred with healpix routines: floor should be used instead*/
  pix_rownum = (long int)floor(pix_rownumd);
  pix_colnum = (long int)floor(pix_colnumd);
  if (pix_rownum >= nside) pix_rownum = nside-1; // this keeps the pixel within the selected region
  if (pix_colnum >= nside) pix_colnum = nside-1;
  if (pix_rownum < 0) pix_rownum = 0; // this keeps the pixel within the selected region
  if (pix_colnum < 0) pix_colnum = 0;
// modified on 2008-01-11

//  printf("^^^^ th = %lf, phi = %lf half_pix_size = %lf, colnum = %li rownum = %li colnumd=%lf rownumd=%lf x=%lf y=%lf reg.=%lf\n", 90-180/PI * th, 180/PI * phi, half_pix_size,pix_colnum, pix_rownum,pix_colnumd,pix_rownumd,x,y,region); // debug stuff



  dbnum = 2.0*cpeds_log2((double)nside);
  byte_num = (long int)round(dbnum);

  /* reducing the pixel number to the number in a particular region */

  pix_num = nside*nside*(long int)region;


  // !!!!!!!! THIS IS TERRIBLY STUPID --- FIX IT !!!!!!!!!!!
  /* setting up the binary representation for finding the pixel x,y position */
/*   long *tab_num = (long*)malloc(byte_num*sizeof(long)); // new int[byte_num]; // commented out on 2008-01-10 to avoid mixing new with malloc */
  long *tab_num = new long[byte_num*(long)sizeof(long)];

/*   commented out on 2008-01-11 when searching for The bug "ang2pix_healpix-bug1", see BUGTRAQ for more info */
/*   nsteps = byte_num/2; */
/*   for (i = 0; i< nsteps; i++) { */
/*     divrow = floor(pix_rownumd / 2.0); tab_num[2*(nsteps-i)-1] = (long)round(pix_rownumd - 2*divrow); pix_rownumd = divrow; */
/*     divcol = floor(pix_colnumd / 2.0); tab_num[2*(nsteps-i)-2] = (long)round(pix_colnumd - 2*divcol); pix_colnumd = divcol; */
/*   } */
/*   commented out on 2008-01-11 when searching for The bug "ang2pix_healpix-bug1", see BUGTRAQ for more info */

// optimized on 2008-01-11
  nsteps = byte_num/2;
  for (i = 0; i< nsteps; i++) { j=2*(nsteps-i);
    tab_num[j-1] = pix_rownum % 2; row = (pix_rownum - tab_num[j-1])/2; pix_rownum = row;
    tab_num[j-2] = pix_colnum % 2; col = (pix_colnum - tab_num[j-2])/2; pix_colnum = col;
  }
// optimized on 2008-01-11


  j = 0;  // sum variable was replaced with j, to decrease the number of variables to allocate on cheap for each call to the routine
/*   for (i=byte_num-1;i>=0;i--) {    sum = sum+pow(2.0,(double)i)*(double)tab_num[byte_num-1-i];  } */
  for (i=byte_num-1;i>=0;i--) {
/*     printf("%li ",tab_num[byte_num-1-i]); // debug stuff */
    j += cpeds_pow_int((long)2,i) * tab_num[byte_num-1-i];
  } // optimized on 2008-01-10 when tracing a bug that causes phi=0 th=-7.5312388177E-01 to return pixnum=2359133 which corresponds to l,b [deg] = 4.6142579202E+01, -1.0445091833E+00 !!!!

/*     printf("\n"); // debug stuff */
/*   free(tab_num); // commented out on 2008-01-10 to avoid mixing new with malloc */
  delete [] tab_num;

  *pix = pix_num + j;
/*   printf("==========>>>>>>>>>>> pix = %li\n\n",*pix); */

/*   if (ordering == 0) { // if ring ordering is required */
/*     sum = 0; */
/*     ringpix_num = cpeds_get_ring_healpix_nest(nside,*pix)-1; */
/*     for (i=0;i<ringpix_num;i++) { sum+=cpeds_get_pixnum_in_ring_healpix(nside,i); } */
/*   } */

// added on 2008-01-11: here a check of the number should be included with the cpeds_pix2ang_healpix routing to check if the
// returned pixel falls close enough (for a given resolution) with the original one. since there are discontinouities in the
// numbering of pixels w.r.t
// their spatial alignments and since the problems often appeared on the boundaries of the regions, it's likely that such test
// could be a good cross check.

//EXPERIMENTAL -- USE OF HEALPIX PROCEDURE
/*   FILE* controlf; */
/*   double check_th,check_phi,check_thH,check_phiH; */
/*   ang2pix_nest(nside,PIsnd-th,phi,&test_pix);  */
/*   if (*pix!=test_pix) {  */
/*     printf("CPEDS_PIXELIZATION_ERROR: nside: %li ang(lb coord) th: %lE phi: %lE, cpeds_ang2pix_healpix: %li, ang2pix_nest: %li\n",nside,PI180inv*th,PI180inv*phi,*pix,test_pix);  */
/*     controlf = fopen("/home/blew/programy/CPEDS/cpeds-pixelization-ang2pix.debug","a");   */
/*     fprintf(controlf,"nside: %li b: %lE l: %lE, cpeds_ang2pix_healpix: %li, ang2pix_nest: %li ",nside,PI180inv*th,PI180inv*phi,*pix,test_pix);  */
/*     fprintf(controlf,"(1) ------------> x = %lf, y = %lf region = %lf ",xNP,yNP,regionNP); */
/*     fprintf(controlf,"(2) x_eq = %lf, y_eq = %lf, x_np = %lf y_np = %lf x_sp = %lf, y_sp = %lf region_np = %lf region_eq = %lf region_sp = %lf ",x_eq,y_eq, x_np, y_np,x_sp, y_sp,region_np,region_eq,region_sp); */
/*     fprintf(controlf,"(3) x = %lf, y=%lf region: %lf ", x,y,region);  */
/*     fprintf(controlf,"^^^^ th = %lf, phi = %lf half_pix_size = %lf, colnum = %li rownum = %li colnumd=%lf rownumd=%lf x=%lf y=%lf reg.=%lf", 90-180/PI * th, 180/PI * phi, half_pix_size,pix_colnum, pix_rownum,pix_colnumd,pix_rownumd,x,y,region); // debug stuff */
/*     pix2ang_nest(nside,*pix,&check_th,&check_phi); */
/*     pix2ang_nest(nside,test_pix,&check_thH,&check_phiH); */
/*     fprintf(controlf," ang: %lE\n",PI180inv*cpeds_ang_n1n2(check_th,check_phi,check_thH,check_phiH)); */

/*     fclose(controlf); */
/*     *pix=test_pix;  */
/*   } */
//EXPERIMENTAL -- USE OF HEALPIX PROCEDURE


  return 0;
}
/**********************************************************************************/
/* returns the approximte size of the pixel in healpix pixelization system for a given nside [rad]. */
double cpeds_pix_size_healpix(long int nside) {
  long int pix_num = 12*nside*nside;

  return sqrt(4*PI/pix_num);
}

//************************************************************************
// returns the total number of rings in the healpix map of a given nside
long int cpeds_get_ring_num_healpix(long int nside) {
  long int num;
  //  if (pix_system == 1 ) {
  num =  2*(2*nside-1)+1; //}
  //  if (pix_system == 2 ) { /* to be done latter  */}
  return  num;
}
/**********************************************************************************/
// returns the number of the ring corresponding to a given pixel in NESTED ordering
// it seems that an ERROR WAS FOUND IN THIS ROUTINE ---- TO BE DEBUGGED !!!!!!!!!!!!!!
long int cpeds_get_ring_healpix_nest (long int nside,long int i) {
  double x,y;
  long int ring, region;
  cpeds_pix2xy(nside,i,&x,&y,1);
  ring = (long int)round((2-y)*(double)nside-1);
  region = cpeds_get_healpix_region(nside,i);
  if ((region >= 4) && (region <8)) {ring = ring + nside;}
  else {
    if (region >= 8) { ring = ring + 2*nside; }
  }

  return ring;
}
/**********************************************************************************/
// returns the number of the ring corresponding to a given pixel in RING ordering
long int cpeds_get_ring_healpix_ring(long int nside,long int i) {
  long int ring, pix;
  pix = 0;    ring = 0;
  while (pix <= i) {
    pix = pix + cpeds_get_pixnum_in_ring_healpix(nside,ring);
    ring++;
  }
  ring--;

  return ring;
}

/**********************************************************************************/
// returns the number of pixels in a given ring
long int cpeds_get_pixnum_in_ring_healpix(long int nside,long int ring) {
  long int pix;
  if (ring < nside) { pix = 4*(ring+1); }
  else {
    if (ring >= nside+(2*nside-1)) { pix  = 4*(cpeds_get_ring_num_healpix(nside)-ring); }
    else { pix = 4*nside; }
  }

  return pix;
}

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
//these procedures for RING NEST converstions are probably not consistent with RING definition of Gorski et. al
// because although the pixels are aligned in the rings the these  pixels in a given ring may not ally with systematically increasing phi coordinate
// but for integrations purposes or FFT this is completly irrelevant - this would become important when the transform was to be performed
// only in some confined part of the sky.
// The problem in this transformation is that ring allied pixels in the map - i.e. the sequences of pixels of length of the subsequent
// numbers of pixels in the subsequent ring numbers are in the same ring - the same lattitude, however the starting pixel in a given
// ring may not be of the lowest phi value - rather it's more less arbitrary  - i.e. the starting longitude in a given ring is at
// the smallest pixel number in nest ordering in this ring that was encountered as the loop iteratively reaches greater values.
//---------------------
// correction from Jan 23, 2006 - now these functions are perfectly consistent with HEALPIX package although terribly slow
//   the floor function doesn't return the argument value in case when the argument is integer already and it should !!!
//   that's why if the argument is integer to within 1e-10 then it uses round function which gives the right value, but if the
//   value is not integer (in this case it can be close to x.5, then floor works ok. In polar zones it's not a problem since
//   there are only distances of type x.5.
// --------------- THIS PROBLEM (and speed)  SHOULD BE FIXED SOMEHOW - or the transformations can be precalculated in a single
// dimentional array
// comment from 2006/03/21: this is the problem of the compilation flags I guess -O2
// and INTRINSICS stuff for skipping recognising the functions by parameetrs in order to increase speed during program execution
// This optimization concerns floor(), ceil(), round() and some other functions like these from the math library. But it works fine as it is now.
void cpeds_conv_nest2ring_healpix(long int nside, double * mapNEST, double * mapRING) {
  long int i,j,ring_num,pix_num;
  double halfpix_size;

  pix_num = 12*nside*nside;
  halfpix_size = 1.0/(double)nside;
  ring_num = cpeds_get_ring_num_healpix(nside);
  long * ringpixtab = new long[ring_num];

  ringpixtab[0] = 0;
  for (i=1;i<ring_num;i++) { ringpixtab[i] = ringpixtab[i-1] + cpeds_get_pixnum_in_ring_healpix(nside,i-1); }

  for (i=0;i<pix_num;i++) {
    j = cpeds_i2j_healpix(nside,i,ringpixtab);
    //j-=1; // less 1 to get the index number as we start indexig from 0
    mapRING[j] = mapNEST[i];
    //printf("**** convN2R: i=%li, ring=%li, j=%li\n\n",i,ring,j);
  }
  //for (i=0;i<pix_num;i++) { tmpmap[i] = map[i]; }
  //delete map;
  delete ringpixtab;

}

void cpeds_conv_nest2ring_healpix(long int nside, cpeds_direction * mapNEST, cpeds_direction * mapRING) {
  long int i,j,ring_num,pix_num;
  double halfpix_size;

  pix_num = 12*nside*nside;
  halfpix_size = 1.0/(double)nside;
  ring_num = cpeds_get_ring_num_healpix(nside);
  long * ringpixtab = new long[ring_num];

  ringpixtab[0] = 0;
  for (i=1;i<ring_num;i++) { ringpixtab[i] = ringpixtab[i-1] + cpeds_get_pixnum_in_ring_healpix(nside,i-1); }

  for (i=0;i<pix_num;i++) {
    j = cpeds_i2j_healpix(nside,i,ringpixtab);
    //j-=1; // less 1 to get the index number as we start indexig from 0
    mapRING[j] = mapNEST[i];
    //printf("**** convN2R: i=%li, ring=%li, j=%li\n\n",i,ring,j);
  }
  //for (i=0;i<pix_num;i++) { tmpmap[i] = map[i]; }
  //delete map;
  delete ringpixtab;

}
/**********************************************************************************/
void cpeds_conv_ring2nest_healpix(long int nside, double * mapRING,double * mapNEST) {
  long int i,j,ring_num,pix_num;
  double halfpix_size;

/*   long int i,j,k,ring,ring_num,nside; */
/*   double * ringpixtab = new double[pix_num]; */
/*   long int *pixinring; */
  pix_num = 12*nside*nside;
  //  nside = (long int)round(sqrt(((double)pix_num)/12));
  halfpix_size = 1.0/(double)nside;
  ring_num = cpeds_get_ring_num_healpix(nside);
  long * ringpixtab = new long[ring_num];

  ringpixtab[0] = 0;
  for (i=1;i<ring_num;i++) { ringpixtab[i] = ringpixtab[i-1] + cpeds_get_pixnum_in_ring_healpix(nside,i-1); }

  for (i=0;i<pix_num;i++) {
    j = cpeds_i2j_healpix(nside,i,ringpixtab);
    mapNEST[i] = mapRING[j];
    //pixinring[ring] = pixinring[ring]+1;
    //printf("**** convR2N: i=%li, ring=%li, j=%li\n\n",i,ring,j);
  }
  //for (i=0;i<pix_num;i++) { tmpmap[i] = map[i]; }
  //delete map;
  delete ringpixtab;

}

void cpeds_conv_ring2nest_healpix(long int nside, cpeds_direction * mapRING, cpeds_direction * mapNEST) {
  long int i,j,ring_num,pix_num;
  double halfpix_size;

/*   long int i,j,k,ring,ring_num,nside; */
/*   double * ringpixtab = new double[pix_num]; */
/*   long int *pixinring; */
  pix_num = 12*nside*nside;
  //  nside = (long int)round(sqrt(((double)pix_num)/12));
  halfpix_size = 1.0/(double)nside;
  ring_num = cpeds_get_ring_num_healpix(nside);
  long * ringpixtab = new long[ring_num];

  ringpixtab[0] = 0;
  for (i=1;i<ring_num;i++) { ringpixtab[i] = ringpixtab[i-1] + cpeds_get_pixnum_in_ring_healpix(nside,i-1); }

  for (i=0;i<pix_num;i++) {
    j = cpeds_i2j_healpix(nside,i,ringpixtab);
    mapNEST[i] = mapRING[j];
    //pixinring[ring] = pixinring[ring]+1;
    //printf("**** convR2N: i=%li, ring=%li, j=%li\n\n",i,ring,j);
  }
  //for (i=0;i<pix_num;i++) { tmpmap[i] = map[i]; }
  //delete map;
  delete ringpixtab;

}
/**********************************************************************************/

// for i'th pixel in nested ordering calculate the number in ring ordering provided the nside value and the array of size ring_num with numberes
// of the piexls opening each ring+1: i.e. ringpixtab[0] = 0 ringpixtab[1] = 4 and so on cumulatively
long cpeds_i2j_healpix(long nside, long i, long* ringpixtab) {
  long int j,ring;
  double x,y,th,phi,X,prev_pix_numD,partreg_pixD,tmp,nsideD=(double)nside;
  long subregion,region,fullreg_pix, partreg_pix, pixinring,prev_pix_num;

  //variable defs.
  ring = cpeds_get_ring_healpix_nest(nside,i); cpeds_pix2xy(nside,i,&x,&y,1);     region = cpeds_get_healpix_region(nside,i);
  pixinring = cpeds_get_pixnum_in_ring_healpix(nside,ring);    subregion = region;
  if (region >= 4) { subregion = region-4; }     if (subregion >= 4) { subregion -= 4; } // subregion = region mod 4 //zone def.
  //printf("region=%li subregion=%li ring=%li, pixinring=%li x=%lE y=%lE\n",region,subregion,ring,pixinring,x,y);
  //-------------
  j = ringpixtab[ring]; // number of all pixels upto the last ring (ring-1)
  //printf("1) j=%li \n",j);

  if ((region < 4) && (y > 1.0)) { // north polar zone
    fullreg_pix = subregion*(long)round((double)pixinring/4.0); // plus number of pixels in this ring in all previous regions
    partreg_pixD = fabs(y-x-2)*(double)nside/2.0; tmp = floor(partreg_pixD);
    partreg_pix = (long int)tmp; // plus number of pixels in this region
    j += fullreg_pix + partreg_pix;
    //printf("fullreg_pix=%li, partreg_pix=%li\n",fullreg_pix, partreg_pix);
  }

  if (((region >=4) && (region < 8)) || ((region < 4) && (y<=1)) || ((region >=8) && (y>=1))) { // equatorial zone def.
    cpeds_pix2ang_healpix(nside,i,&th,&phi,1); X = phi/PI*4.0;//-halfpix_size;
    //prev_pix_numD = X*((double)nside)/2.0;   tmp = floor((double)prev_pix_numD);
    prev_pix_numD = X*nsideD/2.0;   //tmp = trunc(prev_pix_numD);
    //prev_pix_num = (long)trunc((X*(double)nside/2.0)); // this is OK without -O2 - i.e if the function is NOT defined as an intrinsic function for the compiler
    //prev_pix_num = (long)tmp;
    if (fabs(prev_pix_numD - round(prev_pix_numD)) < 1e-10) { prev_pix_num = (long int)round(prev_pix_numD); } else { prev_pix_num = (long int)floor(prev_pix_numD); } // this is done so because the floor function doesn't return the argument value in case when the argument is integer already and it should !!! that's why if the argument is integer to within 1e-10 then it uses round function which gives the right value, but if the value is not integer (in this case it can be close to x.5, then floor works ok. In polar zones it's not a problem since there are only distances of type x.5.
    j += prev_pix_num;  // plus pixels in the equatorial zone with smaller phi.
    //if ((long)prev_pix_numD != prev_pix_num)
    //printf("i=%li prev_pix_num=%li prev_pix_numD=%.20lf floor(prev_pix_num)=%lf  ||||||| X=%lE phi=%lE th=%lE\n",i,prev_pix_num,prev_pix_numD,tmp,X,phi*180/PI,th*180/PI);
  }

  if ((region >= 8) && (y <1.0)) { // south polar zone
    fullreg_pix = subregion*(long)round((double)pixinring/4.0); // plus number of pixels in this ring in all previous regions
    partreg_pixD = fabs(x+y)*(double)nside/2.0;  tmp = floor(partreg_pixD);
    partreg_pix = (long int)tmp; // plus number of pixels in this region
    j += fullreg_pix + partreg_pix;
/*     printf("i=%li, fullreg_pix=%li, partreg_pix=%li, partreg_pixD=%lf\n",i,fullreg_pix, partreg_pix,partreg_pixD); */
/*     printf("region=%li subregion=%li ring=%li, pixinring=%li x=%lE y=%lE\n",region,subregion,ring,pixinring,x,y); */

  }

  return j;
}


// returns the total number of pixels in all rings above a given ring
long cpeds_get_pix_num_above_ring_healpix(long nside, long ring) {
  long i,pix_num_above = 0;

//  ring--;
  for (i=0;i<ring;i++) { pix_num_above += cpeds_get_pixnum_in_ring_healpix(nside,i); }
  return pix_num_above;
}

long cpeds_get_healpix_pix_num(long nside) {
  return 12*nside*nside;
}
/**********************************************************************************/
long* cpeds_get_pixnum_in_ring_healpix(long int nside, long* N) {
	*N=cpeds_get_ring_num_healpix(nside);
	long *t = new long[*N];

	for (int i = 0; i < *N; i++) {
		t[i]=cpeds_get_pixnum_in_ring_healpix(nside,i);
	}
	return t;
}

/**********************************************************************************/
double* cpeds_get_ring_colatitudes_healpix(long int nside, long* N) {
	long *r=cpeds_get_pixnum_in_ring_healpix(nside,N);
	double *t=new double[*N];
	double th,phi;
	long j=0;
	for (int i = 0; i < *N; i++) {
		cpeds_pix2ang_healpix(nside, j,&th, &phi, 0);
		t[i]=th;
		j+=r[i];
	}
	return t;
}

//double * cpeds_mk_ThetaVals_vect(long nside,














/**********************************************************************************/





























/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/


/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
/*                       Tetrahedron Based Pixlization Scheme (TeBaPix)           */
//Isoangular - roughly - but better to implement isocahedron system by M. Tegmark
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

/* // BASIC DEFINITIONS */

/* // the coordinate system is spherical with theta=0 at x=0,y=0,z=sqrt(3) and this is the  */
/* // top of the tetrahedron */

/* typedef struct { */
/*   double x,y,z,th,phi; */
/* } cpeds_vertex; */

/* typedef struct {   */
/*   cpeds_vertex v[4]; */
/* } cpeds_tetrahedron; */

/* // definition of thetrahedron vertecses */
/* cpeds_tetrahedron cpeds_tetra; */
/* cpeds_tetra.v[0].x = 0; cpeds_tetra.v[0].y = 0; cpeds_tetra.v[0].z = sqrt(3); cpeds_tetra.v[0].th = 0; cpeds_tetra.v[0].phi = 0;  */
/* cpeds_tetra.v[1].x = 0; cpeds_tetra.v[1].y = sqrt(8/3); cpeds_tetra.v[1].z = -1/sqrt(3); cpeds_tetra.v[1].th = 0; cpeds_tetra.v[1].phi = 0;  */
/* cpeds_tetra.v[2].x = -sqrt(2); cpeds_tetra.v[2].y = -sqrt(2/3); cpeds_tetra.v[2].z = -1/sqrt(3); cpeds_tetra.v[2].th = 0; cpeds_tetra.v[2].phi = 0;  */
/* cpeds_tetra.v[3].x = sqrt(2); cpeds_tetra.v[3].y = -sqrt(2/3); cpeds_tetra.v[3].z = -1/sqrt(3); cpeds_tetra.v[3].th = 0; cpeds_tetra.v[3].phi = 0;  */




/* // methods */

/* long int cpeds_get_pix_num_tebapix(long int nfold) { */
/*   return 4*powf(2,(double)(nfold-1)); */
/* } */

/* //cpeds_ */





































































void mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;

  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
/*       ID = (int)fmod(J,2); */
      ID = J % 2;
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }

}




void ang2pix_nest( const long nside, double theta, double phi, long *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off
   * will be consistent for every resolution
   * =======================================================================
   */

  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;

  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  if( theta<0. || theta>pi ) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  if( !setup_done ) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }

  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */

  if( za<=z0 ) { /* equatorial region */

    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */

    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;

    if( ifp==ifm ) face_num = (int)(ifp %4 ) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)(ifp % 4 ); /* (half-)faces 0 to 3 */
    else face_num = (int)(ifm %4 ) + 8;           /* (half-)faces 8 to 11 */

    ix = (int)(jm % ns_max);
    iy = ns_max - (int)(jp % ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */

    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */

    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp );

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);

    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }

  ix_low = (int)(ix % 128);
  ix_hi  =     ix/128;
  iy_low = (int)(iy % 128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow((double)(ns_max/nside),(int)2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*nside*nside); /* in {0, 12*nside**2 - 1} */
}







void pix2ang_nest( long nside, long ipix, double *theta, double *phi) {

  /*
    c=======================================================================
    subroutine pix2ang_nest(nside, ipix, theta, phi)
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (NESTED)
    c     for a parameter nside
    c=======================================================================
  */

      int npix, npface, face_num;
      int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
      int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
      double z, fn, fact1, fact2;
      double piover2=0.5*M_PI;
      int ns_max=8192;

      static int pix2x[1024], pix2y[1024];
      //      common /pix2xy/ pix2x, pix2y

      int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
      //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
      //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
      jrll[0]=2;
      jrll[1]=2;
      jrll[2]=2;
      jrll[3]=2;
      jrll[4]=3;
      jrll[5]=3;
      jrll[6]=3;
      jrll[7]=3;
      jrll[8]=4;
      jrll[9]=4;
      jrll[10]=4;
      jrll[11]=4;
      jpll[0]=1;
      jpll[1]=3;
      jpll[2]=5;
      jpll[3]=7;
      jpll[4]=0;
      jpll[5]=2;
      jpll[6]=4;
      jpll[7]=6;
      jpll[8]=1;
      jpll[9]=3;
      jpll[10]=5;
      jpll[11]=7;


      if( nside<1 || nside>ns_max ) {
        fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
        exit(0);
      }
      npix = 12 * nside*nside;
      if( ipix<0 || ipix>npix-1 ) {
        fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
        exit(0);
      }

      /* initiates the array for the pixel number -> (x,y) mapping */
      if( pix2x[1023]<=0 ) mk_pix2xy(pix2x,pix2y);

      fn = 1.*nside;
      fact1 = 1./(3.*fn*fn);
      fact2 = 2./(3.*fn);
      nl4   = 4*nside;

      //c     finds the face, and the number in the face
      npface = nside*nside;

      face_num = ipix/npface;//  ! face number in {0,11}
      ipf = (int)(ipix % npface);//  ! pixel number in the face {0,npface-1}

      //c     finds the x,y on the face (starting from the lowest corner)
      //c     from the pixel number
      ip_low = (int)(ipf % 1024);//       ! content of the last 10 bits
      ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
      ip_med = (int)(ip_trunc % 1024);//  ! content of the next 10 bits
      ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

      //c     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

      //c     computes the z coordinate on the sphere
      //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
      jr =  jrll[face_num]*nside - jrt - 1;
      //      cout << "face_num=" << face_num << endl;
      //      cout << "jr = " << jr << endl;
      //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
      //      cout << "----------------------------------------------------" << endl;
      nr = nside;//                  ! equatorial region (the most frequent)
      z  = (2*nside-jr)*fact2;
      kshift = (int)((jr - nside) % 2);
      if( jr<nside ) { //then     ! north pole region
         nr = jr;
         z = 1. - nr*nr*fact1;
         kshift = 0;
      }
      else {
        if( jr>3*nside ) {// then ! south pole region
         nr = nl4 - jr;
         z = - 1. + nr*nr*fact1;
         kshift = 0;
        }
      }
      *theta = acos(z);

      //c     computes the phi coordinate on the sphere, in [0,2Pi]
      //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
      jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
      if( jp>nl4 ) jp = jp - nl4;
      if( jp<1 )   jp = jp + nl4;

      *phi = (jp - (kshift+1)*0.5) * (piover2 / nr);

}


void mk_pix2xy(int *pix2x, int *pix2y) {

  /* =======================================================================
   * subroutine mk_pix2xy
   * =======================================================================
   * constructs the array giving x and y in the face from pixel number
   * for the nested (quad-cube like) ordering of pixels
   *
   * the bits corresponding to x and y are interleaved in the pixel number
   * one breaks up the pixel number by even and odd bits
   * =======================================================================
   */

  int i, kpix, jpix, IX, IY, IP, ID;
  for (i = 0; i < 1023; i++) pix2x[i]=0;

  for( kpix=0;kpix<1024;kpix++ ) {
    jpix = kpix;
    IX = 0;
    IY = 0;
    IP = 1 ;//              ! bit position (in x and y)
    while( jpix!=0 ){// ! go through all the bits
      ID = (int)(jpix % 2);//  ! bit value (in kpix), goes in ix
      jpix = jpix/2;
      IX = ID*IP+IX;

      ID = (int)(jpix % 2);//  ! bit value (in kpix), goes in iy
      jpix = jpix/2;
      IY = ID*IP+IY;

      IP = 2*IP;//         ! next bit (in x and y)
    }

    pix2x[kpix] = IX;//     ! in 0,31
    pix2y[kpix] = IY;//     ! in 0,31
  }

}
