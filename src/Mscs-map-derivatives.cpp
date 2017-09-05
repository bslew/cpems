#include "Mscs-map-derivatives.h"

/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* Implementation of the map_derivatives class */
/* **************************************************************************************************** */
/* **************************************************************************************************** */

/* **************************************************************************************************** */
mapDerivative::mapDerivative(bool th, bool phi, bool th2, bool phi2, bool thphi, map_class* map) {

  object_name = "map_derivatives";
  TH = th;
  PHI = phi;
  TH2 = th2;
  PHI2 = phi2;
  THPHI = thphi;
  m = map;

  if (m==NULL) { printf("|%s> * The map given is NULL. EXITING. \n",object_name.c_str()); exit(0); }
  if (!m->mapLoaded()) { printf("|%s> * The map is not loaded. EXITING. \n",object_name.c_str()); exit(0); }
  if (m->maskLoaded()) masked = true; else masked = false;
  if (!m->coordLoaded()) m->set_map_coord(0,0);
  if (!m->isRing()) m->conv_nest2ring(m->map);

  // set the rows number: assuming that the map is in the equilateral pixelization system.
  rownum = m->get_rows_num();
  rownumleo=rownum-1;
  currow=0;
  curpix=0;

  ringbeg=new long[rownum];
  if (m->pix_system == PIXSYS_healpix) { resolve_pix_rings(); }

}

/* **************************************************************************************************** */
mapDerivative::~mapDerivative() {
  delete [] ringbeg;
}

/* **************************************************************************************************** */
void mapDerivative::set_accuracy(long N) {
  ACC = N;
}


/* **************************************************************************************************** */
matrix<double>* mapDerivative::get_allDrow(long N) {

  // check the row number
  if (N<0 || N>=rownum) return NULL;

  if (N==0) { // do the first row
    // get the first row
    currow=0;
    get_row(currow, X1, &Y1, Z1, M1, &S1);
    // get 2nd row
    currow++;
    get_row(currow, X3, &Y3, Z3, M3, &S3);


    //
    // interpolate (multi)mask using linear interpolation
    //
    if (S1<=S3) { 
      //intepolate
      M1int = cpeds_interpolate(X1,M1,S1,X3,S3,"linear");
      M3int = M3;
      // regularize the (multi)mask after interpolation
      correct_mask(M1int,S3);

      // interpolate Z using cubic spline only outside of the mask
      Z1int = cpeds_interpolate(X1,Z1,S1,X3,S3,"cspline_periodic");
      Z3int = M3;

    } 
    else {
      //intepolate
      M3int = cpeds_interpolate(X3,M3,S3,X1,S1,"linear");
      M1int = M1;
      // regularize the (multi)mask after interpolation
      correct_mask(M3int,S1);
    }


    

    //
    // extend the mask
    //

    //
    // calculate derivatives
    //

    // formulate the answer

    // store the current state
  } 
  else {
    if (N < rownumleo) {



    } 
    else { // do the last row

    }
  }
  
}



/* **************************************************************************************************** */

void mapDerivative::get_row(long N, double* X, double* Y, double *Z, double* M, long* S ) {
  long i;
  long fromPix=ringbeg[currow];
  long toPix=ringbeg[currow+1];
  direction n;

  curpixinrow = toPix-fromPix;

  X = new double[curpixinrow];
  Y = new double[curpixinrow];
  Z = new double[curpixinrow];
  M = new double[curpixinrow];

  for (i=fromPix; i<toPix;i++) {
    n=m->get_C(i); X[i]=n.l; 
    Z[i]=m->get_T(i);
    M[i]=m->get_m(i);
  }
  Y=PIsnd-n.b;
  *S=curpixinrow;  
}


/* **************************************************************************************************** */
void mapDerivative::resolve_pix_rings() {
  long i,j;

  // precalculate the cumulative number of pixels ring by ring and store in pix_num_cuml
  ringbeg[0] = 0;  

  for (i=1;i<rownum;i++) { 
    j = cpeds_get_pixnum_in_ring_healpix(m->nside,i-1);
    ringbeg[i] = ringbeg[i-1] + j;
  }
  
}

/* **************************************************************************************************** */
void mapDerivative::correct_mask(double* M, long S) {
  long i;

  for (i=0;i<S;i++) { M[i]=round(M[i]); }

}
