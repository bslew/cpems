#define _ISOC99_SOURCE
#include <features.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <cpgplot.h>
#include "matrix.h"
#include "Mscs-topo-dodecahedron.h"
//#include "cpeds-consts.h"
//#include "cpeds-math.h"

// CONVENTIONS:

// MAP CONVENTION
// the map object passed to this object is assumed to be in ring ordering
// and with seeded coordinates and optionally with a mask

// RADIAN AND DEGREES
// the dodecahedron parameters from outside are given in radians

// COORDINATE SYSTEMS
// the coordinate system throughout this object is consistent with the Mscs package (l,b) however
// the rotation matrixes uses the healpix like CS. thetaphi (TP) for rotations; but the I/O for these methods is (l,b)

// CIRCLE NUMBERING
// the circles are numberes from 0 to 11. 0th is the top, 1-5 upper 5, 6-bottom, 7-11 lower 5
// the matching circles are then: 0-6, 1-7, 2-8, 3-9, 4-10 and 5-11.

/********************************************************************************/
/* dodecahedron::dodecahedron(double a,double s, map_class * T) { // creates a dodec. in "0" reference position */
/*   map = T; */
/*   nside = map->nside; */
/*   pix_num = map->pix_num; */
/*   ring_num = cpeds_get_ring_num_healpix(nside); */
/*   circ_pix_num = 2*((long)round((double)nside*a/22.5)); // must be even to match the points 180 deg apart in phase */

/*   printf("|dodec>  -- nside: %li\n",nside); */
/*   printf("|dodec>  -- ring_num: %li\n",ring_num); */
/*   printf("|dodec>  -- circ_pix_num: %li\n",circ_pix_num); */

/*   initiate_dodecahedron_structures(""); */
  
/*   define_Dfaces0(); // define the faces in zero position in l,b coords. */

/*   Dl=0; Db=0; Dg=0; Da=PI/180*a; Ds=PI/180*s; // initiate dodec. params. */
/*   Dth=PIsnd-Db; */

/*   derive_D0 (Da,Ds); // find the set of points on circles at "0" orientation of dodec */

/* } */
/********************************************************************************/
// amax parameter defines the maximal size of the circle in the search. if it is given then all circles will be probed with the resolution corresponding to the size
// of this circle regardless of their size. if it's negative then the number of pixels per circle is derived dynamically to each circle size a according to 
// the formula defined in the init_dodecahedron_structures()
// rho - parameter scalling the number of pixels for probing the circles.
// rho = 8 - a value estimated for wmap_ilc map with no extra smoothing via FFT analysis of the D(0;90;0;10;0) which implied that around 108 pixels per circle of 
//           size 10deg are sufficient to reconstruct the fluctuations around this circle.
//           for 2deg smoothed map 31 pixels should be enough for circle of size 10 deg. and map ns=256
//           for other map smoothings some other value should be taken.

dodecahedron::dodecahedron(double l, double b, double g, double a, double s, map_class * T, double amax) { // creates a dodec with requested parameters

  nside = T->nside;
  pix_num = T->pix_num;
  ring_num = cpeds_get_ring_num_healpix(nside);
  Drho=initiate_rho(); // rho is initiated once only and should be fitted optimally to the data being tested 

  Dl=l; Db=b; Dg=g; Da=a; Ds=s; if (amax > 0) { Damax=amax; Damax_b = true; } else { Damax = a; Damax_b = false; } // initiate dodec. params.
  
  //circ_pix_num = 2*((long)round((double)nside*a/22.5)); // must be even to match the points 180 deg apart in phase
/*   circ_pix_num = 16*((long)round((double)nside*a/PI)); // must be even to match the points 180 deg apart in phase */

  printf("|dodec>  -- nside: %li\n",nside);
  printf("|dodec>  -- ring_num: %li\n",ring_num);
  initiate_dodecahedron_structures("all",T);
  printf("|dodec>  -- circ_pix_num: %li\n",circ_pix_num);

  define_Dfaces0();

  //Dl=PI/180*l; Db=PI/180*b; Dg=PI/180*g; Da=PI/180*a; Ds=PI/180*s; // initiate dodec. params.
  
  derive_D0(a,s); // find the set of points on circles at "0" orientation of dodec
  //exit(0);
  copyD(D,D0); copyD(Drotg,D); 
  Dl=0; Db=PIsnd; Dg=0; Da=a; Ds=s; // initiate dodec. params.
  Dth=PIsnd-Db;
  derive_D(l,b,g,a,s); //  find the set of points on circles in  requested orientation of dodec

}

/********************************************************************************/
dodecahedron::~dodecahedron() { // creates a dodec. in "0" reference position
  kill_dodecahedron_structures();

  if (btab != NULL) delete [] btab;
  if (thtab != NULL) delete [] thtab;
  if (ltab != NULL) delete [] ltab;
  if (pirtab != NULL) delete [] pirtab;
  if (ringtab != NULL) delete [] ringtab;
  delete D0;
  delete D;
  delete Drotg;
  delete Di;

  for (i=0;i<3;i++) {   delete [] S[i]; }
  delete [] S;
  delete [] w;
  delete [] map; delete [] mask; delete [] msq;
}
/********************************************************************************/
void dodecahedron::kill_dodecahedron_structures() {
  for (i=0;i<12;i++) {
    if (D0[i] != NULL) { delete [] D0[i]; D0[i] = NULL; }
    if (D[i] != NULL)  { delete [] D[i]; D[i] = NULL; }
    if (Drotg[i] != NULL)  { delete [] Drotg[i]; Drotg[i] = NULL; }
    if (Di[i] != NULL) { delete [] Di[i]; Di[i] = NULL; }
  }
}
/********************************************************************************/
void dodecahedron::init_dodecahedron_structures(double a) {
  double loca;
  double cpnd;
  //loca = a;
  //loca = 20*PI/180; // setup as constant circle size to get same amount of pixels for different circles, if this works should be given as parameter to the object (e.g. a_pix_num_max)

  if (Damax_b) loca = Damax; else loca = a;
/*   circ_pix_num = rho*((long)round((double)nside*loca/PI)); // must be even to match the points 180 deg apart in phase;  */

  // this is the fit the measurments on 20 sims. of convergence of variances of S values to Sinf value. 
  //the line was fitted so that 1% acurracy at 1 sigma CL was obtained. for rho=8 ns=256 for all circle sizes
  cpnd = (double)rho/32.0 * (double)nside/256.0 * (loca*180.0/PI * 3.39725 + 76.8495);
  circ_pix_num = (long)round(cpnd); // must be even to match the points 180 deg apart in phase; 
  // if not even then set to closest even number 
  if (circ_pix_num % 2 == 1) { if (cpnd-(double)((circ_pix_num-1)/2) < (double)((circ_pix_num+1)/2)-cpnd) circ_pix_num--; else circ_pix_num++; }
  if (circ_pix_num == 0) circ_pix_num=1; //safety trigger
  printf("  -- setting circ_pix_num to: %li for circle size %lf [deg] and nside: %li\n",circ_pix_num,loca*180/PI, nside);

  //circ_pix_num = 16*((long)round((double)nside*loca/PI)); // must be even to match the points 180 deg apart in phase
  kill_dodecahedron_structures();
  init_dodecahedron_Dstructures();
}
/********************************************************************************/
/* void dodecahedron::init_dodecahedron_structures(double a, long rho) { */
/*   double loca; */
  
/*   //loca = a; */
/*   //loca = 20*PI/180; // setup as constant circle size to get same amount of pixels for different circles, if this works should be given as parameter to the object (e.g. a_pix_num_max) */
/*   if (Damax_b) loca = Damax; else loca = a; */
/*   circ_pix_num = rho*((long)round((double)nside*loca/PI)); // must be even to match the points 180 deg apart in phase; 8 - a value estimated for wmap_ilc map with no extra smoothing via FFT analysis of the D(0;90;0;10;0) which implied that around 108 pixels per circle of size 10deg are sufficient to reconstruct the fluctuations around this circle. */
/*   // for 2deg smoothed map 31 pixels should be enough for circle of size 10 deg. */
/*   // for other map smoothings some other value should be taken. */
/* /\*     circ_pix_num = 16*((long)round((double)nside*loca/PI)); // must be even to match the points 180 deg apart in phase *\/ */

/*   init_dodecahedron_structures2(); */
/* } */

/********************************************************************************/
void dodecahedron::init_dodecahedron_Dstructures() {
  for (i=0;i<12;i++) {
    if (D0[i] == NULL) D0[i] = new direction[circ_pix_num];
    if (D[i] == NULL) D[i] = new direction[circ_pix_num]; 
    if (Drotg[i] == NULL) Drotg[i] = new direction[circ_pix_num]; 
    if (Di[i] == NULL) Di[i] = new long[circ_pix_num]; 
  }
}

/********************************************************************************/
long dodecahedron::initiate_rho() { rho = 16; return rho;}
long dodecahedron::initiate_rho(long rho_loc) {
  if (rho_loc == 0) { 
    initiate_rho();  }
  else {
    if (rho_loc == 1) { //  perform some test first on the data to find the best value
      printf("WARNING !! not implemented yet.Will use the default value\n");
      initiate_rho();    } 
  }
  return rho;
}
/********************************************************************************/
void dodecahedron::set_rho(long rho_loc) { rho = rho_loc;  kill_dodecahedron_structures(); init_dodecahedron_structures(Da); }
/********************************************************************************/
// handlers
double dodecahedron::get_Dl() { return Dl; }
double dodecahedron::get_Db() { return Db; }
double dodecahedron::get_Dg() { return Dg; }
// returns the orientation parameter of the dodec given tirections of the two adjacent face centers. the angle between them must be 63....[deg]  as for dodec
double dodecahedron::get_Dg(double l1, double b1, double l2, double b2) { 
  direction n1,n2;
  // find the dodec params we need
  n1.l=l1; n1.b=b1; n2.l=l2; n2.b=b2;
  n1=RyRz(n1,-PIsnd+b1,-l1);  n2=RyRz(n2,-PIsnd+b1,-l1);  
  return n2.l;
}
double dodecahedron::get_Da() { return Da; }
double dodecahedron::get_Ds() { return Ds; }

/********************************************************************************/
/* void dodecahedron::calculate_circ_pix_num(long rho, double a) { */
/* } */
/********************************************************************************/
void dodecahedron::initiate_dodecahedron_structures(string what, map_class * T) {
  direction tmp;
  double b;
  long k;

  // create arrays for circles in 0 position, arbitrary position and in pixel number representation
  D0 = new direction*[12];
  D = new direction*[12];
  Drotg = new direction*[12];
  Di = new long*[12];
  w = new double[7];
  for (k=0;k<12;k++) {    D0[k] = NULL; D[k] = NULL; Drotg[k] = NULL; Di[k] = NULL;  }

  Clifford_rot = 36.0*PI/180;
  Clifford_rot2 = 2.0*Clifford_rot;

  init_dodecahedron_structures(Da);
  //zero statistics arrays
  S = new double*[3];
  for (j=0;j<3;j++) { 
    S[j] = new double[8];
    for (i=0;i<8;i++) { S[j][i] = 0;  } //K[i] = 0; }
  }

  // initiate help arrays
  btab = new double[ring_num]; // array with all theta values in the map 
  thtab = new double[ring_num]; // array with all theta values in the map 
  ltab = new double[pix_num]; // array with all l values in the map
  pirtab = new long[ring_num]; // array with number of pixels in each ring in the map
  ringtab = new long[ring_num]; // array with indexes of pixels that start a given ring
  map = new double[pix_num];
  msq = new double[pix_num];
  mask = new double[pix_num];

  tmp = T->get_C(0);       
  b = tmp.b; k=0; btab[k] = b; thtab[k] = PIsnd-b; pirtab[k] = 4; ringtab[0] = 0;
  for (i=0;i<pix_num;i++) { 
    tmp = T->get_C(i);    ltab[i] = tmp.l; // copy all l vals. from map
    //printf("i=%li ltab[i]=%lf  \n",i,180/PI*ltab[i]);
    if (b != tmp.b) {
      //printf("k=%li btab[k]=%lf  pirtab=%li ringtab=%li\n",k,btab[k],pirtab[k],ringtab[k]);
      k++; b=tmp.b; btab[k]=b; thtab[k]=PIsnd-b;// copy all b vals. from map
      pirtab[k] = cpeds_get_pixnum_in_ring_healpix(nside,k); // derive pix.numers in rings
      ringtab[k] = ringtab[k-1] + pirtab[k-1];
    }
    //import the map
/*     map[i] = T->get_T(i)/sqrt(T->varianceT); */
    map[i] = T->get_T(i)/T->meanT;
/*     map[i] = (double)i; // DEBUG thing */
    mask[i] = T->get_m(i); // import the mask; if no mask was given this will make a transparent mask
    map[i] *= mask[i]; // map - mask merge - in case the mask wasn't given in the command line
    msq[i] = map[i]*map[i]; // calculate the T^2*mask 
/*     printf("map[i]=%lf msq[i]=%lf, mask[i]=%lf\n",map[i],msq[i],mask[i]); */
  }

  T->kill_map();
  //exit(0);

}
/********************************************************************************/
void dodecahedron::define_Dfaces0() {  // definition of dodecahedron faces centers
  Dface0[0].l = 0;  Dface0[0].b = PIsnd; // 0'th face
  for (i=1;i<6;i++) { Dface0[i].l = (double)(i-1)*72.0*PI/180;  Dface0[i].b = PIsnd-63.434948822922*PI/180;  } // 1-5 faces
  Dface0[6].l = 0;  Dface0[6].b = -PIsnd; // 6'th face
  for (i=7;i<12;i++) { Dface0[i].l = ((double)((i-5)*72%360)+36.0)*PI/180;  Dface0[i].b = PIsnd-116.565051177*PI/180;  } // 7'th -- 11'th faces

  for (i=0;i<12;i++) {
    if (i<6) printf("D0faces: nr:%li, l: %lf b: %lf\n",i,Dface0[i].l*180/PI,Dface0[i].b*180/PI);
    if (i>=6) printf("D0faces: nr:%li, l: %lf b: %lf ang: %lf\n",i,Dface0[i].l*180/PI,Dface0[i].b*180/PI,180/PI*cpeds_ang_n1n2(PIsnd-Dface0[i].b,Dface0[i].l,PIsnd-Dface0[i-6].b,Dface0[i-6].l));
  }
}

/********************************************************************************/
// shift points on circles in arrays by requested phase p.
void dodecahedron::match_circles(direction **D, double p) {
  long j;
  long shift,ien;
  direction * tmp = new direction[circ_pix_num];

  // copy circles with: b --> -b, l --> l+180
  for (j=6;j<12;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i].b = -D[j-6][i].b; if (D[j-6][i].l < PI) D[j][i].l = PI + D[j-6][i].l; else D[j][i].l = D[j-6][i].l - PI; }

  // phase shift 180 + phase shift p
/*   shift = (long)(round((double)circ_pix_num*p/twoPI+(double)(circ_pix_num/2))); */
  shift = (long)(round((double)circ_pix_num*(p+PI)/twoPI));
  if (shift < 0) {shift = circ_pix_num + shift; }
  if (shift > circ_pix_num) { shift=shift-circ_pix_num; }
  ien = circ_pix_num-shift; 
  for (j=6;j<12;j++) {
    for (i=0;i<circ_pix_num;i++) { tmp[i] = D[j][i]; } // copy the circle
    for (i=0;i<ien;i++) { D[j][i] = tmp[i+shift]; }    for (i=ien;i<circ_pix_num;i++) { D[j][i] = tmp[i-ien]; } // make shift
  }
  delete [] tmp;
/*   for (i=0;i<circ_pix_num;i++) { printf("%lf,%lf ",D[3][i].l,D[3][i].b); } printf("\n"); // copy the circle */
/*   for (i=0;i<circ_pix_num;i++) { printf("%lf,%lf ",D[9][i].l,D[9][i].b); } printf("\n"); // copy the circle */

}

/********************************************************************************/
// shift points on circles in arrays by requested phase p.
void dodecahedron::shift_circlesi(double p) {
  long j;
  long shift,ien;
  double * tmp = new double[circ_pix_num];


  // phase shift p
  shift = (long)(round((double)circ_pix_num*p/twoPI));
/*   printf("shift by: %li\n",shift); */
  if (shift < 0) {shift = circ_pix_num + shift; }
  if (shift > circ_pix_num) { shift=shift-circ_pix_num; }
/*   printf("shift by: %li\n",shift); */
  ien = circ_pix_num-shift; 
  for (j=6;j<12;j++) {
    for (i=0;i<circ_pix_num;i++) { tmp[i] = Di[j][i]; } // copy the circle
    for (i=0;i<ien;i++) { Di[j][i] = tmp[i+shift]; }    for (i=ien;i<circ_pix_num;i++) { Di[j][i] = tmp[i-ien]; } // make shift
  }
  delete [] tmp;

}


/********************************************************************************/
void dodecahedron::derive_D0(double a, double s) {  //find the set of points on circles at "0" orientation of dodec
  double phi,dphi;
  kill_dodecahedron_structures();
  init_dodecahedron_structures(a);

  dphi = twoPI/(double)circ_pix_num;


  // top circle: C_0 = C(th=a,phi)
  phi = 0;
  for (i=0;i<circ_pix_num;i++) { D0[0][i].l = phi;    D0[0][i].b = PIsnd-a; phi+=dphi;  }//printf("---------a:%lE, l1:%lE b1:%lE\n",a,D0[0][i].l,D0[0][i].b);}

  // upper 1st circle: C_1 = Rz(Dface0[1].l)Rx(90-Dface0[1].b) C_0
  //for (i=0;i<circ_pix_num;i++) { D0[1][i] = RzRx(D0[0][i],Dface0[1].l,PIsnd-Dface0[1].b); }
  for (i=0;i<circ_pix_num;i++) { D0[1][i] = RzRy(D0[0][i],Dface0[1].l,PIsnd-Dface0[1].b); }

  // upper i=2,3,4,5 circles: C_i = Rz(Dface0[i].l) C_1
  for (j=2;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D0[j][i] = Rz(D0[1][i],Dface0[j].l); }

  // bottom circle: C_6 = Rz(s) C(pi-a,phi)
/*   phi = 0; */
/*   for (i=0;i<circ_pix_num;i++) { D0[6][i].l = phi+s;    D0[6][i].b = -PIsnd+a; phi+=dphi; } // this is temporary */

  // lower 1st circle:  C_7 = Rz(Dface0[1].l) Rx(90-Dface0[1].b) Rz(s) C_0(a,-phi)
/*   for (i=0;i<circ_pix_num;i++) { D0[7][i] = RzRyRz(D0[0][(circ_pix_num-i)%circ_pix_num],Dface0[7].l,PIsnd-Dface0[7].b,-s); } */

  // lower i=8,9,10,11 circle: C_i = Rz(Dface0[i].l C_7(th,phi)
/*   for (j=8;j<12;j++)   */
/*     for (i=0;i<circ_pix_num;i++) { D0[j][i] = Rz(D0[7][i],Dface0[j-6].l); } */

  //match_circles(D0,36*PI/180); // shift the arrays of the matched circles (6-11) by PI
  match_circles(D0,s); // shift the arrays of the matched circles (6-11) by PI


/*   for (j=0;j<12;j++) { */
/*     printf("D0: circle: %li \n",j); */
/*     for (i=0;i<circ_pix_num;i++) {      */
/*       if (j<6) printf("l:%lf b: %lf\n",180/PI*D0[j][i].l, 180/PI*D0[j][i].b);  */
/*       if (j>=6 && i<circ_pix_num-1) printf("l:%lf b: %lf ang:%lf ang2: %lf\n",180/PI*D0[j][i].l, 180/PI*D0[j][i].b,180/PI*cpeds_ang_n1n2(PIsnd-D0[j][i].b,D0[j][i].l,PIsnd-D0[j-6][i].b,D0[j-6][i].l),180/PI*cpeds_ang_n1n2(PIsnd-D0[j][i].b,D0[j][i].l,PIsnd-D0[j-6][i+1].b,D0[j-6][i+1].l));} */
/*   } */

  //Da=a; Ds=s;
  //derive_D(0,0,0,Da,Ds);
}
/********************************************************************************/
//find the set of points on circles int requested orientation of dodec given dodec params.
// D(Db,Dl,Dg,Da,Ds) = Rz(Dl)Rx(PIsnd-Db)Rz(Dg) D0(0,0,0,Da,Ds)
void dodecahedron::derive_D(double l, double b, double g, double a,double s) {  
/*   derive_D0(a,s); */
/*   copyD(D,D0); */
/*   RzRyRzD(D,l,PIsnd-b,g); */
 
/*   if (Ds != s || Da != a) { derive_D0(a,s); copyD(D,D0); RzD(D,g); copyD(Drotg,D); RzRyD(D,l,PIsnd-b); } */
/*   else { */
/*     if (Dg != g) { copyD(D,D0); RzD(D,g); copyD(Drotg,D); RzRyD(D,l,PIsnd-b); } */
/*     else { */
/*       if (Db != b) { copyD(D,Drotg); RzRyD(D,l,PIsnd-b); } */
/*       else { */
/* 	if ( Dl != l) { RzD(D,l-Dl); } */
/*       } */
/*     } */
/*   } */


// current implementation
  if (Ds != s || Da != a || Drho != rho) { derive_D0(a,s); copyD(D,D0); RzRyRzD(D,l,PIsnd-b,g); }
  else {
    if (Dg != g || Db != b) { copyD(D,D0); RzRyRzD(D,l,PIsnd-b,g); }
    else {
      if ( Dl != l) { RzD(D,l-Dl); }
    }
  }

// experimental implementation begin -- not finished
/*   if (Da != a || Drho != rho) { derive_D0(a,s); copyD(D,D0); RzRyRzD(D,l,PIsnd-b,g); } */
/*   else { */
/*     if (Dg != g || Db != b) { copyD(D,D0); RzRyRzD(D,l,PIsnd-b,g); } */
/*     else { */
/*       if ( Dl != l) { RzD(D,l-Dl); } */
/*     } */
/*   } */
// experimental implementation end

  Dl=l; Db=b; Dg=g; Da=a; Ds=s; Drho = rho;
  Dth = PIsnd-Db;

  derive_Di();
}
/********************************************************************************/
//find the set of points on circles in requested orientation of dodec given by two of it's face centers
void dodecahedron::derive_D(double l1, double b1, double l2, double b2, double a,double s) {
  double l,b,g;
  direction n1,n2;

  // find the dodec params we need
  n1.l=l1; n1.b=b1; n2.l=l2; n2.b=b2;
/*   printf("----------->n1 %lE %lE\n",180/PI*n1.l,180/PI*n1.b); */
/*   printf("----------->n2 %lE %lE\n",180/PI*n2.l,180/PI*n2.b); */
  n1=RyRz(n1,-PIsnd+b1,-l1);  n2=RyRz(n2,-PIsnd+b1,-l1);
/*   printf("----------->n1 %lE %lE\n",180/PI*n1.l,180/PI*n1.b); */
/*   printf("----------->n2 %lE %lE\n",180/PI*n2.l,180/PI*n2.b); */
  
  l=l1; b=b1; g=n2.l; 

  // give them to the other method
  derive_D(l,b,g,a,s);
}
/********************************************************************************/
// updates the set of D projected onto map pixelization directions
void dodecahedron::derive_Di() {  
  double b,l,th;
  long ring,pix,ring2,pix2;
  double ang,ang2;
/*   double * subltab; */

  for (j=0;j<12;j++)  
    for (i=0;i<circ_pix_num;i++) { 
      l = D[j][i].l;    b = D[j][i].b; th=PIsnd-b;
      ring = cpeds_find_value(th,thtab,ring_num,0,ring_num); if (th <= thtab[ring] && ring > 0) ring2 = ring-1; else ring2 = ring+1;
      if (ring2 == ring_num) ring2--; 
/*       printf("1------> looking for b: %lf getting %lf \n",b,btab[ring]); */
/*       printf("2------> looking for b: %lf getting %lf \n",b,btab[ring2]); */

      pix = cpeds_find_value(l,ltab,pix_num,ringtab[ring],pirtab[ring]); pix2 = cpeds_find_value(l,ltab,pix_num,ringtab[ring2],pirtab[ring2]); 
/*       printf("1------> looking for l: %lf getting %lf \n",l,ltab[pix]); */
/*       printf("2------> looking for l: %lf getting %lf \n",l,ltab[pix2]); */
      ang = cpeds_ang_n1n2(PIsnd-b,l,PIsnd-btab[ring],ltab[pix]); ang2 = cpeds_ang_n1n2(PIsnd-b,l,PIsnd-btab[ring2],ltab[pix2]);
      if (ang < ang2) Di[j][i] = pix; else Di[j][i] = pix2;
/*       printf("*** Di: circ: %li, pix %li, num %li num2 %li ang %lf ang2 %lf Di %li, Ti %lE\n",j,i,pix,pix2,180/PI*ang, 180/PI*ang2,Di[j][i],map[Di[j][i]]); */
    }

}

/********************************************************************************/
// calculates S statistics on map T with dodec. in requested orientation defined by Di (D)
// the orientation given to this method is in RADIANS !!!
double ** dodecahedron::Sstat(double l, double b, double g, double a,double s) {
  derive_D(l,b,g,a,s);
  return Sstat_cur();
}


/********************************************************************************/
// calculates S statistics on map T with dodec. in current orientation defined by Di (D)

double ** dodecahedron::Sstat_cur() {
  double T1, T2, T1T2, v1,v2,v1tot,v2tot;
  double pixoncirc1, pixoncirc2, pixoncirc ,pixoncirctot1, pixoncirctot2, pixoncirctot;
  long k,p1,p2;

  for (k=0;k<3;k++) {
/*     if (k==1) { match_circles(D,Clifford_rot+Ds); derive_Di(); } */
/*     if (k==2) { match_circles(D,-Clifford_rot+Ds); derive_Di(); } */


    // this implementation is much faster and works identically with the above one.
    // caveat is: after such shifts there D and Di do not match. to make them match again use match_dircles and derive_Di routine which does D --> Di
    if (k==1) { shift_circlesi( Clifford_rot); }
    if (k==2) { shift_circlesi(-Clifford_rot); shift_circlesi(-Clifford_rot);} // this is not the same as -Clifford_rot2 ! 

    S[k][7]=0; v1tot=v2tot=0; pixoncirctot=pixoncirctot1=pixoncirctot2=0;
    for (j=0;j<6;j++) {
      S[k][j] = 0;  v1=v2=0;  pixoncirc=pixoncirc1=pixoncirc2=0;
      for (i=0;i<circ_pix_num;i++) { 
	//printf("T1: %li, T2: %li\n",Di[j][i],Di[j+6][i]); 

	p1=Di[j][i]; p2 = Di[j+6][i];	T1T2=map[p1]*map[p2]; 

/* 	S[k][j] += T1T2;  v1+=msq[p1]; v2+=msq[p2]; pixoncirc1 += mask[p1]; pixoncirc2 += mask[p2]; pixoncirc += mask[p1]*mask[p2];  */
	S[k][j] += T1T2;  v1+=msq[p1]*mask[p2]; v2+=msq[p2]*mask[p1]; pixoncirc += mask[p1]*mask[p2];
/* 	printf("------> v1 %lf v2 %lf pixoncirc1 %lf pixoncirc2 %lf pixoncirc %lf S:%lf msq[p1]:%lf msq[p2]:%lf map[p1=%li]:%lf\n",v1,v2,pixoncirc1,pixoncirc2,pixoncirc, S[k][j],msq[p1],msq[p2],p1,map[p1]); */
      }
      S[k][7]+=S[k][j]; v1tot+=v1; v2tot+=v2; //pixoncirctot1 += pixoncirc1; pixoncirctot2 += pixoncirc2; pixoncirctot += pixoncirc;

      // calculation of S for individual circle pairs
      //pixoncirc1--; pixoncirc2--;
/*       if (pixoncirc1 > 0 && pixoncirc2 > 0 && pixoncirc > 0) { v1/=pixoncirc1; v2/=pixoncirc2; S[k][j] = 2*S[k][j]/(v1 + v2)/pixoncirc; } else S[k][j] = 0; */
      if (S[k][j] != 0) { S[k][j] = 2*S[k][j]/(v1 + v2); } 
/*       printf("--- v1 %lf v2 %lf pixoncirc1 %lf pixoncirc2 %lf pixoncirc %lf S:%lf\n",v1,v2,pixoncirc1,pixoncirc2,pixoncirc, S[k][j]); */
      w[j] = pixoncirc;  //  corrected for weightning of uncomplete or missing circles
    }

    // calculation of S as simple average from the six circles -- this is no longer supported since it can be easily calculated from the stored data afterwards
/*     for (j=0;j<6;j++) { S[k][6] += S[k][j]; } S[k][6]/=6.0;       */

    // calculation of S as weight average from the six circles
    S[k][6]=0; w[6]=0;
    for (j=0;j<6;j++) { S[k][6] += S[k][j]*w[j]; w[6]+=w[j]; } S[k][6]/=w[6];      // S stat with correction for the proper weightning proportional to number of pixels in circles due to sky cut

/*     printf("-------> ptot1 %lf ptot2 %lf ptot %lf\n",pixoncirctot1,pixoncirctot2,pixoncirctot); */
    // calculation of S as correlation over all set of matching pixels from all 6 circles
    //pixoncirctot1--; pixoncirctot2--;
/*    if (pixoncirctot > 0 && pixoncirctot1 > 0 && pixoncirctot2 > 0) { v1tot/=pixoncirctot1; v2tot/=pixoncirctot2; S[k][7]=2*S[k][7]/(v1tot+v2tot)/pixoncirctot; } else S[k][7] = 0; */
    if (S[k][7] != 0) { S[k][7]=2*S[k][7]/(v1tot+v2tot); } 
  }

  return S;
}


void dodecahedron::save_dodecahedron_circles(string fname) {
  FILE * tmpf;
  filenamestr tmpch;
  long j,i;
  direction tmp;

  // saving increasing values along circles - for testing purposes
  sprintf(tmpch,"%s.test",fname.c_str());
  tmpf = fopen(tmpch,"w");

  for (j=0;j<12;j++) {  
    for (i=0;i<circ_pix_num;i++) { 
      fprintf(tmpf,"%lE %lE 0.2 %li\n",180/PI*D[j][i].l,180/PI*D[j][i].b,i+((j % 6)+1)*circ_pix_num); 
    }
/*     printf("%lf %lf\n",180/PI*Dface0[j].l,180/PI*Dface0[j].b); */
  }
  fclose(tmpf); 

  // saving face centers
  sprintf(tmpch,"%s.faces",fname.c_str());
  tmpf = fopen(tmpch,"w");

  for (j=0;j<12;j++) {  
    tmp = RzRyRz(Dface0[j],Dl,PIsnd-Db,Dg);
    fprintf(tmpf,"%lf %lf\n",180/PI*tmp.l,180/PI*tmp.b); 
    printf("dodecahedron face %li [Dl=%lf Db=%lf Dg=%lf Da=%lf]: (l,b) = (%lf, %lf)\n",j,180/PI*Dl,180/PI*Db,180/PI*Dg,180/PI*Da,180/PI*tmp.l,180/PI*tmp.b); 

  }
  fclose(tmpf); 


  // saving temperatures along circles
  sprintf(tmpch,"%s.temp",fname.c_str());
  tmpf = fopen(tmpch,"w");

  for (i=0;i<circ_pix_num;i++) { 
    fprintf(tmpf,"%li ",i); 
    for (j=0;j<6;j++)  {  fprintf(tmpf,"%lE ",map[Di[j][i]]*mask[Di[j][i]]);     }
    for (j=6;j<12;j++) {  fprintf(tmpf,"%lE ",map[Di[j][i]]*mask[Di[j][i]]);     }
    fprintf(tmpf,"\n"); 
  }
  fclose(tmpf);

  // saving directions along circles
  sprintf(tmpch,"%s.dirs",fname.c_str());
  tmpf = fopen(tmpch,"w");

  for (i=0;i<circ_pix_num;i++) { 
    fprintf(tmpf,"%li ",i); 
    for (j=0;j<6;j++)  {  fprintf(tmpf,"%lE %lE ",180/PI*D[j][i].l,180/PI*D[j][i].b);     }
    for (j=6;j<12;j++) {  fprintf(tmpf,"%lE %lE ",180/PI*D[j][i].l,180/PI*D[j][i].b);     }
    fprintf(tmpf,"\n"); 
  }
  fclose(tmpf);

}
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
// HELP METHODS

// rotations methods: these work in (th, phi) CS. NOT (l,b) 
// so input is in th,phi and output in th,phi CS
direction dodecahedron::Rz(direction p,double Az) {
  direction pp;
  pp.l = (p.l+Az);  pp.b = p.b;
  cpeds_check_phi(&(pp.l));
  return pp;
}

direction dodecahedron::Rx(direction p,double Ax) {
  direction pp;
  double x,y,z,xx,yy,zz,th,phi,sinAx,cosAx;

  th=PIsnd-p.b; phi = p.l; // change from (l,b) CS to (th,phi) CS
  x = cpeds_sph2cart(0,th,p.l);  y = cpeds_sph2cart(1,th,p.l);   z = cpeds_sph2cart(2,th,p.l);   sinAx=sin(Ax); cosAx=cos(Ax);

  xx=x;
  yy=y*cosAx+z*sinAx;
  zz=-y*sinAx+z*cosAx;

  th = cpeds_cart2sph(0,xx,yy,zz); phi = cpeds_cart2sph(1,xx,yy,zz);
  cpeds_check_thphi(&th,&phi);
  pp.b=PIsnd-th; pp.l=phi; // change from (th,phi) CS to (l,b) CS
  return pp;
}


direction dodecahedron::Ry(direction p,double Ay) {
  direction pp;
  double x,y,z,xx,yy,zz,th,phi,sinAy,cosAy;

  th=PIsnd-p.b; phi = p.l; // change from (l,b) CS to (th,phi) CS
  x = cpeds_sph2cart(0,th,p.l);  y = cpeds_sph2cart(1,th,p.l);   z = cpeds_sph2cart(2,th,p.l);   sinAy=sin(Ay); cosAy=cos(Ay);
/*   printf("xyz: %lf %lf %lf\n",x,y,z); */

  xx=x*cosAy+z*sinAy;
  yy=y;
  zz=-x*sinAy+z*cosAy;

/*   printf("xx yy zz: %lf %lf %lf\n",xx,yy,zz); */
  th = cpeds_cart2sph(0,xx,yy,zz); phi = cpeds_cart2sph(1,xx,yy,zz);
/*   printf("thth phiphi: %lf %lf\n",th,phi); */
  cpeds_check_thphi(&th,&phi);
/*   printf("thth phiphi after check: %lf %lf\n",th,phi); */
  pp.b=PIsnd-th; pp.l=phi; // change from (th,phi) CS to (l,b) CS
  return pp;
}


direction dodecahedron::RzRx(direction p,double Az, double Ax) {
  direction pp;
  pp=Rx(p,Ax);  pp=Rz(pp,Az);
  return pp;
}

direction dodecahedron::RzRy(direction p,double Az, double Ay) {
  direction pp;
  pp=Ry(p,Ay);  
/*     printf("%lf %lf\n",180/PI*Dface0[j].l,180/PI*Dface0[j].b); */
  pp=Rz(pp,Az);
/*     printf("%lf %lf\n",180/PI*Dface0[j].l,180/PI*Dface0[j].b); */
  return pp;
}

direction dodecahedron::RyRz(direction p,double Ay, double Az) {
  direction pp;
  pp=Rz(p,Az);  
/*     printf("%lf %lf\n",180/PI*Dface0[j].l,180/PI*Dface0[j].b); */
  pp=Ry(pp,Ay);
/*     printf("%lf %lf\n",180/PI*Dface0[j].l,180/PI*Dface0[j].b); */
  return pp;
}

direction dodecahedron::RzRxRz(direction p, double Azz,double Ax, double Az) {
  direction pp;
  pp=Rz(p,Az); pp=Rx(pp,Ax);  pp=Rz(pp,Azz);
  return pp;
}

direction dodecahedron::RzRyRz(direction p, double Azz,double Ay, double Az) {
  direction pp;
  pp=Rz(p,Az); pp=Ry(pp,Ay);  pp=Rz(pp,Azz);
  return pp;
}

/********************************************************************************/
// copies the data from D2 to D1 structures. assumes the space is initielized
void dodecahedron::copyD(direction **D1, direction **D2) {
  for (j=0;j<12;j++)  
    for (i=0;i<circ_pix_num;i++) { D1[j][i] = D2[j][i]; }
 }
/********************************************************************************/
// rotates whole dodecahedron around z axix by Az
 void dodecahedron::RzD(direction **D, double Az) {
  for (j=0;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i] = Rz(D[j][i],Az); }
  match_circles(D,Ds); // shift the arrays of the matched circles (6-11) by PI
}
/********************************************************************************/
// rotates whole dodecahedron around z x and z axix by Azz, Ax, Az angles
void dodecahedron::RzRxRzD(direction **D, double Azz, double Ax, double Az) {
  for (j=0;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i] = RzRxRz(D[j][i],Azz,Ax,Az); }
  match_circles(D,Ds); // shift the arrays of the matched circles (6-11) by PI
}
/********************************************************************************/
// rotates whole dodecahedron around z y and z axix by Azz, Ay, Az angles
void dodecahedron::RzRyRzD(direction **D, double Azz, double Ay, double Az) {
  for (j=0;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i] = RzRyRz(D[j][i],Azz,Ay,Az); }
  match_circles(D,Ds); // shift the arrays of the matched circles (6-11) by PI
}
/********************************************************************************/
// rotates whole dodecahedron around z y and z axix by Azz, Ay, Az angles
void dodecahedron::RyRzD(direction **D, double Ay, double Az) {
  for (j=0;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i] = RyRz(D[j][i],Ay,Az); }
  match_circles(D,Ds); // shift the arrays of the matched circles (6-11) by PI
}
/********************************************************************************/
// rotates whole dodecahedron around z y and z axix by Azz, Ay, Az angles
void dodecahedron::RzRyD(direction **D, double Az, double Ay) {
  for (j=0;j<6;j++)  
    for (i=0;i<circ_pix_num;i++) { D[j][i] = RzRy(D[j][i],Az,Ay); }
  match_circles(D,Ds); // shift the arrays of the matched circles (6-11) by PI
}
