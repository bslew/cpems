#include <stdlib.h>
#include <math.h>
#include "cpeds-consts.h"
// the spherical coordinate system is with theta = 0 at the north pole

double CPEDS_sph2cart(int coord, double th, double phi) { // angles passed in radians
  // coord - 0 - x, 1 - y, 2 - z
  double x,y,z;

  //  ath = dr*ath; aphi=dr*aphi; o = dr*o;

  if (coord == 0) { // x
    x = sin(th)*cos(phi);
    return x;
  }

  if (coord == 1) { // y
    y = sin(th)*sin(phi);
    return y;
  }
   
  if (coord == 2) { // z
    z = cos(th);
    return z;
  }
  return -1;
}

double CPEDS_cart2sph(int coord, double x, double y,double z) {
  // coord - 0 - theta, 1 - phi
  double th=0,phi=0;
  double dr,rd;//,PI;
  //PI = 3.141592654;
  dr=PI/180; rd=180/PI;
  //  ath = dr*ath; aphi=dr*aphi; o = dr*o;

  //printf("%f %f %f\n",x,y,z);
  if (coord == 0) {
    th = rd*acos(z/sqrt(x*x+y*y+z*z)); return th;}

  if (coord == 1) {
    if (y > 0) {
      if (x > 0) {phi = rd*atan(y/x); }
      if (x < 0) {phi = rd*atan(y/x) + 180; } 
      if (x == 0) { phi = 90; }}
    if (y < 0) {
      if (x > 0) {phi = rd*atan(y/x) + 360; }
      if (x < 0) {phi = rd*atan(y/x) + 180; } 
      if (x == 0) { phi = 270; }}
    if (y == 0) {
      if (x > 0) {phi = 0; }
      if (x < 0) {phi = 180; }}
    return phi; }
  return -1;
}

//****************************************************************************************************************
//calculates an angle between the two direction on the sky in radians
double CPEDS_ang_n1n2(double th1, double phi1,double th2,double phi2) {
  double x1,y1,z1,x2,y2,z2;
  double ang;

  x1 = CPEDS_sph2cart(0,th1,phi1);   y1 = CPEDS_sph2cart(1,th1,phi1);   z1 = CPEDS_sph2cart(2,th1,phi1);
  x2 = CPEDS_sph2cart(0,th2,phi2);   y2 = CPEDS_sph2cart(1,th2,phi2);   z2 = CPEDS_sph2cart(2,th2,phi2);

  ang = acos( (x1*x2 + y1*y2 + z1*z2)/( sqrt((x1*x1 + y1*y1 + z1*z1) * (x2*x2 + y2*y2 + z2*z2)) ));

  return ang;
}
