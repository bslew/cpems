//#include <math.h>
//#include "/home/blew/programy/CPEDS/cpeds-defs.h"

double CPEDS_cart2sph(int coord, double x, double y,double z) {
  // coord - 0 - theta, 1 - phi
  double th,phi,dr,rd,PI;
  PI = 3.141592654;
  dr=DEG2RAD; rd=RAD2DEG;
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

}
