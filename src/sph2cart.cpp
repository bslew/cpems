//#include <math.h>
//#include "/home/blew/programy/CPEDS/cpeds-defs.h"

// the spherical coordinate system is with theta = 0 at the north pole

double CPEDS_sph2cart(int coord, double th, double phi) { // angles passed in [deg]s
  // coord - 0 - x, 1 - y, 2 - z
  double x,y,z,dr,rd,PI;
  PI = 3.141592654;

  //  ath = dr*ath; aphi=dr*aphi; o = dr*o;

  if (coord == 0) { // x
    x = sin(PI/180*th)*cos(PI/180*phi);
    return x;
  }

  if (coord == 1) { // y
    y = sin(PI/180*th)*sin(PI/180*phi);
    return y;
  }
   
  if (coord == 2) { // z
    z = cos(PI/180*th);
    return z;
  }

}
