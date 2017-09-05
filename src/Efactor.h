/* standard library for CPEDS routines */

//#include "cpeds-defs.h"

//**************************************************************************************
//  E(z) factor

double EE(double Wr0, double Wm0, double Wl0, double z, double w) {
  double tmp;
  
  tmp = Wr0*pow((1+z),4)+Wm0*pow((1+z),3)+Wl0*pow((1+z),3*(w+1))+(1-Wr0-Wm0-Wl0)*pow((1+z),2);
//  if (tmp >= 0 ) { return sqrt(tmp);  } 
//  else { return -sqrt(-tmp); }
//  return tmp;
  return sqrt(tmp);
}



//**************************************************************************************
