/* standard library for CPEDS routines
sigma_p - comoving distance to the object at given redshift z
in given cosmological model described by omegas and eq. od state with 
factor w0 and H0.

unit: [Mpc] 

H0 in km/s/Mpc
*/
#include "Efactor.h"
#include "cpeds.h"


double o_p(double z, double Wr0, double Wm0, double Wl0, double w0, double H0) {
double sigma_p;
  ZMINP = z;
  z=1.0/NF; z = skok(z);
  sigma_p = 0.0;

  do {
    z = z+ skok(z);
    
    sigma_p = sigma_p + skok(z)/EE(Wr0,Wm0,Wl0,z,w0);
    
//  cout << model[22][i] << "  " << z << "\n";
//    if ((1 - model[22][i]) >= 2/NP ) { Tlss = model[16][i]; tlss = model[2][i]; zlss = model[0][i];}
//    if ((1 - model[23][i]) >= 2/NP ) { Tgal = model[16][i]; tgal = model[2][i]; zgal = model[0][i];}
//    pisz_raport(z,i,2);
  } while (z+0.1*skok(z) <= ZMINP );    
  sigma_p =  c/(H0*1000.0)*sigma_p;
  return sigma_p;

}

