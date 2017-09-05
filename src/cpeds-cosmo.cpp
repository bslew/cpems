#include "cpeds-cosmo.h"
#include "Mscs-function.h"


//**************************************************************************************//
double cpeds_cmbfMax(double T) {
  return 2.821439*CPEDS_kB*T/CPEDS_h;
}
//**************************************************************************************//
double cpeds_Wtot0(double Wr0, double Wm0, double Wl0) {
  return Wr0+Wm0+Wl0;
}
//**************************************************************************************//
double cpeds_Wtot(double Wr0, double Wm0, double Wl0, double z, double w) {
  return cpeds_Wr(Wr0,Wm0,Wl0,z,w)+cpeds_Wm(Wr0,Wm0,Wl0,z,w)+cpeds_Wl(Wr0,Wm0,Wl0,z,w);
}

//**************************************************************************************//
void cpeds_set_points_number_per_logz(long np) {
  CPEDS_NP=(double)np;
}

//**************************************************************************************//
double cpeds_deltaz(double z) {
  double dz=0;
  if (z > 0)   dz=((double)CPEDS_ZINF/pow(10.0,floor(log10((double)CPEDS_ZINF/z)) ))/(double)CPEDS_NP; 
  if (z <= 0)  dz=1.0/(double)CPEDS_NP;
  return dz;
}

//**************************************************************************************//
double cpeds_Efactor(double Wr0, double Wm0, double Wl0, double z, double w) {
  double EE;
  double opz=z+1.0;
  double opz2=opz*opz;
  double opz3=opz2*opz;
  double opz4=opz2*opz2;

  EE=sqrt(Wr0*opz4 + Wm0*opz3 + Wl0*pow(opz,3.0*(1.0+w)) + (1.0-cpeds_Wtot0(Wr0,Wm0,Wl0))*opz2);
  return EE;
}
// TO BE IMPLEMENTED THE E FACTOR FOR GENERAL QUINESSENCE MODELS WITH w = w(z)
// E(z) = sqrt ( ... + Wde*exp(3 int_0^z (1+w(z))/(1+z) dz ) + .... )

//**************************************************************************************//
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w) {
  double zz,d,dz;
  
  zz=CPEDS_ZINF;  d=0.0;
  if (z==0.0) z=CPEDS_ZMIN;
  
  do {
    dz=cpeds_deltaz(zz);
    zz=zz-dz;
    d+=dz / cpeds_Efactor(Wr0,Wm0,Wl0,zz+0.5*dz,w);
  } while (zz-0.1*dz >= z);

  return d*CPEDS_c/1.0e8;  
}

/***************************************************************************************/
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w, double** Z, double** dp, long* N) {
  double zz,d,dz;
  long n=0;
  
  zz=CPEDS_ZINF;  d=0.0;
  if (z==0.0) z=CPEDS_ZMIN;
  do {
	  n++;
	  dz=cpeds_deltaz(zz);
	  zz=zz-dz;
  } while (zz-0.1*dz >= z);	  
  
  *Z = new double[n];
  *dp = new double[n];
  *N=n;
//  printf("n: %li",n);
  
  zz=CPEDS_ZINF;  d=0.0;
  if (z==0.0) z=CPEDS_ZMIN;
  
  long i=0;
  do {
    dz=cpeds_deltaz(zz);
    zz=zz-dz;
    d+=dz / cpeds_Efactor(Wr0,Wm0,Wl0,zz+0.5*dz,w);
    (*Z)[i]=zz+0.5*dz;
    (*dp)[i]=d;
//    printf("i: %li\n",i);
    i++;
  } while (zz-0.1*dz >= z);

  double coH100=CPEDS_c/1.0e8;
  for (long i = 0; i < n; i++) {
	(*dp)[i]*=coH100;
  }

  return d*CPEDS_c/1.0e8;  
}


//**************************************************************************************//
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w, double h) {
  return cpeds_dp(Wr0,Wm0,Wl0,z,w)/h;
}

//**************************************************************************************//
/* cpeds_queue<cpedsDirection*>* cpeds_age_of_universe(double Wr0, double Wm0,double Wl0,double z, double w) { */
double cpeds_age_of_universe(double Wr0, double Wm0,double Wl0,double z, double w) {
  double dz;
  mscsFunction t("age vs z");
  double zz,Z;

  if (z<0) { // integrate future
    zz=z;
    printf("integrating into future\n");
    do {
      dz=cpeds_deltaz(zz);
      t.newPoint(zz, 1.0 / ( (1+zz)*cpeds_Efactor(Wr0,Wm0,Wl0,zz,w) )); 
      zz+=dz;
    } while (zz <= -CPEDS_ZMIN);
  }

  zz=double(CPEDS_ZMIN);
  if (z>zz) zz=z;
  do {
    dz=cpeds_deltaz(zz);
    t.newPoint(zz,1.0 / ( (1+zz)*cpeds_Efactor(Wr0,Wm0,Wl0,zz,w) )); 
    zz+=dz;
  } while (zz <= CPEDS_ZINF);
  
  return t.integrate()/(100000.0/CPEDS_MPC) / GYR;
}

//**************************************************************************************

double cpeds_age_of_universe(double Wr0, double Wm0,double Wl0,double z, double w, double h) {
  return cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w)/h;
}


//**************************************************************************************
/* // redshift corresponding to given time in [Gyr] TO BE IMPLEMENTED */
/* double cpeds_redshift_at_time(double Wr0, double Wm0,double Wl0,double t, double w, double h) { */
/*   return cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w)/h; */
/* } */


//**************************************************************************************
double cpeds_lbtd(double Wr0, double Wm0,double Wl0,double z, double w) {
//	printf("%lf\n",cpeds_age_of_universe(Wr0,Wm0,Wl0,0.0,w));
//	printf("%lf\n",cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w));
  return CPEDS_c*GYR/CPEDS_MPC/1e3* (cpeds_age_of_universe(Wr0,Wm0,Wl0,0.0,w)-cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w));
}

//**************************************************************************************
double cpeds_lbtd(double Wr0, double Wm0,double Wl0,double z, double w, double h) {
  return cpeds_lbtd(Wr0,Wm0,Wl0,z,w)/h;
}

//**************************************************************************************
double cpeds_comoving_distance(double Wr0, double Wm0,double Wl0,double z, double w) {
  return cpeds_dp(Wr0,Wm0,Wl0,0.0,w) - cpeds_dp(Wr0,Wm0,Wl0,z,w);
/*   return cpeds_dp(Wr0,Wm0,Wl0,z,w); */
}
//**************************************************************************************
double cpeds_comoving_distance(double Wr0, double Wm0,double Wl0,double z, double w, double h) {
  return cpeds_comoving_distance(Wr0,Wm0,Wl0,z,w)/h;
}

//**************************************************************************************
double cpeds_comoving_distance_rec(double Wr0, double Wm0,double Wl0, double w, double f) {
  return cpeds_dp(Wr0,Wm0,Wl0,0.0,w) - f*cpeds_dp(Wr0,Wm0,Wl0,CPEDS_ZREC,w);
}

//**************************************************************************************
double cpeds_comoving_distance_rec(double Wr0, double Wm0,double Wl0, double w, double h, double f) {
  return cpeds_comoving_distance_rec(Wr0,Wm0,Wl0,w,f)/h;
}

//**************************************************************************************

double cpeds_acoustic_horizon_multipole(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h) {
  return cpeds_angular_diameter_distance(Wr0,Wm0,Wl0,z,w,h)/(cpeds_sound_horizon(Wr0,Wb0,Wm0,Wl0,z,w,h)/(1.0+z));
}

//**************************************************************************************
double cpeds_acoustic_horizon_angle(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h) {
  return PI180inv/cpeds_acoustic_horizon_multipole(Wr0,Wb0,Wm0,Wl0,z,w,h);
}

//**************************************************************************************
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0) {
  double Wk0 = cpeds_Wk0(Wr0,Wm0,Wl0);
  double Rc = CPEDS_c/1.0e8/sqrt(fabs(Wk0));
  if (Wk0 > 0) Rc=-Rc;
  return Rc;
}

//**************************************************************************************
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0, double h) {
  return cpeds_curvature_radius(Wr0,Wm0,Wl0)/h;
}
//**************************************************************************************
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0, double z, double w, double h) {
  return cpeds_curvature_radius(cpeds_Wr(Wr0,Wm0,Wl0,z,w),cpeds_Wm(Wr0,Wm0,Wl0,z,w),cpeds_Wl(Wr0,Wm0,Wl0,z,w))/h;
}
//**************************************************************************************
double cpeds_Wr(double Wr0, double Wm0, double Wl0, double z, double w) {
  double opz=1.0+z;
  double e2=cpeds_Efactor(Wr0,Wm0,Wl0,z,w);
  e2*=e2;
  return Wr0*opz*opz*opz*opz/e2;
}
//**************************************************************************************
double cpeds_Wm(double Wr0, double Wm0, double Wl0, double z, double w) {
  double opz=1.0+z;
  double e2=cpeds_Efactor(Wr0,Wm0,Wl0,z,w);
  e2*=e2;
  return Wm0*opz*opz*opz/e2;
}
//**************************************************************************************
double cpeds_Wb(double Wr0, double Wb0, double Wl0, double z, double w) {
  return cpeds_Wm(Wr0,Wb0,Wl0,z,w);
}
//**************************************************************************************
double cpeds_Wl(double Wr0, double Wm0, double Wl0, double z, double w) {
  double opz=1.0+z;
  double e2=cpeds_Efactor(Wr0,Wm0,Wl0,z,w);
  e2*=e2;

  return Wl0*pow(opz,3.0*(1.0+w))/e2;
}

//**************************************************************************************
double cpeds_Wk0(double Wr0, double Wm0,double Wl0) {
  return 1.0-Wr0-Wm0-Wl0;
}
//**************************************************************************************
double cpeds_Wk(double Wr0, double Wm0,double Wl0, double z, double w) {
  return cpeds_Wk0(cpeds_Wr(Wr0,Wm0,Wl0,z,w),cpeds_Wm(Wr0,Wm0,Wl0,z,w),cpeds_Wl(Wr0,Wm0,Wl0,z,w));
}

//**************************************************************************************
double cpeds_angular_diameter_distance(double Wr0, double Wm0,double Wl0, double z, double w, double h) {
  double dA,X;
  double Wk0=cpeds_Wk0(Wr0,Wm0,Wl0);
  double Rc;

//  X = cpeds_lbtd(Wr0,Wm0,Wl0,z,w,h);
  X = cpeds_comoving_distance(Wr0,Wm0,Wl0,z,w,h);

  if (Wk0 == 0.0) { // flat case
    dA = X; 
  }
  else {
    Rc=cpeds_curvature_radius(Wr0,Wm0,Wl0,h);

    if (Rc < 0.0) { // negative curvature
      dA = -Rc * sinh( -X/Rc ) ; }
    else { // positive curvature
      dA = Rc * sin( X/Rc ); }
  }

  dA/=(1.0+z);
  return dA;
}
/***************************************************************************************/
double cpeds_get_dA_from_ang(double Wr0, double Wb0, double Wm0, double  Wl0, double* Z, double w0, double h, double size) {
	double zz;
	double z;
//	mscsFunction t;
//	
//	zz=double(CPEDS_ZMIN);
//	if (z>zz) zz=z;
//	do {
//		dz=cpeds_deltaz(zz);
//		t.newPoint(zz,1.0 / ( (1+zz)*cpeds_Efactor(Wr0,Wm0,Wl0,zz,w) )); 
//		zz+=dz;
//	} while (zz <= CPEDS_ZINF);

	
}
//**************************************************************************************
double cpeds_sound_horizon(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h) {
  double rhoC=cpeds_rhoC(Wr0,Wm0,Wl0,z,w,h);
  return cpeds_dp(Wr0,Wm0,Wl0,z,w,h)*cpeds_sound_speed(cpeds_Wr(Wr0,Wm0,Wl0,z,w)*rhoC,cpeds_Wb(Wr0,Wb0,Wl0,z,w)*rhoC);
}

//**************************************************************************************
// sound horizon [Gpc]
// assumes that the sound speed is simply c/Sqrt(3) 

/* double cpeds_sound_horizon(double Wr0, double Wm0,double Wl0, double z, double w, double h) { */
/*   return cpeds_sound_horizon(Wr0,Wm0,Wl0,z,w)/h; */
/* } */

//**************************************************************************************
double cpeds_luminosity_distance(double Wr0, double Wm0,double Wl0, double z, double w, double h) {
  return (1.0+z)*(1.0+z)*cpeds_angular_diameter_distance(Wr0,Wm0,Wl0,z,w,h);
}

//**************************************************************************************
double cpeds_hubble(double Wr0, double Wm0,double Wl0, double z, double w, double h) {
  return h*100*cpeds_Efactor(Wr0,Wm0,Wl0,z,w);
}
//**************************************************************************************
double cpeds_rhoC(double Wr0, double Wm0,double Wl0, double z, double w) {
	double H2=1.0e5/CPEDS_MPC;
	
	if (z!=0) {
		H2*=cpeds_Efactor(Wr0,Wm0,Wl0,z,w);
		H2*=H2;
	}
	else {
		H2*=H2;
	}
	return 3*H2 / ( 8*PI*CPEDS_G);
}
//**************************************************************************************
double cpeds_rhoC(double Wr0, double Wm0,double Wl0, double z, double w, double h) {
  return cpeds_rhoC(Wr0,Wm0,Wl0,z,w)*h*h;
}
/***************************************************************************************/
double cpeds_rhoC0(double h) {
	double H2=1.0e5/CPEDS_MPC*h;
	H2*=H2;
	return 3*H2 / ( 8*PI*CPEDS_G);	
}

//**************************************************************************************
double cpeds_sound_speed(double rhor, double rhob) {
  return 1.0/sqrt(3.0*(1+3.0*rhob/(4.0*rhor)));
}



//**************************************************************************************
double cpeds_fX_e(double T, double eta) {
	double x,P;
	double acc=1e-10;
	
    P = 4.0*sqrt(2)*zetaR3/sqrt(PI) * eta * pow(CPEDS_kB*T/CPEDS_m_e/CPEDS_c/CPEDS_c,1.5) * exp(cpeds_BB(T)*CPEDS_e/CPEDS_kB/T); // / r_w(m_base(T),T*e/kB/gb); // byla poprawka na wielowzbudzeniowy gaz...
    if (P < 1.0e+200 )  {x = (sqrt(1.0+4.0*P)-1.0)/(2.0*P); } else {x = 0;}
    if (x < acc) {x = 0;}

    return x;
}


//**************************************************************************************
double cpeds_BB(double T) {
  return Ry*CPEDS_c*CPEDS_h/CPEDS_e/pow(cpeds_m_base(T),2);
}
/***************************************************************************************/
double cpeds_m_base (double T) {
  double tmp;

  if (gb*CPEDS_kB*T < Ry*CPEDS_h*CPEDS_c) { 
    tmp = floor(sqrt(Ry*CPEDS_h*CPEDS_c/(Ry*CPEDS_h*CPEDS_c-gb*CPEDS_kB*T)));
    if (tmp > cpeds_neff(T)) { return cpeds_neff(T); } else { return tmp; }}
  else return cpeds_neff(T);
}
/***************************************************************************************/
double cpeds_neff(double T) {  return 1.0; }

/***************************************************************************************/
double cpeds_calculate_comptonY(double Wb0, double h, double mueinv) {
	double rhoC0=cpeds_rhoC0(h);
	return CPEDS_sigma_T* CPEDS_kB * rhoC0 * Wb0 / (mueinv * CPEDS_m_p * CPEDS_m_e * CPEDS_c*CPEDS_c) * CPEDS_MPC;
}
/***************************************************************************************/
double cpeds_calculate_meanMolecularWeight(double Xp) {
	return 4.0/(8.0 - 5.0*(1.0-Xp));
}
/***************************************************************************************/
double cpeds_convert_gadget_internalEnergy_to_temp(double u, double gamma, double Xp) {
	double T;
	//
	// convert internal energy into SI units
	// 
	u*= GADGET2_UnitEnergy_in_cgs/GADGET2_UnitMass_in_g * 1e-4; // conversion from (km/s)^2 to (m/s)^2
	//
	// convert internal energy into temperature of the gas in SI units
	// 
	double MeanWeight=cpeds_calculate_meanMolecularWeight(Xp)*CPEDS_m_p;
	 
	T= MeanWeight / CPEDS_kB * (gamma-1) * u;
	return T;
}

//**************************************************************************************
//! deceleration parameter q: 
/*! \f$ q = -\frac{\ddot a a}{\dot a^2} \f$*/
/*!!!!!!! UNDER DEVELOPMENT !!!!!!!!1 */
/* double cpeds_deceleration_parameter_q(double _Wr0, double _Wm0, double _Wl0, double _z, double _w0) { */

/*   a0/a = 1+z; */
/*   da/a = H dt; */
/*   ln a = int H dt; */
  
/* } */


//**************************************************************************************
//! scale factor a
/*!!!!!!! UNDER DEVELOPMENT !!!!!!!!1 */
/*! \f$ q = -\frac{\ddot a a}{\dot a^2} \f$*/

/* double cpeds_scale_factor_a(double _Wr0, double _Wm0, double _Wl0, double _z, double _w0) { */
/* } */
