#include <math.h>
#include <time.h>
#include "cpeds-consts.h"

long cpeds_seed_offset=0;

const double PI=3.1415926535897932;
//const double PIsnd=1.57079632679490; // PI/2w
const double PIsnd=PI/2.0; // PI/2w
//const double twoPI=6.28318530717959; // 2*PI
const double twoPI=2.0*PI; // 2*PI
const double fourPI=12.56637061435917; // 4*PI
const double twooverPI=0.636619772367581; // 2/PI
const double sqrtPI=1.77245385090552;
//const double PI180=0.0174532925199433; // PI/180
const double PI180=PI/180.0; // PI/180
//const double PI180inv=57.295779513082321; // 180/PI
const double PI180inv=180.0/PI; // 180/PI
const double CPEDS_sr2deg2=3282.806350011744; // 180.0/PI*180.0/PI;
const double ANGSECRAD=4.84813681109535993e-6;
const double CPEDS_MPC=3.0856775814914E+22;	// definicja megaparseka [m] // AU=149 597 870 700 m
const double CPEDS_G=6.67384E-11;	// definicja stalej grawitacyjnej [N*m^2*kg^-2] or [m^3/kg/s^2]
const double CPEDS_Gcgs=6.67384E-8;//=6.67384E-8;	// definicja stalej grawitacyjnej [cm^3/g/s^2]
const double CPEDS_h=6.62606876E-34;  // stala Plancka [J*s]
const double CPEDS_kB=1.3806503E-23; // sta�a Boltzmanna [J/K]
const double CPEDS_epsilon0=8.854187817E-12; // sta�a przenikalnosci elektrycznej pr�ni
const double CPEDS_c=299792458;	// definicja predkosci swiatla [m/s]
const double CPEDS_e= 1.602176462E-19;	// ladunek elektronu [C]
const double CPEDS_m_ee=510998.902;   // masa elektronu w [ev]
const double CPEDS_m_p=1.6726e-27;    //masa protonu w [kg]
const double CPEDS_0K=-273.15;    // absolute zero temperature [C]

const double CPEDS_m_e= CPEDS_m_ee*CPEDS_e/CPEDS_c/CPEDS_c; // masa elektronu w [kg]
const double CPEDS_r_e= CPEDS_e/(4*PI*CPEDS_epsilon0*CPEDS_m_ee); // klasyczny promie� elektronu [m]
//const double sigma_T= 8*PI/3*pow(r_e,2); // przekr�j czynny na rozpraszanie Tomsona
const double CPEDS_sigma_T= 8*PI/3*pow(CPEDS_r_e,2); // przekr�j czynny na rozpraszanie Tomsona
const double CPEDS_R=8.3144621;    // gas constant [J/moll/K]
const double GYR= 1.0E+9*365*24*3600;	//definicja miliarda lat w sekundach
const double calc_th= 1.0E-100;       // granica dokladnosci obliczen numerycznych
const double Ry= 10973731.56854; // stala rydberga [m^-1] dla nieskonczonej masy jadra
//const double RyH= m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
const double RyH= CPEDS_m_e*CPEDS_m_p/(CPEDS_m_e + CPEDS_m_p) * pow(CPEDS_e,4)/(8*CPEDS_c*pow(CPEDS_epsilon0,2)*pow(CPEDS_h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
const double zetaR3= 1.202; // funkcja zeta Riemanna trzeciego rzedu (zeta(3)) zeta(n) = sum_x=1^inf 1/(x^n)
const double gb= 1.0;//2.7; // 2.7 to wspolczynnik dla bozonow we wzorze na energie srednie <E> = rho(p)/n(p)


const double CPEDS_Jy=1e-26; //[J/m^2]=[J/m^2/Hz/s] // Jansky definition








// COSMOLOGICAL STUFF FOR cpeds_cosmo

const double  H0= 71.0;  			// definicja stalek hubbla [km/s/Mpc]
const double  S0= 1.0;			// definicja czynnika skali
const double  T_CBR0= 2.726;			// definicja temperatury CBR [K]
const double  T_CMB0_L1_WMAP= 3.346;      // definicja dipola kosmicznego
const double  T_CMB0_L1_WMAP_l= 263.85;  // definicja dipola kosmicznego - dlugosc galaktyczna
const double  T_CMB0_L1_WMAP_b= 48.25;  // definicja dipola kosmicznego - szerokosc galaktyczna
const double  rho_r0= 0.26038;			// definicja gestosci promieniowania dzisiaj [eV*cm^-3]
const double  w0= -1.0;		// definicja wspolczynnika rownania stanu prozni
//double  rho_c0  3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
const double  rho_c0= 3*pow(H0*1000/CPEDS_MPC,2)/(8*PI*CPEDS_G); // dzisiejsza gestosc krytyczna
//double Wr0=0.0000422*pow(H0/100,2); // ten wzor jest dziwny - pochodzi z peeblsa
//double Wr0 = rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
const double Wr0= rho_r0/rho_c0*CPEDS_e*1.0E+6/CPEDS_c/CPEDS_c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
const double  Wdm0=0.27;
const double  Wb0=0.044;
const double  Wl0=0.74;

//double  Wm0=Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
const double  Wm0= Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
//double  Wtot0=Wr0+Wb0+Wdm0+Wl0;
const double  Wtot0= Wr0+Wb0+Wdm0+Wl0;
//double  Rc = c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
const double  Rc= CPEDS_c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
//double  Wk0 = 1.0-Wtot0;
const double  Wk0= 1.0-Wtot0;

const double  Nniu= 3.0;
//double  eta = 2.68E-8*Wb0*pow(H0,2) * 1.0E-4; eta  6.1E-10;
const double eta= 6.5E-10;
//double  B  Ry*h*c; // energia wi�zania atomu wodoru [eV]
const double  B= Ry*CPEDS_h*CPEDS_c; // energia wi�zania atomu wodoru [eV]
//double  n_r0 = 421.8 * pow((T_CBR0/2.75),3) * 1.0E+6;
const double n_r0=410.4 * 1.0E+6;// g�sto�� liczbowa foton�w [m^-3]
const double  Xe_threshold= 0.1; // zakladany prog jonizacji - stosunek gestosci liczbowych wolnych elektronow do wszystkich barionow

const double CPEDS_hydrogen_mass_fraction= 0.76;  /* mass fraction of hydrogen */
const double CPEDS_helium_mass_fraction= 0.24;  /* mass fraction of hydrogen */


  // sigma_rec jest sta�� cz�ci� przekroju czynnego na rekombinacje usrednionego po makswellowskim rozkladzie
  // pr�dko�ci elektron�w - dlatego w og�lno�ci zalezy on od temperatury i dlatego trzeba go ostatecznie liczyc w p�tli

//double  a0  pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
const double  a0= pow(CPEDS_h/(2*PI),2)/CPEDS_m_ee/CPEDS_e*CPEDS_c*CPEDS_c/exp(2.0); // h^2/(m_ee*exp(2))
  //  sigma_rec  (pow(2,10)*pow(PI,1.5))/3/exp(4)*pow(e*e/h*2*PI/c,3)*a0*B*e/(h/2/PI);
const double  sigma_rec_padi= 1.4E-13*1.0E-6/CPEDS_c; // wyprowadzenie tego wzoru nie jest takie proste [m^2]
const double  sigma_bf= 7.9e-22; // [m^2]

//double  r_gal0 = 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
const double  r_gal0= 10*(100/H0)*1.0E-3*CPEDS_MPC; // rozmiary galaktyki
//double  n_gal0 = 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
const double  n_gal0= 0.02 * pow(H0/100,3) * pow(CPEDS_MPC,-3); // przestrzenna gestosc liczbowa galaktyk

const double CPEDS_ZREC=1088.0; // redshift of the decoupling (doesn not depend much on the cosmological parameters): peack of the visibility function: taken from astro-ph/0302209

// integration limits defs. stuff
double CPEDS_ZINF=1e5; // infinite redshift
const double CPEDS_ZMIN=1e-5; // minimal redshift for past
double CPEDS_NP= 2000.0; 	// number of points probing an order of z range for redshifts > 0


// usual astronomical definitions
const double CPEDS_SOLAR_MASS=1.98892e30; // solar mass [kg]
const double CPEDS_AU=1.4959787069098932e+11; // astronomical unit [m]
const double CPEDS_SEC_IN_YR=31557600; // seconds in year

const double CPEDS_ECLIPTIC_NPOLE_GAL_L=96.38397;
const double CPEDS_ECLIPTIC_NPOLE_GAL_B=29.81144;
const double CPEDS_EQUATOR_NPOLE_GAL_L=122.93192;
const double CPEDS_EQUATOR_NPOLE_GAL_B=27.12825;

const double CPEDS_ECLIPTIC_NPOLE_EQ_RA=270.0; // right ascension of the ecliptic north pole [deg] (barycentric, for epoch 2000)
const double CPEDS_ECLIPTIC_NPOLE_EQ_DEC=66.56071; // declination of the ecliptic north pole [deg] (barycentric, for epoch 2000)
const double CPEDS_ECLIPTIC_SPOLE_EQ_RA=90.0; // right ascension of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
const double CPEDS_ECLIPTIC_SPOLE_EQ_DEC=-66.56071; // declination of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
const double CPEDS_TAI_UTC=37.0;  // This is the total count of leap seconds. 37 since Jan 1st 2017. The actual value can be found in http://maia.usno.navy.mil/ser7/ser7.dat
const double CPEDS_TT_TAI_OFFSET=32.184; // definition of the time difference between terrestial time and the international atomic time [s]
const double CPEDS_TT_TAI=32.184; // 32.184s = time difference between terrestial time and the international atomic time [s]
// 
// MATHEMATICAL STUFF
//


const int cpeds_bicubic_interpolation_ccoef_wt_d[16][16]= {
		{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
		{-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0},
		{2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1},
		{0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1},
		{-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
		{9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2},
		{-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2},
		{2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
		{-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1},
		{4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1}
}; 

//
// GADGET-2 STUFF
//
const double GADGET2_UnitMass_in_g= 1e13*CPEDS_SOLAR_MASS;         /*  code mass unit in g/h = 10^10 Msol_kg / h */
const double GADGET2_UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h = kpc, the actual unit in the code is kpc/h */
const double GADGET2_UnitVelocity_in_cm_per_s= 1.0e5; // = km /s
const double GADGET2_UnitTime_in_s= GADGET2_UnitLength_in_cm / GADGET2_UnitVelocity_in_cm_per_s;
const double GADGET2_UnitEnergy_in_cgs= GADGET2_UnitMass_in_g * pow(GADGET2_UnitLength_in_cm,2) / pow(GADGET2_UnitTime_in_s,2); // J/h
const double GADGET2_UnitDensity_in_cgs= GADGET2_UnitMass_in_g/ pow(GADGET2_UnitLength_in_cm,3);
