#ifndef CPEDS_CONSTS
#define CPEDS_CONSTS

//#define const double PI=3.14159265358979;
extern long cpeds_seed_offset;

extern const double PI;//=3.14159265358979;
extern const double PIsnd;//=1.57079632679490; // PI/2w
extern const double twoPI;//=6.28318530717959; // 2*PI
extern const double fourPI;//=6.28318530717959; // 2*PI
extern const double twooverPI;//=6.28318530717959; // 2*PI
extern const double sqrtPI;
extern const double PI180;
extern const double PI180inv;
extern const double ANGSECRAD;
extern const double CPEDS_sr2deg2;
extern const double CPEDS_MPC;//=3.08566758074E+22;	// definicja megaparseka [m]
extern const double CPEDS_G;//=6.6726E-11;	// definicja stalej grawitacyjnej [N*m^2*kg^-2] or [m^3/kg/s^2]
extern const double CPEDS_Gcgs;//=6.67384E-8;	// definicja stalej grawitacyjnej [cm^3/g/s^2]
extern const double CPEDS_h;//=6.62606876E-34;  // stala Plancka [J*s]
extern const double CPEDS_kB;//=1.3806503E-23; // sta�a Boltzmanna [J/K]
extern const double CPEDS_epsilon0;//=8.854187817E-12; // sta�a przenikalnosci elektrycznej pr�ni
extern const double CPEDS_c;//=299792458;	// definicja predkosci swiatla [m/s]
extern const double CPEDS_e;//= 1.602176462E-19;	// ladunek elektronu [C]
extern const double CPEDS_m_ee;//=510998.902;   // masa elektronu w [ev]
extern const double CPEDS_m_p;//=1.6726e-27;    //masa protonu w [kg]
extern const double CPEDS_0K;    // absolute zero temperature [C]

extern const double CPEDS_m_e;//= m_ee*e/c/c; // masa elektronu w [kg]
extern const double CPEDS_r_e;//= e/(4*PI*epsilon0*m_ee); // klasyczny promie� elektronu [m]
//extern const double sigma_T;//= 8*PI/3*pow(r_e,2); // przekr�j czynny na rozpraszanie Tomsona
extern const double CPEDS_sigma_T;//= 8*PI/3*pow(r_e,2); // przekr�j czynny na rozpraszanie Tomsona
extern const double CPEDS_R;    // gas constant [J/moll/K]
extern const double GYR;//= 1.0E+9*365*24*3600;	//definicja miliarda lat w sekundach 
extern const double calc_th;//= 1.0E-100;       // granica dokladnosci obliczen numerycznych
extern const double Ry;//= 10973731.56854; // stala rydberga [m^-1] dla nieskonczonej masy jadra
//extern const double RyH;//= m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
extern const double RyH;//= m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
extern const double zetaR3;//= 1.202; // funkcja zeta Riemanna trzeciego rzedu (zeta(3)) zeta(n) ;//= sum_x;//=1^inf 1/(x^n)
extern const double gb;//= 1.0;//2.7; // 2.7 to wspolczynnik dla bozonow we wzorze na energie srednie <E> ;//= rho(p)/n(p)

extern const double CPEDS_Jy;//=1e-26 [J/m^2] // Jansky definition

// DEFINICJA MODELU STANDARDOWEGO - DOMYSLNEGO

extern const double  H0;//= 71.0;  			// definicja stalek hubbla [km/s/Mpc]
extern const double  S0;//= 1.0;			// definicja czynnika skali
extern const double  T_CBR0;//= 2.725;			// definicja temperatury CBR [K]
extern const double  T_CMB0_L1_WMAP;//= 3.346;      // definicja dipola kosmicznego
extern const double  T_CMB0_L1_WMAP_l;//= 263.85;  // definicja dipola kosmicznego - dlugosc galaktyczna
extern const double  T_CMB0_L1_WMAP_b;//= 48.25;  // definicja dipola kosmicznego - szerokosc galaktyczna
extern const double  rho_r0;//= 0.26038;			// definicja gestosci promieniowania dzisiaj [eV*cm^-3]
extern const double  w0;//= -1.0;		// definicja wspolczynnika rownania stanu prozni
//double  rho_c0  3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
extern const double  rho_c0;//= 3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
//double Wr0;//=0.0000422*pow(H0/100,2); // ten wzor jest dziwny - pochodzi z peeblsa
//double Wr0 ;//= rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
extern const double Wr0;//= rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
extern const double  Wdm0;//=0.27; 
extern const double  Wb0;//=0.044;
extern const double  Wl0;//=0.74;

//double  Wm0;//=Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
extern const double  Wm0;//= Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
//double  Wtot0;//=Wr0+Wb0+Wdm0+Wl0;  
extern const double  Wtot0;//= Wr0+Wb0+Wdm0+Wl0;  
//double  Rc ;//= c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
extern const double  Rc;//= c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
//double  Wk0 ;//= 1.0-Wtot0;
extern const double  Wk0;//= 1.0-Wtot0;
  
extern const double  Nniu;//= 3.0;
//double  eta ;//= 2.68E-8*Wb0*pow(H0,2) * 1.0E-4; eta  6.1E-10;
extern const double eta;//= 6.5E-10;
//double  B  Ry*h*c; // energia wi�zania atomu wodoru [eV]
extern const double  B;//= Ry*h*c; // energia wi�zania atomu wodoru [eV]
//double  n_r0 ;//= 421.8 * pow((T_CBR0/2.75),3) * 1.0E+6; 
extern const double n_r0;//=410.4 * 1.0E+6;// g�sto�� liczbowa foton�w [m^-3]
extern const double  Xe_threshold;//= 0.1; // zakladany prog jonizacji - stosunek gestosci liczbowych wolnych elektronow do wszystkich barionow 

extern const double CPEDS_hydrogen_mass_fraction;  /* mass fraction of hydrogen */
extern const double CPEDS_helium_mass_fraction;  /* mass fraction of hydrogen */


  // sigma_rec jest sta�� cz�ci� przekroju czynnego na rekombinacje usrednionego po makswellowskim rozkladzie 
  // pr�dko�ci elektron�w - dlatego w og�lno�ci zalezy on od temperatury i dlatego trzeba go ostatecznie liczyc w p�tli
  
//double  a0  pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
extern const double  a0;//= pow(h/(2*PI),2)/m_ee/e*c*c/exp(2.0); // h^2/(m_ee*exp(2))
  //  sigma_rec  (pow(2,10)*pow(PI,1.5))/3/exp(4)*pow(e*e/h*2*PI/c,3)*a0*B*e/(h/2/PI);
extern const double  sigma_rec_padi;//= 1.4E-13*1.0E-6/c; // wyprowadzenie tego wzoru nie jest takie proste [m^2]
extern const double  sigma_bf;//= 7.9e-22; // [m^2]

//double  r_gal0 ;//= 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
extern const double  r_gal0;//= 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
//double  n_gal0 ;//= 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
extern const double  n_gal0;//= 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
 
extern const double CPEDS_ZREC; // redshift of the recombination (doesn not depend much on the cosmological parameters)

// integration limits defs. stuff 

extern double CPEDS_ZINF; // infinite redshift
extern const double CPEDS_ZMIN; // infinite redshift
extern  double CPEDS_NP; 	// number of points probing an order of z range for redshifts > 0


// usual astronomical definitions
extern const double CPEDS_SOLAR_MASS; // solar mass [kg]
extern const double CPEDS_AU; // astronomical unit [m]
extern const double CPEDS_SEC_IN_YR; // number of seconds in year = 365.25 * 24 * 3600

extern const double CPEDS_ECLIPTIC_NPOLE_GAL_L; // galactic longitude of the ecliptic north pole [deg]
extern const double CPEDS_ECLIPTIC_NPOLE_GAL_B; // galactic latitude of the ecliptic north pole [deg]
extern const double CPEDS_EQUATOR_NPOLE_GAL_L; // galactic longitude of the earth equator north pole [deg]
extern const double CPEDS_EQUATOR_NPOLE_GAL_B; // galactic latitude of the earth equator north pole [deg]

extern const double CPEDS_ECLIPTIC_NPOLE_EQ_RA; // right ascension of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
extern const double CPEDS_ECLIPTIC_NPOLE_EQ_DEC; // declination of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
extern const double CPEDS_ECLIPTIC_SPOLE_EQ_RA; // right ascension of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
extern const double CPEDS_ECLIPTIC_SPOLE_EQ_DEC; // declination of the ecliptic north pole [deg] (barycentric, for epoch 2000) [deg]
extern const double CPEDS_TAI_UTC; // DeltaAT: the actual time difference between UTC and international atomic clock time as for 1 July 2012 - total count of leap seconds in UTC,  [s]
extern const double CPEDS_TT_TAI_OFFSET; // the time difference between terrestial time and the international atomic time [s]
extern const double CPEDS_TT_TAI; // 32.184s = time difference between terrestial time and the international atomic time [s]

#define CPEDS_JD2000 2451545.0
#define CPEDS_B1950 2433282.4235


// 
// MATHEMATICAL STUFF
//

extern const int cpeds_bicubic_interpolation_ccoef_wt_d[16][16];




//
// GADGET-2 STUFF
//

extern const double GADGET2_UnitMass_in_g;         /*  code mass unit in g/h */
extern const double GADGET2_UnitLength_in_cm;   /*  code length unit in cm/h */
extern const double GADGET2_UnitVelocity_in_cm_per_s;
extern const double GADGET2_UnitTime_in_s;
extern const double GADGET2_UnitEnergy_in_cgs;
extern const double GADGET2_UnitDensity_in_cgs;




#endif
