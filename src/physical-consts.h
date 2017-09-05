#define PI 3.141592654
#define MPC 3.08566758074E+22	// definicja megaparseka [m]
#define G 6.6726E-11	// definicja stalej grawitacyjnej [N*m^2*kg^-2]
#define h 6.62606876E-34  // stala Plancka [J*s]
#define kB 1.3806503E-23 // sta³a Boltzmanna [J/K]
#define epsilon0 8.854187817E-12 // sta³a przenikalnosci elektrycznej pró¿ni
#define c 299792458	// definicja predkosci swiatla [m/s]
#define e  1.602176462E-19	// ladunek elektronu [C]
#define m_ee 510998.902   // masa elektronu w [ev]
#define m_p 1.6726e-27    //masa protonu w [kg]

#define m_e  m_ee*e/c/c; // masa elektronu w [kg]
#define r_e  e/(4*PI*epsilon0*m_ee); // klasyczny promieñ elektronu [m]
#define sigma_T  8*PI/3*pow(r_e,2); // przekrój czynny na rozpraszanie Tomsona
#define sigma_T  8*PI/3*pow(r_e,2); // przekrój czynny na rozpraszanie Tomsona
#define GYR  1.0E+9*365*24*3600;	//definicja miliarda lat w sekundach 
#define calc_th  1.0E-100;       // granica dokladnosci obliczen numerycznych
#define Ry  10973731.56854; // stala rydberga [m^-1] dla nieskonczonej masy jadra
#define RyH  m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
#define RyH  m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
#define zetaR3  1.202; // funkcja zeta Riemanna trzeciego rzedu (zeta(3)) zeta(n) = sum_x=1^inf 1/(x^n)
#define gb  1.0;//2.7; // 2.7 to wspolczynnik dla bozonow we wzorze na energie srednie <E> = rho(p)/n(p)

// DEFINICJA MODELU STANDARDOWEGO - DOMYSLNEGO

#define  H0  71.0;  			// definicja stalek hubbla [km/s/Mpc]
#define  S0  1.0;			// definicja czynnika skali
#define  T_CBR0  2.725;			// definicja temperatury CBR [K]
#define  T_CMB0_L1_WMAP  3.346;      // definicja dipola kosmicznego
#define  T_CMB0_L1_WMAP_l  263.85;  // definicja dipola kosmicznego - dlugosc galaktyczna
#define  T_CMB0_L1_WMAP_b  48.25;  // definicja dipola kosmicznego - szerokosc galaktyczna
#define  rho_r0  0.26038;			// definicja gestosci promieniowania dzisiaj [eV*cm^-3]
#define  w0  -1.0;		// definicja wspolczynnika rownania stanu prozni
//double  rho_c0  3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
#define  rho_c0  3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
//double Wr0=0.0000422*pow(H0/100,2); // ten wzor jest dziwny - pochodzi z peeblsa
//double Wr0 = rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
#define Wr0  rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
#define  Wdm 0.27; 
#define  Wb0 0.044;
#define  Wl 0.74;

//double  Wm0=Wb0+Wdm0; // dopóki sie nie wprowadzie roznego rownanie stanu dla barionów i DM to mozna je termodynamicznie traktowac tak samo i l±cznie
#define  Wm0  Wb0+Wdm0; // dopóki sie nie wprowadzie roznego rownanie stanu dla barionów i DM to mozna je termodynamicznie traktowac tak samo i l±cznie
//double  Wtot0=Wr0+Wb0+Wdm0+Wl0;  
#define  Wtot0  Wr0+Wb0+Wdm0+Wl0;  
//double  Rc = c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
#define  Rc  c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
//double  Wk0 = 1.0-Wtot0;
#define  Wk0  1.0-Wtot0;
  
#define  Nniu  3.0;
//double  eta = 2.68E-8*Wb0*pow(H0,2) * 1.0E-4; eta  6.1E-10;
#define eta  6.5E-10;
//double  B  Ry*h*c; // energia wi±zania atomu wodoru [eV]
#define  B  Ry*h*c; // energia wi±zania atomu wodoru [eV]
//double  n_r0 = 421.8 * pow((T_CBR0/2.75),3) * 1.0E+6; 
#define n_r0 410.4 * 1.0E+6;// gêsto¶æ liczbowa fotonów [m^-3]
#define  Xe_threshold  0.1; // zakladany prog jonizacji - stosunek gestosci liczbowych wolnych elektronow do wszystkich barionow 
  
  // sigma_rec jest sta³± czê¶ci± przekroju czynnego na rekombinacje usrednionego po makswellowskim rozkladzie 
  // prêdko¶ci elektronów - dlatego w ogólno¶ci zalezy on od temperatury i dlatego trzeba go ostatecznie liczyc w pêtli
  
//double  a0  pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
#define  a0  pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
  //  sigma_rec  (pow(2,10)*pow(PI,1.5))/3/exp(4)*pow(e*e/h*2*PI/c,3)*a0*B*e/(h/2/PI);
#define  sigma_rec_padi  1.4E-13*1.0E-6/c; // wyprowadzenie tego wzoru nie jest takie proste [m^2]
#define  sigma_bf  7.9e-22; // [m^2]

//double  r_gal0 = 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
#define  r_gal0  10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
//double  n_gal0 = 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
#define  n_gal0  0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
 
