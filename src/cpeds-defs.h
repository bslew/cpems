#include <math.h>

double ZMINP = 0.001;// minimalne z do ktorego ma byc calkowana przeszlosc
double ZMINF = 0.0;	// wartosc z do ktorej ma byc calkowana przyszlosc
double ZMAX = 1.0E+4;	// wartosc z od ktorego ma byc calkowana przeszlosc

#define ZINF 1.0E+15// wartosc z odpowiadajaca z nieskonczonemu od ktorej zaczyna sie calkowanie
#define NP 2000	// ilosc punktow przypadajacych na rzad wielkosci w przeszlosci
#define NF 1000	// ilosc punktow przypadajacych na odcinek 0..-1 w przyszlosci
#define M 35		// liczba parametrow - do definiowania wielkosci tablicy
#define MM 17000	// liczba cel w talicy - powinna byc wieksza od NP*log(ZMAX/ZMINP)+NF i 
#define NAZWA_KATALOGU "./"		// nazwa katalogu wyjsciowego
#define NAZWA_PLIKU "model1"		// nazwa pliku wyjsciowego


#define PI 3.141592654
double DEG2RAD = PI/180;
double RAD2DEG = 180/PI;
#define MPC 3.08566758074E+22	// definicja megaparseka [m]
#define G 6.6726E-11	// definicja stalej grawitacyjnej [N*m^2*kg^-2]
#define h 6.62606876E-34  // stala Plancka [J*s]
#define kB 1.3806503E-23 // sta�a Boltzmanna [J/K]
#define epsilon0 8.854187817E-12 // sta�a przenikalnosci elektrycznej pr�ni
#define c 299792458	// definicja predkosci swiatla [m/s]
#define e  1.602176462E-19	// ladunek elektronu [C]
#define m_ee 510998.902   // masa elektronu w [ev]
#define m_p 1.6726e-27    //masa protonu w [kg]

double m_e = m_ee*e/c/c; // masa elektronu w [kg]
double r_e = e/(4*PI*epsilon0*m_ee); // klasyczny promie� elektronu [m]
//double sigma_T = 8*PI/3*pow(r_e,2); // przekr�j czynny na rozpraszanie Tomsona
#define sigma_T = 8*PI/3*pow(r_e,2); // przekr�j czynny na rozpraszanie Tomsona
double GYR = 1.0E+9*365*24*3600;	//definicja miliarda lat w sekundach 
double calc_th = 1.0E-100;       // granica dokladnosci obliczen numerycznych
double Ry = 10973731.56854; // stala rydberga [m^-1] dla nieskonczonej masy jadra
//double RyH = m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
#define RyH = m_e*m_p/(m_e + m_p) * pow(e,4)/(8*c*pow(epsilon0,2)*pow(h,3)); // stala Rydberga ze poprawka na skonczona mase jadra
double zetaR3 = 1.202; // funkcja zeta Riemanna trzeciego rzedu (zeta(3)) zeta(n) = sum_x=1^inf 1/(x^n)
double gb = 1.0;//2.7; // 2.7 to wspolczynnik dla bozonow we wzorze na energie srednie <E> = rho(p)/n(p)

// DEFINICJA MODELU STANDARDOWEGO - DOMYSLNEGO

double  H0 = 71.0;  			// definicja stalek hubbla [km/s/Mpc]
double  S0 = 1.0;			// definicja czynnika skali
double  T_CBR0 = 2.7255;			// definicja temperatury CBR [K]
double  T_CMB0_L1_WMAP = 3.346;      // definicja dipola kosmicznego
double  T_CMB0_L1_WMAP_l = 263.85;  // definicja dipola kosmicznego - dlugosc galaktyczna
double  T_CMB0_L1_WMAP_b = 48.25;  // definicja dipola kosmicznego - szerokosc galaktyczna
double  rho_r0 = 0.26038;			// definicja gestosci promieniowania dzisiaj [eV*cm^-3]
double  w0 = -1.0;		// definicja wspolczynnika rownania stanu prozni
//double  rho_c0 = 3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
#define  rho_c0 = 3*pow(H0*1000/MPC,2)/(8*PI*G); // dzisiejsza gestosc krytyczna
//double Wr0=0.0000422*pow(H0/100,2); // ten wzor jest dziwny - pochodzi z peeblsa
//double Wr0 = rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
#define Wr0 = rho_r0/rho_c0*e*1.0E+6/c/c;	// definicja danego modelu jesli nie bedzie podany z komendy lini
double  Wdm0=.27; 
double  Wb0=0.044;
double  Wl0=.74;

//double  Wm0=Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
#define  Wm0 = Wb0+Wdm0; // dop�ki sie nie wprowadzie roznego rownanie stanu dla barion�w i DM to mozna je termodynamicznie traktowac tak samo i l�cznie
//double  Wtot0=Wr0+Wb0+Wdm0+Wl0;  
#define  Wtot0 = Wr0+Wb0+Wdm0+Wl0;  
//double  Rc = c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
#define  Rc = c/(H0*1e6*sqrt(fabs(1.0-Wtot0))); // promien krzywizny wszechswiata [Gpc]
//double  Wk0 = 1.0-Wtot0;
#define  Wk0 = 1.0-Wtot0;
  
double  Nniu = 3.0;
//double  eta = 2.68E-8*Wb0*pow(H0,2) * 1.0E-4; eta = 6.1E-10;
double eta = 6.5E-10;
//double  B = Ry*h*c; // energia wi�zania atomu wodoru [eV]
#define  B = Ry*h*c; // energia wi�zania atomu wodoru [eV]
//double  n_r0 = 421.8 * pow((T_CBR0/2.75),3) * 1.0E+6; 
double n_r0= 410.4 * 1.0E+6;// g�sto�� liczbowa foton�w [m^-3]
double  Xe_threshold = 0.1; // zakladany prog jonizacji - stosunek gestosci liczbowych wolnych elektronow do wszystkich barionow 
  
  // sigma_rec jest sta�� cz�ci� przekroju czynnego na rekombinacje usrednionego po makswellowskim rozkladzie 
  // pr�dko�ci elektron�w - dlatego w og�lno�ci zalezy on od temperatury i dlatego trzeba go ostatecznie liczyc w p�tli
  
//double  a0 = pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
#define  a0 = pow(h/(2*PI),2)/m_ee/e*c*c/exp(2); // h^2/(m_ee*exp(2))
  //  sigma_rec = (pow(2,10)*pow(PI,1.5))/3/exp(4)*pow(e*e/h*2*PI/c,3)*a0*B*e/(h/2/PI);
double  sigma_rec_padi = 1.4E-13*1.0E-6/c; // wyprowadzenie tego wzoru nie jest takie proste [m^2]
double  sigma_bf = 7.9e-22; // [m^2]

//double  r_gal0 = 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
#define  r_gal0 = 10*(100/H0)*1.0E-3*MPC; // rozmiary galaktyki
//double  n_gal0 = 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
#define  n_gal0 = 0.02 * pow(H0/100,3) * pow(MPC,-3); // przestrzenna gestosc liczbowa galaktyk
 
