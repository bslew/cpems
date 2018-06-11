/*!
 \file cosmocalc program
 \brief cosmological parameters calculator
 \details 



 TODO: add calculation of flux for planets (Jupiter, Moon, Venus, Mars) for given time with given beam and at given frequency,
 assuming planet temperature, black body radiance, and classical Cassegrain beam size
 
 \date 2009/05/27 01:53:56 
 \author Bartosz Lew
 */

#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <iostream>
#include <tclap/CmdLine.h>
#include <QtCore/QStringList>
#include <QtCore/QString>
#include <QtCore/QFileInfo>
#include "cpeds-cosmo.h"
#include "cpeds-function-cosmo.h"
#include "Mscs-function3dregc.h"
#include "cpeds-direction.h"
#include "cpeds-direction_set.h"
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "cpeds-list.h"
#include  "cpeds-msgs.h"
#include "cpeds-point_set.h"
#include "gsl/gsl_sf_bessel.h"



#include <libnova/dynamical_time.h>

//#include "Mscs-alms.h"

#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


double _ns, _A, _k0, _s0scale;
string _tf, _refraction, _gravPot, _refrMoonAtDusk;
bool _sigma0, _sigmaR, _notfnorm, _BBKS_transferFunction, _TSZfreqFn,  _RTres;
double _Wr0, _a0, _Wcdm0, _Wm0, _Wb0, _Wl0, _w0, _h, _rmin, _rmax, _dr, _fmin,
		_fmax, _X, _Xrec, _Tcmb0, _angdA, _y0, _mueinv1, _TSZfreq, _gravSoft,
		_freq, _nu2l, _l2nu, _kSZbetall;
double _z, _size, _zmin, _zmax, _K2KeV, _taper;
bool _visibility, _vsa, _vsz, _vst, _thA, _lA, _Rc, _Wtot0, _Wtot, _Wk, _agea,
		_agez, _dp, _lbtd, _lbtd_rec, _dA, _dL, _dac_rec, _dac, _lac, _thac,
		_cs, _ang, _q, _a, _Tcmb, _cmbfMax, _dAvsz, _Xvsz, _lbtdvsz, _partMass;
bool _Wl, _Wr, _Wm, _Wb, _Wl0fit2flat, _rhoC, _boxesInField, _Wvsz, _lightCone,
		_printGadgetUnits, _a0H0inv, _WkRTopHat, _Efact;
bool _JD, _doJD2cal, _dou2cal;
double _JD2cal, _JDut1, _u2cal;
string _cal2JD;

int _verbosity;
bool _script_mode, _refrfn, _testCoor, _dontSaveRun, _testVtotRT32;
long _fN, _NP, _cells, _tfcol;
long _rsize;
double _Lmin, _TinkerDeltaRho;
string _outputFile, _lbdeg2dms, _inCoordFmt;

bool _radec2000toGal;
bool _radec2000toNow,_nutate, _aberrate;

void parseOptions(int argc, char** argv);
//-----------------------------------

void say(cpedsMsgs& msgs, string msg, string key, string value);
double caclulate_sigma0(double R, cpedsMsgs& msgs,
		mscsFunction *transfer_function = NULL);
void caclulate_sigmaR(cpedsMsgs msgs);
void doRefraction(cpedsMsgs& msgs);
void doMoonRefractionAtDusk(cpedsMsgs& msgs);
void doModel4();
void doTestVtotRT32();


/*!
 \brief calculates beam power pattern for the Cassagrain telescopes with parabolic illumination
 \details 
 @param d - primary mirror aperture [m]
 @param ds - secondary mirror aperture [m]
 @param lambda - wavelength [m]
 @param taper - illumination power at the edge of the aperture relative to the central value [db] (eg. -12 db)
 @param thetaMax [deg]
 @param dTheta [deg]
 
 @return

 \date Apr 25, 2014, 12:01:48 AM
 \author Bartosz Lew
 */
mscsFunction calculate_beam_power_pattern(double d, double ds, double lambda,
		double taper, double thetaMax, double dTheta);

void calculateGravitationalPotential(string infile);

/*!
 \brief calculate nonrelativistic frequency depencence of TSZ effect
 \details 
 @param f - frequenct [GHz]
 @return

This is the f(x) factor which when multiplied by comptonization parameter gives 
the relative change of thermodynamic temperature of CMB

 \date Oct 16, 2013, 10:39:06 PM
 \author Bartosz Lew
 */
double getSZfreqFact(double f);
double getSZintensityFact(double f);
/*!
 \brief calculate nonrelativistic frequency depencence of TSZ effect
 \details 
 @param fmin - [GHz]
 @param fmax - [GHz]
 @param n - number of points in the output function
 
 @return

 \date Oct 16, 2013, 10:39:59 PM
 \author Bartosz Lew
 */
mscsFunction calculate_TSZfreqFn(double fmin, double fmax, long n);
/*!
	\brief calculate nonrelativistic frequency depencence of TSZ effect intensity 
	\details 
    @param fmin - [GHz]
    @param fmax - [GHz]
    @param n - number of points in the output function
	@return DeltaI = deltaT*g(nu) * (x^4 e^x)/(e^x-1)^2 * 2 (kB Tcmb0)^3 / (hc)^2
	
	where deltaT = (T-Tcmb)/Tcmb 
	and DeltaI = I-I0 = I-B(nu,Tcmb)

	\date Jan 2, 2016, 2:59:48 PM
	\author Bartosz Lew
*/
mscsFunction calculate_TSZintensityFn(double fmin, double fmax, long n);

/*!
	\brief  claculate kSZ intensity dependence
	\details 
    @param fmin - [GHz]
    @param fmax - [GHz]
    @param n - number of points in the output function
	@return DeltaI = B_nu(Tcmb) * x e^x/(e^x-1) [W/m^s/Hz/sr]

	\date Jan 2, 2016, 5:06:21 PM
	\author Bartosz Lew
*/
mscsFunction calculate_deltaToT_to_intnsity_conversionFn(double fmin, double fmax, long n);

double HMSstrToAng(string HHMMSS,string separator=":");
double DMSstrToAng(string HHMMSS,string separator=":");

double doCal2JDconversion(string cal2JD);



int main(int argc, char **argv) {
	cpedsMsgs msgs("cosmocalc");
	double tmp;
	parseOptions(argc, argv);

	if (_radec2000toGal) {
		FILE* f=stdin;
		double ra,dec;
//		char rac[100],decc[100];
		string rastr,decstr;
		for (std::string line; std::getline(std::cin, line);) {
			QString str=line.c_str();
			str=str.simplified();
			QStringList qsl=str.split(" ");
			rastr=qsl[0].toStdString();
			decstr=qsl[1].toStdString();

			ra=0;
			dec=0;
			if (_inCoordFmt=="HMSdms") {
				ra=HMSstrToAng(rastr,":");
				dec=DMSstrToAng(decstr,":");
			}
			if (_inCoordFmt=="deg") {
				QString ras=rastr.c_str();
				QString decs=decstr.c_str();
				ra=ras.toDouble();
				dec=decs.toDouble();
			}
//			printf("%lf %lf\n",ra,dec);
			
			DirectionRaDec radec(ra*PI180,dec*PI180,JD2000);
			cpedsDirection lb=radec.toGal();
			lb*=PI180inv;
			printf("%lE %lE\n",lb.lon(),lb.lat());
		}
//		fclose(f);
		exit(0);
	}
	
	if (_radec2000toNow) {
		FILE* f=stdin;
		double ra,dec;
//		char rac[100],decc[100];
		string rastr,decstr;
		for (std::string line; std::getline(std::cin, line);) {
			QString str=line.c_str();
			str=str.simplified();
			QStringList qsl=str.split(" ");
			rastr=qsl[0].toStdString();
			decstr=qsl[1].toStdString();

			ra=0;
			dec=0;
			if (_inCoordFmt=="HMSdms") {
				ra=HMSstrToAng(rastr,":");
				dec=DMSstrToAng(decstr,":");
			}
			if (_inCoordFmt=="deg") {
				QString ras=rastr.c_str();
				QString decs=decstr.c_str();
				ra=ras.toDouble();
				dec=decs.toDouble();
			}
//			printf("%lf %lf\n",ra,dec);
			double JDut1=cpeds_julian_time();
			if (_JDut1!=0) JDut1=_JDut1;
			DirectionRaDec radec(ra*PI180,dec*PI180,JD2000);
			radec.toJD(JDut1);
			if (_nutate) radec.nutate(0);
			if (_aberrate) radec.aberrate(radec.epoch());
			radec.print_direction("radec now");
			radec.unaberrate(radec.epoch());
			radec.nutate(-1);
			radec.print_direction("unaberrated, mean: radec now");
			
			double rt4_ra,rt4_dec,rt4_ra2,rt4_dec2;
			rt4_ra=ra;
			rt4_dec=dec;
			rt4_ra=RT4_PrecessionRA(JD2000,JDut1,ra,dec);
			rt4_dec=RT4_PrecessionDEC(JD2000,JDut1,ra,dec);
			printf("RT4 precessed ra: %lf dec: %lf\n",rt4_ra, rt4_dec);
			rt4_ra2=rt4_ra; //deg
			rt4_dec2=rt4_dec; //deg
			double dpsi,deps,eObliquity;
			cpeds_rt4_nutation(JDut1,&rt4_ra2,&rt4_dec2,&dpsi,&deps,&eObliquity);
			printf("RT4 nutation, dphi [\"]: %lf deps [\"]: %lf, eObliquity: %lf\n",dpsi,deps,eObliquity);
			
			printf("RT4 precessed + nutation, ra: %lf dec: %lf\n",rt4_ra2, rt4_dec2);

			cpeds_rt4_nutation_and_aberration(JDut1,&rt4_ra,&rt4_dec);

			printf("RT4 precessed + nutation + aberration, ra: %lf dec: %lf\n",rt4_ra, rt4_dec);
			printf("ln_get_jde 	( 	JD2000=%.14lf	): %.14lf\n",JD2000,ln_get_jde 	( 	JD2000	));
			//			printf("%lE %lE\n",lb.lon(),lb.lat());
		}
//		fclose(f);
		exit(0);
	}
	
	
	
	if (_dontSaveRun == false) {
		msgs.setSaveRunWriteMode('a');
		msgs.saveThisRun(argc, argv, ".runlog");
	}
	
	if (_z > 10000) CPEDS_ZINF = 1e10;
	
	cpeds_set_points_number_per_logz(_NP);
	
	if (_script_mode) {
		msgs.warning(
				"We are in the transition stage. Not all messages obey the script_mode settings or verbocity settings",
				High);
	}
	
	if (_verbosity > 0) {
		printf("COSMOLOGICAL PARAMETERS USED:\n");
		printf("Wr0: %lE\n", _Wr0);
		printf("Wb0: %lE\n", _Wb0);
		printf("Wcdm0: %lE\n", _Wcdm0);
		printf("Wm0: %lE\n", _Wm0);
		printf("Wl0: %lE\n", _Wl0);
		printf("Wtot0: %lE\n", _Wl0 + _Wm0 + _Wr0);
		printf("w: %lE\n", _w0);
		printf("h: %lE\n", _h);
		printf("z: %lE\n", _z);
		
		printf("\n\n");
	}
	
	if (_dp) {
		printf("particle horizon at redshift %lE [Gpc]: %lE\n", _z,
				cpeds_dp(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
	}
	if (_lbtd) {
		double lbtd = cpeds_lbtd(_Wr0, _Wm0, _Wl0, _z, _w0, _h);
		printf("lookback time distance to redshift %lE [Gpc]: %lE\n", _z, lbtd);
		printf("lookback time distance to redshift %lE [Gly]: %lE\n", _z,
				lbtd / (CPEDS_c * GYR / CPEDS_MPC / 1e3));
	}
	if (_X) {
		printf("comoving distance to redshift %lE  [Gpc]: %lE\n", _z,
				cpeds_comoving_distance(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
	}
	if (_Xrec) {
		printf("comoving distance to recombination  [Gpc]: %lE\n",
				cpeds_comoving_distance_rec(_Wr0, _Wm0, _Wl0, _w0, _h, 1.0));
	}
	if (_Rc) {
		printf("cutvature radius at redshift %lE [Gpc]: %lE\n", _z,
				cpeds_curvature_radius(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
	}
	if (_dL) {
		printf("luminosity distance to redshift %lE [Gpc]: %lE\n", _z,
				cpeds_luminosity_distance(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
	}
//	if (_dL) { say(msgs,"luminosity distance to redshift %lE [Gpc]: %lE\n",_z, cpeds_luminosity_distance(_Wr0,_Wm0,_Wl0,_z,_w0,_h)); }
	if (_dA) {
		printf("angular size distance to redshift %lE [Mpc]: %lE\n", _z,
				1000.0
						* cpeds_angular_diameter_distance(_Wr0, _Wm0, _Wl0, _z,
								_w0, _h));
	}
	if (_dac) {
		printf("sound horizon at redshift %lE [Gpc]: %lE\n", _z,
				cpeds_sound_horizon(_Wr0, _Wb0, _Wm0, _Wl0, _z, _w0, _h));
	}
	if (_lac) {
		printf("acoustic horizon multipole at redshift %lE : %lE\n", _z,
				cpeds_acoustic_horizon_multipole(_Wr0, _Wb0, _Wm0, _Wl0, _z,
						_w0, _h));
	}
	if (_thac) {
		printf("acoustic horizon angle at redshift %lE : %lE\n", _z,
				cpeds_acoustic_horizon_angle(_Wr0, _Wb0, _Wm0, _Wl0, _z, _w0,
						_h));
	}
	if (_cs) {
		printf("sound speed in c at  redshift %lE : %lE\n", _z,
				cpeds_sound_speed(
						cpeds_Wr(_Wr0, _Wm0, _Wl0, _z, _w0)
								* cpeds_rhoC(_Wr0, _Wm0, _Wl0, _z, _w0),
						cpeds_Wb(_Wr0, _Wb0, _Wl0, _z, _w0)
								* cpeds_rhoC(_Wr0, _Wm0, _Wl0, _z, _w0)));
	}
	if (_ang) {
		printf(
				"angular size [deg] of an object of size %lE [Mpc] at redshift %lE : %lE\n",
				_size, _z,
				2.0 * PI180inv
						* atan(
								_size / 2
										/ (1000.0
												* cpeds_angular_diameter_distance(
														_Wr0, _Wm0, _Wl0, _z,
														_w0, _h))));
		//	  printf("angular size [deg] of an object of size %lE [Mpc] at redshift %lE : %lE\n",_size,_z, PI180inv*_size/(1000.0*cpeds_angular_diameter_distance(_Wr0,_Wm0,_Wl0,_z,_w0,_h))); 
	}
	if (_Tcmb) {
		printf("CMB temperature at redshift %lE : %lE (%lE eV)\n", _z,
				_Tcmb0 * (1 + _z), _Tcmb0 * (1 + _z) * CPEDS_kB / CPEDS_e);
	}
	if (_rhoC) {
		printf("Critical density at redshift %lE [kg/m^3]: %lE\n", _z,
				cpeds_rhoC(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
		printf("Critical density at redshift %lE [g/cm^3]: %lE\n", _z,
				1.0e-3 * cpeds_rhoC(_Wr0, _Wm0, _Wl0, _z, _w0, _h));
		printf("Critical density at redshift %lE [proton/m^3]: %lE\n", _z,
				cpeds_rhoC(_Wr0, _Wm0, _Wl0, _z, _w0, _h) / CPEDS_m_p);
	}
	
	//  if (_angdA>0) { double dA=cpeds_get_dA_from_ang(_Wr0,_Wb0,_Wm0,_Wl0,&_z,_w0,_h,_size);  printf("angular size distance to object of size %lE located at redshift %lE is [Mpc]: %lE\n",_size,_z,dA ); } 
	
	if (_cmbfMax) {
		printf(
				"frequency at which CMB spectral energy density peaks at redshift %lf is [GHz] : %lE\n",
				_z, cpeds_cmbfMax(_Tcmb0 * (1 + _z)) * 1e-9);
		printf("the corresponding wavelenght is [mm] : %lE\n",
				CPEDS_c / cpeds_cmbfMax(_Tcmb0 * (1 + _z)) * 1e3);
	}
	
	if (_a) {
		cpedsFunctionCosmo at("a vs t");
		at.make_a_vs_t(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h, _NP);
		at.save(_outputFile + "--a_vs_t.txt");
	}
	
	if (_a0H0inv) {
		printf(
				"co-moving length scale that enters horizon today (a0 H0)^-1 [Mpc h^-1] : %lE\n",
				CPEDS_c / (_a0 * 100000));
	}
	
	if (_q) {
		/* printf("q at redshift %lE : %lE\n",_z, cpeds_deceleration_parameter_q(_Wr0,_Wm0,_Wl0,_z,_w0));  */
		cpedsFunctionCosmo qt("q vs t");
		cpedsFunctionCosmo *qvsz = new cpedsFunctionCosmo("q vs z");
		qt.make_q_vs_t(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h, _NP, qvsz);
		qt.save(_outputFile + "--q_vs_t.txt");
		qvsz->save(_outputFile + "--q_vs_z.txt");
	}
	
	if (_dAvsz) {
		cpedsFunctionCosmo dAvsz("dA vs z");
		//	    cpedsFunctionCosmo dAvst("dA vs t");
		dAvsz.make_dA_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h);
		dAvsz.save(_outputFile + "--dA_vs_z.txt");
		//	    dAvst.save(_outputFile+"--dA_vs_t.txt");
	}
	if (_Xvsz) {
		cpedsFunctionCosmo dAvsz("dA vs z");
		cpedsFunctionCosmo Xvsz("X vs z");
		dAvsz.make_dA_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h, 0, &Xvsz);
		Xvsz.save(_outputFile + "--X_vs_z.txt");
		//	    dAvst.save(_outputFile+"--dA_vs_t.txt");
	}
	if (_lbtdvsz) {
		cpedsFunctionCosmo lbtdvsz("lbtd vs z");
		lbtdvsz.make_lbtd_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h, _NP);
		lbtdvsz.save(_outputFile + "--lbtd_vs_z.txt");
		//	    dAvst.save(_outputFile+"--dA_vs_t.txt");
	}
	
	if (_Wvsz) {
		cpedsFunctionCosmo Wxvsz("Omega_x vs z");
		Wxvsz.make_Wr_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0);
		Wxvsz.save(_outputFile + "--Wr_vs_z.txt");
		Wxvsz.make_Wm_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0);
		Wxvsz.save(_outputFile + "--Wm_vs_z.txt");
		Wxvsz.make_Wl_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0);
		Wxvsz.save(_outputFile + "--Wl_vs_z.txt");
		Wxvsz.make_Wtot_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0);
		Wxvsz.save(_outputFile + "--Wtot_vs_z.txt");
		Wxvsz.make_Wk_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0);
		Wxvsz.save(_outputFile + "--Wk_vs_z.txt");
		
	}
	
	if (_visibility) {
		//    cpedsFunctionCosmo visibt("visibility function vs t");
		//    cpedsFunctionCosmo *visibz=new cpedsFunctionCosmo("visibility function vs z");
		//    visibt.make_visibility_vs_t(_zmin,_zmax,_Wr0,_Wm0,_Wl0,_w0,_h,_NP,visibz);
		//    visibt.save(_outputFile+"--visibility_vs_t.txt");
		//    visibz->save(_outputFile+"--visibility_vs_z.txt");
	}
	
	if (_Wr) {
		printf("Omega_r at redshift %lE : %lE\n", _z,
				cpeds_Wr(_Wr0, _Wm0, _Wl0, _z, _w0));
	}
	if (_Wm) {
		printf("Omega_m at redshift %lE : %lE\n", _z,
				cpeds_Wm(_Wr0, _Wm0, _Wl0, _z, _w0));
		printf("Omega_cdm at redshift %lE : %lE\n", _z,
				cpeds_Wm(_Wr0, _Wm0, _Wl0, _z, _w0) * _Wcdm0 / _Wm0);
	}
	if (_Wb) {
		printf("Omega_b at redshift %lE : %lE\n", _z,
				cpeds_Wm(_Wr0, _Wm0, _Wl0, _z, _w0) * _Wb0 / _Wm0);
	}
	if (_Wl) {
		printf("Omega_l at redshift %lE : %lE\n", _z,
				cpeds_Wl(_Wr0, _Wm0, _Wl0, _z, _w0));
	}
	if (_Wk) {
		printf("Omega_k at redshift %lE : %lE\n", _z,
				cpeds_Wk(_Wr0, _Wm0, _Wl0, _z, _w0));
	}
	if (_Wtot) {
		printf("Omega_tot at redshift %lE : %lE\n", _z,
				cpeds_Wtot(_Wr0, _Wm0, _Wl0, _z, _w0));
	}
	
	if (_refraction != "") {
		doRefraction(msgs);
	}
	if (_refrMoonAtDusk!= "") {
		doMoonRefractionAtDusk(msgs);
	}
	if (_testCoor) {
		doModel4();
	}
	if (_testVtotRT32) {
		doTestVtotRT32();
	}
	
	if (_agez) {
		tmp = cpeds_age_of_universe(_Wr0, _Wm0, _Wl0, _z, _w0, _h);
		if (_script_mode) say(msgs, "", "age at redshift", msgs.toStr(tmp));
		else {
			if (tmp > 1.0) {
				printf("age at redshift %lE [Gyr] : %lE\n", _z, tmp);
			} else {
				if (tmp > 1e-3) {
					printf("age at redshift %lE [Myr] : %lE\n", _z, tmp * 1e3);
				} else {
					if (tmp > 1e-6) {
						printf("age at redshift %lE [kyr] : %lE\n", _z,
								tmp * 1e6);
					} else {
						if (tmp > 1e-9) {
							printf("age at redshift %lE [yr] : %lE\n", _z,
									tmp * 1e9);
						} else {
							//						if (tmp > 1e-9*365*24*3600) {				
							printf("age at redshift %lE [s] : %lE\n", _z,
									tmp * 1e9 * 365 * 24 * 3600);
							//						}
							
						}
					}
				}
			}
		}
	}
	
	if (_K2KeV != 0) {
		printf("gas temperature for T=%lE [K] in eV is: %lE\n", _K2KeV,
				cpeds_K2eV(_K2KeV));
	}
	if (_Efact) {
		printf("E(z=%lE): %lE\n", _z, cpeds_Efactor(_Wr0,_Wm0,_Wl0,_z,_w0));
	}
	
	if (_partMass) {
		double mb, mcdm;
		double rhoC = cpeds_rhoC0(_h);
		mb = _Wb0 * rhoC / (1.0e10 * CPEDS_SOLAR_MASS) * pow(CPEDS_MPC, 3)
				* pow(_Lmin / _cells, 3);
		mcdm = _Wcdm0 * rhoC / (1.0e10 * CPEDS_SOLAR_MASS) * pow(CPEDS_MPC, 3)
				* pow(_Lmin / _cells, 3);
		printf("mb [10^10 Msol/h]: %lE\n", mb * _h);
		printf("mcdm [10^10 Msol/h]: %lE\n", mcdm * _h);
		printf("mb [10^10 Msol]: %lE\n", mb);
		printf("mcdm [10^10 Msol]: %lE\n", mcdm);
	}
	
	if (_lbdeg2dms != "") {
		cpedsDirection n;
		QString s = _lbdeg2dms.c_str();
		QStringList qsl = s.split(',');
		n.lon() = qsl[0].toDouble() * PI180;
		n.lat() = qsl[1].toDouble() * PI180;
		n.print_direction();
	}
	
	if (_RTres) {
		printf("aperture [m]: %lf\n", _size);
		printf("frequency [GHz]: %lf\n", _freq);
		double res = cpeds_calculateAngularResolution(_size, _freq);
		printf("radio telescope angular resolution [deg]: %lf\n", res);
		printf("radio telescope angular resolution [arcmin]: %lf\n", res * 60);
		printf("radio telescope angular resolution [arcsec]: %lf\n",
				res * 3600);
		printf("radio telescope beam solid angle [sr]: %lE\n",
				pow(res * PI180, 2));
		
		double lambda = CPEDS_c / (_freq * 1e9);
		double size_sub = _size / 10;
		double taper = _taper;
		double thetaMax = res * 5;
		mscsFunction Pb = calculate_beam_power_pattern(_size, size_sub, lambda,
				taper, thetaMax, thetaMax / 1000);
		string fname = "beam_power_pattern";
		fname += "-size_" + msgs.toStrf(_size, 0) + "m-freq_"
				+ msgs.toStrf(_freq, 3) + "GHz-taper_" + msgs.toStrf(_taper, 3)
				+ "db";
//		printf("beta: %lf\n",beta);
		Pb.save(fname);
		Pb -= 0.5;
		mscsVector<double> fwhp = Pb.findRoot();
		if (fwhp.size() == 2) {
			printf("derived FWHP for given taper [deg]: %lf\n",
					fwhp[1] - fwhp[0]);
		}
	}
	
	if (_sigma0) {
		caclulate_sigma0(_s0scale, msgs, NULL);
	}
	if (_sigmaR) {
		caclulate_sigmaR(msgs);
	}
	
	if (_y0) {
		double rhoC0 = cpeds_rhoC0(_h);
		double y0 = CPEDS_sigma_T * CPEDS_kB * rhoC0 * _Wb0
				/ (_mueinv1 * CPEDS_m_p * CPEDS_m_e * CPEDS_c * CPEDS_c)
				* CPEDS_MPC;
		double y0ksz=y0*CPEDS_m_e * CPEDS_c * CPEDS_c / CPEDS_kB;
		printf("y0 [K^-1 MPC^-1]: %lE\n", y0);
		printf("y0ksz [MPC^-1]: %lE\n", y0ksz);
		
		
	}
	
	if (_TSZfreq > 0) {
		double x = CPEDS_h * _TSZfreq * 1.0e9 / (CPEDS_kB * _Tcmb0);
		double fnu = x * cosh(x / 2) / sinh(x / 2) - 4;
		double gnu = x * x * x * x * exp(x) / pow(exp(x) - 1, 2);
		double I0=2.0*pow(CPEDS_kB*_Tcmb0,3)/pow(CPEDS_h*CPEDS_c,2)/1.0e-26; // Jy/sr
		printf("TSZ frequency factor fnu for f=%lf GHz (x=%lf, lambda=%lf mm): %lE\n",_TSZfreq, x, CPEDS_c / (_TSZfreq * 1.0e9) * 1e3, fnu);
		printf("TSZ intensity factor gnu for f=%lf GHz (x=%lf, lambda=%lf mm): %lE\n",_TSZfreq, x, CPEDS_c / (_TSZfreq * 1.0e9) * 1e3, gnu);
		printf("TSZ intensity factor I0 [Jy/sr]: %lE\n", I0);
		printf("TSZ intensity factor I0 [Jy/arcmin2]: %lE\n", I0/pow(PI180inv*60.0,2));
		printf("Delta I/ (deltaT) = I0*gnu [Jy/sr]: %lE\n", I0*gnu);
		printf("TSZfrequency factor fnu*gnu for f=%lf GHz (x=%lf, lambda=%lf mm): %lE\n",
				_TSZfreq, x, CPEDS_c / (_TSZfreq * 1.0e9) * 1e3, fnu * gnu);
		double yTSZ=1e-4;
		printf("(Delta T)_TSZ/(Delta I)_TSZ [mK/mJy/sr]:%lE \n",	_Tcmb0/gnu/I0);
		if (_TSZfreqFn) {
			mscsFunction f = calculate_TSZfreqFn(_fmin, _fmax, _NP);
			f.save("TSZfrequency-dependence.txt");
			// calculate temperature change duet to TSZ for compton y-parmeter y=1e-4
			f*=yTSZ*_Tcmb0*1000;
			f.save("TSZ-DeltaT-y_1e-4-frequency-dependence.txt");			
			f.clearFunction();
			// intensity change due to TSZ
			f = calculate_TSZintensityFn(_fmin, _fmax, _NP);
			f*=yTSZ;
			f/=CPEDS_Jy; // convert to Jy/sr
			f.save("TSZ-DeltaI-y_1e-4-frequency-dependence.txt");			
			f.clearFunction();
			
			f.mkPlank(_fmin*1e9,_fmax*1e9,(_fmax-_fmin)*1e9/_NP,_Tcmb0);
			f.scaleX(1e-9); // convert to GHz
			f*=double(CPEDS_c/fourPI); // convert to intensity [W/m^2/Hz/sr]
			f/=CPEDS_Jy;
			f.save("CMBintensity-frequency-dependence.txt");
			f.derivative(true);
			f.save("CMBintensity_derivative-frequency-dependence.txt");
			f.clearFunction();
//		}
//		if (_kSZfreqFn) {
			//
			// kinetic effect
			//
			// Delta T / T =~ - ykSZ (1+z) int betall rho/rho_av dX
			// betall = 1/500, rho/rho_av = 500
			double betall=_kSZbetall;
			double yKSZ = cpeds_calculate_comptonY(_Wb0,_h,_mueinv1) * CPEDS_m_e*CPEDS_c*CPEDS_c/CPEDS_kB; // Mpc^-1
			double dX = 2; // Mpc
			double drhorhoav=500;
			yKSZ=0.01;
			printf("y^kSZ: %lE\n",yKSZ);
//			double dTkSZ=-yKSZ*(1+_z) * betall * drhorhoav * dX * _Tcmb0; // K
			double DTkSZ=-yKSZ*betall*_Tcmb0;
			printf("(Delta T)_kSZ [mK]: %lE\n",DTkSZ*1000);
			f.mkConst(_fmin,_fmax,(_fmax-_fmin)/_NP,DTkSZ*1000);
			f.save("kSZ-DeltaT-beta_"+msgs.toStrf(betall,4)+".txt");
			f.clearFunction();
			// calculate kSZ intensity change frequency dependence
			f=calculate_deltaToT_to_intnsity_conversionFn(_fmin,_fmax,_NP);
			f.save("deltaT2DeltaIconversion.txt");
			f*=DTkSZ/_Tcmb0;
			f.save("kSZ-DeltaI-beta_"+msgs.toStrf(betall,4)+".txt");
			f.derivative(true);
			f.save("kSZ-DeltaI-deriv-beta_"+msgs.toStrf(betall,4)+".txt");
		}
	}
	
	if (_printGadgetUnits) {
		printf("GADGET2_UnitMass_in_g: %lE\n", GADGET2_UnitMass_in_g);
		printf("GADGET2_UnitLength_in_cm: %lE\n", GADGET2_UnitLength_in_cm);
		printf("GADGET2_UnitVelocity_in_cm_per_s: %lE\n",
				GADGET2_UnitVelocity_in_cm_per_s);
		printf("GADGET2_UnitTime_in_s: %lE\n", GADGET2_UnitTime_in_s);
		printf("GADGET2_UnitEnergy_in_cgs: %lE\n", GADGET2_UnitEnergy_in_cgs);
		printf("GADGET2_UnitDensity_in_cgs: %lE\n", GADGET2_UnitDensity_in_cgs);
	}
	
	if (_BBKS_transferFunction) {
		cpedsFunctionCosmo bbks("bbks");
		bbks.make_BBKS_transferFunction(1e-4, 10, 0.01, _Wb0, _Wtot0, _h);
		bbks.save("BBKS.tf");
	}
	
	if (_lightCone) {
		
	}
	
	if (_boxesInField) {
		
		double z = _zmax;
		cpedsFunctionCosmo dAvsz("dA vs z");
		cpedsFunctionCosmo Xvsz("X vs z");
		msgs.say("deriving X vs z", High);
		dAvsz.make_dA_vs_z(_zmin, _zmax, _Wr0, _Wm0, _Wl0, _w0, _h, 0, &Xvsz);
		//	  // convert to Gpc
		//	  Xvsz/=double(1000);
		msgs.say("sorting X vs z", High);
		Xvsz.sortFunctionArgAscending();
		//	  Xvsz.invert();
		double *Z = Xvsz.extractArguments();
		double *chi = Xvsz.extractValues();
		long N = Xvsz.pointsCount();
		
		msgs.say("numer of cells along dimension: " + msgs.toStr(_cells), High);
		double D0; // comoving box size of the simulation
		double X0; // comoving distance to the center of the zeroth box
		double Di; // i'th comoving box size
		double Xi; // i'th comoving distance to the center of the i'th box
		
		mscsFunction X, D;
		D0 = (1 + _zmax) * _angdA * PI180
				* cpeds_angular_diameter_distance(_Wr0, _Wm0, _Wl0, _zmax, _w0,
						_h) * 1000;
		X0 = cpeds_comoving_distance(_Wr0, _Wm0, _Wl0, _zmax, _w0, _h) * 1000;
		X.newPoint(_zmax, X0);
		D.newPoint(_zmax, D0);
		
		long iz = 0;
		do {
			iz++;
			printf(
					"%li) z: %lf, X [Mpc]: %lf, L [Mpc]: %lf    || cell size[Mpc]: %lf,  k1=1/L[Mpc^-1]: %lf  kmax[Mpc^-1]: %lf\n",
					iz, z, X0, D0, D0 / _cells, 1.0 / D0, _cells / D0);
			if (D0 > _Lmin) {
				Xi = (2 * X0 - D0) / (2 * X0 + D0) * X0;
				Di = (2 * X0 - D0) / (2 * X0 + D0) * D0;
				//		  z=Xvsz.f(double(Xi),&iz);
				z = Z[cpeds_find_value(Xi, chi, N, 0, N)];
				X.newPoint(z, Xi);
				D.newPoint(z, Di);
				X0 = Xi;
				D0 = Di;
			} else {
				Xi = X0 - Di;
				z = Z[cpeds_find_value(Xi, chi, N, 0, N)];
				X.newPoint(z, Xi);
				D.newPoint(z, Di);
				X0 = Xi;
				D0 = Di;
			}
		} while (z > Z[0]);
		
		long n = X.pointsCount();
		matrix<double> m(n, 3);
		for (long i = 0; i < n; i++) {
			m(i, 0) = X.getX(i);
			m(i, 1) = X.f(i);
			m(i, 2) = D.f(i);
		}
		
		cpeds_matrix_save(m, "simulation-boxes");
		
		delete[] Z;
		delete[] chi;
	}
	
	if (_gravPot != "") {
		calculateGravitationalPotential(_gravPot);
	}
	
	if (_l2nu > 0) {
		printf("frequency [GHz]: %f\n", CPEDS_c / (_l2nu * 1.0e-2) / 1.0e9);
	}
	if (_nu2l > 0) {
		printf("lambda [cm]: %f\n", CPEDS_c / (_nu2l * 1.0e9) * 100);
	}
	
	if (_JD) {
		double JD_UTC=cpeds_julian_time();
		printf("JD: %.15lf (%s)\n",JD_UTC,cpeds_JDToYMDhms(JD_UTC).c_str());
	}
	
	if (_doJD2cal) {
		printf("%s\n",cpeds_JDToYMDhms(_JD2cal).c_str());
		
		rt4_time_date dt=cpeds_RT4_cal_date(_JD2cal);
		printf("%i-%i-%i\n",dt.year,dt.month,dt.day);
		printf("JD=%.14f\n",cpeds_rt4_JulianDay(dt.year,dt.month,dt.day));
//		printf("JD=%.14f\n",JulianDay(2000,0,0));
		
	}
	if (_cal2JD!="") {
		printf("%.15lf\n",doCal2JDconversion(_cal2JD));
	}
	
	if (_dou2cal) {
		printf("%s\n",cpeds_unixToYMDhms(_u2cal).c_str());
		
//		rt4_time_date dt=cpeds_RT4_cal_date(_JD2cal);
//		printf("%i-%i-%i\n",dt.year,dt.month,dt.day);
//		printf("JD=%.14f\n",cpeds_rt4_JulianDay(dt.year,dt.month,dt.day));
	}

	return 0;
}

void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd(
				"cosmological calculator\n calculates vatious cosmological distances and other useful quantities from well known formulae",
				' ', "");
		
		// 
		// Define arguments
		//

		SwitchArg dontSaveRun("", "dontSaveRun",
				"whether or not save how we were called. It may be important to disable logging in case of very frequency calls to this program "
						"from eg python script.", false);
		cmd.add(dontSaveRun);
		
		// verbosity
		ValueArg<int> verbosity("v", "verb", "level of vebocity. The lower the quieter (minimal 0). Useful for using in scripts", false, 1, "int");
		cmd.add(verbosity);
		SwitchArg script_mode("", "script",	"putputs in script mode. Does not make any unit conversions to human readable form. Useful for scripts", false); cmd.add(script_mode);
		
		// distances
		SwitchArg dp("", "dp", "calculates particle horizon upto given redshift", false); cmd.add(dp);
		
		SwitchArg lbtd("", "lbtd", "calculates look-back time distance upto given redshift from today",	false);	cmd.add(lbtd);
		SwitchArg X("", "X","calculates comoving distance (conformal distance) upto given redshift from today",	false); cmd.add(X);
		SwitchArg Xrec("", "Xrec","calculates  comoving distance (conformal distance) distance  @ recombination",false); cmd.add(Xrec);
		SwitchArg dA("", "dA","calculates angular diameter distance to requested redshift",	false);	cmd.add(dA);
		SwitchArg dL("", "dL","calculates luminosity distance to requested redshift", false); cmd.add(dL);
		SwitchArg dac_rec("", "dac_rec","calculates sound horizon at decoupling", false); cmd.add(dac_rec);
		SwitchArg dac("", "dac", "calculates sound horizon", false); cmd.add(dac);
		SwitchArg Efact("", "Efact", "calculates Efactor", false); cmd.add(Efact);
		
		// time
		SwitchArg agez("", "agez","calculates age of the Universe for a given redshift", false); cmd.add(agez);
		SwitchArg agea("", "agea", "calculates age of the Universe for a given scale factor",false); cmd.add(agea);
		
		// dependent params
		SwitchArg Wtot0("", "Wtot0", "calculates Omega_tot0", false); cmd.add(Wtot0);
		SwitchArg Wtot("", "Wtot", "calculates Omega_tot for given redshift", false); cmd.add(Wtot);
		SwitchArg Wk("", "Wk", "calculates Omega_k for given redshift", false);	cmd.add(Wk);
		SwitchArg Wr("", "Wr", "calculates Omega_r for given redshift", false);	cmd.add(Wr);
		SwitchArg Wm("", "Wm", "calculates Omega_m for given redshift", false);	cmd.add(Wm);
		SwitchArg Wb("", "Wb", "calculates Omega_b for given redshift", false);	cmd.add(Wb);
		SwitchArg Wl("", "Wl", "calculates Omega_de for given redshift", false); cmd.add(Wl);
		SwitchArg Rc("", "Rc", "calculates curavture radius for given redshift", false); cmd.add(Rc);
		SwitchArg lac("", "lac", "acoustic horizon multipole", false); cmd.add(lac);
		SwitchArg thac("", "thac", "acoustic horizon angle at redshift z", false); cmd.add(thac);
		SwitchArg cs("", "cs", "sound speed at redshift z in units of c", false); cmd.add(cs);
		ValueArg<double> size("", "size", "physical size of an object of given size in Mpc at given redshift z ( default: 100 Mpc )", false, -1, "double");	cmd.add(size);
		SwitchArg ang("", "ang", "calculate angular size of an object (in deg) at redshift z of size given by size", false); cmd.add(ang);
		ValueArg<double> angdA("", "angdA",	"angle [rad] for calculation the corresponding dA and z for provided size",	false, -1, "double"); cmd.add(angdA);
		SwitchArg dAvsz("", "dAvsz", "derive the angular size distance [Mpc] vs redshift ", false);	cmd.add(dAvsz);
		SwitchArg Xvsz("", "Xvsz",	"derive the comoving distance [Mpc] vs redshift ", false);	cmd.add(Xvsz);
		SwitchArg lbtdvsz("", "lbtdvsz","derive the comoving distance [Mpc] vs redshift ", false); cmd.add(lbtdvsz);
		SwitchArg Wvsz("", "Wvsz", "derive the Omega_rml vs redshift ", false);	cmd.add(Wvsz);
		
		SwitchArg Tcmb("", "Tcmb", "calculate CMB temperature at redshift z",false); cmd.add(Tcmb);
		ValueArg<double> Tcmb0("", "Tcmb0",	"current CMB temperature  ( default: 2.726 K )", false, 2.726, "double"); cmd.add(Tcmb0);
		SwitchArg cmbfMax("", "cmbfMax","calculate the frequency and corresponding wavelenght at which CMB energy density spectrum peaks",false); cmd.add(cmbfMax);
		SwitchArg rhoC("", "rhoC", "calculates the critical density for given cosmology [kg/m^3]", false); cmd.add(rhoC);
		SwitchArg radec2000toGal("", "radec2000toGal","convert stdin ra2000,dec2000 coordinates to l b",false); cmd.add(radec2000toGal);
		SwitchArg radec2000toNow("", "radec2000toNow","convert stdin ra2000,dec2000 coordinates to ra,dec at current UTC time epoch",false); cmd.add(radec2000toNow);
		SwitchArg nutate("", "nutate","radec2000toNow sub-option: include nutation in conversion",false); cmd.add(nutate);
		SwitchArg aberrate("", "aberrate","radec2000toNow sub-option: include annual aberration in conversion",false); cmd.add(aberrate);
		
		ValueArg<string> inCoordFmt("", "inCoordFmt","sub-option for --radec2000toGal. Defines input coordinates format (default: deg)."
				"Other possible values are: HMSdms for HH:MM:SS.S dd:mm:ss.s",false,"deg","string"); cmd.add(inCoordFmt);

		
		ValueArg<double> K2KeV("", "K2kev",	"prints temperature in KeV given temperature in Kelvins for electron gas", false, 0, "double");	cmd.add(K2KeV);
		
		SwitchArg vsz("", "vsz","calculates requested quantities in function of redshift in ranges given by zmin zmax",	false);	cmd.add(vsz);
		SwitchArg vsa("", "vsa","calculates requested quantities in function of scale factor in ranges given by amin amax",	false);	cmd.add(vsa);
		SwitchArg vst("", "vst","calculates requested quantities in function of cosmological time in ranges given by tmin tmax",false); cmd.add(vst);
		SwitchArg q("q", "deceleration","calculates the deceleration parameter as a function of time [Gyr] and z in redshift range from zmin to zmax > zmin",false); cmd.add(q);
		SwitchArg a("a", "scaleFactor",	"calculates the scale factor as a function of time [Gyr] in redshift range from zmin to zmax > zmin",false);cmd.add(a);
		SwitchArg visibility("", "visibility",	"calculates the visibility function as a function of time [Gyr] and z in redshift range from zmin to zmax > zmin",	false); cmd.add(visibility);
		SwitchArg sigma0("", "sigma0", "calculates current sqrt(variance) (sigma0(R)) of the density perturbations contrast smoothed with a filter of choice (gaussian by default"
				"but can also be top-hat --WkRTopHat) at a requested scale\n"
						"For this to work you need to supply:\n"
						"- the transfer function file with tf parameter\n"
						"- the smoothing scale with s0scale parameter\n"
						"- and the primordial power spectrum parameters with Pk* parameters",
				false);
		cmd.add(sigma0);
		SwitchArg sigmaR("", "sigmaR","same as sigma0 but calculates as a function of scale R [Mpc h^-1]. This option also calculates Tinker mass function and Press-Schechter mass function"
						"for the requested cosmology."
						"The range of calculation is given by options\n"
						"- --rmin, --rmax, --dr\n", false);		cmd.add(sigmaR);
		ValueArg<double> s0scale("", "s0scale",	"comoving scale for smoothing [Mpc h^-1]", false, 0, "double");		cmd.add(s0scale);
		ValueArg<string> tf("", "tf","file with transfer function to use with arguments in units of k/h",	false, "", "string");		cmd.add(tf);
		ValueArg<long> tfcol("", "tfcol",
				"column in transfer function file to read. Useful in case of multiple column transfer function file as in case for CAMB\n"
						"in camb: \n"
						"col0 = k/h\n"
						"col1 = CDM tf\n"
						"col2 = baryon\n"
						"col3 = phoron\n"
						"col4 = massless neutrino\n"
						"col5 = massive neutrino\n"
						"col6 = total\n"
						"(default: 1)", false, 1, "long");		cmd.add(tfcol);
		SwitchArg notfnorm("", "notfnorm",	"do not normalize the input transfer function", false);		cmd.add(notfnorm);
		ValueArg<double> PkA("", "PkA",	"power spectrum amplitude - the amplitude of the comoving curvature perturbation Delta_R^2 WMAP9 = 2.41e-9 ",false, 2.41e-9, "double");	cmd.add(PkA);
		ValueArg<double> Pkk0("", "Pkk0",
				"power spectrum pivot scale: 0.002 Mpc^-1 for wmap9", false,
				0.002, "double");
		cmd.add(Pkk0);
		ValueArg<double> Pkns("", "Pkns",
				"power spectrum spectral index (0.972 from WMAP9)", false,
				0.963, "double");
		cmd.add(Pkns);
		SwitchArg BBKS_transferFunction("", "BBKStf",
				"calculate BBKS transfer function", false);
		cmd.add(BBKS_transferFunction);
		
		// ranges of dumps
		
		ValueArg<double> amin("", "amin",
				"minimal scale factor for tabular computations", false, 0,
				"double");
		cmd.add(amin);
		ValueArg<double> amax("", "amax",
				"maximal scale factor for tabular computations", false, 0,
				"double");
		cmd.add(amax);
		ValueArg<double> tmin("", "tmin",
				"minimal time for tabular computations", false, 0, "double");
		cmd.add(tmin);
		ValueArg<double> tmax("", "tmax",
				"maximal time for tabular computations", false, 0, "double");
		cmd.add(tmax);
		ValueArg<double> zmin("", "zmin",
				"minimal redshift for tabular computations", false, 0,
				"double");
		cmd.add(zmin);
		ValueArg<double> zmax("", "zmax",
				"maximal redshift for tabular computations", false, 0,
				"double");
		cmd.add(zmax);
		ValueArg<double> z("z", "z", "redshift ( default: 0 )", false, 0,
				"double");
		cmd.add(z);
		
		// cosmological stuff budget parameters
		ValueArg<double> Wr0("", "Wr0",
				"current radiation density ( default: 4.5e-5/h^2 )", false,
				4.5e-5, "double");
		cmd.add(Wr0);
		ValueArg<double> Wb0("", "Wb0",
				"current baryon matter density ( default: 0.045 )", false,
				0.044, "double");
		cmd.add(Wb0);
		ValueArg<double> Wcdm0("", "Wcdm0",
				"current dark matter density (CDM) ( default: 0.237 )", false,
				0.237, "double");
		cmd.add(Wcdm0);
		ValueArg<double> Wm0("", "Wm0",
				"current matter density (CDM+baryons). If this is set then Wcdm0 is set according to the value of Wb0 so that Wm0=Wb0+Wcdm0( default: 0.282 )",
				false, 0.282, "double");
		cmd.add(Wm0);
		ValueArg<double> Wl0("", "Wl0",
				"current dark energy density ( default: 0.7 )", false, 0.719,
				"double");
		cmd.add(Wl0);
		SwitchArg Wl0fit2flat("", "Wl0fit2flat",
				"current DE density will be set as the Wtot=1", false);
		cmd.add(Wl0fit2flat);
		//	ValueArg<double> a0("", "a0", "current expansion factor ( default: 1 )", false,1,"double");	cmd.add( a0);
		
		// DE EOS
		ValueArg<double> w0("", "w0", "EOS for DE ( default: -1 )", false, -1,
				"double");
		cmd.add(w0);
		
		ValueArg<double> h("", "h",
				"reduced current Hubble constant ( default: 0.7 )", false, 0.71,
				"double");
		cmd.add(h);
		ValueArg<double> a0("", "a0", "scale factor today ( default: 1 )",
				false, 1, "double");
		cmd.add(a0);
		
		// ranges of lookback time distance and density of probing real space tf.
		ValueArg<double> rmin("", "rmin",
				"conformal minimal distance  ( default: 13000 [Mpc] )", false,
				13000, "double");
		cmd.add(rmin);
		ValueArg<double> rmax("", "rmax",
				"conformal maximal distance  ( default: 14500 [Mpc] )", false,
				14500, "double");
		cmd.add(rmax);
		ValueArg<double> dr("", "dr",
				"step with which to probe the ranges of rmin/max.  ( default: 7 [Mpc] )\n"
						"In case of --sigmaR option this is dlog10(r).", false,
				7, "double");
		cmd.add(dr);
		
		ValueArg<double> fmin("", "fmin",
				"minimal fraction of lbtd_SLS to alternatively define the rmax  ( default: 0.01  )",
				false, 0.01, "double");
		cmd.add(fmin);
		ValueArg<double> fmax("", "fmax",
				"maximal fraction of lbtd_SLS to alternatively define the rmin  ( default: 2  )",
				false, 2.0, "double");
		cmd.add(fmax);
		ValueArg<long> fN("", "fN",
				"step with which to probe the ranges of fmin/max  ( default: 100  )",
				false, 100, "long");
		cmd.add(fN);
		
		ValueArg<long> NP("", "NP",
				"number of points per log10 interval of the redshift space used for integrations  ( default: 200000  )",
				false, 200000, "long");
		cmd.add(NP);
		
		ValueArg<string> output_file("o", "output", "output file name prefix",
				false, "out", "string");
		cmd.add(output_file);
		ValueArg<string> lbdeg2dms("", "lbdeg2dms",
				"print direction in hms, dms format for input direction given in deg (eg. 12.12,13.13)",
				false, "", "string");
		cmd.add(lbdeg2dms);
		
		SwitchArg lightCone("", "lightCone",
				"calculates the simulation boxes sizes and redshifts to fill the comoving volume up to requested redshift within requested field of view and with requested thickness of a slice. The important parameters to be set are: angdA, min and zmax",
				false);
		cmd.add(lightCone);
		SwitchArg boxesInField("", "boxesInField",
				"calculates the simulation boxes sizes and redshifts to fill the comoving volume up to requested redshift within requested field of view. The important parameters to be set are: angdA, min and zmax",
				false);
		cmd.add(boxesInField);
		ValueArg<long> cells("", "cells",
				"number of cells in simulation along one dimention (only relevant if used with boxesInField option) (default: 512)",
				false, 512, "long");
		cmd.add(cells);
		ValueArg<double> Lmin("", "Lmin",
				"Minimal comoving size of the simulation box [Mpc] (default: 50)",
				false, 50, "double");
		cmd.add(Lmin);
		SwitchArg y0("", "y0",
				"calculates Compton y0-parameter and KSZ tau_0 parameter "
				"for the requested Hubble constant and baryon density Wb0",
				false);
		cmd.add(y0);
		ValueArg<double> mueinv1("", "mueinv1",
				"number of electrons per proton (default: 1.136)", false, 1.136,
				"double");
		cmd.add(mueinv1);
		SwitchArg partMass("", "partMass",
				"calculates the particle mass for the current cosmology in units of 10^10 Msol/h for baryons and CDM."
						"Combine this option with --Lmin and --cells.", false);
		cmd.add(partMass);
		SwitchArg a0H0inv("", "a0H0inv",
				"co-moving length scale that enters horizon today [Mpc h^-1]",
				false);
		cmd.add(a0H0inv);
		SwitchArg WkRTopHat("", "WkRTopHat",
				"trigger using top hat kernel instead of default gaussian for sigma_0 calculation",
				false);
		cmd.add(WkRTopHat);
		ValueArg<double> TinkerDeltaRho("", "TinkerDeltaRho",
				"Tinker mass function overdensity (default: 200)", false, 200,
				"double");
		cmd.add(TinkerDeltaRho);
		
		ValueArg<string> refraction("", "refr",
				"calculates zenith distance in space from observed zenith distance in atmosphere."
						"Accepts a string in format: ZD[deg],T[C],P[mbar],H[percent],lambda[cm],alt above sea[m],latitude[deg],Tlapse[K/m],accuracy."
						"(default: '' eg: '45,20,1013,100,1,113,54,0.0065,1e-8')",
				false, "", "string");
		cmd.add(refraction);
		SwitchArg refrfn("", "refrfn",
				"calculates various refraction dependencies", false);
		cmd.add(refrfn);
		ValueArg<string> refrMoonAtDusk("", "refrMoonAtDusk",
				"calculates Moon shape deformed by refraction at rise."
				"Accepts a string in format: ZDstart [deg], ZDend [deg], ZDstep [deg], Moon angular size [deg], T[C],P[mbar],H[percent],lambda[cm],alt above sea[m],latitude[deg],Tlapse[K/m],accuracy"
				"(default: '' eg: '80,92,0.1,0.5,20,1013,100,1,113,54,0.0065,1e-8')",
				 false, "","string");
		cmd.add(refrMoonAtDusk);
		SwitchArg testCoor("", "testCoor",
				"calculates correction table Model4 of the meteofix control system",
				false);
		cmd.add(testCoor);
		SwitchArg testVtotRT32("", "testVtotRT32",
				"calculates RT32 radial velocity wrt a fixed direction", false);
		cmd.add(testVtotRT32);
		SwitchArg printGadgetUnits("", "printGadgetUnits",
				"print Gadget Units defined in CPEDS", false);
		cmd.add(printGadgetUnits);
		ValueArg<double> TSZfreq("", "TSZfreq",
				"calculate SZ frequency factor [GHz]."
						"Additionally if fmin and fmax options are set, then this option saves the frequency relation"
						"for the requensted frequency range with --NP points.  (default: 0 )",
				false, 0, "double");
		cmd.add(TSZfreq);
		ValueArg<double> l2nu("", "l2nu",
				"wavelength [cm] to frequnecy [GHz] converter (default: -1 - not calculated)",
				false, -1, "double");
		cmd.add(l2nu);
		ValueArg<double> nu2l("", "nu2l",
				"frequnecy [GHz] to wavelength [cm] converter (default: -1 - not calculated)",
				false, -1, "double");
		cmd.add(nu2l);
		
		ValueArg<string> gravPot("", "gravPot",	"calculate gravitational potential for a set of points"
						"provided in a file with this option. The file should contain distribution of x y z mass."
						"Optionally use with --gravPotAt option to calculate the grativational potential at specified location. (this is not implemented yet)."
						"The provided file name has also a potential to generate dataset if special file name is provided."
						"Recognised special names are: mkUniBall.", false, "", "string");	cmd.add(gravPot);
		ValueArg<double> gravSoft("", "gravSoft","gravitational potential softening"
						"0 - no softening, -1 - use mean average interparticle separation as L/N, or specified as any other value"
						"(default: 0 )", false, 0, "double"); cmd.add(gravSoft);
		
		ValueArg<double> freq("", "freq",
				"frequency in [GHz].for various calculation, eg for calculation of cassegratin radio telescope angular resolution"
						"  (default: 30 )", false, 30, "double");	cmd.add(freq);
		ValueArg<double> taper("", "taper",	"taper in power for RT illumination at the edge of aperture used for --RTres calculation"
						" [db] (default: -12 )", false, -12, "double");	cmd.add(taper);
		
		SwitchArg RTres("", "RTres","calculate radio telescope angular resolution for given aperture and frequency. "
				"To specify aperture in [meters] use --size option"
						"and to specify frequency use --freq", false);	cmd.add(RTres);
		ValueArg<double> JD2cal("", "JD2cal","print date/time for given JD_UTC (default: dont do anything)", false, 0,"double");		cmd.add(JD2cal);
		ValueArg<string> cal2JD("", "cal2JD","print JD_UTC given date/time. Eg. 2016-05-13-17-40-01.234 (default: dont do anything)", false, "","string");		cmd.add(cal2JD);
		ValueArg<double> u2cal("", "u2cal","print date/time for given unix time (default: don't do anything). If set to -1"
				"then the current time is taken", false, -1,"double");		cmd.add(u2cal);
		SwitchArg JD("", "JD",	"print current julian date/time", false); cmd.add(JD);
		ValueArg<double> JDut1("", "JDut1","provides JD_ut1 time for --radec2000toNow test (default: current time)", false, 0,"double");		cmd.add(JDut1);
		
		ValueArg<double> kSZbetall("", "kSZbetall","kSZ effect option: "
				"LOS cluster velocity (positive for receeding clusters) in units of c (default: 1/300)", false, 1.0/300,"double");
		cmd.add(kSZbetall);
		
		//direction conversions
		/* ValueArg<double> ra("", "ra", "right ascension ( default: 0.0 )", false,0.0,"double");	cmd.add( ra ); */
		/* ValueArg<double> dec("", "dec", "declination ( default: 0.0 )", false,0.0,"double");	cmd.add( dec ); */

		//
		// Parse the command line.
		//
		cmd.parse(argc, argv);
		
		//
		// Set variables
		//
		
		/* 	vector<string>	v4; */
		/* 	v4 = cmf.getValue(); if (v4.size() > 0) { _cmf_size = (long)v4.size(); for (i=0;i<(long)v4.size();i++) { _cmf[i] = v4[i]; } } */
		/* 	_cmf_U = cmf_U.getValue(); */
		/* 	_cmf_r = cmf_r.getValue(); */
		/* 	_cmf_r0 = cmf_r0.getValue(); */

		_dp = dp.getValue();
		_lbtd = lbtd.getValue();
		_X = X.getValue();
		_Xrec = Xrec.getValue();
		_dA = dA.getValue();
		_dL = dL.getValue();
		_dac_rec = dac_rec.getValue();
		_dac = dac.getValue();
		
		_agez = agez.getValue();
		_agea = agea.getValue();
		
		_Wtot0 = Wtot0.getValue();
		_a0 = a0.getValue();
		_Wtot = Wtot.getValue();
		_Wk = Wk.getValue();
		_Rc = Rc.getValue();
		_lac = lac.getValue();
		_thac = thac.getValue();
		_cs = cs.getValue();
		_ang = ang.getValue();
		_Tcmb = Tcmb.getValue();
		_Tcmb0 = Tcmb0.getValue();
		_cmbfMax = cmbfMax.getValue();
		_rhoC = rhoC.getValue();
		_size = size.getValue();
		_angdA = angdA.getValue();
		_dAvsz = dAvsz.getValue();
		_Xvsz = Xvsz.getValue();
		_lbtdvsz = lbtdvsz.getValue();
		_vsz = vsz.getValue();
		_vsa = vsa.getValue();
		_vst = vst.getValue();
		_visibility = visibility.getValue();
		_Wvsz = Wvsz.getValue();
		
		_boxesInField = boxesInField.getValue();
		_lightCone = lightCone.getValue();
		_cells = cells.getValue();
		_Lmin = Lmin.getValue();
		
		_q = q.getValue();
		_a = a.getValue();
		_Wr0 = Wr0.getValue();
		_Wb0 = Wb0.getValue();
		_Wcdm0 = Wcdm0.getValue();
		if (Wm0.isSet()) {
			_Wm0 = Wm0.getValue();
			_Wcdm0 = _Wm0 - _Wb0;
		} else _Wm0 = _Wb0 + _Wcdm0;
		if (Wm0.isSet() and Wl0.isSet()) _Wtot0 = _Wm0 + _Wl0 + _Wr0;
		_Wl0 = Wl0.getValue();
		_Wl0fit2flat = Wl0fit2flat.getValue();
		if (_Wl0fit2flat) {
			_Wl0 = 1.0 - _Wm0 - _Wr0;
		}
		
		_Wr = Wr.getValue();
		_Wm = Wm.getValue();
		_Wb = Wb.getValue();
		_Wl = Wl.getValue();
		_w0 = w0.getValue();
		_h = h.getValue();
		if (!Wr0.isSet()) {
			_Wr0 /= _h * _h;
		}
		_z = z.getValue();
		_zmax = zmax.getValue();
		_zmin = zmin.getValue();
		
		_rmin = rmin.getValue();
		_rmax = rmax.getValue();
		_dr = dr.getValue();
		
		_fmin = fmin.getValue();
		_fmax = fmax.getValue();
		_fN = fN.getValue();
		
		_NP = NP.getValue();
		
		_notfnorm = notfnorm.getValue();
		_BBKS_transferFunction = BBKS_transferFunction.getValue();
		_sigma0 = sigma0.getValue();
		_sigmaR = sigmaR.getValue();
		_s0scale = s0scale.getValue();
		_tf = tf.getValue();
		_tfcol = tfcol.getValue();
		_A = PkA.getValue();
		_ns = Pkns.getValue();
		_k0 = Pkk0.getValue();
		
		_outputFile = output_file.getValue();
		_verbosity = verbosity.getValue();
		_script_mode = script_mode.getValue();
		
		_K2KeV = K2KeV.getValue();
		_y0 = y0.getValue();
		_mueinv1 = mueinv1.getValue();
		
		_refraction = refraction.getValue();
		_refrfn = refrfn.getValue();
		_refrMoonAtDusk = refrMoonAtDusk.getValue();
		_testCoor = testCoor.getValue();
		_testVtotRT32 = testVtotRT32.getValue();
		
		_lbdeg2dms = lbdeg2dms.getValue();
		_printGadgetUnits = printGadgetUnits.getValue();
		_a0H0inv = a0H0inv.getValue();
		_WkRTopHat = WkRTopHat.getValue();
		_partMass = partMass.getValue();
		_TinkerDeltaRho = TinkerDeltaRho.getValue();
		
		_dontSaveRun = dontSaveRun.getValue();
		
		_TSZfreqFn = false;
		_TSZfreq = TSZfreq.getValue();
		if (fmin.isSet() and fmax.isSet()) _TSZfreqFn = true;
//		_kSZfreqFn = false;
		
		_gravPot = gravPot.getValue();
		_gravSoft = gravSoft.getValue();
		
		_freq = freq.getValue();
		_RTres = RTres.getValue();
		_taper = taper.getValue();
		_l2nu = l2nu.getValue();
		_nu2l = nu2l.getValue();
		_kSZbetall=kSZbetall.getValue();
		
		
		_JD=JD.getValue();
		_JDut1=JDut1.getValue();
		_JD2cal=JD2cal.getValue();
		_u2cal=u2cal.getValue();
		if (u2cal.isSet()) _dou2cal=true; else _dou2cal=false;
		_doJD2cal=JD2cal.isSet();
		_cal2JD=cal2JD.getValue();
		
		_radec2000toGal=radec2000toGal.getValue();
		_radec2000toNow=radec2000toNow.getValue();
		_nutate=nutate.getValue();
		_aberrate=aberrate.getValue();
		_inCoordFmt=inCoordFmt.getValue();
		
		_Efact=Efact.getValue();
		
		
		/* 	_save_overplot_as = save_overplot_as.getValue(); */
	} catch (ArgException& e) {
		cout << "ERROR: " << e.error() << " " << e.argId() << endl;
	}
}

void say(cpedsMsgs& msgs, string msg, string key, string value) {
	if (_verbosity == 0) {
		printf("%s : %s\n", key.c_str(), value.c_str());
	} else {
		msgs.say(msg, High);
	}
}
/***************************************************************************************/
void doRefraction(cpedsMsgs& msgs) {
	msgs.say("Calculating refraction", High);
	
	double T, P, H, alt, lat, lambda, Tlapse, acc, ZDobs;
	QString str = _refraction.c_str();
	QStringList qsl = str.split(",");
	if (qsl.size() != 9) {
		msgs.criticalError("Wrong number of parameters given. see help", High);
	}
	for (long i = 0; i < qsl.size(); i++) {
		switch (i) {
			case 0:
				ZDobs = qsl[i].toDouble();
				break;
			case 1:
				T = qsl[i].toDouble();
				break;
			case 2:
				P = qsl[i].toDouble();
				break;
			case 3:
				H = qsl[i].toDouble();
				break;
			case 4:
				lambda = qsl[i].toDouble();
				break;
			case 5:
				alt = qsl[i].toDouble();
				break;
			case 6:
				lat = qsl[i].toDouble();
				break;
			case 7:
				Tlapse = qsl[i].toDouble();
				break;
			case 8:
				acc = qsl[i].toDouble();
				break;
				
			default:
				break;
		}
	}
	
	double ZDspace = cpeds_refraction(ZDobs, alt, T, P, H, lambda, lat, Tlapse,
			acc);
	msgs.say("Refraction parameters:", High);
	msgs.say("T [C]: %lf", T, Medium);
	msgs.say("P [mbar]: %lf", P, Medium);
	msgs.say("H [\%]: %lf", H, Medium);
	msgs.say("lambda [cm]: %lf", lambda, Medium);
	msgs.say("alt [m]: %lf", alt, Medium);
	msgs.say("latitude [deg]: %lf", lat, Medium);
	msgs.say("Tlapse [K/m]: %lf", Tlapse, Medium);
	msgs.say("acc: %lE", acc, Medium);
	msgs.say("ZDspace [deg]: %.6lf, ZDobs [deg]: %.6lf, ref [deg]: %lE",
			ZDspace, ZDobs, ZDspace - ZDobs, Low);
	
	if (_refrfn) {
		mscsFunction ref, refKBpos, refKBneg, diff;
		double h;
		for (long zd = 0; zd <= 90; zd++) {
			ZDobs = zd;
			ref.newPoint(ZDobs,
					cpeds_refraction(ZDobs, alt, T, P, H, lambda, lat, Tlapse,
							acc) - ZDobs);
			h = double(90) - ZDobs;
			refKBpos.newPoint(ZDobs,
					(double(90) - cpeds_refractionKB(h * PI180, 1) * PI180inv)
							- ZDobs);
			refKBneg.newPoint(ZDobs,
					(double(90) - cpeds_refractionKB(h * PI180, -1) * PI180inv)
							- ZDobs);
		}
		ref.save("refr-sla.ZD");
		refKBpos.save("refrKBpos.ZD");
		refKBneg.save("refrKBneg.ZD");
		diff = ref - refKBpos.absoluteValue();
		diff *= double(1000);
		string fname = str.replace(',', '_').toStdString();
		diff.save("diff--" + fname + "_less_RT32KBneg.mdeg");
	}
	
}
/***************************************************************************************/
void doMoonRefractionAtDusk(cpedsMsgs& msgs) {
	msgs.say("Simulating Moon refraction", High);
	// collect run parameters
	
	
	double T, P, H, alt, lat, lambda, Tlapse, acc, ZDobs_start, ZDobs_end, ZDobs_step, MoonSize;
	double ZDobs, ZDspaceSet,ZDobsSet;
	double minLon,maxLon,lowEdgeElevObs,highEdgeElevObs;
	double lowEdgeElevTrue,highEdgeElevTrue;
	double extentHorizObs,extentHorizTrue, extentVertObs, extentVertTrue;

	QString str = _refrMoonAtDusk.c_str();
	QStringList qsl = str.split(",");
	if (qsl.size() != 12) {
		printf("Nargs: %li\n",qsl.size());
		msgs.criticalError("Wrong number of parameters given. see help", High);
	}
	for (long i = 0; i < qsl.size(); i++) {
		switch (i) {
			case 0:
				ZDobs_start = qsl[i].toDouble();
				break;
			case 1:
				ZDobs_end = qsl[i].toDouble();
				break;
			case 2:
				ZDobs_step= qsl[i].toDouble();
				break;
			case 3:
				MoonSize = qsl[i].toDouble();
				break;
			case 4:
				T = qsl[i].toDouble();
				break;
			case 5:
				P = qsl[i].toDouble();
				break;
			case 6:
				H = qsl[i].toDouble();
				break;
			case 7:
				lambda = qsl[i].toDouble();
				break;
			case 8:
				alt = qsl[i].toDouble();
				break;
			case 9:
				lat = qsl[i].toDouble();
				break;
			case 10:
				Tlapse = qsl[i].toDouble();
				break;
			case 11:
				acc = qsl[i].toDouble();
				break;
				
			default:
				break;
		}
	}

	msgs.say("Caclulation parameters:", High);
	msgs.say("ZDobs start [deg]: %lf", ZDobs_start, Medium);
	msgs.say("ZDobs end [deg]: %lf", ZDobs_end, Medium);
	msgs.say("ZDobs step [deg]: %lf", ZDobs_step, Medium);
	msgs.say("Moon angular diameter [deg]: %lf", MoonSize, Medium);
	msgs.say("T [C]: %lf", T, Medium);
	msgs.say("P [mbar]: %lf", P, Medium);
	msgs.say("H [\%]: %lf", H, Medium);
	msgs.say("lambda [cm]: %lf", lambda, Medium);
	msgs.say("alt [m]: %lf", alt, Medium);
	msgs.say("latitude [deg]: %lf", lat, Medium);
	msgs.say("Tlapse [K/m]: %lf", Tlapse, Medium);
	msgs.say("acc: %lE", acc, Medium);

	double ZDspace = cpeds_refraction(ZDobs_start, alt, T, P, H, lambda, lat, Tlapse,acc);
//	msgs.say("ZDspace [deg]: %.6lf, ZDobs [deg]: %.6lf, ref [deg]: %lE",
//			ZDspace, ZDobs_start, ZDspace - ZDobs_start, Low);

	//
	// calculate true2observed ZD conversion function
	//
	
	mscsFunction obsZD; // x - space ZD [deg] y - observed ZD [deg]
	double h;
	ZDobs=0;
	while (ZDobs <=95) {
		obsZD.newPoint(ZDobs, cpeds_refraction(ZDobs, alt, T, P, H, lambda, lat, Tlapse, acc) );
//		h = double(90) - ZDobs;
//		refKBpos.newPoint(ZDobs,
//				(double(90) - cpeds_refractionKB(h * PI180, 1) * PI180inv)
//						- ZDobs);
//		refKBneg.newPoint(ZDobs,
//				(double(90) - cpeds_refractionKB(h * PI180, -1) * PI180inv)
//						- ZDobs);
		ZDobs+=0.1;
	}
	obsZD.invert();
	obsZD.save("spaceZD-obsZD");
	
	//
	// generate moon at zenith
	//
	cpedsDirectionSet dsMoonTrue, dsMoonObs, dsMoonSpaceIni;
	long Mpts=1000;
	double t;
	for (long i = 0; i < Mpts; i++) {
		cpedsDirection d(twoPI/Mpts*i,(90.0-MoonSize/2)*PI180); // "l,b" coords
		dsMoonTrue.append(d);
	}

	// refract moon at zenith -- no need I think

	
	dsMoonSpaceIni=dsMoonTrue;

	ZDspaceSet=ZDobs_start;
	while (ZDspaceSet<ZDobs_end) {
		
		
		// rotate space Moon to initial ZDspace
		for (long i = 0; i < Mpts; i++) {
			cpedsDirection d=dsMoonSpaceIni[i];
			d.Rx(ZDspaceSet*PI180);
			d*=PI180inv; // convert to deg
			dsMoonTrue[i]=d;
		}
		dsMoonObs=dsMoonTrue;
		
		
		// refract Moon at initial ZDobs (this will cause that the initial ZD will be slightly different from the requested one)
		for (long i = 0; i < Mpts; i++) {
			cpedsDirection d=dsMoonTrue[i];
			ZDspace=90.0-d.lat();
			ZDobs = obsZD.finter(ZDspace,"linear");
			dsMoonObs[i].setLat(90.0-ZDobs); // set observed elevation		
		}
		
		// save Moon space
		dsMoonTrue.save("MoonSpace");
		dsMoonTrue.getRanges(minLon,maxLon,lowEdgeElevTrue,highEdgeElevTrue);
		
		// save Moon obs
		dsMoonObs.save("MoonObs");
		dsMoonObs.getRanges(minLon,maxLon,lowEdgeElevObs,highEdgeElevObs);
		extentHorizObs=maxLon-minLon;
		extentVertObs=highEdgeElevObs-lowEdgeElevObs;
		printf("minLon obs [deg]: %lf\n",minLon);
		printf("maxLon obs [deg]: %lf\n",maxLon);
		printf("horizontal extent obs [deg]: %lf\n",extentHorizObs*cos(ZDobs*PI180));
		printf("vertical extent obs [deg]: %lf\n",extentVertObs);
		printf("horiz/vert obs: %lf\n",extentHorizObs/extentVertObs);

		// save Moon without atmosphere at elevation of the observed Moon in case with the atmosphere
		ZDobsSet=obsZD.finter(ZDspaceSet,"linear");
		for (long i = 0; i < Mpts; i++) {
			cpedsDirection d=dsMoonSpaceIni[i];
			d.Rx(ZDobsSet*PI180);
			d*=PI180inv; // convert to deg
			dsMoonTrue[i]=d; 
		}
		dsMoonTrue.save("MoonSpaceShifted");
		printf("elev obs [deg]: %lf\n",90.0-ZDobsSet);
		printf("elev true [deg]: %lf\n",90.0-ZDspaceSet);

		// measure distances
		dsMoonTrue.getRanges(minLon,maxLon,lowEdgeElevTrue,highEdgeElevTrue);
		extentHorizTrue=maxLon-minLon;
		extentVertTrue=highEdgeElevTrue-lowEdgeElevTrue;
		printf("horizontal extent true [deg]: %lf\n",extentHorizTrue);
		printf("vertical extent true [deg]: %lf\n",extentVertTrue);
		printf("horiz/vert true: %lf\n",extentHorizTrue/extentVertTrue);
		printf("dist to lower edge from centre: %lf\n",90.0-ZDobsSet-lowEdgeElevObs);
		printf("dist to upper edge from centre: %lf\n",highEdgeElevObs-(90.0-ZDobsSet));
		
		
		exit(0);
	}
	
//	//
//	// calculate Moon deformation due to refraction
//	//
//	
//	ZDobs=ZDobs_start;
//		mscsFunction ref, refKBpos, refKBneg, diff;
//		double h;
//		for (long zd = 0; zd < 90; zd++) {
//			ZDobs_start = zd;
//			ref.newPoint(ZDobs_start,
//					cpeds_refraction(ZDobs_start, alt, T, P, H, lambda, lat, Tlapse,
//							acc) - ZDobs_start);
//			h = double(90) - ZDobs_start;
//			refKBpos.newPoint(ZDobs_start,
//					(double(90) - cpeds_refractionKB(h * PI180, 1) * PI180inv)
//							- ZDobs_start);
//			refKBneg.newPoint(ZDobs_start,
//					(double(90) - cpeds_refractionKB(h * PI180, -1) * PI180inv)
//							- ZDobs_start);
//		}
//		ref.save("refr-sla.ZD");
//		refKBpos.save("refrKBpos.ZD");
//		refKBneg.save("refrKBneg.ZD");
//		diff = ref - refKBpos.absoluteValue();
//		diff *= double(1000);
//		string fname = str.replace(',', '_').toStdString();
//		diff.save("diff--" + fname + "_less_RT32KBneg.mdeg");
	
}

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
/*!
 \brief this is a function returning the corrections in A,h for requested A,h which are used in the control system meteofix
 \details 
 @param
 @return

 \date Jun 5, 2013, 6:57:13 PM
 \author Bartosz Lew
 */
double arcsin(double x) {
	if (x > 1.) x = 1;
	if (x < -1.) x = -1.;
	return asin(x);
}

void Model4(double AZ, double ZD, double *dAZ, double *dZD) {
	
	//const double p[13] = {
	//-0.032367, -0.049374, 0.051222, -0.054595, -0.000161, 0.002850,
	//0.001370, 0.020659, -0.012753, 0.008571, 0.004812, -0.006485, 0.028568 };
	
	const double p[13] = { -3.2367E-02, -4.9374E-02, 5.1222E-02, -5.4595E-02,
			-1.6113E-04, 2.8503E-03, 1.3701E-03, 2.0659E-02, -1.2753E-02,
			8.5710E-03, 4.8125E-03, -6.4854E-03, 2.8568E-02 };
	
	double pi180 = M_PI / 180.;
	
	AZ = AZ * pi180;
	double Alt = (90. - ZD) * M_PI / 180.;
	
	double xi = p[4] * pi180;
	double zeta = p[5] * pi180;
	double sigma = p[2] * pi180;
	double beta = p[3] * pi180;
	double sh = sin(Alt);
	double ch = cos(Alt);
	
	double sinxi = sin(xi);
	double sinzeta = sin(zeta);
	double se = sqrt(sinxi * sinxi + sinzeta * sinzeta);
	
	if (ch < 0.) {
		se = -se;
	}
	
	double ce = sqrt(1. - se * se);
	double alfa = atan2(sin(zeta), sin(xi));
	double AT = atan2(ch * sin(alfa - AZ), sh * se - ch * ce * cos(alfa - AZ))
			- atan2(sin(alfa), -ce * cos(alfa));
	
	double hT = asin(ce * sh + se * ch * cos(alfa - AZ));
	double AT_Az = AT - AZ;
	if (ch < 0.) {
		AT_Az = M_PI - AT_Az;
		hT = M_PI - hT;
	}
	
	*dAZ = arcsin((sin(sigma) * sin(hT) + sin(beta)) / (cos(hT) * cos(sigma)));
	
	double hb = atan2(sin(hT) * cos(sigma) + cos(hT) * sin(sigma) * sin(*dAZ),
			cos(hT) * cos(*dAZ));
	
	*dAZ = fmod((*dAZ + AT_Az) * 180. / M_PI, 360.) + p[0] + p[8] * sin(2. * AZ)
			+ p[9] * cos(2. * AZ) + p[11] * sin(3. * AZ) * sin(Alt)
			+ p[12] * cos(Alt) * cos(AZ / 4.);
	
	*dZD = (hb - Alt) * 180. / M_PI + p[1] + p[6] * ch + p[7] * sh
			+ p[10] * sin(2. * AZ);
	*dZD = -(*dZD);
}
void doModel4() {
	double aa, hh;
	double da, dzd;
	matrix<double> mdA(90, 360);
	matrix<double> mdZD(90, 360);
	for (long A = -180; A < 180; A++) {
		for (long h = 0; h < 90; h++) {
			aa = A;
			hh = h;
			Model4(aa, 90.0 - hh, &da, &dzd);
			mdA(h, A + 180) = da;
			mdZD(h, A + 180) = dzd;
		}
		
	}
	cpeds_matrix_save(mdA, "mdA");
	cpeds_matrix_save(mdZD, "mdZD");
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
void caclulate_sigmaR(cpedsMsgs msgs) {
	double r = _rmin;
	double logr = log10(r);
	mscsFunction sigmaR("sigmaR");
	mscsFunction sigmaM("sigmaM");
	mscsFunction dlnsMinvdM("dlnsMinvdM");
	cpedsFunctionCosmo tinkerFact, tinker("tinker_mass_fn"), tinkerCuml;
	cpedsFunctionCosmo PSfact, PS("PS_mass_fn"), PScuml;
	
	//
	// load transfer function
	//
	
	mscsFunction tf("tf");
	tf.load(_tf, true, 0, _tfcol);
	if (!_notfnorm) tf /= tf.Y(0); // normalize to one at large scales
	// convert argumets from units h/Mpc to Mpc^-1
	tf.scaleX(_h);
	
	cpedsMsgs msgsLoc(msgs);
	msgsLoc.setVerbosity(Zero);
	
	double M, rhom = cpeds_rhoC0(1.0) * _Wm0; // kg/m^3 h^2
	while (r < _rmax) {
		sigmaR.newPoint(r, caclulate_sigma0(r, msgsLoc, &tf));
		M = 4 * PI / 3 * r * r * r * rhom * CPEDS_MPC * CPEDS_MPC * CPEDS_MPC; // kg h^-1
		sigmaM.newPoint(M, caclulate_sigma0(r, msgsLoc, &tf));
		msgs.say("R [Mpc h^-1]: %lf", r, Low);
		logr += _dr;
		r = exp10(logr);
//		printf("R: %lE\n",r);
	}
	double deltaRho = _TinkerDeltaRho;
	tinkerFact = tinker.make_tinker_mass_function(sigmaM, _Wm0, deltaRho);
	tinkerCuml = tinker.integrateX_Xmax();
	PSfact = PS.make_PS_mass_function(sigmaM, _Wm0);
	PScuml = PS.integrateX_Xmax();
	
	sigmaM.scaleX(1.0 / (CPEDS_SOLAR_MASS));
	
	string outf;
	string outfTinker;
//	QFileInfo qfi(_tf.c_str());
//	outf=qfi.baseName().toStdString();
	outf = "mfn";
	outfTinker = "-Drho_" + msgs.toStrf(deltaRho, 0);
	
	sigmaR.save(outf + ".sigmaR");
	sigmaM.save(outf + ".sigmaM");
	
	tinker.save(outf + outfTinker + ".tinker");
	tinker.logX(10);
	tinker.logY(10);
	tinker.save(outf + outfTinker + ".tinkerLog10");
	tinkerCuml.save(outf + outfTinker + ".tinkerCuml");
	tinkerFact.save(outf + outfTinker + ".tinkerFact");
	tinkerFact.logX(10);
	tinkerFact.logY(10);
	tinkerFact.save(outf + outfTinker + ".tinkerFactlog10");
	
	PS.save(outf + ".PS");
	PS.logX(10);
	PS.logY(10);
	PS.save(outf + ".PSlog10");
	PScuml.save(outf + ".PScuml");
	PSfact.save(outf + ".PSfact");
	PSfact.logX(10);
	PSfact.logY(10);
	PSfact.save(outf + ".PSfactlog10");
	
	msgs.setLogFileName(outf + ".sigmaR.params");
	msgs.loggingOn();
	msgs.say("Cosmological parameters used for this sigma(R) relation", High);
	msgs.say("Wb0: %lf", _Wb0, Low);
	msgs.say("Wcdm0: %lf", _Wcdm0, Low);
	msgs.say("Wl0: %lf", _Wl0, Low);
	msgs.say("h: %lf", _h, Low);
	msgs.say("Delta_R^2: %lE", _A, Low);
	msgs.say("k0 [Mpc^-1]: %lf", _k0, Low);
	msgs.say("ns: %lf", _ns, Low);
	msgs.say("transfer function file: " + _tf, Low);
	msgs.say("transfer function column: %li", _tfcol, Low);
	msgs.say("window function top-hat (if false then gaussian is used): ",
			_WkRTopHat, Low);
	msgs.say("Rmin [Mpc h^-1]: %lf", _rmin, Low);
	msgs.say("Rmax [Mpc h^-1]: %lf", _rmax, Low);
	msgs.say("dr [Mpc h^-1]: %lf", _dr, Low);
	msgs.say("Files", High);
	msgs.say(".sigmaR - sigma(R) relation with R[Mpc h^-1]", Medium);
	msgs.say(
			".sigmaM - sigma(M) relation with M[Msol h^-1] where M is the average mass within sphere of radius R",
			Medium);
	msgs.say(".tinker - tinker mass function dn/dM [Mpc^-3 h^3 Msol^-1 h^1]",
			Medium);
	msgs.say(
			".tinkerCuml - cumulative tinker mass function int_Mmin^M dn/dM [Mpc^-3 h^3]",
			Medium);
	msgs.say(
			".tinkerFact - tinker mass function dn/dM * M^2/rho_m (consistent with Tinker paper function)",
			Medium);
	msgs.say(".tinkerLog10 - tinker mass function log10(dn/dM) vs log10(M) ",
			Medium);
	msgs.say(
			".tinkerFactlog10 - tinker mass function log10(dn/dM * M^2/rho_m) vs log10(M) (consistent with Tinker paper function)",
			Medium);
	msgs.say(
			".PS - Press-Schechter mass function dn/dM [Mpc^-3 h^3 Msol^-1 h^1]",
			Medium);
	msgs.say(
			".PSlog10 - Press-Schechter mass function log10(dn/dM) vs log10(M) ",
			Medium);
	msgs.say(
			".PScuml - cumulative Press-Schechter mass function int_Mmin^M dn/dM [Mpc^-3 h^3]",
			Medium);
	msgs.say(
			".PSfact - Press-Schechter mass function dn/dM * M^2/rho_m (consistent with Tinker paper function)",
			Medium);
	msgs.say(
			".PSfactlog10 - Press-Schechter mass function log10(dn/dM * M^2/rho_m) vs log10(M) (consistent with Tinker paper function)",
			Medium);
	msgs.say(
			".Pk - matter power spectrum for the provided transfer function: P(k_h) [h^3/Mpc^3], k_h [h/Mpc]",
			Medium);
	
	msgs.loggingOff();
}

/***************************************************************************************/
double caclulate_sigma0(double R, cpedsMsgs& msgs,
		mscsFunction *transfer_function) {
	mscsFunction tmp("tmp", msgs.getVerbosity());
	mscsFunction Pk("Pk", msgs.getVerbosity());
	mscsWindowFunction WkR("gauss", msgs.getVerbosity());
	
	cpedsFunctionCosmo Pk_h; // Pk in units of (h/Mpc)^3, k: [h/Mpc]
	
	// tf^2 part
	mscsFunction tf("tf");
	if (transfer_function == NULL) {
		tf.load(_tf, true, 0, _tfcol);
		if (!_notfnorm) tf /= tf.Y(0); // normalize to one at large scales
		tf.save("cosmocalc.sigma0.tf");
		Pk_h.make_Pk(tf, _A, _k0, _ns, _h);
		Pk_h.save("cosmocalc.sigma0.Pk");
		// convert argumets from units h/Mpc to Mpc^-1
		tf.scaleX(_h);
		//	  mscsFunction tmp2=tf.interpolate(tf.getx(0),false,"cspline");
		//	  tmp2.save("tf.inter.out");
	} else {
		tf = (*transfer_function);
	}
	
	// 4/25 * (k/a0/H0)^4/k = 4/25 * k^3/(a0*H0)^4 part    for CDM model
	double* k = tf.extractArguments();
	tmp.importFunction(k, k, tf.pointsCount());
	delete k; //[Mpc]
	Pk = tmp;
	WkR = tmp;
	tmp = tmp * tmp * tmp;
	tmp *= double(4.0 / 25.0);
	tmp /= double(pow(_a0 * _h * 100000.0 / CPEDS_c, 4.0)); // H [m/s/Mpc]
			
	if (_z != 0) {
		msgs.warning("sigma(R) calculation was not tested on z!=0. Be warned.",
				Top);
	}
	// LCDM model linear growth correction factor
	double Wtot = cpeds_Wtot(_Wr0, _Wm0, _Wl0, _z, _w0);
	if (fabs(Wtot - 1.0) < 0.01) { // for flat case use approximation formula
		if (_Wl0 != 0) {
			double lcdmFact, g, Wm = cpeds_Wm(_Wr0, _Wm0, _Wl0, _z, _w0), Wl =
					cpeds_Wl(_Wr0, _Wm0, _Wl0, _z, _w0);
			msgs.say("Wm: %lE", Wm, Medium);
			msgs.say("Wl: %lE", Wl, Medium);
			g = 5.0 / 2.0 * Wm
					/ (pow(Wm, 4.0 / 7.0) - Wl
							+ (1.0 + Wm / 2) * (1.0 + Wl / 70));
			msgs.say("g: %lE", g, Medium);
			//			printf("(): %lE\n",( 1.0/70.0 + 209.0*W/140 - W*W/140 + pow(W,4.0/7.0)));
			//			printf("(): %lE\n",( 1.0/70.0 + 209.0*W/140 - W*W/140 ));
			//			printf("(): %lE\n",pow(W,4.0/7.0));
			lcdmFact = pow(g / Wm, 2);
			msgs.say(
					"LCDM wrt standard CDM model sigma8 boost factor (g/W)^2 (Carrol 1992 p.515, eq.29): %lE",
					lcdmFact, Medium);
			tmp *= lcdmFact;
		} else {
			msgs.say(
					"no need for extra growth suppresion factors in case of flat model without cosmological constant",
					Medium);
		}
	} else { // non-flat cases require integration
		msgs.warning("NON-FLAT MODELS ARE NOT CORRECTLY IMPLEMENTED YET.", Top);
		double lcdmFact, g, Wm = cpeds_Wm(_Wr0, _Wm0, _Wl0, _z, _w0), Wl =
				cpeds_Wl(_Wr0, _Wm0, _Wl0, _z, _w0);
		msgs.say("Wm: %lE", Wm, Medium);
		msgs.say("Wl: %lE", Wl, Medium);
		
		g = 5.0 / 2.0 * Wm / (1.0 + Wm / 2 + pow(Wm, 4.0 / 7.0));
		msgs.say("g: %lE", g, Medium);
		
		/*
		 //Aug 2, 2013 12:58:45 AM - DEBUG BEGIN accorgin to eq.6.10
		 
		 cpedsFunctionCosmo aH;
		 double ast=1e-4;
		 double aen=1;
		 double da=ast;
		 double a=ast;
		 double z;
		 while (a<aen) {
		 z=1.0/a-1;
		 aH.newPoint(a,pow(a*cpeds_hubble(_Wr0,_Wm0,_Wl0,z,_w0,_h),-3));
		 a+=da;
		 }
		 g=5.0/2.0 * Wm*Wm* pow(100*_h,3)*aH.integrate();
		 
		 //DEBUG - END
		 */

		/*
		 //Aug 2, 2013 12:59:13 AM - DEBUG BEGIN according to eq.29 and 9 from Carroll 1992
		 
		 cpedsFunctionCosmo dadtau;
		 double ast=1e-4;
		 double aen=1;
		 double da=ast;
		 double a=ast;
		 while (a<aen) {
		 dadtau.newPoint(a,1.0+Wm*(1.0/a-1)+Wl*(a*a-1));
		 a+=da;
		 }
		 dadtau.power(-1.5);
		 g=5.0/2.0 * Wm * dadtau.integrate();
		 printf("g: %lE\n",g);
		 
		 //DEBUG - END
		 */
//			
		lcdmFact = pow(g / Wm, 2);
		msgs.say(
				"non-flat LCDM wrt standard CDM model sigma8 boost factor (g/Wm)^2 (Carrol 1992 p.515, eq.29): %lE\n",
				lcdmFact, Medium);
		
		tmp *= lcdmFact;
	}
//		// flat LCDM model correction factor
//		if (_Wl0!=0) {
//			double lcdmFact,g,W=cpeds_Wm(_Wr0,_Wm0,_Wl0,_z,_w0);
//			printf("Wm: %lE\n",W);
//			g=5.0/2.0*W / ( 1.0/70.0 + 209.0*W/140 - W*W/140 + pow(W,4.0/7.0));
//			printf("g: %lE\n",g);
////			printf("(): %lE\n",( 1.0/70.0 + 209.0*W/140 - W*W/140 + pow(W,4.0/7.0)));
////			printf("(): %lE\n",( 1.0/70.0 + 209.0*W/140 - W*W/140 ));
////			printf("(): %lE\n",pow(W,4.0/7.0));
//			lcdmFact=pow(g/W,2);
//			printf("LCDM wrt standard CDM model sigma8 boost factor (g/Wtot)^2 (Carrol 1992, also in textbook Liddle & Lyth p.144): %lE\n",lcdmFact);
//			tmp*=lcdmFact;
//			if (_Wl0<0) {
//				printf("WARNING not sure that this formula supports negative Wl\n");
//			}
//		}
//		else {
//			double Wtot=cpeds_Wtot(_Wr0,_Wm0,_Wl0,_z,_w0);
//			if (fabs(Wtot-1.0)>0.01) { // assume open CDM model
//				double ocdmFact,g,W=cpeds_Wm(_Wr0,_Wm0,_Wl0,_z,_w0);
//				printf("Wm: %lE\n",W);
//				g=5.0/2.0*W / ( 1.0 + W/2 + pow(W,4.0/7.0));
//				printf("g: %lE\n",g);
//				ocdmFact=pow(g/W,2);
//				printf("OCDM wrt standard CDM model sigma8 boost factor (g/Wtot)^2 (Carrol 1992, also in textbook Liddle & Lyth p.148): %lE\n",ocdmFact);
//				tmp*=ocdmFact;
//			}
//		}
	
	// primordial power spectrum part
	Pk /= double(_k0);
	Pk = Pk.power(_ns - 1);
	Pk *= _A;
	
	//
	// smoothing kernel part
	//
// this works good
	if (_WkRTopHat) WkR.mkTopHatKernel(R / _h, 1);
	else WkR.mkGaussianKernel(R / _h);
	
	/* sigma0^2(R) = int P_g(k) * W(kR)^2 dk/k
	 * where 
	 * P_g(k) = 4/25 * (k/a0/H0)^4 T^2(k) P_R(k)
	 * where P_R(k) - is the power spectrum of the co-moving curvature perturbation
	 */
	tmp = tmp * Pk * WkR * WkR * tf * tf;
	tmp.interpolate(tmp.getx(0), true, "linear");
	if (transfer_function == NULL) tmp.save("cosmocalc.sigma0.intfn");
	double sigma0 = sqrt(tmp.integrate());
	
	msgs.say(
			"sigma0 for the scale %lE [Mpc h^-1] (for the age of the provided transfer function) is: %lE",
			R, sigma0, Medium);
	msgs.say("assumed Delta_R^2(k): %lE\n", _A, Medium);
	msgs.say("assumed k0 [Mpc^-1]: %lE\n", _k0, Medium);
	msgs.say("assumed ns: %lE\n", _ns, Medium);
	
	return sigma0;
}
/***************************************************************************************/
double getSZfreqFact(double f) {
	double x = CPEDS_h * f * 1.0e9 / (CPEDS_kB * _Tcmb0);
	return x * cosh(x / 2) / sinh(x / 2) - 4;
}
/***************************************************************************************/
//double getSZintensityFact(double f) {
//	double x = CPEDS_h * f * 1.0e9 / (CPEDS_kB * _Tcmb0);
//	return x * cosh(x / 2) / sinh(x / 2) - 4;
//}
/***************************************************************************************/
mscsFunction calculate_TSZfreqFn(double fmin, double fmax, long n) {
	mscsFunction fn;
	double df = (fmax - fmin) / n;
	double f = fmin;
	for (unsigned long i = 0; i < n; i++) {
		fn.newPoint(f, getSZfreqFact(f));
		f += df;
	}
	return fn;
}
/***************************************************************************************/
mscsFunction calculate_TSZintensityFn(double fmin, double fmax, long n) {
	mscsFunction fn;
	double df = (fmax - fmin) / n;
	double f = fmin;
	for (unsigned long i = 0; i < n; i++) {
		double x = CPEDS_h * f * 1.0e9 / (CPEDS_kB * _Tcmb0);
		double gnu = x * x * x * x * exp(x) / pow(exp(x) - 1, 2);

		fn.newPoint(f, 2.0*pow(CPEDS_kB*_Tcmb0,3)/pow(CPEDS_h*CPEDS_c,2)*getSZfreqFact(f)*gnu); // / cpeds_black_body_radiance_Bnu(f*1e9,_Tcmb0));
		f += df;
	}
	return fn;
}
/***************************************************************************************/
mscsFunction calculate_deltaToT_to_intnsity_conversionFn(double fmin, double fmax, long n) {
	mscsFunction fn;
	double df = (fmax - fmin) / n;
	double f = fmin;
	for (unsigned long i = 0; i < n; i++) {
		double x = CPEDS_h * f * 1.0e9 / (CPEDS_kB * _Tcmb0);
		double gnu = x * exp(x) / (exp(x) - 1);

		fn.newPoint(f, 1.0/cpeds_radiance_to_temp_conversionFactor(f*1e9,_Tcmb0));
		f += df;
	}
	return fn;
}

/***************************************************************************************/
void calculateGravitationalPotential(string infile) {
	cpedsPointSet3D ps;
	cpedsPointSet3D gravPot;
	cpedsPoint3D p, p2;
	double V;
	double d;
	
	if (infile == "mkUniBall") {
		cpedsRNG rns("uniform");
		rns.setMinMax(-1, 1);
		double x, y, z, M;
		long n = 0;
		while (n < _NP) {
			x = rns.getRN();
			y = rns.getRN();
			z = rns.getRN();
			M = 1;
			if (sqrt(x * x + y * y + z * z) <= 1) {
				ps.append(cpedsPoint3D(x, y, z), M);
				n++;
			}
		}
		ps.save(infile);
	} else
		if (infile == "mk1Dclub") {
			cpedsRNG rns("uniform");
			cpedsRNG rnsG("gaussian");
			rns.setMinMax(-10, 0);
			rns.setMeanVariance(0, 1);
			double x, y, z, M;
			long n = 0;
			while (n < _NP) {
				x = rns.getRN();
				y = 0;
				z = 0;
				M = 1;
				ps.append(cpedsPoint3D(x, y, z), M);
				x = rnsG.getRN();
				ps.append(cpedsPoint3D(x, y, z), M);
				n++;
			}
			ps.save(infile);
		} else {
			ps.load(infile);
		}
	
	gravPot = ps;
	double gravSoft = _gravSoft;
	
	if (_gravSoft == -1) {
		double xmin, xmax, ymin, ymax, zmin, zmax;
		ps.getRanges(xmin, xmax, ymin, ymax, zmin, zmax);
		gravSoft = sqrt(
				pow((xmax - xmin), 2) + pow((ymax - ymin), 2)
						+ pow((zmax - zmin), 2)) / _NP;
		printf("gravitational softening: %lE\n", gravSoft);
	}
	
	gravPot = ps.calculateGravitationalPotential(gravSoft);
//	for (unsigned long i = 0; i < ps.size(); i++) {
//		V=0;
//		p=ps[i];
//		for (unsigned long j = 0; j < ps.size(); j++) {
//			if (i!=j) {
////				printf("i: %li, j: %li\n",i,j);
//				p2=ps[j];
//				d=p.dist(p2);
//				V-=ps.val(j)/(d+gravSoft);
//			}
//		}
//		V*=CPEDS_G;
//		gravPot.values()[i]=V;
//	}
	
	gravPot.save(infile + ".V");
//	gravPot2.save(infile+".V2");
	ps.massCenter().print_point("mass center");
	
	mscsFunction3dregc gravPotGrid;
	gravPotGrid.calculateGravitationalPotential(ps, 0.01, gravSoft);
	gravPotGrid.saveHDF5("gravPotGrid.hdf5", "V");
	gravPotGrid.savetxtlin("gravPotGrid.txt", true);
	gravPotGrid.getMinValueCell().print_point(
			"minimum gravitational potential location");
	printf("minimal grav pot: %lE\n", gravPotGrid.getMinValue());
	
}
/***************************************************************************************/
mscsFunction calculate_beam_power_pattern(double d, double ds, double lambda,
		double taper, double thetaMax, double dTheta) {
	cpedsFunctionCosmo f;
	return f.calculateBeamPowerPattern(d, ds, lambda, taper, thetaMax, dTheta);
	
//	double th=-thetaMax;
//	double z,q,u,J1z,J2z, J1kz,J2kz;
//	double kappa=ds/d;
//	double k2=kappa*kappa;
//	double U0betaz,U0betakkzk,U;
//	double beta=1.0-pow(10,taper/20);
//	printf("beta: %lf\n",beta);
//	u=(k2 * (2.0-k2*beta))/(2.0-beta);
//	
//	mscsFunction Uth;
//	
//	while (th<thetaMax) {
//		q=sin(th*PI180)/lambda;
//		z=PI*d*q;
//		J1z=gsl_sf_bessel_J1(z);
//		J2z=gsl_sf_bessel_Jn(2,z);
//		J1kz=gsl_sf_bessel_J1(kappa*z);
//		J2kz =gsl_sf_bessel_Jn(2,kappa*z);
//		U0betaz=4.0/(2.0-beta)/z * ( (1.0-beta)*J1z + 2.0*beta*J2z/z );
//		U0betakkzk=4.0/(2.0-beta*k2)/(z*kappa) * ( (1.0-beta*k2)*J1kz + 2.0*beta*k2*J2kz/(z*kappa) );
//		
//		U=(U0betaz - u*U0betakkzk)/(1.0-u);
//		
//		Uth.newPoint(th,U*U);
//		
//		th+=dTheta;
//	}
//	return Uth;
}
/***************************************************************************************/
void doTestVtotRT32() {
	DirectionRaDec radec;
	double rah, ram, ras, decd, decm, decs;
	double JD = cpeds_julian_time(2014, 1, 5, 14.0);
	double JD0 = CPEDS_JD2000;
	rah = 10;
	ram = 10;
	ras = 10;
	decd = 10;
	decm = 10;
	decs = 10;
	
	radec.set(cpeds_HMSToAng(rah, ram, ras), cpeds_DMSToAng(decd, decm, decs));
	radec.setEpochYr(2000);
	radec *= PI180;
	radec.toJD(JD, true);
	
	double aLMST, Vsun, Vobs, Vdop, ra, dec;
	ra = radec.ra();
	dec = radec.dec();
	cpeds_vdrt32_(&JD, &ra, &dec, &aLMST, &Vsun, &Vobs, &Vdop);
	
	printf("JD: %.10lf\n", JD);
	radec.print_direction("requested direction", true);
	printf("aLMST [deg]: %.10lf\n", aLMST * PI180inv);
	printf("Vsun [km/s]: %.10lf\n", Vsun);
	printf("Vobs [km/s]: %.10lf\n", Vobs);
	printf("Vdop [km/s]: %.10lf\n", Vdop);
	
	//  czestotliwosc w pasmie spektrografu na jakiej pojawi sie linia o podanym vlsr 
	double fset = 61.0;
	
	//  syntezer RS przy spektrografie
	double lo_0 = 489.0;
	
	//  syntezer sc
	double lo_1 = 4300.0;
	
	double f0 = 6667; // MHz
	radec.set(36.7671 * PI180, 61.8728 * PI180);
	radec.setEpochYr(2000);
	double Vtot = radec.VradialLSR(cpeds_julian_time());
	double vlsr_peak = -46.4;
	lo_1 = f0 * (((Vtot - vlsr_peak) / (CPEDS_c * 1.0e-3)) + 1.0) - fset - lo_0;
	
	lo_1 /= 1000;
//	   lo_2/=lo_2/1000;
	radec.print_direction("src dir");
	printf("JD: %lE\n", cpeds_julian_time());
	printf("Vtot [km/s]: %lE\n", Vtot);
	printf("vlsr_peak [km/s]: %lE\n", vlsr_peak);
//	   printf("losl [GHz]: %lE\n",Tsys_observation_params.sllof);
	printf("losc [GHz]: %lE\n", lo_1);
	
}
/***************************************************************************************/
double DMSstrToAng(string DDMMSS,string separator) {
	QString s=DDMMSS.c_str();
	QStringList qsl=s.split(separator.c_str());
	return cpeds_DMSToAng(qsl[0].toDouble(),qsl[1].toDouble(),qsl[2].toDouble());
}
/***************************************************************************************/
double HMSstrToAng(string HHMMSS,string separator) {
	QString s=HHMMSS.c_str();
	QStringList qsl=s.split(separator.c_str());
	return cpeds_HMSToAng(qsl[0].toDouble(),qsl[1].toDouble(),qsl[2].toDouble());
}
/***************************************************************************************/
double doCal2JDconversion(string cal2JD) {
	QString str=cal2JD.c_str();
	QStringList qsl=str.split("-");
	return cpeds_julian_time(qsl[0].toInt(),qsl[1].toInt(),qsl[2].toInt(),qsl[3].toDouble()+qsl[4].toDouble()/60+qsl[5].toDouble()/3600);
}
/***************************************************************************************/
