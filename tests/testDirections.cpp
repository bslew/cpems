/*!
 * \file testDirections.cpp
 *
 *  Project: cpeds
 *  Created on: Dec 5, 2012 1:49:44 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <tclap/CmdLine.h>
#include <QtCore/QStringList>
#include <QtCore/QString>
#include "cpeds-msgs.h"
//#include "Mscs-function.h"
#include "matrix.h"
#include "cpedsMoonDirection.h"
#include "cpedsPlanetDirection.h"
#include "cpeds-direction.h"
#include "Mscs-function.h"
#include "novas.h"

//#include "rt32com"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString, _track, _planet;
bool _calcMoon, _calcPlanet, _utime2JD, _novaVSnovas, _novaReversabilityTest, _testRTlstCalculation;
bool _calcSolstice;
double _lon, _lat, _altitude, _P, _T, _JD, _da, _dt, _timeLimit, _ut;
double _ra,_dec;

void parseOptions(int argc, char** argv);
string getProgramVersionString();
void calculateMoon(cpedsMsgs& msgs);
void calculateNextSolstice(cpedsMsgs& msgs);
void calculatePlanet(string planetName, cpedsMsgs& msgs);
void novaVSnovas(cpedsMsgs& msgs);
void novaReversabilityTest(cpedsMsgs& msgs);
double RTcontrol_local_sidereal_time(double UT, double JD,double longitude );
double RTcontrol_JulianDay(int tye,int tmo,int tda);
void testRTlstCalculation();
void testConversionDeg2HMSDMS(double ra,double dec, double jdutc, cpedsMsgs& msgs);


int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".testDirections",false,"",High);
	msgs.setSaveRunWriteMode('w');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	if (_ra!=-1) { testConversionDeg2HMSDMS(_ra,_dec, 0,msgs); }
	
	if (_novaVSnovas) novaVSnovas(msgs);
	if (_novaReversabilityTest) novaReversabilityTest(msgs);
	if (_testRTlstCalculation) testRTlstCalculation();
	
	if (_utime2JD) {
//		time_t rawtime=time(0);
//		double ss=(double)rawtime+3600.0;
//		printf("%lf\n",cpeds_timeSec_to_julian_time(rawtime));
		printf("%lf\n",cpeds_timeSec_to_julian_time(_ut));
		exit(0);
	}
	
	if (_calcMoon) {		calculateMoon(msgs);			}
	if (_calcSolstice) {		calculateNextSolstice(msgs);			}
	if (_planet!="none") {		calculatePlanet(_planet, msgs);			}

	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("testDirections\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program calculates the coordinates of bodies in the Solar System, and controls the RT32 telescope to track them in the sky.\n\n. "
				"\n"
				"example usage: \n"
				"testDirections --moon - to check the current position of the moon \n"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
		ValueArg<double> lon("", "lon", "longitude of the observatory location (-180,180) [deg]", false, 18.56406,"double");	cmd.add( lon );
		ValueArg<double> lat("", "lat", "geodetic latitude of the observatory location (-90,90) [deg]", false, 53.09546,"double");	cmd.add( lat );
		ValueArg<double> altitude("", "alt", "altitude of the observatory over the ellipsoid (133) [m]", false, 133.61,"double");	cmd.add( altitude );
		ValueArg<double> pressure("P", "press", "current atmospheric pressure (1012) [mbar]", false, 1012.0,"double");	cmd.add( pressure );
		ValueArg<double> temp("T", "temp", "current atmospheric temperature (20) [Celsius]", false, 20.0,"double");	cmd.add( temp );
		ValueArg<double> JD("", "JD", "JD for which to calculate the stuff (default: 0 = current UT time) [JD]", false, 0,"double");	cmd.add( JD );
		ValueArg<double> da("", "da", "days after JD for which to calculate the stuff (default: 0 = current JD setting) [JD]", false, 0,"double");	cmd.add( da );
		ValueArg<double> dt("", "dt", "time step in calculating moon positions (default: 24 h) [h]", false, 24,"double");	cmd.add( dt );
		ValueArg<double> timeLimit("", "period", "time period over which to calculate the ephemerides (default: 30 JD) [JD]", false, 30,"double");	cmd.add( timeLimit );
		ValueArg<double> ut("", "ut", "unix UTC time to be converted to JD UTC date (default:  1.35577627082611608505E+09) ", false, 1.35577627082611608505E+09,"double");	cmd.add( ut );
//		
		SwitchArg solstice("","calcSolstice", "calculate next solstice (default: false)", false);	cmd.add( solstice );
		SwitchArg moon("","moon", "show moon direction in the sky and for the following month (default: false)", false);	cmd.add( moon );
//		SwitchArg track("","track", "send RT32 to the selected object (default: false)", false);	cmd.add( moon );
		SwitchArg utime2JD("","utime2JD", "convert unix time to JD. Use --ut option to set the unix time (default: false)", false);	cmd.add( utime2JD );
		SwitchArg novaVSnovas("","novaVSnovas", "compare coordinates conversion radec->ah between the two libraries (default: false)", false);	cmd.add( novaVSnovas );
		SwitchArg novaReversabilityTest("","novaRev", "test reversability of nova directions conversion ra,dec->ah->ra,dec (default: false)", false);	cmd.add( novaReversabilityTest );
		SwitchArg testRTlstCalculation("","testRTlst", "test continuity of lst calculation in RT control system (default: false)", false);	cmd.add( testRTlstCalculation );
		
		ValueArg<double> ra("", "ra", "right ascension for testing conversion from deg to hms [deg] [0,360) (default: -1 = not used ) ", false, -1.0,"double");	cmd.add( ra );
		ValueArg<double> dec("", "dec", "declination for testing conversion from deg to dms [deg] [-90,90] (default: -100 = not used ) ", false, -100.0,"double");	cmd.add( dec );
//		
		
//		ValueArg<string> planet("","planet","send RT32 to the selected object ",false,"","string"); cmd.add(outfile);
//		ValueArg<string> outdir("","odir","output directory name",false,"./","string"); cmd.add(outdir);
//		
		std::vector<string> allowedStr;
		allowedStr.push_back("none");
		allowedStr.push_back("Moon");
		allowedStr.push_back("Moondsc");
		allowedStr.push_back("Mercury");
		allowedStr.push_back("Venus");
		allowedStr.push_back("Mars");
		allowedStr.push_back("Jupiter");
		allowedStr.push_back("Saturn");
		allowedStr.push_back("Uranus");
		allowedStr.push_back("Neptune");
		allowedStr.push_back("Pluto");
		allowedStr.push_back("Sun");
//		
		ValuesConstraint<string>* allowedStrNew;
		allowedStrNew = new ValuesConstraint<string>(allowedStr);
//		ValueArg<string> track("", "traj", "track","send RT32 to the selected object", false,"none", allowedStrNew);	cmd.add( track );
		ValueArg<string> planet("", "planet", "planet name (default: none)", false,"none", allowedStrNew);	cmd.add( planet );
		allowedStr.clear();

		//     UnlabeledValueArg<string> input_file("input_file","file to plot","nofile","string");	cmd.add( input_file );
		//     ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		//     ValueArg<string> proj("p","proj","projection to use [mall - Mollweide/ aitoff etc] (default: mall)",false,"mall","string"); cmd.add(proj);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_lon = lon.getValue();
		_lat = lat.getValue();
		_altitude=altitude.getValue();
		_P=pressure.getValue();
		_T=temp.getValue();
		_calcMoon= moon.getValue();
		_JD = JD.getValue();
		_da=da.getValue();
		_dt=dt.getValue();
		_timeLimit=timeLimit.getValue();
		_planet=planet.getValue();
		_calcSolstice=solstice.getValue();
		_ut=ut.getValue();
		_utime2JD=utime2JD.getValue();
		_novaVSnovas=novaVSnovas.getValue();
		_testRTlstCalculation=testRTlstCalculation.getValue();
		
		
		_novaReversabilityTest=novaReversabilityTest.getValue();
		
		_ra=ra.getValue();
		_dec=dec.getValue();
//		_nLissajous=nLissajous.getValue();
//		_trackHoriz = horizTrack.getValue();
//		
//		_outfile = outfile.getValue();
//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}

string getProgramVersionString() {
	string rev;
#ifdef GENERATE_CVS_VERSION_STRING
	rev="$Revision: 1.1 $";
	QString qrev=rev.c_str();
	QStringList qrevl=qrev.split(' ');
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif

	return rev;
}



void calculateMoon(cpedsMsgs& msgs) {
	cpedsDirection obs(_lon,_lat, _altitude);
	cpedsMoonDirection moon(obs,_altitude, _P, _T);
	obs*=PI180;
	obs.print_direction("my location:",true,0);
	double dscSize;
	msgs.say("the moon is currently at:",High);
	if (_JD==0) _JD=cpeds_julian_time();
	msgs.say("JD [UTC]: %.16lf",_JD,High);
	double JDutc=_JD+_da;
	DirectionRaDec moonNow;
//		moonNow=moon.nowAtRaDec(_JD+_da, 2000);
//		exit(0);
	moonNow=moon.nowAtRaDec(JDutc, JDutc,true);
	moonNow.toJD(JDutc);
	moonNow.print_direction("radec (geocentric)");
	moonNow=moon.nowAtRaDec(JDutc, JDutc,false);
	moonNow.print_direction("radec (topocentric)");
	moonNow.toJD(CPEDS_JD2000);
	moonNow.print_direction("radec (topocentric)");
	moonNow.toJD(JDutc);
	DirectionAh moonNowAh=moonNow.toAh(obs,JDutc,0,0,0,_P,_T,false);
	moonNowAh.check().print_direction("Ah (topocentric, no refraction, UT1=UTC, no polar motion)");
	moonNowAh=moonNow.toAh(obs,JDutc,0,0,0,_P,_T,true);
	moonNowAh.check().print_direction("Ah (topocentric, with refraction, UT1=UTC, no polar motion)");
	msgs.say("distance: %lf [km]",moon.distance(_JD)/1000,High);

	msgs.say("position angle of the bright limb: %lf", moon.brightLimbPA(JDutc)*PI180inv,High);
	msgs.say("illuminated disk fraction: %lf",moon.illuminatedDiskFraction(JDutc),High);
	msgs.say("phase angle: %lf",moon.opticalPhase(JDutc)*PI180inv,High);

	DirectionRaDec dsc=moon.getDarkSideCenter(JDutc,&dscSize);
	dsc.print_direction("dark side disk center radec (topocentric, no refraction, UT1=UTC, no polar motion):");
	moonNow.toJD(CPEDS_JD2000);
	moonNow.check().print_direction("radec (topocentric, no refraction, UT1=UTC, no polar motion):");
	moonNow.toJD(JDutc);
	moonNow.check().print_direction("radec (topocentric, no refraction, UT1=UTC, no polar motion):");
	moonNowAh=dsc.toAh(obs,JDutc,0,0,0,_P,_T,true);
	moonNowAh.check().print_direction("dark side disk center Ah (topocentric, with refraction, UT1=UTC, no polar motion)");
	msgs.say("dark side disk size [']: %lf\n",double(dscSize*60),High);
	
	// make calculation for the next 30 days
	double jdtmp,JD=cpeds_julian_time();
	jdtmp=JD;
	double period=_timeLimit;
	double dt=_dt/24;
	double timeLimit=JD+period;
	long N=period/dt+1;
	long i=0;
	while (JD<=timeLimit) {			i++;			JD+=dt;		} // first pass to learn how big matrix to allocate
	N=i;
	matrix<double> m(N,13);
	i=0; JD=jdtmp;
	while (JD<=timeLimit) {
		moonNow=moon.nowAtRaDec(JD,JD);
		moonNow.toJD(CPEDS_JD2000);
		moonNow*=double(PI180inv);
//		moon*=double(PI180inv);
		m(i,0)=JD;
		m(i,1)=moonNow.lon();
		m(i,2)=moonNow.lat();
		m(i,3)=moon.opticalPhase(JD)*PI180inv;
		m(i,4)=moon.illuminatedDiskFraction(JD);
		m(i,5)=moon.brightLimbPA(JD)*PI180inv;
		m(i,6)=moon.distance(JD)/1000;
		m(i,7)=moon.angularSize(JD)*PI180inv;
		moonNow=moon.getDarkSideCenter(JD,&dscSize);
		moonNowAh=moonNow.toAh(obs,JD,0,0,0,_P,_T,true);
		m(i,8)=dscSize*60;
		moonNow.toJD(CPEDS_JD2000);
		moonNow*=PI180inv;
		moonNowAh*=PI180inv;
		m(i,9)=moonNow.ra();
		m(i,10)=moonNow.dec();
		m(i,11)=moonNowAh.A();
		m(i,12)=moonNowAh.h();
			
		i++;
		JD+=dt;
//			msgs.say("calculating position, step: %li", i,Medium);
	}
	cpeds_matrix_save(m,"moon-data.txt");
	FILE* f=fopen("moon-data.readme","w");
	fprintf(f,"Columns description:\n");
	fprintf(f,"0 - JD [utc]\n");
	fprintf(f,"1 - moon center ra [deg] (J2000, topocentric)\n");
	fprintf(f,"2 - moon center dec [deg] (J2000, topocentric)\n");
	fprintf(f,"3 - phase [deg]\n");
	fprintf(f,"4 - illuminated disk fraction\n");
	fprintf(f,"5 - bright limb mid point position angle [deg] (north cusp, eastwards)\n");
	fprintf(f,"6 - distance Earth Center- Moon center [km]\n");
	fprintf(f,"7 - angular size [deg]\n");
	fprintf(f,"8 - dark side disk maximal angular size [']\n");
	fprintf(f,"9 - dark side center ra [deg] (J2000, topocentric)\n");
	fprintf(f,"10 - dark side center dec [deg] (J2000, topocentric)\n");
	fprintf(f,"11 - dark side center A [deg] (for JD, topocentric)\n");
	fprintf(f,"12 - dark side center h [deg] (for JD, topocentric)\n");
	fclose(f);
//	fn.save("moon-radec.txt");

}
void calculatePlanet(string planetName, cpedsMsgs& msgs) {
	cpedsDirection obs(_lon,_lat, _altitude);
	cpedsPlanetDirection planet(planetName,obs,_altitude, _P, _T);
	obs*=PI180;
	obs.print_direction("my location:",true,0);
	double dscSize;
	msgs.say("the planet %s is currently at:", planetName,High);
	if (_JD==0) _JD=cpeds_julian_time();
	double JDutc=_JD+_da;
	long yyyy,mm,dd,HH,MM;
	double SS;
	cpeds_JDToYMDhms(JDutc,&yyyy,&mm,&dd,&HH,&MM,&SS);
	std::cout << "JDutc: " << JDutc 
			<< ", " << yyyy << "-";
	std::cout.width(2);
	std::cout.fill('0');
	std::cout << mm << "-" << dd << " " << HH << ":" << MM << ":";
	std::cout.setf(ios::showpoint);
	std::cout << SS << " UTC\n";
	
	DirectionRaDec nowAt;
	nowAt=planet.nowAtRaDec(JDutc, JDutc,true);
	nowAt.check().print_direction("radec (geocentric)");
	nowAt=planet.nowAtRaDec(JDutc, JDutc,false);
	nowAt.print_direction("radec (topocentric)");
	nowAt.toJD(CPEDS_JD2000);
	nowAt.print_direction("radec (topocentric)");
	nowAt.toJD(JDutc);
	DirectionAh nowAtAh=nowAt.toAh(obs,JDutc,0,0,0,_P,_T,false);
	nowAtAh.check().print_direction("Ah (topocentric, no refraction, UT1=UTC, no polar motion)");
	nowAtAh=nowAt.toAh(obs,JDutc,0,0,0,_P,_T,true);
	nowAtAh.check().print_direction("Ah (topocentric, with optical refraction, UT1=UTC, no polar motion)");
	
	if (planetName == "Sun") {
		std::cout << "ang.size [deg]: " << planet.angularSize()*PI180inv << "\n";
	}

	// make calculation for the next days
	double jdtmp,JD=_JD;
	if (_JD==0) JD=cpeds_julian_time();
	jdtmp=JD;
	double period=_timeLimit;
	double dt=_dt/24;
	double timeLimit=JD+period;
	long N=period/dt+1;
	long i=0;
	while (JD<=timeLimit) {			i++;			JD+=dt;		} // first pass to learn how big matrix to allocate
	N=i;
	matrix<double> m(N,13);
	i=0; JD=jdtmp;
	while (JD<=timeLimit) {
		nowAt=planet.nowAtRaDec(JD,JD,0,true);
//		nowAt.toJD(CPEDS_JD2000);
		nowAt*=double(PI180inv);
		m(i,0)=JD;
		m(i,1)=nowAt.lon();
		m(i,2)=nowAt.lat();
		m(i,3)=planet.getDistance()/CPEDS_AU;
		nowAt=planet.nowAtRaDec(JD,JD);
		nowAtAh=nowAt.toAh(obs,JD,0,0,0,_P,_T,true);
		nowAtAh*=PI180inv;
		m(i,4)=nowAtAh.A();
		m(i,5)=nowAtAh.h();
			
		i++;
		JD+=dt;
//			msgs.say("calculating position, step: %li", i,Medium);
	}
	cpeds_matrix_save(m,planetName+"-data.txt");
	string fname=planetName+"-data.readme";
	FILE* f=fopen(fname.c_str(),"w");
	fprintf(f,"Columns description:\n");
	fprintf(f,"0 - JD [utc]\n");
	fprintf(f,"1 - planet center ra [deg] (J2000, topocentric)\n");
	fprintf(f,"2 - planet center dec [deg] (J2000, topocentric)\n");
	fprintf(f,"3 - distance to the planet [AU]\n");
	fprintf(f,"4 - planet center A [deg] (for JD, topocentric, with optical refrection)\n");
	fprintf(f,"5 - planet center h [deg] (for JD, topocentric, with optical refrection)\n");
	fclose(f);
//	fn.save("moon-radec.txt");
	
}

void calculateNextSolstice(cpedsMsgs& msgs) {
	cpedsDirection obs(_lon,_lat, _altitude);
	cpedsPlanetDirection planet("Sun",obs,_altitude, _P, _T);
	obs*=PI180;
	obs.print_direction("my location:",true,0);
	double dscSize;
	if (_JD==0) _JD=cpeds_julian_time();
	double JDutc=_JD+_da;
	double JD=JDutc;
	DirectionRaDec nowAt;
	
	// make calculation for the next days
	double period=_timeLimit;
	double dt=_dt/24;
	double timeLimit=JDutc+period;
	long N=period/dt+1;
	long i=0;
//	while (JD<=timeLimit) {			i++;			JD+=dt;		} // first pass to learn how big matrix to allocate
//	N=i;
//	matrix<double> m(N,4);
	FILE* f=fopen("Solstice.txt","w");
	i=0; JD=JDutc;
	mscsFunction JDdec;
	while (JD<=timeLimit) {
		nowAt=planet.nowAtRaDec(JD,JD,0,true);
		nowAt*=double(PI180inv);
//		m(i,0)=JD;
//		m(i,1)=nowAt.lon()*PI180inv;
//		m(i,2)=nowAt.lat()*PI180inv;
//		m(i,3)=planet.getDistance()/CPEDS_AU;
//		nowAtAh=nowAt.toAh(obs,JD,0,0,0,_P,_T,false);
//		m(i,4)=nowAtAh.A()*PI180inv;
//		m(i,5)=nowAtAh.h()*PI180inv;
		JDdec.newPoint(JD,nowAt.lat()*PI180inv);
		string tstr=cpeds_JDToYMDhms(JD,"%li-%02li-%02li %02li:%02li:%06.3lf");
		fprintf(f,"%s %.15lf %.10lf %.10lf\n",tstr.c_str(),JD,nowAt.lon(),nowAt.lat());
			
		i++;
		JD+=dt;
	}
	fclose(f);
	JD=JDdec.getX(JDdec.getMinValueIdx());
	string tstr=cpeds_JDToYMDhms(JD,"%li-%02li-%02li %02li:%02li:%06.3lf");
	printf("(min dec) Solstice on: %s\n", tstr.c_str() );
	printf("JD: %.15lf\n", JD );
	printf("dt: %.15lf [s]\n", _dt*3600 );
	JD=JDdec.getX(JDdec.getMaxValueIdx());
	tstr=cpeds_JDToYMDhms(JD,"%li-%02li-%02li %02li:%02li:%06.3lf");
	printf("(max dec) Solstice on: %s\n", tstr.c_str() );
	printf("JD: %.15lf\n", JD );
	printf("dt: %.15lf [s]\n", _dt*3600 );
}

void novaVSnovas(cpedsMsgs& msgs) {
	cpedsDirection obs(18.0,54.0,133);
	DirectionRaDec radec(250,5);
	DirectionAh Ah_nova, Ah_novas;
	double JD=2456302.871628;
	
	radec*=PI180;
	obs*=PI180;
	obs.print_direction("observer");
	radec.print_direction("radec");
	printf("JD: %.10lf\n",JD);
	Ah_nova=radec.toAh(obs,JD,false,false);
	Ah_nova.check();
	Ah_nova.print_direction("nova");
	Ah_novas=radec.toAh(obs,JD,0,0,0,0,0,false);
	Ah_novas.check();
	Ah_novas.print_direction("novas");
	printf("nova - novas\n");
	printf("dA [\"]: %lf\n",(Ah_nova.lon()-Ah_novas.lon())*3600*PI180inv);
	printf("dh [\"]: %lf\n",(Ah_nova.lat()-Ah_novas.lat())*3600*PI180inv);

	printf("second time scale errors impact on coordinates transformation\n");
	double sec=double(1)/3600/24;
	JD=2456302.871628+sec;
	printf("JD: %.10lf\n",JD);
	Ah_nova=radec.toAh(obs,JD,false,false);
	Ah_nova.check();
	Ah_nova.print_direction("nova");
	Ah_novas=radec.toAh(obs,JD,0,0,0,0,0,false);
	Ah_novas.check();
	Ah_novas.print_direction("novas");
	printf("nova - novas\n");
	printf("dA [\"]: %lf\n",(Ah_nova.lon()-Ah_novas.lon())*3600*PI180inv);
	printf("dh [\"]: %lf\n",(Ah_nova.lat()-Ah_novas.lat())*3600*PI180inv);
	
	
	printf("testing sidereal time calculations [seconds]\n");
//	JD=2456302.871628;
	JD=2456302.6;
//	JD=CPEDS_JD2000;
	printf("JD : %.10lf\n",JD);
	double NsecPerDay=24*3600;
	double Jday=floor(JD);
	double JDut=(JD-Jday);
	double UT=(JDut+0.5)*24*3600;
//	double UT=(JDut-0.5)*24*3600;
//	if (UT>=NsecPerDay) { UT=UT-NsecPerDay;  }
	printf("UT [s]: %.10lf (%lf h)\n",UT,UT/3600);
	printf("Jday [JD]: %.10lf\n",Jday);
	printf("JDut [JD]: %.10lf\n",JDut);
	
	double lonDt=obs.lon()*PI180inv/15.0*3600.0;
	double lstRT=RTcontrol_local_sidereal_time(UT,Jday,lonDt);
	printf("RT lst: %.15lf\n",lstRT);
	double mlst;
	double alst;
	sidereal_time(Jday,JDut,0,1,0,0,&alst);
	alst=alst*3600+lonDt;
	printf("novas apparent st method 0: %.15lf\n",alst);
	sidereal_time(Jday,JDut,0,0,0,0,&mlst);
	mlst=mlst*3600+lonDt;
	printf("novas mean st method 0: %.15lf\n",mlst);
	printf("diff novas mean-RTlst 0: %.15lf\n",mlst-lstRT);
	printf("diff novas app-RTlst 0: %.15lf\n",alst-lstRT);
	
	
//	sidereal_time(JD,0,0,1,0,0,&gst);
//	printf("novas apparent st method 0, split full-0: %.15lf\n",gst*3600);
//	sidereal_time(JD,0,0,0,0,0,&gst);
//	printf("novas mean st method 0, split full-0: %.15lf\n",gst*3600);
//	sidereal_time(Jday,JDut,0,1,1,0,&gst);
//	printf("novas apparent st method 1: %.15lf\n",gst*3600);
//	sidereal_time(Jday,JDut,0,0,1,0,&gst);
//	printf("novas mean st method 1: %.15lf\n",gst*3600);
	
	printf("julian day comparizon :\n");
	JD=cpeds_julian_time(2000,1,1,12);
	printf("cpeds: %.10lf\n",JD);
//	julian_date(2000,1,1,0);
	
	printf("comparizon with wiki definition:\n");
	double d=JD-CPEDS_JD2000;
	printf("wiki mean st : %.15lf\n",fmod(18.697374558 + 24.06570982441908 * d,24.0)*3600+lonDt);
	
	double x;
//	double dx=1.0/12;
	mscsFunction lstDiff,applst,meanlst,RTlst;
	double hour;
	for (long yy = 1980; yy <= 2013; yy++) {
		for (long i = 1; i <= 12; i++) {
			hour=21;
			JD=cpeds_julian_time(yy,i,1,hour);
			Jday=floor(JD);
			JDut=(JD-Jday);
//			UT=(JDut+0.5)*24*3600;
//			if (UT>=NsecPerDay) { UT=UT-NsecPerDay;  }
			UT=hour*3600;
			printf("-----------------\n",JD);
			printf("JD : %.10lf\n",JD);
			printf("JD RT : %.10lf  (month count from 0)\n",RTcontrol_JulianDay(yy,i,0));
//			printf("JD RT : %.10lf  (month count from 1)\n",RTcontrol_JulianDay(yy,i,1));
			printf("UT [s]: %.10lf\n",UT);
			printf("Jday [JD]: %.10lf\n",Jday);
			printf("JDut [JD]: %.10lf\n",JDut);
			sidereal_time(Jday,JDut,0,0,0,0,&mlst);
			mlst=mlst*3600+lonDt;
			sidereal_time(Jday,JDut,0,1,0,0,&alst);
			alst=alst*3600+lonDt;
			d=alst-RTcontrol_local_sidereal_time(UT,Jday,lonDt);
			x=yy+double(i)/12;
			RTlst.newPoint(x,RTcontrol_local_sidereal_time(UT,Jday,lonDt));
			applst.newPoint(x,alst);
			meanlst.newPoint(x,mlst);
					
			printf("x: %lf, JD: %.10lf, diff(app-RTctrl): %lf, diff(app-mean): %.10lf\n",x,JD,d,(alst-mlst));
			
		}
	}
	lstDiff=applst-RTlst;
	lstDiff.save("lstDiff-novasApp-RTcontrol");
	lstDiff=meanlst-RTlst;
	lstDiff.save("lstDiff-novasMean-RTcontrol");
	applst.save("applst");
	meanlst.save("meanlst");
	lstDiff=applst-meanlst;
	lstDiff.save("lstDiff-app-mean");
	
}


double RTcontrol_local_sidereal_time(double UT, double JD,double longitude )  
{
	double d;
	d = JD - 2451545.5E0 + UT / 86400.0;
//	d = JD - 2451545.0E0 + UT / 86400.0;
	d = fmod( 1.55811454642E0 + fmod(d + d, 2.0E0) + d * (0.547581870159E-2 +
			 d * (1.61549E-15 -d * 1.473E-24) ), 2.0E0);
	return 	fmod( 43200.0 * fmod( d + 2.0, 2.0E0) + longitude, 86400.0);
}
void novaReversabilityTest(cpedsMsgs& msgs) {
	cpedsDirection obs(18.0,54.0,133);
	DirectionAh ah(150,50);
	DirectionRaDec radec;
	
	ah*=PI180;
	obs*=PI180;
	radec=ah.toRaDec(obs,CPEDS_JD2000,false,false);
	ah.print_direction("ah");
	radec.print_direction("radec");
	ah=radec.toAh(obs,CPEDS_JD2000,false,false);
	ah.print_direction("ah2");
}

double RTcontrol_JulianDay(int tye,int tmo,int tda)    
{
	long i,j;
	j = tye + (tmo - 9) / 7;
	i = tmo + 10 - ((tmo + 9) / 12) * 12;
	return (double)( 1471086l + 365 * j + (j + 1000004l) / 4 + 30 * i + 
		i / 2 + ( (i / 7) * i ) % 2 + i / 12 + tda + (j + 1000400l) / 
		400 - (j + 1000100l) / 100 + 7502 );
};
void testRTlstCalculation() {

	double yy=2013;
	double mm=2;
	double dd=10;
	double h;
	double m=0;
	double s=0;
	double JD;
	double JDay;
	double UT;
	double lon=double(18.0)/15*3600;
	double alst,rtlst;
	for (long Y = 2000; Y <= 2013; Y++) {
		yy=Y;
		for (long M = 1; M <= 12; M++) {
			mm=M;
			for (long h = 0; h < 24; h++) {
				JDay=RTcontrol_JulianDay(yy,mm,dd);
				JD=cpeds_julian_time(yy,mm,dd,h);
				//		Jday=JulianDay(yy,mm,dd);
				UT=double(h)*3600+m*60+s;
				sidereal_time(JD,0,0,1,0,0,&alst);		alst=alst*3600+lon; if (alst>86400) alst-=86400;
				rtlst=RTcontrol_local_sidereal_time(UT,JDay,lon);
				printf("yy:%lf mm: %lf dd: %lf, h: %li m: %lf s: %lf, JD: %.10lf, UT: %lf, "
						"lst: %.10lf alstNOVAS: %.10lf, diff: %lf"
						"\n", 
						yy,mm,dd,h,m,s,JDay,UT, rtlst, alst, alst-rtlst);
				printf("  JDcpeds: %.10lf\n", JD);
				printf("--------------------\n");
			}
		}
	}


}
/***************************************************************************************/
void testConversionDeg2HMSDMS(double ra,double dec, double jdutc,cpedsMsgs& msgs) {
	QString cmd,Sradec;
	char radec[1000];

	double H,M,S,d,m,s;
	
	cpeds_angToHMS(ra,&H,&M,&S);
	cpeds_angToDMS(dec,&d,&m,&s);

	sprintf(radec,"%ih%im%fs %id%im%fs",int(H),int(M),S,int(d),int(m),s);
	Sradec=QString(radec);
	cmd="ps "+Sradec;
	
	msgs.say("--->sending to RT: "+cmd.toStdString(),High);

}











