/*!
 * \file cpedsMoonDirection.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: blew
 */

#include "cpedsMoonDirection.h"
#include <libnova/lunar.h>
#include <libnova/transform.h>
#include <libnova/sidereal_time.h>
#include <libnova/aberration.h>
#include "cpeds-consts.h"
#include "cpeds-project.h"
#include "cpeds-direction_set.h"
#include "novas.h"
#include "eph_manager.h"
#include "Mscs-function.h"

cpedsMoonDirection::cpedsMoonDirection(cpedsDirection obsLoc, double altitude, double pressure, double temperature) {
	_COproperties.observer=obsLoc;
	_COproperties.observerAltitude=altitude;
	_COproperties.pressure=pressure;
	_COproperties.temperature=temperature;
	
//	char *tmpch = new char[_COproperties.ephemeridesFile.length()];
	char *tmpch = new char[2000];
	sprintf(tmpch,"%s/external_packages/novas/JPLEPH",getenv("CPEDS_DIR"));
	_COproperties.ephemeridesFile=tmpch;
//	strcpy(tmpch,_COproperties.ephemeridesFile.c_str());
//	printf("JDstart: %lf, JDend: %lf\n",JDbeg, JDend);
//	exit(0);
	//	_COproperties.ephemeridesFile="/home/blew/programy/CPEDS/external_packages/novas3.0/de421/JPLEPH";
	_COproperties.distance=0;
	_COproperties.phase=0;
	_COproperties.equatorialRadius=1738140; // m
	_COproperties.polarRadius=1735970; // m
	_COproperties.COname="Moon";
	_COproperties.novas_obj_type=0; 
	_COproperties.novas_body_number=11;  // for the moon;
	_COproperties.bodyTemperature=235;
//	_COproperties.bodyTemperature=197; // mean lunation brightness temperature at 3.4 mm; 1968ApJ...151L.149E
	

	delete [] tmpch;

}

cpedsMoonDirection::cpedsMoonDirection(const cpedsMoonDirection& parent) : cpedsDirection(parent) {
	_COproperties=parent._COproperties;
	
}

cpedsMoonDirection::~cpedsMoonDirection() {

}


DirectionRaDec cpedsMoonDirection::nowAtRaDec(double JD, double epochJD, double ut1_utc, bool geocentric) {

	
	if (JD==0) JD=cpeds_julian_time();
	if (epochJD==2000) epochJD=CPEDS_JD2000;
//	printf("JD: %lf\n",JD);
	
	
	cat_entry dummy_star;
	object moon;
	short int error;
//	double pos[3], vel[3];
	double ra,dec,dis, jd_tt;
	jd_tt=JD+(CPEDS_TAI_UTC+CPEDS_TT_TAI_OFFSET)/86400.0;
//	printf("jd_tt: %15.8lf\n",jd_tt);
//	double jd[2];
//	jd[0]=long(JD);
//	jd[1]=JD-long(JD);

	on_surface obs;
	obs.longitude=_COproperties.observer.lon();
	obs.latitude=_COproperties.observer.lat();
	obs.height=_COproperties.observerAltitude;
	obs.pressure=_COproperties.pressure;
	obs.temperature=_COproperties.temperature;

//	printf("lon: %lf\n",obs.longitude);
//	printf("lat: %lf\n",obs.latitude);
//	printf("H: %lf\n",obs.height);
//	printf("P: %lf\n",obs.pressure);
//	printf("T: %lf\n",obs.temperature);

	double JDbeg, JDend;
	short int DEnum;
//	short int error=ephem_open(tmpch,&JDbeg, &JDend, &DEnum);
	error=ephem_open(_COproperties.ephemeridesFile.c_str(),&JDbeg, &JDend, &DEnum);
	if (error!=0) { 
		printf("cannot load the file: %s\n",_COproperties.ephemeridesFile.c_str());
		printf("the ephem_open error ws: %i\n",error);
		printf("make sure to set CPEDS_DIR environment variable to point to the source directory of the CPEMS package\n");
		exit(-1);
	}
#ifdef DEBUG
	printf("opening ephemerides file: %s\n",_COproperties.ephemeridesFile.c_str());
#endif

	make_cat_entry ("DUMMY","xxx",0,0.0,0.0,0.0,0.0,0.0,0.0,&dummy_star);
	error = make_object(_COproperties.novas_obj_type,_COproperties.novas_body_number,_COproperties.COname.c_str(),&dummy_star,&moon);
	if (error!=0) {
		printf("make object error: %i \n",error);		
	}



//	printf("JD:%lf\n",jd_tt);
	if (geocentric)
		error = app_planet(jd_tt,&moon,0,&ra,&dec,&dis);
	else
		error = topo_planet(jd_tt,&moon,CPEDS_TAI_UTC + CPEDS_TT_TAI_OFFSET - ut1_utc,&obs,0,&ra,&dec,&dis);
	if (error!=0) {
		printf("app/topo_planet error: %i \n",error);		
	}
	ephem_close();
#ifdef DEBUG
	printf("closing ephemerides file: %s\n",_COproperties.ephemeridesFile.c_str());
#endif

//	printf("ra [deg] %lf, dec [deg]: %lf topo planet error: %i \n",ra*15,dec,error);

	
	DirectionRaDec radec(ra*15.0*PI180,dec*PI180,JD);
	radec.toJD(epochJD);
//	lon()=radec.ra()*15.0*PI180;
//	lat()=radec.dec()*PI180; // BLcomment (Jan 7, 2015, 1:40:50 PM): commented as it doesn't seem to make sense.
	lon()=radec.ra();
	lat()=radec.dec();
	_COproperties.JD=JD;
	_COproperties.epochJD=epochJD;
	_COproperties.distance=dis*CPEDS_AU;


	
	
	
	
	
//#ifdef WITH_NOVAS

	
//	if (JD==0) JD=cpeds_julian_time();
//	if (epochJD==2000) epochJD=CPEDS_JD2000;
//	struct ln_equ_posn equ;
//	ln_get_lunar_equ_coords_prec(JD, &equ,0);
//	equ.ra*=PI180;
//	equ.dec*=PI180;
//	DirectionRaDec radec(equ.ra,equ.dec,JD);
//	radec.print_direction("debug");
//	exit(0);
//	radec.toJD(epochJD);
//	lon()=radec.ra();
//	lat()=radec.dec();
//	_COproperties.JD=JD;
//	_COproperties.epochJD=epochJD;
//#endif
	
	return radec;
}


/***************************************************************************************/
cpedsMoonDirection& cpedsMoonDirection::operator=(const cpedsMoonDirection& rhs) {
	cpedsDirection::operator=(rhs);
	_COproperties=rhs._COproperties;
	return (*this);
}
/***************************************************************************************/
double cpedsMoonDirection::illuminatedDiskFraction(double JD) {
	if (JD==0) JD=cpeds_julian_time();
	return ln_get_lunar_disk(JD);
}
double cpedsMoonDirection::opticalPhase(double JD) {
	if (JD==0) JD=cpeds_julian_time();
	_COproperties.phase=ln_get_lunar_phase(JD)*PI180;
	return _COproperties.phase;
}
double cpedsMoonDirection::brightLimbPA(double JD) {
	if (JD==0) JD=cpeds_julian_time();
	return ln_get_lunar_bright_limb(JD)*PI180;
}
/***************************************************************************************/
double cpedsMoonDirection::distance(double JD) {
	if (JD==0) JD=cpeds_julian_time();
	_COproperties.distance=ln_get_lunar_earth_dist(JD)*1000;
	return _COproperties.distance;
}
/***************************************************************************************/
double cpedsMoonDirection::angularSize(double JD, int which) {
	double avRadii=0;
	
	if (which==0)
		avRadii=0.5*(_COproperties.equatorialRadius+_COproperties.polarRadius);
	if (which==-1)
		avRadii=_COproperties.polarRadius;
	if (which==1)
		avRadii=_COproperties.equatorialRadius;
	
//	printf("avRadii: %lf [km]\n",avRadii/1000);
//	printf("distance: %lf [km]\n",distance(JD)/1000);
	double ang=atan(avRadii/distance(JD))*2.0;
	_COproperties.angularSize=ang;
	return ang;
}
/***************************************************************************************/
DirectionRaDec cpedsMoonDirection::getDarkSideCenter(double JD, double* dsCtrSize) {
	cpedsDirectionSet ds;
	cpedsDirection n0=nowAtRaDec(JD,JD);
	cpedsDirection n=n0;
	n.lat()+=angularSize(JD)/2.0;
#ifdef DEBUG
	n0.print_direction("n0");
	n.print_direction("n");
#endif
	ds.append(n);
//	ds.append(n0);
	cpedsProject proj(ds,"stere");
	cpedsPointSet3D ps;
	ps=proj.projectOnPlane(n0*PI180inv);
#ifdef DEBUG
	ps.print();
#endif
	
//	ds=proj.projectOnSphere(n0*PI180inv);
//	ds.print();
//	exit(0);

#ifdef DEBUG
	n0.print_direction("n0");
#endif
	cpedsPoint3D bsl; // bright side limb midpoint
	cpedsPoint3D dsl; // dark side limb midpoint
	cpedsPoint3D term; // terminator midpoint
	cpedsPoint3D dsctr; // dark side center point
	
	// calculate position of the bright side limb midpoint
	bsl=ps[0];
	bsl.Rz(brightLimbPA(JD)); 
	// calculate position of the dark side limb midpoint
	dsl=bsl;
	dsl.Rz(PI); 
	
	// calculate the position of the terminator midpoint
	term=ps[0];
//	term.print_point("term0");
	term.Rz(-PIsnd);
//	term.print_point("termPI2");
	term.Ry(-opticalPhase(JD));
//	term.print_point("termPI2Ry");
	term.Rz(brightLimbPA(JD)-PIsnd);
//	term.print_point("termPI2RyPA");
#ifdef DEBUG
	bsl.print_point("bsl");
	dsl.print_point("dsl");
	term.print_point("term");
#endif
	
	dsctr.set(dsl.x()-(dsl.x()-term.x())/2.0, dsl.y()-(dsl.y()-term.y())/2.0, dsl.y()-(dsl.y()-term.y())/2.0);
	
	proj.points().clear();
	proj.points().append(bsl);
	proj.points().append(dsl);
	proj.points().append(term);
	proj.points().append(dsctr);
//	proj.points().append(ps[1]);
#ifdef DEBUG
	proj.points().save("moon-proj-dirs.txt");
#endif
	ds=proj.projectOnSphere(n0*PI180inv);
#ifdef DEBUG
	ds.print();
#endif
	double bsl2term=cpeds_ang_n1n2(ds[0].get_direction(),ds[2].get_direction());
	double dsl2term=cpeds_ang_n1n2(ds[1].get_direction(),ds[2].get_direction());
	double bsl2dsl=cpeds_ang_n1n2(ds[0].get_direction(),ds[1].get_direction());
//	double bsl20=cpeds_ang_n1n2(ds[0].get_direction(),ds[4].get_direction());
//	double dsl20=cpeds_ang_n1n2(ds[1].get_direction(),ds[4].get_direction());
	
#ifdef DEBUG
	printf("ang bsl2dsl : %lf\n",bsl2dsl*PI180inv);
	printf("ang bsl2term: %lf\n",bsl2term*PI180inv);
	printf("ang dsl2term: %lf\n",dsl2term*PI180inv);
	printf("ang bsl2term+dsl2term: %lf\n",(bsl2term+dsl2term)*PI180inv);
	printf("ang err: %lf\n",(bsl2term+dsl2term-angularSize(JD))*PI180inv);
//	printf("ang bsl20+dsl20: %lf\n",(bsl20+dsl20)*PI180inv);
	printf("moon angular size: %lf\n",angularSize(JD)*PI180inv);
//	printf("moon min angular size: %lf\n",angularSize(JD,-1)*PI180inv);
//	printf("moon max angular size: %lf\n",angularSize(JD,1)*PI180inv);
#endif
	DirectionRaDec dsc; // dark side center
	dsc=ds[3];
//	dsc.setEpoch(_COproperties.epochJD);
	
	if (dsCtrSize!=NULL) (*dsCtrSize)=dsl2term*PI180inv;
#ifdef DEBUG
	dsc.print_direction("dsc");
#endif
	
	dsc.setEpoch(JD);
	dsc.check();
	return dsc;
}

/***************************************************************************************/
double cpedsMoonDirection::fluxPerBeam(double freq, double fwhm, double T) {
	if (T==0) T=_COproperties.bodyTemperature;
	double I=cpeds_black_body_radiance_Bnu(freq, T);
	double thetaMax=angularSize()/2;
	double dx=0.001*fwhm;
	double beamMax=5*fwhm;
	mscsFunction beam,tmp;
	beam.mkGauss(0,beamMax,dx,1,0,cpeds_fwhm2sigma(fwhm));
	tmp.mkSin(0,beamMax,dx,twoPI,0,1);
//	beam.save("beam");
	beam*=tmp;
	if (thetaMax>beamMax) thetaMax=beamMax;
//	double sourceSolidAngle=PI*pow(angSize,2)/4.0;
//	double beamSolidAngle=1.133*fwhm*fwhm;
	double flux = twoPI*I*beam.integrate(0,thetaMax); // 2pi*int I(th) * beam(th) sin(th) dth; where I(th) = 1 ; th< angSize/2 and 0 otherwise
	printf("planet: %s\n",_COproperties.COname.c_str());
	printf("freq: %lE\n",freq);
	printf("T: %lE\n",T);
	printf("I: %lE\n",I);
	printf("fwhm: %lE\n",fwhm*PI180inv);
	printf("angSize: %lE\n",angularSize()*PI180inv);
	printf("thetaMax: %lE\n",thetaMax*PI180inv);
	printf("flux: %lE\n\n",flux);
	if (flux<0) flux=-1;
	return flux;

}
