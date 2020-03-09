/*!
 * \file cpedsPlanetDirection.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: blew
 */

#include "cpedsPlanetDirection.h"
#include "Mscs-function.h"

cpedsPlanetDirection::cpedsPlanetDirection(string planetName, cpedsDirection obsLoc, double altitude, double pressure, double temperature) : cpedsMoonDirection(obsLoc,altitude,pressure,temperature) {
	_COproperties.COname=planetName;
	_COproperties.novas_obj_type=0; 
	_COproperties.novas_body_number=5;  // for the Jupiter
//	if (planetName=="Mercury") _COproperties.novas_body_number=1;
//	if (planetName=="Venus") _COproperties.novas_body_number=2;
//	if (planetName=="Earth") _COproperties.novas_body_number=3;
//	if (planetName=="Mars") _COproperties.novas_body_number=4;
//	if (planetName=="Jupiter") _COproperties.novas_body_number=5;
//	if (planetName=="Saturn") _COproperties.novas_body_number=6;
//	if (planetName=="Uranus") _COproperties.novas_body_number=7;
//	if (planetName=="Neptune") _COproperties.novas_body_number=8;
//	if (planetName=="Pluto") _COproperties.novas_body_number=9;
//	if (planetName=="Sun") _COproperties.novas_body_number=10;

	setPlanetProperties(planetName);
}

cpedsPlanetDirection::~cpedsPlanetDirection() {
}


cpedsPlanetDirection::cpedsPlanetDirection(const cpedsPlanetDirection& parent) {
	this->operator=(parent);
}
cpedsPlanetDirection::cpedsPlanetDirection(const cpedsMoonDirection& parent) : cpedsMoonDirection(parent) {
	
}

//virtual DirectionRaDec cpedsPlanetDirection::nowAtRaDec(double JD, double epochJD, double ut1_utc, bool geocentric) {
//}

const cpedsPlanetDirection& cpedsPlanetDirection::operator=(const cpedsPlanetDirection& rhs) {
	cpedsMoonDirection::operator=(rhs);
	return (*this);
}
const cpedsPlanetDirection& cpedsPlanetDirection::operator=(const cpedsMoonDirection& rhs) {
	cpedsMoonDirection::operator=(rhs);
	return (*this);	
}
/***************************************************************************************/
double cpedsPlanetDirection::angularSize(double JD, int which) {
	double ang=-1, r=0;
	if (JD==0) {
		if (which==-1) r=_COproperties.polarRadius;
		if (which==0) r=(_COproperties.polarRadius+_COproperties.equatorialRadius)/2;
		if (which==1) r=_COproperties.equatorialRadius;
		if (getDistance()!=0)
			ang=2.0*atan(r/getDistance());
	}
	return ang;
}
/***************************************************************************************/
void cpedsPlanetDirection::setPlanetProperties(string planetName) {
	if (planetName=="Mercury") {
		_COproperties.novas_body_number=1;
		_COproperties.equatorialRadius=2439700;
		_COproperties.polarRadius=2439700;
		_COproperties.bodyTemperature=-1;
	}
	if (planetName=="Venus") {
		_COproperties.novas_body_number=2;
		_COproperties.equatorialRadius=6051800;
		_COproperties.polarRadius=6051800;
		_COproperties.bodyTemperature=296; // 1968ApJ...151L.149E
		
	}
	if (planetName=="Earth") {
		_COproperties.novas_body_number=3;
		_COproperties.equatorialRadius=-1;
		_COproperties.polarRadius=-1;
	}
	if (planetName=="Mars") {
		_COproperties.novas_body_number=4;
		_COproperties.equatorialRadius=3396200;
		_COproperties.polarRadius=3376200;
		_COproperties.bodyTemperature=230; // Rough estimate http://www.astronomycafe.net/qadir/q2681.html
		_COproperties.bodyTemperature=163; // 1968ApJ...151L.149E
	}
	if (planetName=="Jupiter") {
		_COproperties.novas_body_number=5;
		_COproperties.equatorialRadius=71492000;
		_COproperties.polarRadius=66854000;
		_COproperties.bodyTemperature=140; // 1968ApJ...151L.149E
	}
	if (planetName=="Saturn") {
		_COproperties.novas_body_number=6;
		_COproperties.equatorialRadius=60268000;
		_COproperties.polarRadius=54364000;
		_COproperties.bodyTemperature=130; // 1968ApJ...151L.149E
	}
	if (planetName=="Uranus") {
		_COproperties.novas_body_number=7;
		_COproperties.equatorialRadius=25559000;
		_COproperties.polarRadius=24973000;
		_COproperties.bodyTemperature=-1; 
	}
	if (planetName=="Neptune") {
		_COproperties.novas_body_number=8;
		_COproperties.equatorialRadius=24764000;
		_COproperties.polarRadius=24341000;
		_COproperties.bodyTemperature=-1; 
	}
	if (planetName=="Pluto") {
		_COproperties.novas_body_number=9;
		_COproperties.equatorialRadius=1153000;
		_COproperties.polarRadius=1153000;
		_COproperties.bodyTemperature=-1; 
	}
	if (planetName=="Sun") {
		_COproperties.novas_body_number=10;
		_COproperties.equatorialRadius=696342000;
		_COproperties.polarRadius=696342000;
		_COproperties.bodyTemperature=6600; // 1968ApJ...151L.149E
	}

	_COproperties.angularSize=angularSize();
}
/***************************************************************************************/
//double cpedsPlanetDirection::fluxPerBeam(double freq, double fwhm, double T) {
//	if (T==0) T=_COproperties.bodyTemperature;
//	double I=cpeds_black_body_radiance_Bnu(freq, T);
//	double thetaMax=angularSize()/2;
//	double dx=0.001*fwhm;
//	double beamMax=5*fwhm;
//	mscsFunction beam,tmp;
//	beam.mkGauss(0,beamMax,dx,1,0,cpeds_fwhm2sigma(fwhm));
//	tmp.mkSin(0,beamMax,dx,twoPI,0,1);
////	beam.save("beam");
//	beam*=tmp;
//	if (thetaMax>beamMax) thetaMax=fwhm;
////	double sourceSolidAngle=PI*pow(angSize,2)/4.0;
////	double beamSolidAngle=1.133*fwhm*fwhm;
//	double flux = twoPI*I*beam.integrate(0,thetaMax); // 2pi*int I(th) * beam(th) sin(th) dth; where I(th) = 1 ; th< angSize/2 and 0 otherwise
//	printf("planet: %s\n",_COproperties.COname.c_str());
//	printf("freq: %lE\n",freq);
//	printf("T: %lE\n",T);
//	printf("I: %lE\n",I);
//	printf("fwhm: %lE\n",fwhm*PI180inv);
//	printf("angSize: %lE\n",angularSize()*PI180inv);
//	printf("thetaMax: %lE\n",thetaMax*PI180inv);
//	printf("flux: %lE\n\n",flux);
//	if (flux<0) flux=-1;
//	return flux;
//}
