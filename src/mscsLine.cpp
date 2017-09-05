/*!
 * \file mscsLine.cpp
 *
 *  Created on: Jan 16, 2012
 *      Author: blew
 */

#include "mscsLine.h"
#include "cpeds-direction.h"

mscsLine::mscsLine(string name) : mscsObject(name) {
	// TODO Auto-generated constructor stub
	setInitialValues();
	initialize_abc();
}
/***************************************************************************************/
mscsLine::mscsLine(cpedsPoint3D p1, cpedsPoint3D p2, string name) : mscsObject(name) {
	_params.x1=p1.x();
	_params.y1=p1.y();
	_params.z1=p1.z();
	_params.x2=p2.x();
	_params.y2=p2.y();
	_params.z2=p2.z();
	initialize_abc();
}
/***************************************************************************************/
mscsLine::mscsLine(double x1, double y1, double z1, double x2, double y2, double z2, double zmin, double zmax, string name) : mscsObject(name) {
	if (zmin==zmax and zmin==0) {
		zmin=z1;
		zmax=z2;
	}
	_params.x1=x1;
	_params.y1=y1;
	_params.z1=z1;
	_params.x2=x2;
	_params.y2=y2;
	_params.z2=z2;
	initialize_abc();
}
/***************************************************************************************/
mscsLine::~mscsLine() {
	// TODO Auto-generated destructor stub
}
/***************************************************************************************/
void mscsLine::setInitialValues() {
	_params.x1=0;
	_params.y1=0;
	_params.z1=0;
	_params.x2=0;
	_params.y2=0;
	_params.z2=0;	
}


/***************************************************************************************/
mscsLine::mscsLine(const mscsLine & parent) : mscsObject(parent) {	*this=parent; }
/***************************************************************************************/
cpedsPoint3D mscsLine::getValue_zparam(double z) {
	cpedsPoint3D p;
	p.setX(_params.a*(z-_params.z1)/_params.c+_params.x1);
	p.setY(_params.b*(z-_params.z1)/_params.c+_params.y1);
	p.setZ(z);
	return p;
}
/***************************************************************************************/
cpedsPoint3D mscsLine::getValue(double t) {
	cpedsPoint3D p;
	p.setX(_params.a*t+_params.x1);
	p.setY(_params.b*t+_params.y1);
	p.setZ(_params.c*t+_params.z1);
	return p;
}
/***************************************************************************************/
mscsLine & mscsLine::translate(const cpedsPoint3D & v) {
//	printf("translating line by: xyz: %lf %lf %lf\n",v.x(),v.y(),v.z());
	_params.x1+=v.x();
	_params.y1+=v.y();
	_params.z1+=v.z();
	_params.x2+=v.x();
	_params.y2+=v.y();
	_params.z2+=v.z();
	initialize_abc();
//	printf("translated line: \n");
//	getP1().print_point("p1");
//	getP2().print_point("p2");
	return *this;
}

/***************************************************************************************/
cpedsPointSet3D mscsLine::tabulateZ(double zmin, double zmax, double dz) {
	cpedsPointSet3D ps;
	cpedsPoint3D p;
	double z=zmin;
	long i=0;
	while (z<zmax) {
		ps.append(getValue_zparam(z));
		z=zmin+i*dz;
		i++;
	}
	return ps;
}
/***************************************************************************************/
double mscsLine::angle(mscsLine & l) const {
	double ang;
	cpedsPoint3D p1,p2;
	p1=getP2()-getP1();
//	p1.print_point("diff");
	p2=l.getP2()-l.getP1();
	return p1.angle(p2);
}
/***************************************************************************************/
double mscsLine::dot(mscsLine & l) const {
	double ang;
	cpedsPoint3D p1,p2;
	p1=getP1()-getP2();
	p2=l.getP1()-l.getP2();
	return p1.dot(p2);
}


/***************************************************************************************/
mscsLine & mscsLine::Rx(double ax) {
//	cpedsPoint3D p1ref(getP1());
//	cpedsPoint3D p2=getP2();
//	
//	// translate line;
////	translate(-getP1());
//	// rotate
//	p2=getP2();
//	p2.Rx(ax);
//	// reset line in a new orientation
//	_params.x1=0;
//	_params.y1=0;	
//	_params.z1=0;
//	_params.x2=p2.x();
//	_params.y2=p2.y();	
//	_params.z2=p2.z();
//	// translate line back
////	translate(p1ref);

	cpedsPoint3D p;
	p=getP1();
	p.Rz(ax);
	_params.x1=p.x();
	_params.y1=p.y();	
	_params.z1=p.z();
	p=getP2();
	p.Rz(ax);
	_params.x2=p.x();
	_params.y2=p.y();	
	_params.z2=p.z();
	initialize_abc();

	
	return *this;
}
/***************************************************************************************/
mscsLine & mscsLine::Ry(double ay) {
//	cpedsPoint3D p1ref(getP1());
//	cpedsPoint3D p2;
//	
//	// translate line;
////	translate(-getP1());
//	// rotate
//	p2=getP2();
//	p2.Ry(ay);
//	// reset line in a new orientation
//	_params.x1=0;
//	_params.y1=0;	
//	_params.z1=0;
//	_params.x2=p2.x();
//	_params.y2=p2.y();	
//	_params.z2=p2.z();
//	// translate line back
////	translate(p1ref);
	

	cpedsPoint3D p;
	p=getP1();
	p.Ry(ay);
	_params.x1=p.x();
	_params.y1=p.y();	
	_params.z1=p.z();
	p=getP2();
	p.Ry(ay);
	_params.x2=p.x();
	_params.y2=p.y();	
	_params.z2=p.z();
	initialize_abc();

	
	return *this;
}
/***************************************************************************************/
mscsLine & mscsLine::Rz(double az) {
//	cpedsPoint3D p1ref(getP1());
//	cpedsPoint3D p2;
	
	// translate line;
//	translate(-getP1());
//	printLine("Rz, after translation");
	// rotate
//	p2=getP2();
//	p2.Rz(az);
//	// reset line in a new orientation
//	_params.x1=0;
//	_params.y1=0;	
//	_params.z1=0;
//	_params.x2=p2.x();
//	_params.y2=p2.y();	
//	_params.z2=p2.z();

	//	printLine("Rz, after rotation");
		// translate line back
	//	translate(p1ref);
	//	printLine("Rz end, after translation");
		
	
	cpedsPoint3D p;
	p=getP1();
	p.Rz(az);
	_params.x1=p.x();
	_params.y1=p.y();	
	_params.z1=p.z();
	p=getP2();
	p.Rz(az);
	_params.x2=p.x();
	_params.y2=p.y();	
	_params.z2=p.z();
	initialize_abc();
	
	return *this;
}

/***************************************************************************************/
mscsLine& mscsLine::operator=(const mscsLine& rhs) {
	_params=rhs._params;
	return *this;	
}
/***************************************************************************************/
cpedsPointSet3D mscsLine::tabulateN(long  N) {
	double tmin=0;
	double tmax=1;
	double dt=double(1.0)/N;
	cpedsPointSet3D ps;
	for (long i = 0; i < N; i++) {
		ps.append(getValue(i*dt));
	}
	return ps;
}
/***************************************************************************************/
mscsLine mscsLine::cross(mscsLine & l) const {
	cpedsPoint3D p1=getP2()-getP1();
	cpedsPoint3D p2=l.getP2()-l.getP1();
//	cpedsPoint3D p0=getP1();
	return mscsLine( 0,0,0, p1.y()*p2.z()-p1.z()*p2.y(),  p1.z()*p2.x()-p1.x()*p2.z(),  p1.x()*p2.y()-p1.y()*p2.x() );
}
/***************************************************************************************/
double mscsLine::angleNoVect(mscsLine & l) const {
	double ang=angle(l);
	if (ang>PIsnd) { ang=PI-ang; }
	return ang;
}
/***************************************************************************************/
double mscsLine::x_zparam(double z) const { return _params.a*(z-_params.z1)/_params.c+_params.x1; }
/***************************************************************************************/
double mscsLine::y_zparam(double z) const { return _params.b*(z-_params.z1)/_params.c+_params.y1; }

/***************************************************************************************/
void mscsLine::initialize_abc() {
	_params.a=_params.x2-_params.x1;
	_params.b=_params.y2-_params.y1;
	_params.c=_params.z2-_params.z1;
}
/***************************************************************************************/
void mscsLine::printLine(string info) const{ 
	msgs->say("Line Info: "+info,High);
	getP1().print_point(getName()+" P1:");
	getP2().print_point(getName()+" P2:");
}

