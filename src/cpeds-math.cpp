#define _ISOC99_SOURCE


#include <features.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
//#include <libnova/transform.h>
#include <libnova/julian_day.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
//#include "/home/blew/programy/CPEDS/cpeds-defs.h"
//#include <ncarg/ngmath.h>
#include "sys/stat.h"
#include "cpeds-math.h"
#include "cpeds-templates.h"
//#include "matrix.h"
#include <fftw3.h>
#include "../external_packages/novas/novas.h"

#ifndef _NO_NAMESPACE
/* using namespace LiDIA; */
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif







//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_sf_legendre.h>
//#include <gsl/gsl_sf_result.h>
//double PI = 3.141592654;
// the spherical coordinate system is with theta = 0 at the north pole

//using namespace LiDIA;


// some usefull declarations
/* double* cpeds_generate_gaussian_distribuant_function(long int size, double s, double x0); */
/* double* cpeds_random_uniform_numbers(double min, double max, int pix_num); */
/* double cpeds_random_uniform_number(double min, double max); */
/* double cpeds_random_gauss_number(double m,double s, long int pix_num); */
/* double* cpeds_random_gauss_numbers(double m, double s,long int pix_num, int method); */
/* long int  cpeds_xy2num(long int x, long int y); // y can only by 0 or 1 -- this is conversion from linear table to rectangular table of size nx2 */
/* long int cpeds_count_numbers_in_bin(long int k, double * t, double bin, double x); */
//***************************************************************************************************
//using namespace std;

double cpeds_sph2cart(int coord, double th, double phi) { // angles passed in radians
	// coord - 0 - x, 1 - y, 2 - z
	double x,y,z;
	
	//  ath = dr*ath; aphi=dr*aphi; o = dr*o;
	
	if (coord == 0) { // x
		x = sin(th)*cos(phi);
		return x;
	}
	
	if (coord == 1) { // y
		y = sin(th)*sin(phi);
		return y;
	}
	
	if (coord == 2) { // z
		z = cos(th);
		return z;
	}
	return -1;
}

double cpeds_cart2sph(int coord, double x, double y,double z) {
	// coord - 0 - theta, 1 - phi
	double th=0,phi=0;
	//double dr,rd;//,PI;
	//PI = 3.141592654;
	//dr=PI/180; rd=180/PI;
	//  ath = dr*ath; aphi=dr*aphi; o = dr*o;
	
	//printf("%f %f %f\n",x,y,z);
	if (coord == 0) {
		th = acos(z/sqrt(x*x+y*y+z*z)); return th;}
	
	if (coord == 1) {
		if (y > 0) {
			if (x > 0) {phi = atan(y/x); }
			if (x < 0) {phi = atan(y/x) + PI; }
			if (x == 0) { phi = PIsnd; }}
		else
			if (y < 0) {
				if (x > 0) {phi = atan(y/x) + twoPI; }
				if (x < 0) {phi = atan(y/x) + PI; }
				if (x == 0) { phi = PI+PIsnd; }}
			else {
				if (y == 0) {
					if (x >= 0) {phi = 0; }
					if (x < 0) {phi = PI; }}
			}
		return phi; }
	return -1;
}


/****************************************************************************************************************/
//returns a cos(angle) between the two directions given in radians on the sky in radians
double cpeds_cosang_n1n2(double th1, double phi1,double th2,double phi2) {
	double x1,y1,z1,x2,y2,z2;
	double cosang;
	
	x1 = cpeds_sph2cart(0,th1,phi1);   y1 = cpeds_sph2cart(1,th1,phi1);   z1 = cpeds_sph2cart(2,th1,phi1);
	x2 = cpeds_sph2cart(0,th2,phi2);   y2 = cpeds_sph2cart(1,th2,phi2);   z2 = cpeds_sph2cart(2,th2,phi2);
	
	cosang = x1*x2 + y1*y2 + z1*z2;
	
	if (cosang < -1.0) cosang=-1.0;
	else {
		if (cosang > 1.0) cosang=1.0;
	}
	
	return cosang;
}

/****************************************************************************************************************/
//returns an angle between the two directions given in radians on the sky in radians
double cpeds_ang_n1n2(double th1, double phi1,double th2,double phi2) {
	double x1,y1,z1,x2,y2,z2;
	double ang,arg;
	
	x1 = cpeds_sph2cart(0,th1,phi1);   y1 = cpeds_sph2cart(1,th1,phi1);   z1 = cpeds_sph2cart(2,th1,phi1);
	x2 = cpeds_sph2cart(0,th2,phi2);   y2 = cpeds_sph2cart(1,th2,phi2);   z2 = cpeds_sph2cart(2,th2,phi2);
	//ang = acos( (x1*x2 + y1*y2 + z1*z2)/( sqrt((x1*x1 + y1*y1 + z1*z1) * (x2*x2 + y2*y2 + z2*z2)) ));
	arg=(x1*x2 + y1*y2 + z1*z2);
	if (arg < -1.0) arg=-1.0;
	else {
		if (arg > 1.0) arg=1.0;
	}
	ang = acos(arg);
	
	/*   printf("|cpeds_ang_n1n2> x1,y1,z1: %.20lE, %.20lE, %.20lE\n",x1,y1,z1); */
	/*   printf("|cpeds_ang_n1n2> x2,y2,z2: %.20lE, %.20lE, %.20lE\n",x2,y2,z2); */
	/*   printf("|cpeds_ang_n1n2> ang = %.20lE ang[deg] = %lf\n",ang,ang*PI180inv); */
	/*   printf("|cpeds_ang_n1n2> phi1, th1= %.20lE, %.20lE\n",phi1,th1); */
	/*   printf("|cpeds_ang_n1n2> phi2, th2= %.20lE, %.20lE\n",phi2,th2); */
	return ang;
}
/****************************************************************************************************************/
double cpeds_ang_n1n2(cpeds_direction n1, cpeds_direction n2) {
	return cpeds_ang_n1n2(PIsnd-n1.b, n1.l, PIsnd-n2.b, n2.l);
}

/****************************************************************************************************************/
/****************************************************************************************************************/
/* Cpeds Rotations */
// angles are given in radians

//! Rotates the point p about the Ox axix by Ax angle given in radians
/*! The coordinate convention is that: */
/*! Ox: (1,0,0) */
/*! Oy: (0,1,0) */
/*! Oz: (0,0,1) */
cpeds_direction cpeds_Rx(cpeds_direction p,double Ax) {
	cpeds_direction pp;
	double x,y,z,xx,yy,zz,th,phi,sinAx,cosAx;
	
	th=PIsnd-p.b; phi = p.l; // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0,th,p.l);  y = cpeds_sph2cart(1,th,p.l);   z = cpeds_sph2cart(2,th,p.l);   sinAx=sin(Ax); cosAx=cos(Ax);
	
	xx=x;
	yy=y*cosAx-z*sinAx;
	zz=y*sinAx+z*cosAx;
	
	th = cpeds_cart2sph(0,xx,yy,zz); phi = cpeds_cart2sph(1,xx,yy,zz);
	cpeds_check_thphi(&th,&phi);
	pp.b=PIsnd-th; pp.l=phi; // change from (th,phi) CS to (l,b) CS
	return pp;
}
void cpeds_Rx(double x, double y, double z,double Ax, double* xr, double* yr, double* zr) {
	double xx,yy,zz,sinAx,cosAx;
	sinAx=sin(Ax); cosAx=cos(Ax);
	xx=x;
	yy=y*cosAx-z*sinAx;
	zz=y*sinAx+z*cosAx;
	*xr=xx;
	*yr=yy;
	*zr=zz;
}

/****************************************************************************************************************/
//! Rotates the point p about the Oy axix by Ay angle given in radians
cpeds_direction cpeds_Ry(cpeds_direction p,double Ay) {
	cpeds_direction pp;
	double x,y,z,xx,yy,zz,th,phi,sinAy,cosAy;
	
	th=PIsnd-p.b; phi = p.l; // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0,th,p.l);  y = cpeds_sph2cart(1,th,p.l);   z = cpeds_sph2cart(2,th,p.l);   sinAy=sin(Ay); cosAy=cos(Ay);
	
	xx=x*cosAy+z*sinAy;
	yy=y;
	zz=-x*sinAy+z*cosAy;
	
	th = cpeds_cart2sph(0,xx,yy,zz); phi = cpeds_cart2sph(1,xx,yy,zz);
	cpeds_check_thphi(&th,&phi);
	pp.b=PIsnd-th; pp.l=phi; // change from (th,phi) CS to (l,b) CS
	return pp;
}
void cpeds_Ry(double x, double y, double z,double Ay, double* xr, double* yr, double* zr) {
	double xx,yy,zz,sinAy,cosAy;
	sinAy=sin(Ay); cosAy=cos(Ay);
	xx=x*cosAy+z*sinAy;
	yy=y;
	zz=-x*sinAy+z*cosAy;
	*xr=xx;
	*yr=yy;
	*zr=zz;
}

/****************************************************************************************************************/
//! Rotates the point p about the Oz axix by Az angle given in radians
cpeds_direction cpeds_Rz(cpeds_direction p,double Az) {
	cpeds_direction pp;
	double x,y,z,xx,yy,zz,th,phi,sinAz,cosAz;
	
	th=PIsnd-p.b; phi = p.l; // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0,th,p.l);  y = cpeds_sph2cart(1,th,p.l);   z = cpeds_sph2cart(2,th,p.l);   sinAz=sin(Az); cosAz=cos(Az);
	
	xx=x*cosAz-y*sinAz;
	yy=x*sinAz+y*cosAz;
	zz=z;
	
	th = cpeds_cart2sph(0,xx,yy,zz); phi = cpeds_cart2sph(1,xx,yy,zz);
	cpeds_check_thphi(&th,&phi);
	pp.b=PIsnd-th; pp.l=phi; // change from (th,phi) CS to (l,b) CS
	return pp;
}
void cpeds_Rz(double x, double y, double z,double Az, double* xr, double* yr, double* zr) {
	double xx,yy,zz,sinAz,cosAz;
	sinAz=sin(Az); cosAz=cos(Az);
	xx=x*cosAz-y*sinAz;
	yy=x*sinAz+y*cosAz;
	zz=z;
	*xr=xx;
	*yr=yy;
	*zr=zz;
}

/***************************************************************************************/
void cpeds_cross_product_3d(double x1, double y1, double z1,double x2, double y2, double z2, double* x3, double* y3, double* z3) {
	*x3 = y1*z2-y2*z1;
	*y3 = z1*x2-z2*x1;
	*z3 = x1*y2-x2*y1;
}

double cpeds_dot_product_3d(double x1, double y1, double z1,double x2, double y2, double z2) {
	return x1*x2+y1*y2+z1*z2;
}
double cpeds_dot_product_3d(cpeds_direction n1, cpeds_direction n2) {
	double x1=cpeds_sph2cart(0, PIsnd-n1.b, n1.l);
	double y1=cpeds_sph2cart(1, PIsnd-n1.b, n1.l);
	double z1=cpeds_sph2cart(2, PIsnd-n1.b, n1.l);
	double x2=cpeds_sph2cart(0, PIsnd-n2.b, n2.l);
	double y2=cpeds_sph2cart(1, PIsnd-n2.b, n2.l);
	double z2=cpeds_sph2cart(2, PIsnd-n2.b, n2.l);
	return cpeds_dot_product_3d(x1,y1,z1,x2,y2,z2);
}
/***************************************************************************************/
void cpeds_swap(double& d1, double& d2) {
	double tmp=d1;
	d1=d2;
	d2=tmp;
}
/***************************************************************************************/
void cpeds_swap(int& d1, int& d2) {
	int tmp=d1;
	d1=d2;
	d2=tmp;
}
/***************************************************************************************/
void cpeds_swap(long& d1, long& d2) {
	long tmp=d1;
	d1=d2;
	d2=tmp;
}

/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

cpeds_direction cpeds_horizontalToEquatorialFirst(double A, double h, double th) {
	double ha, ha1, ha2, d;
	double acc=1e-10;
	cpeds_direction p;
	
	d=asinl( sinl(th)*sinl(h) + cosl(th)*cosl(h)*cosl(A) );
	
	if ( fabs(PIsnd-d) < acc || fabs(PIsnd+d) < acc ) ha=0;
	else {
		ha=atan2l( -cosl(h) * sinl(A), sinl(h) * cosl(th) -	cosl(h) * cosl(A) * sinl(th) );
		
		
		//		ha1=-cos(h)*sin(A)/cos(d); // this is the sine of the angle
		//		ha2=(sin(h)*cos(th) - cos(h)*cos(A)*sin(th))/cos(d); // this is the cosine of the angle
		//		if (fabs(ha1)<0) ha1=0; // removing -0
		//		if (fabs(ha2)<0) ha2=0; // removing -0
		//		//    printf("sin ha = %lE, cos ha = %lE ",ha1,ha2);
		//		
		//		if (ha1==0 && ha2==0) { ha=0; }
		//		else {
		//			if (ha1==0 && ha2 >0) { ha=0; }
		//			else {
		//				if (ha1==0 && ha2 < 0) { ha=PI; }
		//				else {
		//					if (ha1>0 && ha2==0) { ha=PIsnd; }
		//					else {
		//						if (ha1<0 && ha2==0) { ha=PI+PIsnd; }
		//						else {
		//							if (ha1 >= 0 && ha2 >= 0) { ha=asin(ha1); } // first quadrant
		//							else {
		//								if (ha1 >= 0 && ha2 < 0) { ha=acos(ha2); } // second qudrant
		//								else {
		//									if (ha1 < 0 && ha2 < 0) { ha=PI-asin(ha1); } // third quadrant
		//									else { ha=twoPI-asin(-ha1); } // fourth quadrant
		//								}
		//							}
		//						}
		//					}
		//				}
		//			}
		//		}
	}
	
	cpeds_check_bl(&d,&ha);
	p.l=ha; p.b=d;
	return p;
	
}
/****************************************************************************************************************/
cpeds_direction cpeds_horizontalToEquatorialFirst(cpeds_direction n, double th) {
	return cpeds_horizontalToEquatorialFirst(n.l,n.b,th);
}
/****************************************************************************************************************/
cpeds_direction cpeds_equatorialFirstToHorizontal(double ha, double d, double th) {
	double A, A1,A2, h;
	double acc=1e-10;
	cpeds_direction p;
	
	if (sin(th)*sin(d) + cos(th)*cos(d)*cos(ha) > 1 || sin(th)*sin(d) + cos(th)*cos(d)*cos(ha) < -1) {    printf("error in cpeds: cpeds_equatorialFirstToHorizontal\n"); exit(0); } // DEBUG STUFF
	h=asin( sin(th)*sin(d) + cos(th)*cos(d)*cos(ha) );
	
	if ( fabs(PIsnd-h) < acc || fabs(PIsnd+h) < acc ) A=0;
	else {
		A1=-cos(d)*sin(ha)/cos(h); // this is the sine of the angle
		A2=(sin(d)*cos(th) - cos(d)*cos(ha)*sin(th))/cos(h); // this is the cossine of the angle
		if (fabs(A1)<acc) A1=+0.0; // removing -0
		if (fabs(A2)<acc) A2=+0.0; // removing -0
		
		if (A1==0 && A2==0) { A=0; }
		else {
			if (A1==0 && A2>0) { A=0; }
			else { if (A1==0 && A2<0) { A=PI; }
			else {
				if (A1>0 && A2==0) { A=PIsnd; }
				else {
					if (A1<0 && A2==0) { A=PI+PIsnd; }
					else {
						/* 	      A=atan2(A1,A2); */
						/* 	      if (A<0) A=PI-A; */
						if (A1 > 0 && A2 > 0) { A=asin(A1); } // first quadrant
						else {
							if (A1 > 0 && A2 < 0) { A=acos(A2); } // second qudrant
							else {
								if (A1 < 0 && A2 < 0) { /* if (A1<-0.5) A=PI+acos(-A2); else */ A=PI-asin(A1); } // third quadrant
								else { A=twoPI-asin(-A1); } // fourth quadrant
							}
						}
						
					}
				}
			}
			}
		}
	}
	
	/*   printf("\n\n A1 = %.20lf,  A2 = %.12lf A=%.20lf\n\n",PI180inv*A1,PI180inv*A2,PI180inv*A); */
	
	cpeds_check_bl(&h,&A);
	
	p.l=A; p.b=h;
	return p;
}
/****************************************************************************************************************/
cpeds_direction cpeds_equatorialFirstToHorizontal(cpeds_direction n, double th) {
	return cpeds_equatorialFirstToHorizontal(n.l,n.b,th);
}

/****************************************************************************************************************/
/****************************************************************************************************************/
cpeds_direction cpeds_equatorialFirstToEquatorialSecond(cpeds_direction n, double lst) {
	
	n.l=lst-n.l;
	cpeds_check_bl(&n.b,&n.l);
	return n;
}

/****************************************************************************************************************/
cpeds_direction cpeds_equatorialSecondToEquatorialFirst(cpeds_direction n, double lst) {
	return cpeds_equatorialFirstToEquatorialSecond(n,lst);
}

/****************************************************************************************************************/
/****************************************************************************************************************/

double cpeds_local_sidereal_time2(double localJulianTime) {
	double d;
	d = localJulianTime - 2451545.5E0;
	d = fmod( 1.55811454642 + fmod(d + d, 2.0) + d * (0.547581870159e-2 + d * (1.61549e-15 -d * 1.473e-24) ), 2.0);
	printf("lst in sec %20.10lf\n",fmod( 43200.0 * fmod( d + 2.0, 2.0), 86400.0));
	printf("lst day frac %20.10lf\n",fmod( 43200.0 * fmod( d + 2.0, 2.0), 86400.0)/86400.0);
	return twoPI*fmod( 43200.0 * fmod( d + 2.0, 2.0), 86400.0)/86400.0;
}

double cpeds_local_sidereal_time(double localJulianTime, double lon) {
	double JD,S,T0,T,UT;
	double loct;
	double JD0;
	JD0=localJulianTime-lon/360.0; // convert to UTC time first;
	JD=floor(JD0);
	loct=JD0-JD; // local solar time at lon=0
	if (loct < 0.5) { JD=JD-0.5; UT=(loct+0.5)*24.0; } else { JD=JD+0.5; UT=(loct-0.5)*24.0; } // this is the Julian Day corresponding to the localJulianTime
	S=JD-2451545.0;
	T=S/36525.0;
	T0=6.697374558 + (2400.051336*T)+(0.000025862*T*T); // printf("T0=%lf\n",T0);
	T0=fmod(T0,24.0);
	if (T0<0) T0+=24.0; // printf("T0=%lf\n",T0);
	loct=UT*1.002737909; // printf("locst=%lf\n",loct); // convert the local solar time to local star time
	T0+=loct; // printf("T0=%lf\n",T0);
	T0=fmod(T0,24.0);  // printf("T0=%lf\n",T0); // this is local star time in hours
	T0+=lon/360.0*24.0;  // printf("T0=%lf\n",T0); // add the effect of the longitude
	return T0/24.0*twoPI; // this is the local star time // this calculation can be improved, I know :)
}
/***************************************************************************************/
double cpeds_local_sidereal_time_novas(double jd_ut1, double ut1_utc, double DeltaAT,double longitude, int LSTtype ) {
	double tt_ut1=DeltaAT+CPEDS_TT_TAI-ut1_utc;
	double lst=0; // hours
//	sidereal_time(jd_ut1,0,tt_ut1,1,1,0,&alst);		
	sidereal_time(jd_ut1,0,tt_ut1,LSTtype,1,0,&lst);		 // this is mean LST which should be used when nutation in enabled converting coordinates
	lst=lst*3600+longitude; 
	if (lst>86400) lst-=86400;
	return lst;
}
/****************************************************************************************************************/
double cpeds_julian_time(long year, long month, long day, double hour) {
	return cpeds_julian_local_time(year,month,day,hour,0);
}
/***************************************************************************************/
double cpeds_timeSec_to_julian_time(double s) {
	double secFrac=s-long(s);
	s=long(s);
	time_t rawtime=(time_t)s;
	//	printf("rawtime: %.10lf\n",double(rawtime));
	struct tm * timeinfo;
	timeinfo=gmtime(&rawtime);
	double hour=timeinfo->tm_hour+double(timeinfo->tm_min)/60.0+double(timeinfo->tm_sec+secFrac)/3600.0;
	
	//	struct timeval tv;
	
	
	return cpeds_julian_time(timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday,hour);
	
}

/****************************************************************************************************************/
double cpeds_julian_local_time(long year, long month, long day, double hour, double longitude) {
	
	long jd12h;
	
	double tjd;
	
	jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
			+ ((long) month - 14L) / 12L) / 4L
			+ 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
			/ 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
					/ 100L) / 4L;
	tjd = (double) jd12h - 0.5 + hour / 24.0;
	
	
	tjd+=longitude/360;
	
	return tjd;
	
};

/****************************************************************************************************************/
void cpeds_JDToYMDhms(double JD,long* year, long* month, long* day, long* hour, long* minute, double* sec) {
	struct ln_date date;
	ln_get_date(JD,&date);
	if (year!=NULL) *year=date.years;
	if (month!=NULL) *month=date.months;
	if (day!=NULL) *day=date.days;
	if (hour!=NULL) *hour=date.hours;
	if (minute!=NULL) *minute=date.minutes;
	if (sec!=NULL) *sec=date.seconds;
}
/***************************************************************************************/
string cpeds_JDToYMDhms(double JD, string fmt) {
	if (JD==-1e9) {
		JD=cpeds_julian_time();
	}
	long YYYY,MM,DD,h,m;
	double s;
	cpeds_JDToYMDhms(JD,&YYYY,&MM,&DD,&h,&m,&s);
	//	if (s==60) {
	//		JD=cpeds_julian_time(YYYY,MM,DD,h,m,s);
	//		cpeds_JDToYMDhms(JD,&YYYY,&MM,&DD,&h,&m,&s);
	//	}
	char buff[100];
	bzero(buff,100);
	sprintf(buff,fmt.c_str(),YYYY,MM,DD,h,m,s);
	string str=buff;
	//	printf("cpeds_JDToYMDhms: DEBUG: |%s|\n",str.c_str());
	return str;
}
/***************************************************************************************/
long cpeds_JDToYear(double JD) {
	if (JD==-1e9) {
		JD=cpeds_julian_time();
	}
	long YYYY,MM,DD,h,m;
	double s;
	cpeds_JDToYMDhms(JD,&YYYY,&MM,&DD,&h,&m,&s);
	return YYYY;
}
//double cpeds_JDToYearFrac(double JD=-1e9) {
//	if (JD==-1e9) {
//		JD=cpeds_julian_time();
//	}
//	long YYYY,MM,DD,h,m;
//	double s;
//	cpeds_JDToYMDhms(JD,&YYYY,&MM,&DD,&h,&m,&s);
//	double yr=double(YYYY)+double(MM)/12+(double(DD)+double(h)+double(m)+s/60)/60))/24)/365.25+
//	return YYYY;	
//}
/****************************************************************************************************************/
double cpeds_julian_time() {
	
	struct timeval  tv;
	struct timezone tz;
	struct tm      *t;
	
	gettimeofday(&tv, &tz);
	t = gmtime(&tv.tv_sec);
	
	
	//	time_t now=time(NULL);
	//	tm *t;
	//	t=gmtime(&now);
	//  printf("y: %li, m: %li, d: %li, h: %li, m:%li, s: %li\n",long(t->tm_year),long(t->tm_mon), long(t->tm_mday),long(t->tm_hour),long(t->tm_min),long(t->tm_sec));
	return cpeds_julian_time(long(t->tm_year)+1900,long(t->tm_mon)+1, long(t->tm_mday),double(t->tm_hour)+double(t->tm_min)/60.0+double(t->tm_sec+tv.tv_usec/1e6)/3600.0);
}
/****************************************************************************************************************/
double cpeds_modified_julian_time(double JD) {
	return JD-2400000.5;
}

/****************************************************************************************************************/
void cpeds_JDToLocalZonalTime(double JD,long timeZone, double longitude, bool dst, long* YY, long* MM, long* DD, long* h, long* m, double* s) {
	// convert from local time to UTC time
	JD-=longitude/360.0;
	// convert to zonal time
	JD+=double(timeZone)/24.0;
	// correct the output for DST
	if (dst) JD+=1.0/24.0; // TODO: is this also valid for western longitudes ?
	// convert to date
	cpeds_JDToYMDhms(JD,YY,MM,DD,h,m,s);
	
}

/****************************************************************************************************************/
double cpeds_ZonalTimeToLocalJD(long YY, long MM, long DD, long h, long m, double s, long timeZone, double longitude, bool dst) {
	// convert to JD civil (zonal) time. This gives to the routine wrong time for the moment but below the time is corrected when it is converted 
	// to linear JD scale 
	double JD=cpeds_julian_local_time(YY,MM,DD,double(h)+double(m)/60.0+s/3600.0,0); //+double(timeZone)/24.0;
	/* double JD=cpeds_julian_local_time(YY,MM,DD,double(h)+double(m)/60.0+s/3600.0,0);  */
	
	// convert to winter civil time if needed
	if (dst) JD-=1.0/24.0; 
	
	// convert to UTC time
	JD-=double(timeZone)/24.0;  /* BLcomment (Nov 23, 2011, 1:24:59 PM): we need this here to correct for the fact that for 
  	  	  	  	  	  	  	  	  cpeds_julian_local_time conversion done above (which require time given for lon=0
  	  	  	  	  	  	  	  	  actually we gave local zonal time. So we must correct for this fact here.
	 */
	// BLmodification (Sep 22, 2011, 10:04:32 PM): uncommented. This should be here.
	
	// convert to local time in JD (i.e. derived for the geographical location of the observer)
	JD+=longitude/360.0;
	return JD;
}

/****************************************************************************************************************/
void cpeds_angToDMS(double ang, double *d, double *m, double *s, double acc) {
	*d=trunc(ang);
	*m=abs(trunc( (ang-(*d))*60.0 ));
	if (*m==60) { *d+=1; *m=0; }
	*s=abs( (abs((ang-(*d))*60.0) -(*m)) * 60.0 );
	*s=round(*s*acc)/acc;
	if (*s==60) { 
		*m+=1; *s=0; 
		if (*m==60) { *d+=1; *m=0; }
		if (*d>=360) *d-=360;
	}
	if (*d>=360) *d-=360;
	if (*d==0) { 
		if (*m==0) {
			if (ang<0) *s=fabs(*s);
		}
		else {
			if (ang<0) { *m=-fabs(*m); *s=fabs(*s); }		  
		}
	}
}
/***************************************************************************************/
void cpeds_angToHMS(double ang, double *h, double *m, double *s) {
	cpeds_angToDMS(ang/15.0,h,m,s);
}
double cpeds_DMSToAng(double d, double m, double s) {
	double deg;
	if (d==0) {
		if (m==0) {
			deg=s/3600.0;
			return deg;
		}
		if (m<0) deg=m/60.0-fabs(s)/3600.0; else deg=m/60.0+fabs(s)/3600.0; 
		return deg;
	}
	if (d<0) deg=d-fabs(m)/60.0-fabs(s)/3600.0; else deg=d+fabs(m)/60.0+fabs(s)/3600.0;
	return deg;
}
double cpeds_HMSToAng(double h, double m, double s) {
	return 15*(h+m/60.0+s/3600.0);
}
/***************************************************************************************/

double cpeds_zenithStarAzimuthAtRising(double lat) { return acos(tan(lat)); }

double cpeds_angleToHorizonAtRising(double A, double lat) { return acos( sin(lat) / sin( acos(cos(A)*cos(lat))) ); }
/***************************************************************************************/
double cpeds_nonplanar_atmospheric_layer_path(double zd, double Ri, double hi, double R0) {
	double Rj=Ri+hi;
	double r=0;
	
	if (zd==0) {
		r=hi;
	}
	else {
		zd*=PI180;
		double sinzd2=sin(zd)*sin(zd);
		r=sqrt( Ri*Ri+Rj*Rj-2.0*(R0*R0+sqrt(Ri*Ri/sinzd2-R0*R0)*sqrt(Rj*Rj/sinzd2-R0*R0))*sinzd2 );
	}
	return r;
}

/***************************************************************************************/
double cpeds_air_mass(double elev, bool fit) {
	double z=(90.0-elev)*PI180;
	double secz=1.0/cos(z);
	double airMass=0;
	double minElev=4;
	if (fit) {
		if (elev>=minElev) {
			airMass= -0.0045 + 1.00672*secz - 0.002234*secz*secz - 0.0006247*secz*secz*secz;		
		}
		else {
			return cpeds_extrapolate_linear(elev,minElev,cpeds_air_mass(minElev),minElev+0.1,cpeds_air_mass(minElev+0.1));
		}
	}
	else {
		airMass=0;
	}
	return airMass;
}
/***************************************************************************************/
double cpeds_K2eV(double temp) {
	return temp*CPEDS_kB/CPEDS_e;
}
double cpeds_K2eV() {
	return CPEDS_kB/CPEDS_e;	
}
double cpeds_K2keV() {
	return CPEDS_kB/CPEDS_e/1.0e3;
}
double cpeds_K2keV(double temp) {
	return temp*CPEDS_kB/CPEDS_e/1.0e3;
}
/***************************************************************************************/
double cpeds_dew_point(double T, double RH, double alpha, double beta, double lambda) {
	double Tdew;
	double tmp=log(RH/100) + beta*T/(lambda+T);
	Tdew=(lambda*tmp) / (beta - tmp);
	return Tdew;
}
/***************************************************************************************/
double cpeds_h2o_RH(double T, double Tdew, double alpha, double beta, double lambda) {
	double tmp=100.0*exp( -T*beta/(T+lambda) + Tdew*beta/(Tdew+lambda) );
	if (tmp>100.0) tmp=100;
	return tmp;
}

/***************************************************************************************/
double cpeds_h2o_Psat(double T) {
	double r1, r2, frac;
	double ln_Psat, ln_T;
	
	/*
	 * Log a warning if T is outside the applicable range of
	 * the Murphy-Koop formula.
	 */
	if (T < 123. || T > 332.) {
		printf("H2O_Psat>> temp. out of range\n");
		printf("the range is: [123,332] K. Requested %lE K\n",T);
		exit(-1);
	}
	
	/*
	 * The following formula gives log(Psat), with Psat in Pa.
	 */
	ln_T = log(T);
	r1 = 54.842763 - 6763.22 / T - 4.210   * ln_T + 0.000367 * T;
	r2 = 53.878    - 1331.22 / T - 9.44523 * ln_T + 0.014025 * T;
	ln_Psat = r1 + r2 * tanh(0.0415 * (T - 218.8));
	frac = 1;
	/*
	 * Convert Psat to mbar.
	 */
	return frac * 0.01 * exp(ln_Psat);	
}
/***************************************************************************************/
double cpeds_h2o_RH2vmr(double RH, double P, double T) {
	return 0.01 * RH * cpeds_h2o_Psat(T-CPEDS_0K) / P;
}
/***************************************************************************************/
double cpeds_h2o_RH(double T, double P, double vmr) {
	return vmr*100*P/cpeds_h2o_Psat(T);
}
/***************************************************************************************/

double cpeds_pressure_vs_altitude(double z, double P0, double T0, double z0, double Lr, double g0, double mu) {
	double R=8.31432;
	double P=P0 * pow(T0/(T0+Lr*(z-z0)),(g0 * mu/(R*Lr)));
	return P;
}
/***************************************************************************************/
double cpeds_radiance_to_temp_conversionFactor(double freq, double T0) {
	double x=CPEDS_h*freq/CPEDS_kB/T0;
	
	double factor=CPEDS_c*CPEDS_c/(2.0*freq*freq*CPEDS_kB)*(exp(x)-1.0)*(exp(x)-1.0)/(exp(x)*x*x);
	return factor*1.0e-26; // conversion to K/Jy
}
/***************************************************************************************/
double cpeds_black_body_radiance_Bnu(double freq, double T0) {
	double x=CPEDS_h*freq/CPEDS_kB/T0;
	//	double Bnu=
	return 2.0*CPEDS_h*freq*freq*freq/(CPEDS_c*CPEDS_c)/(exp(x)-1) / CPEDS_Jy;
}
/***************************************************************************************/
double cpeds_TSZEgnu_factor(double freq, double T0) {
	double x=CPEDS_h*freq/CPEDS_kB/T0;
	double gnu=x*cosh(x/2)/sinh(x/2)-4;
	return gnu;
}
/***************************************************************************************/
double cpeds_refraction(double ZDobs, double alt, double T, double P, double H, double lambda, double lat, double Tlapse, double acc) {
	double ref;
	T+=273.15; // convert to K
	ZDobs*=PI180; // convert to rad
	H/=100.0; // convert to range 0-1
	lambda*=1e4; // convert to um
	lat*=PI180; // convert to rad
	slarefro_(&ZDobs, &alt,&T,&P,&H,&lambda,&lat,&Tlapse, &acc,&ref);	
	ref=(ZDobs+ref)*PI180inv;
	return ref;
}

double ctgKB(double x) {
	return 2.908882e-4/tan(x+2.227e-3/(x+0.07679));
}
double cpeds_refractionKB(double h, int iObs) {
	
	if ((h < -0.0145444) || (iObs == 0)) return h;
	if (iObs < 0) {
		if (h < 0.) return h;
		double r = 1/tan(h+7.31/(h*180./M_PI+4.4)*M_PI/180.);
		return h-r/(60.*180.)*M_PI;
	}
	else {
		double r = ctgKB(h);
		return ctgKB(h+r)-1.745e-5*sin(882.*r+0.2269)+h;
	}
}

/***************************************************************************************/
double cpeds_calculateAngularResolution(double aperture, double freq) {
	double D=aperture;
	double lambda=CPEDS_c/(freq*1.0e9);
	double fwhm;
	//	fwhm=1.22*lambda/D * PI180inv; //deg
	
	fwhm=2.0*asin(1.80102*lambda/(PI*D)) * PI180inv;
	
	return fwhm;
}
/***************************************************************************************/

/******************************************************************************************/
/*                                                                                        */
/******************************************************************************************/
double RT4_PrecessionRA(double juda0,double juda,double alpha0,double delta0){
	double 	r2d = 180.0/M_PI; 
	double 	d2r = M_PI/180.0; 
	
	double T0, T;
	double dzetaA,zedA,thetaA;
	double tmpRA;
	T0 = (juda0 - 2451545.0) / 36525.0;
	T = (juda - juda0) / 36525.0;
	dzetaA = (2306.218 + 1.397 * T0 + 0.302 * T + 0.018 * T * T) * T / 3600.0;
	zedA = dzetaA + 0.00022028 * T * T;
	thetaA = (2004.311 - 0.853 * T0 - 0.427 * T - 0.042 * T * T) * T / 3600.0;
	tmpRA = zedA + r2d * atan2( sin( d2r * (alpha0 + dzetaA) ),
			cos(d2r * thetaA) * cos( d2r * (alpha0 + dzetaA) ) -
			tan(d2r * delta0) * sin(d2r * thetaA) );
	if (tmpRA < 0.0)
		tmpRA += 360.0;
	return tmpRA;
}

double RT4_PrecessionDEC(double juda0,double juda,double alpha0,double delta0) {
	double 	r2d = 180.0/M_PI; 
	double 	d2r = M_PI/180.0; 
	double T0, T;
	double dzetaA,thetaA;
	T0 = (juda0 - 2451545.0) / 36525.0;
	T = (juda - juda0) / 36525.0;
	dzetaA = (2306.218 + 1.397 * T0 + 0.302 * T + 0.018 * T * T) * T / 3600.0;
	thetaA = (2004.311 - 0.853 * T0 - 0.427 * T - 0.042 * T * T) * T / 3600.0;
	return r2d * asin( sin(d2r * thetaA) * cos(d2r * delta0) *
			cos( d2r * (alpha0 + dzetaA) ) + cos(d2r * thetaA) *
			sin(d2r * delta0) );
}

rt4_time_date cpeds_RT4_cal_date(double tjd) {
	long int jd, k, m, n;
	rt4_time_date now;
	double djd;
	
	short int month, day, year;
	double hour;
	djd = tjd + 0.5;
	jd = (long int) djd;
	
	hour = fmod (djd,1.0) * 24.0;
	
	k = jd + 68569L;
	n = 4L * k / 146097L;
	
	k = k - (146097L * n + 3L) / 4L;
	m = 4000L * (k + 1L) / 1461001L;
	k = k - 1461L * m / 4L + 31L;
	
	month = (short int) (80L * k / 2447L);
	day = (short int) (k - 2447L * (long int) month / 80L);
	k = (long int) month / 11L;
	
	month = (short int) ((long int) month + 2L - 12L * k);
	year = (short int) (100L * (n - 49L) + m + k);
	
	now.min = (int)( (hour - (int)hour) * 60.0 );
	now.sec = (int)( ( ( (hour - (int)hour) * 60.0) - now.min ) 
			* 60.0 );
	now.usec = ( ( ( ( ( (hour - (int)hour) * 60.0) - 
			now.min ) * 60.0 ) - now.sec ) * 1000000.0 );
	now.hour = (int)(hour);
	now.day = day;
	now.month = month;
	now.year = year;
	
	return now;
}


double cpeds_rt4_JulianDay(int tye,int tmo,int tda) {
	long i,j;
	j = tye + (tmo - 9) / 7;
	i = tmo + 10 - ((tmo + 9) / 12) * 12;
	return (double)( 1471086l + 365 * j + (j + 1000004l) / 4 + 30 * i + 
			i / 2 + ( (i / 7) * i ) % 2 + i / 12 + tda + (j + 1000400l) / 
			400 - (j + 1000100l) / 100 + 7502 );
};

void cpeds_rt4_nutation(double DJ,double *dRA,double *dDEC, double* dPsi, double* dEps,double* eObliquity){
	double RA=*dRA; // deg
	double DEC=*dDEC; //deg
	double pi=3.14159265359;
	double sec2rad=4.848136811095e-6;
	double deg2rad=0.0174532925199;
	double ddj,SL,aM,Om,dl2,sinE,cosE,dpsi,deps;
	double eL,cosL,sinL,dRAa,dDECa;
	double cosRA=cos(deg2rad*RA);
	double sinRA=sin(deg2rad*RA);
	double cosDEC=cos(deg2rad*DEC);
	double tanDEC=tan(deg2rad*DEC);
	double mObliq=0.0;
	
	//*dRA = 0.0;
	//*dDEC = 0.0;
	
	if(cosDEC<1e-10){
		cosDEC=1e-10;
	}
	
	ddj = DJ-2451545.0;
	SL = (280.460+0.9856474*ddj)*deg2rad;		/* mean longitude of Sun */
	aM = (357.528+0.9856003*ddj)*deg2rad;		/* mean anomaly */
	Om =(47.8-0.05295376*(DJ-2453004.5))*deg2rad;	/* Omega */
	
	/* twice longitude of Moon:*/
	dl2 = fmod(218.316+481267.881*ddj/36525.0,180.0)*pi/90.0;
	//	printf("ddj: %lE\n",ddj);
	//	mObliq=(23.439-4.0-7.0*ddj); //deg
	mObliq=(23.439-3.56e-7*ddj); //deg
	//	mObliq=84373.70/3600; // deg
	
	if (eObliquity!=NULL) { *eObliquity=mObliq*3600; }
	sinE = sin(mObliq*deg2rad); /* sin(ecliptic_obliquity) */
	//	sinE = sin((23.439-4.0-7.0*ddj)*deg2rad); /* sin(ecliptic_obliquity) */
	cosE = sqrt(1.0-sinE*sinE);
	
	/* Compute nutation (precision ~0.04") */
	dpsi=(-17.2*sin(Om) - 1.319*sin(SL+SL) + 0.206*sin(2.0*Om) + 0.143*sin(aM) - 0.227*sin(dl2))*sec2rad;
	deps=(+9.203*cos(Om) + 0.574*cos(SL+SL) + 0.098*cos(dl2) - 0.09*cos(2.0*Om))*sec2rad;
	if (dPsi!=NULL) *dPsi=dpsi/sec2rad;
	if (dEps!=NULL) *dEps=deps/sec2rad;
	//	printf("dpsi: %lE, deps: %lE\n",dpsi*PI180inv,deps*PI180inv);
	RA += ( (cosE + sinE*sinRA*tanDEC)*dpsi - cosRA*tanDEC*deps )/deg2rad;
	DEC += ( sinE*cosRA*dpsi + sinRA*deps )/deg2rad;
	
	*dRA=RA; // convert to deg
	*dDEC=DEC; // convert to deg
	
}

/** nutation
 * 
 * Obliczanie poprawek wspolrzednych zrodla na efekty NUTACJI i
 * ABERRACJI rocznej. Na koncu pokazano jak dolozyc ABERRACJE DOBOWA
 *  - trzeba dostarczyc kat godzinny HA, lub lokalny czas gwiazdowy ST 
 * w radianach i obliczyc HA = ST - RA. Ten efekt ma amplitude ok. 0.2"
 * Wspolrzedne zasadniczo powinny byc wczesniej przeprecesowane na 
 * date DJ. Poprawione wspolrzedne to: RAnew=RA+dRA, DECnew=DEC+dDEC.
 *
 * Computes correction to equatorial coordinates due to effects of
 * aberation and nutation
 * Tests have shown that in the years 2005-2056 errors of these
 * approximations are smaller than 0.3"
 * @param DJ   -  Julian day
 * @param dRA  -  (corrected) right ascension (in degrees)
 * @param dDEC -  (corrected) declination (in degrees)
 *!!! HA   - local Hour Angle (in radians) <- add this argument and
 *            uncomment last but one line in this subprogram !!!)
 * @author K. Borkowski
 * @version 23.06.2006 
 * 
 * BLcomment (Jun 13, 2016, 1:52:49 PM): The input for this routine 
 * (dRa and dDEC) is in degrees.
 * The output (dRa and dDEC) is now also in degrees.
 * 
 * 
 */
void cpeds_rt4_nutation_and_aberration(double DJ,double *dRA,double *dDEC){
	double RA=*dRA;
	double DEC=*dDEC;
	double  pi=3.14159265359;
	double sec2rad=4.848136811095e-6;
	double deg2rad=0.0174532925199;
	double ddj,SL,aM,Om,dl2,sinE,cosE,dpsi,deps;
	double eL,cosL,sinL,dRAa,dDECa;
	double sinRA=sin(deg2rad*RA);
	double cosRA=cos(deg2rad*RA);
	double sinDEC=sin(deg2rad*DEC);
	double cosDEC=cos(deg2rad*DEC);
	double tanDEC=tan(deg2rad*DEC);
	double mObliq=0.0;
	
	//*dRA = 0.0;
	//*dDEC = 0.0;
	
	if(cosDEC<1e-10){
		cosDEC=1e-10;
	}
	
	ddj = DJ-2451545.0;
	SL = (280.460+0.9856474*ddj)*deg2rad;		/* mean longitude of Sun */
	aM = (357.528+0.9856003*ddj)*deg2rad;		/* mean anomaly */
	Om =(47.8-0.05295376*(DJ-2453004.5))*deg2rad;	/* Omega */
	
	/* twice longitude of Moon:*/
	dl2 = fmod(218.316+481267.881*ddj/36525.0,180.0)*pi/90.0;
	
	mObliq=(23.439-3.56e-7*ddj); //deg
	//	mObliq=(23.439-4.0E-7*ddj); //deg
	//	mObliq=(23.439-7.0*ddj/36525); //deg
	//	mObliq=8.437370E+04/3600; // deg
	
	sinE = sin(mObliq*deg2rad); /* sin(ecliptic_obliquity) */
	//	sinE = sin((23.439-4.0-7.0*ddj)*deg2rad); /* sin(ecliptic_obliquity) */
	cosE = sqrt(1.0-sinE*sinE);
	
	/* Compute nutation (precision ~0.04") */
	dpsi=(-17.2*sin(Om) - 1.319*sin(SL+SL) + 0.206*sin(2.0*Om) + 0.143*sin(aM) - 0.227*sin(dl2))*sec2rad;
	deps=(+9.203*cos(Om) + 0.574*cos(SL+SL) + 0.098*cos(dl2) - 0.09*cos(2.0*Om))*sec2rad; //deg
	*dRA += ( (cosE + sinE*sinRA*tanDEC)*dpsi - cosRA*tanDEC*deps ) / deg2rad; //deg
	*dDEC += ( sinE*cosRA*dpsi + sinRA*deps ) / deg2rad;
	
	/* Now compute annual aberration of light  (precision ~0.3") */
	eL = SL + (1.915*sin(aM)+0.02*sin(2.0*aM))*deg2rad; /* ecliptic longitude [rad] */
	cosL = cos(eL);
	sinL = sin(eL);
	dRAa = -20.49*sec2rad/cosDEC*(sinL*sinRA + cosL*cosRA*cosE); // [rad]
	dDECa= -20.49*sec2rad*(sinL*sinDEC*cosRA + cosL*(sinE*cosDEC - cosE*sinDEC*sinRA)); // [rad]
	
	/* Add elliptical aberration (11.25*pi/12 = 2.945243112740)*/
	*dRA += ( dRAa - sec2rad*0.341*sin(RA*deg2rad+2.94524311274)/cosDEC ) / deg2rad; // deg
	*dDEC += ( dDECa - sec2rad*(0.341*cos(RA*deg2rad+2.94524311274)*sinDEC + 0.029*cosDEC) ) / deg2rad; // deg
	
	/* Add diurnal aberration. For RT32: 0.32"*cos(phi') = ~9.34e-7 rad
cSupply Hour Angle in radians to cos(HA) and uncomment the following line

	 *dRA += + 9.34e-7/cosDEC*cos(HA);

	 */
	
	
}



/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/

/****************************************************************************************************************/
// remember to initialize the srand(time(0)) for this function
double* cpeds_random_uniform_numbers(double min, double max, long pix_num) {
	long int i;
	double *uvec = new double[pix_num];
	for (i=0;i<pix_num;i++) { uvec[i] = cpeds_random_uniform_number(min,max); }
	return uvec;
}
/****************************************************************************************************************/
// remember to initialize the srand(time(0)) for this function
double cpeds_random_uniform_number(double min, double max) {
	double tmp;
	if (min == max) { return min; } else
		do {    tmp=min+(max-min)*(double)(rand())/(double)(RAND_MAX); }  while ((tmp > max) || (tmp < min));
	
	/*  do {
    tmp=min+max/(RAND_MAX+min)*double(rg.Random()); //printf("%f \n",tmp);
  }  while ((tmp > max) || (tmp < min));
	 */
	return tmp;
}

/****************************************************************************************************************/
gsl_rng* cpeds_random_uniform_numbersGSL_init(long* seed, long seed_offset, const gsl_rng_type* generator) {
	gsl_rng* r;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = generator;
	r = gsl_rng_alloc (T);
#ifdef DEBUG
	printf("*** initiating rng stuff, using seed %li\n",*seed);
#endif
	if (*seed==0) {
		*seed=time(0);
		gsl_rng_set(r,*seed+seed_offset);
#ifdef DEBUG
		printf("	*** will use time based seed %li\n",*seed);
#endif
	}
	else gsl_rng_set(r,*seed+seed_offset);
#ifdef DEBUG
	printf("*** seed+seed_offset is: %li\n",*seed+seed_offset);
#endif
	return r;
}
/****************************************************************************************************************/
/* ! returns array of size pix_num of random uniform numbers drawn from GSL DIEHARD survivor RNG */
/* ! with values from min to max and initial seed taken from the system time offshifted by the seed_offset in case of */
/* ! parallel implementations/ runs whatever,  the fast option is now NOT IMPLEMENTED  -- will be used to indicate  */
/* ! a very freuqent calls to this routine (faster than once per second) to avoid taking the same seeds from the system time. */

double* cpeds_random_uniform_numbersGSL(double min, double max, long pix_num,long seed_offset, bool fast) {
	const gsl_rng_type * T;
	gsl_rng * r;
	
	unsigned long int rng_seed = (unsigned long int)time(NULL)+cpeds_seed_offset;
	double * tmp=new double[pix_num];
	double delta=max-min;
	long i;
	
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,rng_seed+seed_offset);
	
	for (i=0;i<pix_num;i++) { tmp[i]=  gsl_rng_uniform(r) * delta + min; }
	gsl_rng_free(r);
	return tmp;
}
/****************************************************************************************************************/
gsl_rng* cpeds_random_uniform_numbersGSL(double min, double max, long num, double **t, gsl_rng* rng) {
	*t=new double[num];
	double delta=max-min;
	long i;
	
	for (i=0;i<num;i++) { (*t)[i]=  gsl_rng_uniform(rng) * delta + min; }
	/*   gsl_rng * r=gsl_rng_clone(rng); */
	/*   gsl_rng_free(rng); */
	/*   return r; */
	return rng;
}
/****************************************************************************************************************/
double cpeds_random_uniform_numberGSL(double min, double max, gsl_rng* rng) {
	double t;
	double delta=max-min;
	t=  gsl_rng_uniform(rng) * delta + min;
	return t;
}

/****************************************************************************************************************/
double* cpeds_random_uniform_numbersGSL(double min, double max, long num, gsl_rng* rng) {
	double *t=new double[num];
	double delta=max-min;
	long i;
	
	for (i=0;i<num;i++) { t[i]=  gsl_rng_uniform(rng) * delta + min; }
	/*   gsl_rng * r=gsl_rng_clone(rng); */
	/*   gsl_rng_free(rng); */
	return t;
}

/****************************************************************************************************************/
// remember to initialize the srand(time(0)) for this function
// returns a sigle random number from a gaussian distribution of mean m and variance s tabulted with pix_num pixles
// for central limit theorem mehtod pix_num is the number of uniform numbers with which the gaussian number is drawn
double cpeds_random_gauss_number(double m,double s, long int pix_num, int method) {
	double *gCDF;
	double a,x,xp,np,delta,deltamin; // p - wartosci  primowane dotycza rozkladu redshiftow; nie primowane rozkladu plaskiego
	long int i;
	//gsl_rng * r;
	
	// this is totaly old slow and crapy. inverrse CDF method with linear search and no interpolation. god damn  !!!! - this should be removed or developed
	if (method == 1) {
		//srand(time(0));
		gCDF = cpeds_generate_gaussian_distribuant_function(pix_num, s, m);
		x = cpeds_random_uniform_number(0,1);
		
		np=0; xp=0;
		deltamin=2; // for normalized distribution
		for (i=0;i<pix_num;i++) {
			delta=fabs(gCDF[cpeds_xy2num(i,1)]-x);
			if (delta < deltamin) { np = gCDF[cpeds_xy2num(i,0)]; deltamin=delta; xp=gCDF[cpeds_xy2num(i,1)];}
			//printf("*** z=%f n=%f zp=%f np=%f i=%i Nmax=%i ndistrZNint=%f delta=%f deltamin=%f****\n",z,n,zp,np,i,Nmax,ndistrZNint[i][1],delta,deltamin);
		}
		
		// printf("*** z=%f n=%f zp=%f np=%f ****\n",z,n,zp,np);
		//  printf("-----------------------------------------------------\n");
		delete gCDF;
	}
	
	// this method gives a single gaussian variate value based on the central limit theorem. pix_num is the number of points summed for one value.
	if (method == 2) {
		
		/*     if (pix_num == 10)  a = 5.47456*s; else // this is for 10 numbers */
		/*     if (pix_num == 20)  a = 7.74659*s; else // this is for 20 numbers : the requested variance is fine; the skewness is fine but kurtosis is a bit biased: 0.95 @ 10^5 call statistics */
		/*     if (pix_num == 40)  a = 10.9586*s; else // this is for 40 numbers : the kurtosis is still less than 3 */
		/*     if (pix_num == 100) a = 17.3221*s; else // this is for 100 numbers: slower but the accuracy is fine for the kurtosis */
		/*       a=1.73178*sqrt((double)pix_num)*s; // interpolation formula for any arbitrary number of requested Gaussian numbers -- shouldn't this be sqrt(3) ? */
		a=1.732050757*sqrt((double)pix_num)*s; // interpolation formula for any arbitrary number of requested Gaussian numbers
		np=0;
		for (i=0;i<pix_num;i++) { np+=cpeds_random_uniform_number(-a,a); }
		np/=(double)pix_num;    np+=m;
	}
	
	// this is just as in the method 2 but the uniform generator is the gsl default generator
	if (method == 3) {
	}
	
	// this is for the internal use only; for testing purposes only; (allowes to steer the a value to find the proper arange through m value. the real m is 0 here by default)
	if (method == 5) {
		a=m;
		pix_num = 40;    np=0;
		for (i=0;i<pix_num;i++) { np+=cpeds_random_uniform_number(-a,a); }
		np/=(double)pix_num;
	}
	
	return np;
	//return gauss_distr[k][1];
}
/****************************************************************************************************************/
double cpeds_random_gauss_number(double m,double s, long int num, gsl_rng* rng) {
	double a=sqrt(3.0*(double)num)*s;
	double *t=cpeds_random_uniform_numbersGSL(-a,a,num,rng);
	m=cpeds_mean_value(t,num)+m;
	delete [] t;
	return m;
}
/****************************************************************************************************************/
double* cpeds_random_gauss_numbers(double m,double s, long n, long int num, gsl_rng* rng) {
	double *t = new double[n];
	for (long i=0;i<n;i++) { t[i]=cpeds_random_gauss_number(m,s,num, rng); }
	return t;
}
/****************************************************************************************************************/
// returns a pix_num gaussian random numbers  of mean m and variance s tabulted with pix_num pixles
// the return value is the pointer to the table of pix_num double numbers
// the function has the ACCURACY parameter that controls so for the small number of nubers the gaussian CDF was accurate enough
// method  - 1 inversion of the CDF -- slow
// method  = 2 ivnersion of the CDF but search algorithm is improved (division by 2)
// method = 3  uses the uniform gsl generator but performs inverse gaussCDF transformaiton to get gaussian distr.
// method = 4  gsl gaussian distribution
// method = 5 uses the default GSL uniform generator and central limit theorem to make gaussian numbers from 20 uniform numbers
// remember to initialize the srand(time(0)) for this function
// after finding the closest mathing point the linear interpolation is performed for the requested point between the adjacent grid points in iversion process
//double* cpeds_random_gauss_numbers(double m, double s,long int pix_num, int method) { return cpeds_random_gauss_numbers(m,s,pix_num,method,(long)0,false); }
double* cpeds_random_gauss_numbers(double m, double s,long int pix_num, int method, long seed_offset, bool sims_fast) {
	double *gCDF,*gtab;
	double x,xp,delta,deltamin; // p - wartosci  primowane dotycza rozkladu redshiftow; nie primowane rozkladu plaskiego
	long int i,j,np=0, gCDFtab=2*pix_num; // the size of the tabulating array is to be twice as large as the number of points drown from this distribution
	long int ACCURACY = 10000; // minimal number of points that tabulate the gaussian CDF. if you request eg. 1 point then the CDF will be tabulated with this number of points
	long int istart=0,iend=gCDFtab;
	//double accuracy;
	double xi,yi,xf,yf; // interpolation stuff !!! dont confuse the random x with these variables - these are totaly independent and in the
	// linear interpolation the y of the linear function is in fact the random uniform x :)
	//srand(time(0));
	const gsl_rng_type * T;
	gsl_rng * r;
	unsigned long int rng_seed;
	long pix_num_total;
	double npd,a;
	FILE *f;
	
	gtab = new double[pix_num];
	if ((method == 1) || (method == 2 ) || (method == 3 )) { // inverse gaussCDF methods
		if (gCDFtab < ACCURACY) { gCDFtab = ACCURACY; }
		gCDF = cpeds_generate_gaussian_distribuant_function(gCDFtab, s, m);
		
		if (method == 3) {
			gsl_rng_env_setup();
			T = gsl_rng_default;
			r = gsl_rng_alloc (T);
			
			if (sims_fast) {
				//----- patch to avoid same generations if called in the same second TODO: do this in some object, cause this is stupid
				rng_seed=(unsigned long)(time(NULL)/1000);
				system("date +%N |colrm 4 10 > .mk_gauss_map.tmpfile");
				f=fopen(".mk_gauss_map.tmpfile","r"); fscanf(f,"%li\n",&i); fclose(f);
				rng_seed=rng_seed*1000+i;
				//----
			}
			else { rng_seed = (unsigned long int)time(NULL); }
			
			printf("rng seed: %li\n",rng_seed);
			gsl_rng_set(r,rng_seed+seed_offset);
		}
		
		for (j=0;j<pix_num;j++) {
			if ((method == 1) || (method == 2 )) { x = cpeds_random_uniform_number(0,1); }
			if (method == 3) { x = gsl_rng_uniform(r); }
			//printf("wylosowana liczna1: %lf %li %li\n",x,j,pix_num);
			
			if (method == 1) {  // inverse gaussCDF method for system RNG -- slow
				deltamin = 2;
				// this is the old version - ok. but slow
				istart = 0; iend = gCDFtab;
				if (x < 0.5) { istart = 0; iend = (long int)ceil((double)gCDFtab/2)+1; } // optimalization stuff only
				if (x >= 0.5) { istart = (long int)trunc((double)gCDFtab/2)-1; iend = gCDFtab; }  // optimalization stuff only
				for (i=istart;i<iend;i++) {
					delta=fabs(gCDF[cpeds_xy2num(i,1)]-x);
					
					if (delta < deltamin) { np = (long int)i; deltamin=delta; xp=gCDF[cpeds_xy2num(i,1)];} // this xp was kept redundandly
					//printf("i=%li x = %lf, xp = %lf  delta = %lf deltamin = %lf istart = %li, iend = %li, np = %li \n",i,x,xp,delta,deltamin,istart,iend,np);
				}
				// this is the old version - ok. but slow
				// printf("*** z=%f n=%f zp=%f np=%f ****\n",z,n,zp,np);
				//printf("-----------------------------------------------------\n");
			}
			
			
			if ((method == 2) || (method == 3)) { // inverse gaussCDF method for system and GSL RNGs respectively
				istart = 0; iend = gCDFtab;
				np = (long int)round(((double)iend - (double)istart)/2); xp=gCDF[cpeds_xy2num(np,1)];
				//accuracy = fabs(gCDF[cpeds_xy2num(np,1)]-gCDF[cpeds_xy2num(np-1,1)]); // this is the worst point in the middle of the range of the gCDF
				delta=xp-x;
				//printf("***accuracy = %lf, x = %lf, xp = %lf  delta = %lf deltamin = %lf istart = %li, iend = %li, np = %li \n",accuracy, x,xp,delta,deltamin,istart,iend,np);
				//while (fabs(delta) >= accuracy ) {
				while (np+1 != iend ) {
					if (delta < 0) { istart = np; }
					if (delta > 0) { iend = np; }
					np = (int)round((double)istart+((double)iend - (double)istart)/2); xp=gCDF[cpeds_xy2num(np,1)];
					delta = xp - x;
					//accuracy = fabs(gCDF[cpeds_xy2num(np,1)]-gCDF[cpeds_xy2num(np-1,1)]); // this is the worst point in the middle of the range of the gCDF
					
					//printf("accuracy = %lf,x = %lf, xp = %lf  delta = %lf deltamin = %lf istart = %li, iend = %li, np = %li \n",accuracy,x,xp,delta,deltamin,istart,iend,np);
				}
			}
			
			// linear interpolation between the selected points
			if (x > xp) {xi = gCDF[cpeds_xy2num(np,0)]; xf = gCDF[cpeds_xy2num(np+1,0)]; yi = gCDF[cpeds_xy2num(np,1)]; yf = gCDF[cpeds_xy2num(np+1,1)];}
			if (x < xp) {xi = gCDF[cpeds_xy2num(np-1,0)]; xf = gCDF[cpeds_xy2num(np,0)]; yi = gCDF[cpeds_xy2num(np-1,1)]; yf = gCDF[cpeds_xy2num(np,1)];}
			
			gtab[j] = (x-yi)*(xf-xi)/(yf-yi) + xi;
			/*     gtab[j] = gCDF[cpeds_xy2num(np,0)]; // this is if you want without interpolation; in that case uncomment the accuracy stuff as well above */
			//printf("wylosowana liczna2: x = %.9lf gauss x = %.9lf\n",x, gtab[j]);
		}
		//return gauss_distr[k][1]
		if (method == 3) { gsl_rng_free(r); }
		
		delete gCDF;
	}
	if (method == 4) {
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		gsl_rng_set(r,(unsigned long int)(time(NULL))+seed_offset);
		//gsl_rng_default_seed = (unsigned long int)0;
		//printf ("generator type: %s\n", gsl_rng_name (r));
		//printf("seed = %li\n",gsl_rng_default_seed);
		//printf("seed = %.15lE\n",(double)time(NULL));
		for (i = 0; i < pix_num; i++) {
			gtab[i] = gsl_ran_gaussian(r,s) + m;
		}
		gsl_rng_free (r);
	}
	
	if (method == 5) {
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		rng_seed = (unsigned long int)time(NULL);
		printf("rng seed: %li\n",rng_seed);
		gsl_rng_set(r,rng_seed+seed_offset);
		
		pix_num_total = pix_num;
		pix_num = 20;    npd=0; a = sqrt((double)(pix_num)*s); // this is for 20 numbers
		for (i = 0; i < pix_num_total; i++) {
			for (j = 0; j < pix_num; j++) { npd += (gsl_rng_uniform(r)-0.5)*2*a;  }
			npd/=(double)pix_num; npd+=m;
			gtab[i] = npd;
		}
		gsl_rng_free (r);
	}
	
	return gtab;
}
/****************************************************************************************************************/
bool cpeds_positive(double d) {
	if (d>=0) return true;
	return false;
}
/****************************************************************************************************************/
long int  cpeds_xy2num(long int x, long int y) { // y can only by 0 or 1 -- this is conversion from linear table to rectangular table of size nx2
	return 2*x+y;
}
/****************************************************************************************************************/
// generates the requested gaussian cumulative probability distribution tabulated with size pixels with variance s and mean x0 and
// returns the address to the table structure of double numbers
// the CDF is stored in a table twice as big as size so the ordering is: 1 double number - x1, 2 double number - CDF(x1) 3 double number is x2 and so on
double* cpeds_generate_gaussian_distribuant_function(long int size, double s, double x0) {
	double nsigma = 5.5; // this is rather sufficient for npix ~ 2.6e7 pixels so for eg. ns=1024
	double x,dx,sf, Nfactor,A;
	double twos2=2.0*s*s;
	double xlex0;
	//double min,max;
	//char tmpch[100];
	long int  i;
	//FILE *f;
	double* cpeds_gaussCDF = new double[2*size];
	
	//sprintf(tmpch,"gassuian_distribuant_function-%s",file_ending);
	//f= fopen("gCDF","w");
	//  x=random_unmber(x0-nsigma*s,x0+nsigma*s);
	
	//min= x0-nsigma*s; max=x0+nsigma*s;
	
	// calculating gauss normalization factor
	/*   x = x0-nsigma*s; dx = 2*nsigma*s/(double)size; */
	/*   for (i=0;i<size;i++) { */
	/*     x=x+dx; */
	/*     Nfactor = Nfactor+ 1/(sqrt(2*PI)*s)*exp(-(x-x0)*(x-x0)/(2*s*s)); */
	/*   } */
	/*   Nfactor = Nfactor*dx; */
	
	x = x0-nsigma*s; dx = 2*nsigma*s/(double)size;
	cpeds_gaussCDF[cpeds_xy2num(0,0)] = 0;    cpeds_gaussCDF[cpeds_xy2num(0,1)] = 0;
	
	Nfactor = 1;
	A = 1/(sqrt(twoPI)*s*Nfactor);
	for (i=0;i<size;i++) {
		cpeds_gaussCDF[cpeds_xy2num(i,0)] = x;
		xlex0=x-x0;
		cpeds_gaussCDF[cpeds_xy2num(i,1)] = cpeds_gaussCDF[cpeds_xy2num(i-1,1)]+A*exp(-xlex0*xlex0/twos2);
		x=x+dx;
	}
	
	// normalization of the CDF
	sf = cpeds_gaussCDF[cpeds_xy2num(size-1,1)];
	for (i=0;i<size;i++) {
		cpeds_gaussCDF[cpeds_xy2num(i,1)] = cpeds_gaussCDF[cpeds_xy2num(i,1)]/sf;
		//fprintf(f,"%i %.9lf %.9lf\n",i,cpeds_gaussCDF[cpeds_xy2num(i,0)],cpeds_gaussCDF[cpeds_xy2num(i,1)]);
	}
	//fclose(f);
	return cpeds_gaussCDF;
	
}
/****************************************************************************************************************/
void cpeds_random_uniform_directions_on_sphere(long pix_num, double **th, double **phi, double fromPhi, double toPhi, double fromTh, double toTh) {
	double tmp;
	long j;
	if (fromTh > toTh) { tmp = fromTh; fromTh=toTh; toTh=tmp; }
	if (fromPhi > toPhi) { tmp = fromPhi; fromPhi=toPhi; toPhi=tmp; }
	
	*phi = cpeds_random_uniform_numbersGSL(fromPhi,toPhi,pix_num,0,false);
	*th = cpeds_random_uniform_numbersGSL(cos(toTh),cos(fromTh),pix_num,10,false);
	
	for (j=0;j<pix_num;j++) { (*th)[j]=acos((*th)[j]); } //http://mathworld.wolfram.com/SpherePointPicking.html
}

/****************************************************************************************************************/
void cpeds_random_uniform_directions_on_sphere(long pix_num, double **th, double **phi ) {
	*phi = cpeds_random_uniform_numbersGSL(0,twoPI,pix_num,0,false);
	*th = cpeds_random_uniform_numbersGSL(0,1,pix_num,10,false);
	long j;
	
	
	for (j=0;j<pix_num;j++) { (*th)[j]=acos(2*(*th)[j]-1); } //http://mathworld.wolfram.com/SpherePointPicking.html
}

/****************************************************************************************************************/

struct LongDouble {
		unsigned char mant[8];
		unsigned int exp;
};

union Real {
		struct LongDouble sld;
		long double ld;
};

int cpeds_isnan (double x) {
	union Real real;
	real.ld = x;
	if((real.sld.exp &0x7FFF) == 0x7FFF)  return 1;
	return 0;
}

/****************************************************************************************************************/
long cpeds_get_min(long v1, long v2) {
	if (v1 < v2) return v1; else return v2;
}

long cpeds_get_min(long v1, long v2, long* i) {
	if (v1 < v2) { *i=1; return v1;} else { *i=2; return v2; }
}

long cpeds_get_max(long v1, long v2) {
	if (v1 > v2) return v1; else return v2;
}

long cpeds_get_max(long v1, long v2, long* i) {
	if (v1 > v2) { *i=1; return v1;} else { *i=2; return v2; }
}


double cpeds_get_max(double v1, double v2) {
	if (v1 > v2) return v1; else return v2;
}
double cpeds_get_max(double v1, double v2, long* i) {
	if (v1 > v2) { *i=1; return v1; }  else { *i=2; return v2; }
}

double cpeds_get_min(double v1, double v2) {
	if (v1 < v2) return v1; else return v2;
}

double cpeds_get_min(double v1, double v2, long* i) {
	if (v1 < v2) { *i=1; return v1; } else { *i=2; return v2; }
}

/****************************************************************************************************************/
// returns the logaritm with base 2 of x
double cpeds_log2(double x) {
	return log10(x)/log10(2.0);
}

/****************************************************************************************************************/
// returns the factorial value of an argument
double cpeds_factorial(long int n) {
	double fact;
	long int i;
	
	if (n == 0) { fact = 1; }
	else {
		fact = 1;
		for (i = 1;i<=n;i++) { fact = fact*(double)i; }
	}
	
	return fact;
}

// returns the fractional factorial value of an argument: i.e. k*(k+1)*...*n (k<=n)
double cpeds_factorial_frac(long int n, long k) {
	double fact;
	long int i;
	
	if (n == 0) { fact = 1; }
	else {
		fact = 1;
		for (i = k;i<=n;i++) { fact = fact*(double)i; }
	}
	
	return fact;
}

long double cpeds_factorial_ld(long int n) {
	long double fact;
	long int i;
	
	if (n == 0) { fact = 1; }
	else {
		fact = 1;
		for (i = 1;i<=n;i++) { fact = fact*(long double)i; }
	}
	
	return fact;
}

/* bigint cpeds_factorial_lidia(long int n) { */
/*   bigint fact; */
/*   long i; */

/*   if (n == 0) { fact = 1; }  */
/*   else { */
/*     fact = 1; */
/*     for (i = 1;i<=n;i++) { fact = fact*i; } */
/*   } */

/*   return fact; */
/* } */

// returns the fractional factorial value of an argument: i.e. k*(k+1)*...*n (k<=n)
long double cpeds_factorial_frac_ld(long int n, long k) {
	long double fact;
	long int i;
	
	if (n == 0) { fact = 1; }
	else {
		fact = 1;
		for (i = k;i<=n;i++) { fact = fact*(long double)i; }
	}
	
	return fact;
}

/* bigint cpeds_factorial_frac_lidia(long int n, long k) { */
/*   bigint fact; */
/*   long int i; */

/*   if (n == 0) { fact = 1; }  */
/*   else { */
/*     fact = 1; */
/*     for (i = k;i<=n;i++) { fact = fact*i; } */
/*   } */

/*   return fact; */
/* } */

/****************************************************************************************************************/
void cpeds_find_minmax_value(long *a, long int size, long *min, long *max, long int *mini, long int *maxi) {
	long int i,lmini,lmaxi;
	long lmax, lmin;
	
	lmin = a[0];
	lmax = a[0];
	lmini = 0; lmaxi = 0;
	for (i=0; i<size;i++) {
		if ( a[i] < lmin ) { lmin = a[i]; lmini = i; }
		if ( a[i] > lmax ) { lmax = a[i]; lmaxi = i; }
	}
	*mini = lmini; *maxi = lmaxi; *min = lmin; *max = lmax;
}

/****************************************************************************************************************/
void cpeds_find_minmax_value(double *a, long int size, double *min, double *max, long int *mini, long int *maxi) {
	long int i,lmini,lmaxi;
	double lmax, lmin;
	
	lmin = a[0];
	lmax = a[0];
	lmini = 0; lmaxi = 0;
	for (i=0; i<size;i++) {
		if ( a[i] < lmin ) { lmin = a[i]; lmini = i; }
		if ( a[i] > lmax ) { lmax = a[i]; lmaxi = i; }
	}
	if (mini!=NULL) *mini = lmini; 
	if (maxi!=NULL) *maxi = lmaxi;
	*min = lmin; *max = lmax;
	/*   printf("!!!size = %li mini=%li, maxi=  %li, min = %lE, max = %lE\n",size,lmini,lmaxi,lmin,lmax); */
	/*   printf("!!!size = %li mini=%li, maxi=  %li, min = %lE, max = %lE\n",size,*mini,*maxi,*min,*max); */
	
	//*mini = 0; *maxi = 0; *min = 0; *max = 0;
}
void cpeds_find_minmax_value(float *a, long int size, float *min, float *max, long int *mini, long int *maxi) {
	long int i;
	long imin,imax;
	*min = a[0];
	*max = a[0];
	*mini = 0; *maxi = 0;
	for (i=0; i<size;i++) {
		if ( a[i] < (*min) ) { *min = a[i]; imin = i; }
		if ( a[i] > (*max) ) { *max = a[i]; imax = i; }
		//printf("a[%li] = %E\n",i,a[i]);
	}
	if (mini!=NULL) *mini=imin;
	if (maxi!=NULL) *maxi=imax;
	//*mini = 0; *maxi = 0; *min = 0; *max = 0;
}

void cpeds_find_minmax_value(complex<double> *a, long int size, complex<double> *min, complex<double> *max, long int *mini, long int *maxi) {
	long int i;
	double mag;
	double MIN,MAX;
	MIN = abs(a[0]);
	MAX = abs(a[0]);
	*mini = 0; *maxi = 0;
	for (i=0; i<size;i++) {
		mag=abs(a[i]);
		if ( mag < MIN ) { MIN = mag; *mini = i; }
		if ( mag > MAX ) { MAX = mag; *maxi = i; }
	}
	*min=a[*mini];
	*max=a[*maxi];
}
/****************************************************************************************************************/
double cpeds_find_max_value(double *a, long int size, long st, long int *maxi) {
	long i;
	double max=a[st];
	long imax=st;
	for (i=st;i<size;i++) if (a[i] > max) { max=a[i]; imax=i; }
	if (maxi!=NULL) *maxi=imax;
	return max;
}
/****************************************************************************************************************/
long cpeds_find_max_value(long *a, long int size, long st, long int *maxi) {
	long i;
	long max=a[st];
	long imax=st;
	for (i=st;i<size;i++) if (a[i] > max) { max=a[i]; imax=i; }
	if (maxi!=NULL) *maxi=imax;
	return max;
}
/****************************************************************************************************************/
// gets the maximal value in the cpeds_pixel  structure in the first field starting from st index
cpeds_point cpeds_find_max_value(cpeds_point *a, long int size, long st, long int *maxi) {
	long i;
	cpeds_point max=a[st];
	(*maxi)=st;
	for (i=st;i<size;i++) if (a[i].x > max.x) { max=a[i]; (*maxi)=i; }
	return max;
}
/****************************************************************************************************************/
double cpeds_find_min_value(double *a, long int size, long st, long int *mini) {
	long i;
	double min=a[st];
	long imin;
	imin=st;
	for (i=st;i<size;i++) if (a[i] < min) { min=a[i]; imin=i; }
	if (mini!=NULL) *mini=imin;
	return min;
}
/****************************************************************************************************************/
long cpeds_find_min_value(long *a, long int size, long st, long int *mini) {
	long i;
	long min=a[st];
	(*mini)=st;
	for (i=st;i<size;i++) if (a[i] < min) { min=a[i]; (*mini)=i; }
	return min;
}
/****************************************************************************************************************/
// gets the mimimal value in the cpeds_pixel  structure in the first field starting from st index
cpeds_point cpeds_find_min_value(cpeds_point *a, long int size, long st, long int *mini) {
	long i;
	cpeds_point min=a[st];
	(*mini)=st;
	for (i=st;i<size;i++) if (a[i].x < min.x) { min=a[i]; (*mini)=i; }
	return min;
}

/****************************************************************************************************************/
long cpeds_sum(long * tab,long int size, bool deleteTab) {
	long int i;
	long tmp=0;
	for (i=0;i<size;i++) { tmp=tmp+tab[i]; }
	if (deleteTab) delete [] tab;
	return tmp;
}
/****************************************************************************************************************/
double cpeds_sum(double * tab,long int size, bool deleteTab) {
	long int i;
	double tmp=0;
	for (i=0;i<size;i++) { tmp=tmp+tab[i]; }
	if (deleteTab) delete [] tab;
	return tmp;	
}
/* ******************************************************************************************** */
long double cpeds_sum(long double * tab,long int size, bool deleteTab) {
	long int i;
	long double tmp=0;
	for (i=0;i<size;i++) { tmp=tmp+tab[i]; }
	if (deleteTab) delete [] tab;
	return tmp;		
}

/****************************************************************************************************************/
complex<double> cpeds_sum(complex<double> * tab,long int size, bool deleteTab) {
	long int i;
	complex<double> tmp=0;
	for (i=0;i<size;i++) { tmp+=tab[i]; }
	if (deleteTab) delete [] tab;
	return tmp;
}
/****************************************************************************************************************/
double cpeds_mean_value(const double * tab,long size, long startFrom) {
	long int i;
	double tmp=0;
	long endBefore=startFrom+size;
	if (size <= 0) tmp=0;
	else {
		for (i=startFrom;i<endBefore;i++) { tmp+=tab[i]; } tmp/=(double)size;
	}
	return tmp;
}
/****************************************************************************************************************/
double cpeds_median_value(const double * tab,long size) {
	if (size==0) return 0;
	if (size==1) return tab[0];
	
	double* c=cpeds_copy_array(tab,size);
	cpeds_sort_data(size,c,12);
	long ictr=size/2;
	if (size % 2==0)
		return (c[ictr]+c[ictr-1])/2;
	else
		return c[ictr];

	delete [] c;
}

/****************************************************************************************************************/
complex<double> cpeds_mean_value(complex<double> * tab,long int size) {
	long int i;
	complex<double> tmp=0;
	if (size <= 0) tmp=0;
	else {
		for (i=0;i<size;i++) { tmp+=tab[i]; } tmp/=(double)size;
	}
	return tmp;
}
/****************************************************************************************************************/
// central moment -- moment taken about the mean value
double cpeds_central_moment(double * tab,long int size, double moment, long startFrom) {
	double mean = cpeds_mean_value(tab,size,startFrom);
	long int i;
	long endBefore=startFrom+size;
	double tmp=0;
	if (size <= 0) tmp=0;
	else {
		for (i=startFrom;i<endBefore;i++) { tmp += pow(tab[i]-mean,moment); }
		tmp = tmp/(double)size;
	}
	return tmp;
}

/****************************************************************************************************************/
complex<double> cpeds_central_moment(complex<double> * tab,long int size, double moment) {
	complex<double> mean = cpeds_mean_value(tab,size);
	long int i;
	complex<double> tmp=0;
	if (size <= 0) tmp=0;
	else {
		for (i=0;i<size;i++) { tmp += pow(tab[i]-mean,moment); }
		tmp = tmp/(double)size;
	}
	return tmp;
}
/****************************************************************************************************************/
// sigma^2
double cpeds_variance(double * tab, long int size, long startFrom) {
	double tmp;
	long endBefore=startFrom+size;
	if (size<=0) tmp=0;
	else {
		if (size==1) tmp=0;
		else {
			tmp=cpeds_central_moment(tab,size,2,startFrom)*size/(size-1.0);
		}
	}
	return tmp;
}
/****************************************************************************************************************/
complex<double> cpeds_variance(complex<double> * tab, long int size) {
	complex<double> tmp;
	
	if (size<=0) tmp=0;
	else {
		if (size==1) tmp=0;
		else {
			tmp=cpeds_central_moment(tab,size,2)*(double)size/(size-1.0);
		}
	}
	return tmp;
}
/***************************************************************************************/
double cpeds_covariance(const double* d1,const double* d2,long n) {
	double cov=0;
	double m1,m2;
	m1=cpeds_mean_value(d1,n);
	m2=cpeds_mean_value(d2,n);
	for (long i = 0; i < n; i++) {
		cov+=(d1[i]-m1)*(d2[i]-m2);
	}
	cov/=n;
	return cov;
}

/****************************************************************************************************************/
double cpeds_rms(double * tab, long int size) {
	double t,rms=0;
	long i;
	for (i=0;i<size;i++) { rms+=tab[i]*tab[i]; } rms/=(double)(size);
	return sqrt(rms);
}
/****************************************************************************************************************/
complex<double> cpeds_rms(complex<double> * tab, long int size) {
	complex<double> t,rms=0;
	long i;
	for (i=0;i<size;i++) { rms+=tab[i]*conj(tab[i]); } rms/=(double)(size);
	return sqrt(rms);
}
/****************************************************************************************************************/
double cpeds_skewness(double * tab, long int size) {
	double tmp;
	if (size > 1) { tmp = cpeds_central_moment(tab,size,3.0)/pow(cpeds_variance(tab,size),1.5); } else tmp = 0;
	return tmp;
}
/****************************************************************************************************************/
double cpeds_kurtosis(double * tab, long int size) {
	double a,b,tmp;
	if (size > 1) {  a = cpeds_central_moment(tab,size,4.0); b=cpeds_variance(tab,size); b*=b; tmp = a/b; } else tmp=0;
	return tmp;
}
/****************************************************************************************************************/
// this calculates te value in th,phi direction of twodimentional gauss function stretched on sphere centered on
// phic thc with sigma s. since the boudary conditions the founcion do not integrate to one so it's unnormalized
// on sphere. the normalization if performed in separate function from optimalization purposes - since for all
// orientations of the gauss on the sphere the normalization will be the same.
double cpeds_gauss_on_sphere(double thc, double phic, double s, double th, double phi) {
	//double cpeds_gaussN(double m, double s) {
	double xth,xphi, gauss; // the crossing point on th, phi plane
	
	xphi = phic+PI; if (xphi > twoPI) { xphi = xphi-twoPI; }
	//xth = thc+PIsnd; if (xth > PI) { xth = xth-PI; }
	xth = PI-thc; //if (xth > PI) { xth = xth-PI;  }
	
	
	// setting the boundary conditions for the sphere ( the plane 0 < phi < 2PI and 0 < th < PI in fact)
	if ((phic < xphi) && (phi >= xphi)) { phi = phi-twoPI; }
	if ((thc < xth) && (th >= xth)) { th = th-PI; }
	
	if ((phic > xphi) && (phi < xphi)) { phi = phi+twoPI; }
	if ((thc > xth) && (th > xth)) { th = th+PI; }
	
	gauss = pow(sqrt(twoPI)*s,-2)*exp(- (pow(th-thc,2)+pow(phi-phic,2))/(2*s*s) );
	//printf("********* gauss %lE\n",gauss);
	return gauss;
	
}
/****************************************************************************************************************/
/****************************************************************************************************************/
// this calculates the normalization value for the cpeds_gauss_on_sphere function
// the only relevant parameter is the sigma of the function
// and the resolution parameter which defines the pixel size
double cpeds_normalize_gauss_on_sphere(double s, long int pix_num) {
	//long int i;
	double integr;
	
	
	/*   for (i=0;i<pix_num;i++) {  */
	/*     integr = integr + cpeds_gauss_on_sphere(0,0,beam_size,PIsnd-map[i].n.b,map[i].n.l); */
	/*     //printf("-- dphi: %lE sin(i) %lE, integr: %lE dth = %lE\n",dphi,SIN[i],integr, dth); */
	/*   } */
	
	/*   th = 0; integr=0; */
	/*   do { phi=0; */
	/*     do { */
	/*       integr = integr + cpeds_gauss_on_sphere(0,0,s,th,phi); */
	/*       phi=phi+dphi;  */
	/*     } while (phi <= twoPI); */
	/*     th = th+dth;   */
	/*   } while (th <= PI); */
	integr = integr*4*PI/(double)pix_num;
	return integr;
}
/****************************************************************************************************************/
double cpeds_gaussN(double x, double A, double m, double s) {
	double gauss;
	
	//gauss = 1/sqrt(2*PI)/s*exp(-(x-m)*(x-m)/(2*s*s));
	gauss = A/(sqrt(twoPI)*s)*exp(-pow(x-m,2)/(2*s*s));
	return gauss;
}
/****************************************************************************************************************/
// method 1 calculates the chi-squared

double cpeds_test_chisq(long int k, double *data, double *err, double A, double mean0, double sigma0) {
	long int i,mini,maxi;
	double chisq=0,x,min,max,dx;
	
	cpeds_find_minmax_value(data,k,&min,&max,&mini,&maxi); x=min;
	
	dx = (max - min)/(double)k;
	//mean = cpeds_mean_value(t,k);
	for (i=0;i<k;i++) {
		chisq = chisq + pow(data[i]-cpeds_gaussN(x,A,mean0,sigma0),2)/pow(err[i],2);
		//chisq = chisq/sigma0;
		x = x+dx;
	}
	
	return chisq;
	
}

/****************************************************************************************************************/
// calculates the chisq value on the unbinned sample t of size k with mean0 and variance0 of the underlying distribution - method 1
// chisq = sum_i (x_i-mean0)^2/sigma0
// calculates the chisq value on the sample t of size k relying on the frequencies of stuff in the data - method 2
// chisq = sum_i (x_i - EX(Nx))^2/EX(Nx), where EX(Nx) = k*phi(x), where phi(x) is taken from N(mean0,sigma0)
// EX(Nx) - is the expected number of given x occurances
// in the method 2 the function also calculates the C statistics C = 2 Sum_i [ exp_freq - obs_frq*ln(exp_freq) ]
// method 2 calculates the Pearson test AFAIR - and the sizes of bins are choosen as to yield the minimal frequency condition given by min_freq variable
// if it is not meet the bin is not taken into calculations, thus this sometimes might be a little bit lower then the actuall value

double cpeds_Pearson_test(long int k, double *t, double A, double mean0, double sigma0, double bin, double * bin_num_real, double * Cstat) {
	long int i,mini,maxi;
	double mean,s2,Xsq=0,x,min,max,dx,expected,observed;
	double min_freq=5; // minimal frequency in a bin
	double nbin=1; // this is how many times the bin mus have been binned in order to include more than 5 counts
	double C=0;
	
	// sort the table for binning
	
	//dx=(max-min)/cpeds_log2((double)k);
	*bin_num_real=0;
	while (x <= max) {
		dx = bin; //dx=(max-min)/1000;
		nbin = 1; expected = 0;
		do {
			expected = expected + (double)k*cpeds_gaussN(x+(nbin-1)*bin,A,mean0,sigma0)*bin;
			observed = (double)cpeds_count_numbers_in_bin(k,t,dx,x,1);
			//printf("(1)obs: x=%lE dx = %lE sample_size=%li, observed = %lE, expected: %lE,||| diff=%lE, dXsq=%lE\n",x,dx,k,observed, expected,observed-expected,pow(observed-expected,2)/expected);
			if ((expected < min_freq) || (observed < min_freq)) { nbin++; dx = nbin*bin; }
			
		} while  (((expected < min_freq) || (observed < min_freq)) && ((x+dx) < max));
		if ((expected > min_freq) && (observed > min_freq)) {
			Xsq = Xsq + pow(observed-expected,2)/expected;
			C = C + expected-observed*log(expected);
			(*bin_num_real)++;
		}
		//printf("obs: x=%lE dx = %lE,  observed = %lE, expected: %lE,||| diff=%lE, dXsq=%lE\n",x,dx,observed, expected,observed-expected,pow(observed-expected,2)/expected);
		//printf("obs: x=%lE dx = %lE,  observed = %lE, expected: %lE, dXsq=%lE\n",x,dx,observed, expected,pow(observed-expected,2)/expected);
		x=x+dx;
	}
	//printf("real bin number in the sample: %lE\n", *bin_num_real);
	*Cstat = C;
	
	return Xsq;
	
}


/****************************************************************************************************************/
bool cpeds_chisq_accepted(double alph, double degr_of_freedom, double chisq) {
	double Qx;
	
	/*   Px = gsl_cdf_chisq_Qinv(alph,degr_of_freedom);  //Px = gsl_cdf_chisq_P(chisq,degr_of_freedom);    */
	/*   if (Px < chisq)  return false; else return true;  */
	
	Qx = gsl_cdf_chisq_Q(chisq,degr_of_freedom);
	if (Qx < alph) return false; else return true;
}
/****************************************************************************************************************/
bool cpeds_chisq_accepted(double alph, double degr_of_freedom, double chisq, double *Qx) {
	double Px;
	//  printf("chisq: %lE\n",chisq);
	if (chisq<0) { *Qx=-1; return false; }
	/*   Px = gsl_cdf_chisq_Qinv(alph,degr_of_freedom);  //Px = gsl_cdf_chisq_P(chisq,degr_of_freedom);    */
	/*   *Qx = gsl_cdf_chisq_Q(chisq,degr_of_freedom);  */
	/*   if (Px < chisq)  return false; else return true;  */
	*Qx = gsl_cdf_chisq_Q(chisq,degr_of_freedom);
	if (*Qx < alph) return false; else return true;
}
/****************************************************************************************************************/
// calculates the chisq statistical significance for the chisq test result x - the probability of the null hypothesis
int cpeds_print_chisq_confidence(double alph, double degr_of_freedom, double chisq) {
	double one_le_alph, Px,Qx;
	int accepted;
	one_le_alph = 1-alph;
	Qx = gsl_cdf_chisq_Q(chisq,degr_of_freedom);  //Px = gsl_cdf_chisq_P(chisq,degr_of_freedom);
	/* printf("P[X^2 = %lE > chi^2(%.3lE;%.0lf)=%lf] = %lE",chisq,alph,degr_of_freedom,Qx,alph);  */
	if (Qx < alph) { printf("     REJECTED at %lf%% conf. lev.\n",100*one_le_alph); accepted = 0;} else { printf("     ACCEPTED at %lf%% conf. lev.\n",100*one_le_alph); accepted = 1;}
	//printf("\n");
	
	/*   Qx = gsl_cdf_chisq_Q(chisq,degr_of_freedom);  */
	printf("Qx(X^2 = %lf, nu=%lE) = %lE %%\n\n",chisq,degr_of_freedom,100*Qx);
	return accepted;
	
}
/****************************************************************************************************************/
// returns the probability that randomly ganerated data will have a chisq value smaller than chisq (P)
// i.e. that the random data driven from the distribibution fits the model better than the data being tested.
double cpeds_print_chisq_accepted_confidence(double degr_of_freedom, double chisq) {
	double P=0;
	/*   P = gsl_cdf_chisq_P(chisq,degr_of_freedom); */
	/*   printf("P(x > chisq=%lE)_{nu=%.0lf} = %lf  ||| P(x <= chisq=%lE)_{nu=%.0lf} = %lf\n",chisq,degr_of_freedom,1-P,chisq,degr_of_freedom,P); */
	return P;
}

//************************************************************************
double cpeds_calculate_chisq(double * tab1, double * tab2, long int num) {
	long int i;
	double Xsq = 0;
	
	for (i=0;i<num;i++) { Xsq = Xsq + pow((tab1[i]-tab2[i]),2); }
	return Xsq;
}

/****************************************************************************************************************/
double* cpeds_tabulate_function(int f, double A, double s, double m, long int num, double min, double max) {
	double * func  = new double[num];
	long int i;
	double x,dx;
	//double nsigma=5;
	
	
	if (f == 1) { //tabulating gauss function
		x = min; dx = (max-min)/(double)num;
		for (i=0;i<num;i++) { func[i] = A * exp(-(x-m)*(x-m)/(2*s*s)); x = x+dx; //printf("dupa0 = %lE\n",func[0]); exit(0);
		//printf("*** CDF: i=%li x=%lE CDF = %lE\n",i,x,func[i]);
		}
	}
	
	if (f==2) { //tabulating CDFgauss function
		x = min; dx = (max-min)/(double)num; func[0] = 0;//exp(-(x-m)*(x-m)/(2*s*s)); //printf("dupa0 = %lE\n",func[0]); exit(0);
		for (i=1;i<num;i++) {
			x = x+dx; func[i] = func[i-1] + exp(-(x-m)*(x-m)/(2*s*s))*dx;
			//printf("*** CDF: i=%li x=%lE CDF = %lE dx=%lE\n",i,x,func[i-1],dx);
		}
		
		for (i=0;i<num;i++) {       func[i] = A*func[i]/func[num-1];    }
		
		
	}
	
	return func;
}


/****************************************************************************************************************/
// this counts the number of data points in a bin starting at x of width bin

long int cpeds_count_numbers_in_bin(long int k, double * t, double bin, double x, int method) {
	long int i,count=0;
	int stop_flag;
	double bin_start=x,bin_end=x+bin;
	
	
	if ( method == 1) { // don't sort data
		for (i=0;i<k;i++) {
			if ((t[i] >= bin_start) && (t[i] < bin_end)) { count++; }
		}
	}
	
	if (method == 2) { // sort the data first
		cpeds_sort_data(k,t,12); // sort data increasing
		stop_flag = 0; i=0;
		do {
			if (t[i] >= bin_start) { count++; }
			if (t[i] >= bin_end) { count--; stop_flag=1; }
			i++;
		} while (stop_flag == 0);
		
	}
	//printf("x=%lE, bin=%lE, count=%li\n",x,bin,count);
	return count;
}
/****************************************************************************************************************/
// returns a pointer to array  with binned input data t. bin number is specified with num
// the x values of the bins are returned in the xbin table of size num.
// the y values of the bins are returned by a pointer
float * cpeds_bin_data_float(long int k, double * t, long num, float * xbin) {
	float * tmp = new float[num];
	double * xbind = new double[num];
	
	double * tmp2=cpeds_bin_data(k,t,num,xbind);
	for (long i=0;i<k;i++) { tmp[i]=(float)tmp2[i]; xbin[i]=(float)xbind[i]; }
	delete [] xbind;
	delete [] tmp2;
	return tmp;
}
/***************************************************************************************/
double * cpeds_bin_data(long int k, double * t, long num, double * xbin, long align, double from, double to, double geometricMultiplier) {
	long int i,mini,maxi,last;
	double min,max,x,bin;
	double * tmp = new double[num];
	//xbin = new double[num]; // we don't allocate space here now.
	
	if (from==to) {
		cpeds_find_minmax_value(t,k,&min,&max,&mini,&maxi);
	}
	else {
		min=from;
		max=to;
	}
	
	x = min; i=0; 
	if (min==max)
		bin = (1e-10)/(double)num;
	else {
		if (geometricMultiplier!=1.0) {
			double m=geometricMultiplier;
			bin=(m-1)/(powf(m,num)-1) * (max-min);
		}
		else bin = (max-min)/(double)num;
	}
	last=num-1;
	//  printf("************bin = %lf\n",bin);
	
	if (geometricMultiplier==1.0) {
		for (i=0;i<last;i++) {
			tmp[i] = cpeds_count_numbers_in_bin(k,t,bin,x,1);
			xbin[i] = x;
			x=x+bin;
			//	    printf("x: %lE, x+bin: %lE\n",x,x+bin);
		}
		tmp[i] = cpeds_count_numbers_in_bin(k,t,max-x+bin*1e-5,x,1);  xbin[i] = x;
	}
	else {
		for (i=0;i<last;i++) {
			tmp[i] = cpeds_count_numbers_in_bin(k,t,bin,x,1);
			xbin[i] = x;
			x=x+bin;
			bin*=geometricMultiplier;
			//	    printf("x: %lE, x+bin: %lE\n",x,x+bin);
		}
		tmp[i] = cpeds_count_numbers_in_bin(k,t,max-x+bin*1e-5,x,1);  xbin[i] = x;	  
	}
	
	//  printf("last bin from x: %lE\n",x);
	//  printf("last bin to x: %lE\n",max+bin*1e-5);
	//  printf("max: %lE\n",max);
	//  printf("last bin count: %lE\n",tmp[i]);
	//  exit(0);
	
	if (geometricMultiplier==1.0) {
		if (align==0) { cpeds_add_value((xbin[1]-xbin[0])/2,xbin,num); }
		if (align==1) { cpeds_add_value((xbin[1]-xbin[0]),xbin,num); }
	}
	else {
		double * xbins = new double[num];
		for (i = 0; i < num; ++i) {		xbins[i]=xbin[i];	}
		if (align==0) { 
			for (i = 0; i < num-1; ++i) {			xbin[i]=(xbins[i+1]+xbins[i])/2;		}
			xbin[num-1]=(xbins[num-1]+max)/2;
		}
		if (align==1) { 
			for (i = 0; i < num-1; ++i) {			xbin[i]=xbins[i+1];		}
			xbin[num-1]=max;
		}
		delete [] xbins;
		
	}
	
	return tmp;
}
/***************************************************************************************/
double* cpeds_bin_function(const double* xin, double* yin, const double* w, long Nin, long* binSize, long Nbins, double** xout, long* Nout, long imin) {
	printf("Nin: %li\n",Nin);
	long i,j,k,binPointsNum;
	binPointsNum=cpeds_sum(binSize,Nbins);
	long imax=imin+binPointsNum-1;
	long Cls, Cbs;//,cls,cbs;
	bool EQweights;
	double wsum,W=0;
	long Nw;
	printf("Binning the function from index: %li to index: %li, together: %li points\n",imin,imax,imax-imin+1);
	
	if (w==NULL) Nw=0; else Nw=Nin;
	// check safety conditions
	if (imin<0) { printf("negative imin arument: %li.  Not binning. \n",imin);     return NULL; }
	// requirement that binSize should be positive numbers
	//	if (!(binSize>0)) { msgs->error("negative bin size given. Check the binSize vector. Not binning",High);     return *this; }
	// if there's more bins than data
	if (imax>Nin-1) {
		printf("WARNING: The binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning vector\n");
		// fit the length of the binning vector to the actual data size if the supplied vector was too large.
		while (imax>Nin-1) {
			if (binSize[Nbins-1]>1) binSize[Nbins-1]--; else { Nbins--; }
			if (w!=NULL) Nw--;
			binPointsNum=cpeds_sum(binSize,Nbins);
			imax=imin+binPointsNum-1;
		}
		
		printf("Binning the function from index: %li to index: %li, together: %li\n",imin,imax,imax-imin+1);
		if (imin==imax) return yin; // no binning is done
	}
	
	
	// // find out the size of the binned C_b
	Cls = imax-imin+1; // size of original C_l for binning
	Cbs = Nbins; // size of binned C_l: bs
	if (Nw == 0) { EQweights = true; printf("Using same weights for all points: True\n"); } else { EQweights = false; printf("Using same weights for all points: False\n"); }
	
	if (Nw !=0 and Nw != Cls) {   printf("ERROR: Weights table size (%li) ) doesn't match the size of the data to be binned (%li)\n",Nw,Cls); }
	
	
	
	// define the structures needed
	double *Cl = new double[Cls]; // original C_l vector for binning
	double *x = new double[Cls]; // original function arguments in the binning range
	double *Cb = new double[Cbs]; // binned C_b vector
	double Cbtmp; // binned range value
	double *M = new double[Cls]; // binning operator
	double *effl = new double[Cbs]; // effective l-value
	
	
	// prepare the power spectrum vector for binning
	for (j=Cls;j>0;j--) {  Cl[j-1] = yin[j-1+imin]; x[j-1]=xin[j-1+imin]; }
	//deletePoint(j-1+imin); } // copy the power spectra part for binning and remove the points that will be binned
	
	//
	// prepare the binning matrix operator
	//
	
	for (j=0;j<Cls;j++)	M[j] = 0;  // zero the binning matrix
	k=0;
	double kmin,kmax;
	for (i=0;i<Cbs;i++) {
		if (EQweights) W = 1/(double)binSize[i];
		effl[i]=0; wsum=0;
		kmax=k+binSize[i];
		for (j=k;j<kmax;j++) { M[j]=0; }
		for (j=k;j<kmax;j++) {
			if (!EQweights) { W = w[j]; wsum+=W; }
			M[j] = W;
			effl[i]+=x[j]*W;
		}
		if (!EQweights) { for (j=k;j<kmax;j++) { M[j]/=wsum; } effl[i]/=wsum; }
		
		// do the binning
		Cbtmp=0;
		for (j=k;j<kmax;j++) { Cbtmp+= M[j]*Cl[j]; }
		Cb[i]=Cbtmp;
		
		k+=binSize[i];
	}
	
	//
	// rewrite the Cb onto the output arrays
	//
	*Nout=Cbs+imin+(Nin-imax-1);
	double* yout=new double[*Nout];
	*xout=new double[*Nout];
	
	for (long i = 0; i < imin; i++) {
		*xout[i]=xin[i];
		yout[i]=yin[i];
	}
	for (i=0;i<Cbs;i++)  {
		(*xout)[i+imin]=effl[i];
		yout[i+imin]=Cb[i];
	}
	for (long i = imax+1; i < Nin; i++) {
		(*xout)[imin+Cbs+i-imax-1]=xin[i];
		yout[imin+Cbs+i-imax-1]=yin[i];
	}
	
	
	// 
	// clean up
	//
	delete [] effl;
	delete [] x;
	delete [] Cl;
	delete [] Cb;
	delete [] M;
	
	return yout;
}

/***************************************************************************************/
//double* cpeds_bin_data_geo(long k, double* t, long from, long to, double binSize, double gm, cpedsList<double>& xbin ) {
//	double b=binSize;
//	long fsize;
//	if (to!=-1) { fsize=to+1; } else fsize=k;
//	i=from;
//	while (i<fsize) {
//		i+=long(round(b));
//		if (i<fsize) { xbin.append(long(round(b))); } else { xbin.append(fsize-i+long(round(b))); }
//		b*=gm;
//	}
//	cpeds_bin_data(k,t,_from,_bintab,_w);
//
//}

double* cpeds_copy_array(const double* t, long N) { 
	if (N==0) return NULL; 
	double* copy=new double[N]; 
	for (long i=0;i<N;i++) { 
		copy[i]=t[i];  
	}
	return copy; 
}

void cpeds_add_value(double val, double* t, long N) { for (long i=0;i<N;i++) { t[i]+=val;  } }
void cpeds_add_value(double val, fftwComplex* t, long N) { for (long i = 0; i < N; i++) {	t[i][0]+=val;	t[i][1]+=val; } }
void cpeds_sub_value(double val, double* t, long N) { for (long i=0;i<N;i++) { t[i]-=val;  } }
void cpeds_sub_value(double val, fftwComplex* t, long N) { for (long i = 0; i < N; i++) {	t[i][0]-=val;	t[i][1]-=val; } }
void cpeds_mul_value(double val, double* t, long N) { for (long i=0;i<N;i++) { t[i]*=val;  } }
void cpeds_mul_value(double val, fftwComplex* t, long N) { for (long i = 0; i < N; i++) {	t[i][0]*=val;	t[i][1]*=val; } }
void cpeds_divide_value(double val, double* t, long N) { for (long i=0;i<N;i++) { t[i]/=val;  } }
void cpeds_divide_value(double val, fftwComplex* t, long N) { for (long i = 0; i < N; i++) {	t[i][0]/=val;	t[i][1]/=val; } }

/* template <class Type> Type *cpeds_bin_data(double *t,long int k, Type *xbin, long num ) { */
/*   long int i,mini,maxi; */
/*   Type min,max,x,bin; */
/*   Type tmp = new Type[num]; */
/*   //xbin = new float[num]; */

/*   cpeds_find_minmax_value(t,k,&min,&max,&mini,&maxi);  */
/*   x = min; i=0; bin = (max-min)/(Type)num;  */
/*   printf("************bin = %lf\n",bin); */
/*   for (i=0;i<num;i++) { */
/*     tmp[i] = (Type)cpeds_count_numbers_in_bin(k,t,bin,x,1); */
/*     xbin[i] = x; */
/*     //printf("tmp[%li] = %E xbin[i] = %E\n",i,tmp[i],xbin[i]); */
/*     x=x+bin; */
/*   } */
/*   return tmp; */
/* } */




/****************************************************************************************************************/
// sorts the data - 12 - increasing, 21 decreasing

int cpeds_qsort_compare12(const void *a, const void *b) {
	double arg1 = *static_cast<const double*>(a);
	double arg2 = *static_cast<const double*>(b);
	
	if(arg1 < arg2) return -1;
	if(arg1 > arg2) return 1;
	return 0;	
}

int cpeds_qsort_compare21(const void *a, const void *b) {
	double arg1 = *static_cast<const double*>(a);
	double arg2 = *static_cast<const double*>(b);
	
	if(arg1 > arg2) return -1;
	if(arg1 < arg2) return 1;
	return 0;		
}

void cpeds_sort_data(long int k, double * t, int direction) {
	long int i,j,endi=k-1;
	double tmp;
	
	if (direction == 12) { // sort increasing
		
		qsort(t,k,sizeof(double),&cpeds_qsort_compare12);
		/*
    for (i=0;i<endi;i++) {
      tmp=t[i];
      t[i]=cpeds_find_min_value(t,k,i,&j);
      t[j]=tmp;
    }
		 */
	}
	
	if (direction == 21) { // sort decreasing
		qsort(t,k,sizeof(double),&cpeds_qsort_compare21);
		/*
    for (i=0;i<endi;i++) {
      tmp=t[i];
      t[i]=cpeds_find_max_value(t,k,i,&j);
      t[j]=tmp;
    }
		 */
	}
	
}

void cpeds_sort_data(long int k, long * t, int direction) {
	long int i,j,endi=k-1;
	long tmp;
	
	if (direction == 12) { // sort increasing
		for (i=0;i<endi;i++) {
			tmp=t[i];
			t[i]=cpeds_find_min_value(t,k,i,&j);
			t[j]=tmp;
		}
	}
	
	if (direction == 21) { // sort decreasing
		for (i=0;i<endi;i++) {
			tmp=t[i];
			t[i]=cpeds_find_max_value(t,k,i,&j);
			t[j]=tmp;
		}
	}
}

// sorts the structure according to the first field of the structure
void cpeds_sort_data(long int k, cpeds_point * t, int direction) {
	long int i,j,endi=k-1;
	cpeds_point tmp;
	
	if (direction == 12) { // sort increasing
		for (i=0;i<endi;i++) {
			tmp=t[i];
			t[i]=cpeds_find_min_value(t,k,i,&j);
			t[j]=tmp;
		}
	}
	
	if (direction == 21) { // sort decreasing
		for (i=0;i<endi;i++) {
			tmp=t[i];
			t[i]=cpeds_find_max_value(t,k,i,&j);
			t[j]=tmp;
		}
	}
}
/***************************************************************************************/
double cpeds_fwhm2sigma(double fwhm) {
	return fwhm*0.4246609001440095;
}
/***************************************************************************************/
void cpeds_smoothGauss_data2D(double* data, long N, long n1, double sl0, double sl1) {
	// get sigma from fwhm
	double s0 = sl0/(2.0*sqrt(2.0*log(2.0))); // * twoPI; // BLcomment (Nov 11, 2011, 6:53:08 PM): why here is twoPI??
	double s1 = sl1/(2.0*sqrt(2.0*log(2.0))); // * twoPI;
	
	// prepare the data
	long sx=N/n1;
	long sy=n1;
	long syf=sy/2+1; // size of the array in y direction in fourier space
	long Nf=sx*syf; // total size of the array in fourier space
	
	double *in=data;
	fftwComplex* out = new fftwComplex[Nf];
	
	// make fft along first dimention
	fftw_plan p;
	int n;
	int inembed, onembed;
	double kx,ky;
	double scale=1.0/sy;
	long ij;
	double tf;
	double s02o2, s12o2;
	
	s02o2=s0*s0/2;
	n=sx;
	onembed=sx/2+1;
	p=fftw_plan_many_dft_r2c(1,&n,sy,in,NULL,sy,1,out,NULL,sx/2+1,1,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	// convolve
	scale=1.0/sx;
	for (long i = 0; i < sx; ++i) {
		for (long j = 0; j <= sy/2; ++j) {			
			kx=double(i)/sx;			
			tf=exp(-kx*kx*s02o2);
			tf*=scale;
			ij=syf*i+j;
			out[ij][0]*=tf;
			out[ij][1]*=tf;
		}
	}
	//	n=sx/2+1;
	onembed=sx;
	p=fftw_plan_many_dft_c2r(1,&n,sy,out,NULL,sx/2+1,1,in,NULL,sy,1,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	
	//	// save the matrix
	//	//cpeds_save_matrix(in,sx,sy,"field-smx");
	//	
	
	/***************************************************************************************/
	// do fft in y direction
	n=sy;
	inembed=sy;
	onembed=sy/2+1;
	//	p=fftw_plan_dft_r2c_2d(sx,sy, &in[0][0],&out[0][0],FFTW_ESTIMATE); 
	p=fftw_plan_many_dft_r2c(1,&n,sx,in,NULL,1,sy,out,NULL,1,sy/2+1,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	// convolve
	scale=1.0/sy;
	s12o2=s1*s1/2;
	for (long i = 0; i < sx; ++i) {
		for (long j = 0; j <= sy/2; ++j) {			
			ky=double(j)/sy;			
			tf=exp(-ky*ky*s12o2);
			tf*=scale;
			ij=syf*i+j;
			out[ij][0]*=tf;
			out[ij][1]*=tf;
		}
	}
	
	// make inverse fft
	n=sy;
	p=fftw_plan_many_dft_c2r(1,&n,sx,out,NULL,1,sy/2+1,in,NULL,1,sy,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	// clean up
	delete [] out;
}

/***************************************************************************************/
void cpeds_fft_data3D(double* in, fftwComplex* out, long n0, long n1, long n2, int dim) {
	fftw_plan p;
	long Nf=n0*n1*(n2/2+1);
	long N=n0*n1*n2;
	
	int n;
	//	double scale;
	
	switch (dim) {
		case 0:
			n=n0;
			//			scale=1.0/n0;
			p=fftw_plan_many_dft_r2c(1,&n,n1*n2,in,NULL,n1*n2,1,out,NULL,n1*n2,1,FFTW_ESTIMATE); 
			fftw_execute(p);
			fftw_destroy_plan(p);
			
			break;
		case 1:
			n=n1;
			p=fftw_plan_many_dft_r2c(1,&n,n2,in,NULL,n2,1,out,NULL,n2,1,FFTW_ESTIMATE); 
			for (long i = 0; i < n0; i++) {		
				fftw_execute(p);
				cpeds_shift_array(in,N,n2*n1);
				cpeds_shift_array(out,Nf,n2*(n1/2+1));
			}
			fftw_destroy_plan(p);
			
			break;
		case 2:
			n=n2;
			p=fftw_plan_many_dft_r2c(1,&n,n0*n1,in,NULL,1,n2,out,NULL,1,n2/2+1,FFTW_ESTIMATE); 
			fftw_execute(p);
			fftw_destroy_plan(p);
			
			break;
		default:
			printf("cpeds_fft_data3D: wrong dim parameter: should be 0,1 or 2");
			break;
	}
}

/***************************************************************************************/
void cpeds_fft_data3D(fftwComplex* in, double* out, long n0, long n1, long n2, int dim) {
	fftw_plan p;
	long Nf=n0*n1*(n2/2+1);
	long N=n0*n1*n2;
	
	int n;
	//	double scale;
	
	switch (dim) {
		case 0:
			n=n0;
			p=fftw_plan_many_dft_c2r(1,&n,n1*n2,in,NULL,n1*n2,1,out,NULL,n1*n2,1,FFTW_ESTIMATE); 
			fftw_execute(p);
			fftw_destroy_plan(p);
			break;
		case 1:
			n=n1;
			p=fftw_plan_many_dft_c2r(1,&n,n2,in,NULL,n2,1,out,NULL,n2,1,FFTW_ESTIMATE); 
			for (long i = 0; i < n0; i++) {		
				cpeds_shift_array(in,Nf,n2*(n1/2+1),false);
				cpeds_shift_array(out,N,n2*n1,false);
				fftw_execute(p);
			}
			fftw_destroy_plan(p);
			break;
		case 2:
			n=n2;
			p=fftw_plan_many_dft_c2r(1,&n,n0*n1,in,NULL,1,n2/2+1,out,NULL,1,n2,FFTW_ESTIMATE); 
			fftw_execute(p);
			fftw_destroy_plan(p);
			
			break;
		default:
			printf("cpeds_fft_data3D: wrong dim parameter: should be 0,1 or 2");
			break;
	}
	
}
/***************************************************************************************/
void cpeds_convolve_gauss(fftwComplex* data, long n0, long n1, long n2, int dim, double fwhm, bool scale) {
	double scalef,t;
	double K;
	long idx;
	double sigma=fwhm/(2*sqrt(2*log(2.0))); //* twoPI;
	
	switch (dim) {
		case 0:
			if (scale) scalef=1.0/n0; else scalef=1;
			for (long i = 0; i <= n0/2; i++) {
				for (long j = 0; j < n1; j++) {
					for (long k = 0; k < n2; k++) {			
						//						printf("%li %li %li\n",i,j,k);
						K=double(i)/n0;
						t=exp(-K*K*sigma*sigma/2);
						t*=scalef;
						idx=k + j*n2 + i*n1*n2;
						data[idx][0]*=t;
						data[idx][1]*=t;
					}
				}
			}
			break;
			
		case 1:
			if (scale) scalef=1.0/n1; else scalef=1;
			for (long i = 0; i < n0; i++) {
				for (long j = 0; j <= n1/2; j++) {
					for (long k = 0; k < n2; k++) {			
						K=double(j)/n1;
						t=exp(-K*K*sigma*sigma/2);
						t*=scalef;
						idx=k + j*n2 + i*(n1/2+1)*n2;
						data[idx][0]*=t;
						data[idx][1]*=t;
					}
				}
			}
			break;
			
		case 2:
			if (scale) scalef=1.0/n2; else scalef=1;
			for (long i = 0; i < n0; i++) {
				for (long j = 0; j < n1; j++) {
					for (long k = 0; k <= n2/2; k++) {			
						K=double(k)/n2;
						t=exp(-K*K*sigma*sigma/2);
						t*=scalef;
						
						idx=k + j*(n2/2+1) + i*n1*(n2/2+1);
						data[idx][0]*=t;
						data[idx][1]*=t;
					}
				}
			}
			break;
			
		default:
			printf("cpeds_fft_data3D: wrong dim parameter: should be 0,1 or 2\ndim: %i",dim);
			exit(0);
			break;
	}
}

/***************************************************************************************/
void cpeds_convolve_gauss(double* data, long n0, long n1, long n2, int dim, double fwhm, bool scale) {
	long Nf=n0*n1*(n2/2+1);
	fftwComplex* out=new fftwComplex[Nf];
	cpeds_fft_data3D(data,out,n0,n1,n2,dim);
	cpeds_convolve_gauss(out,n0,n1,n2,dim,fwhm,true);
	cpeds_fft_data3D(out,data,n0,n1,n2,dim);	
	delete [] out;
}

/****************************************************************************************************************/
// finds the value closest to the value val in the array t of double values of size ts searching in num values from
// array index starting at start.
// assumes that the array is an array of values SORTED  in rising order.
// the seartch is performed with half division method
// returns the array index of the closest value
long cpeds_find_value(double val,double * t,long ts, long start, long num) {
	long d,m,u,id;
	double dd,du;
	
	d=start; u=d+num-1; m=(long)floor((double)(d+u)/2.0);
	
	while (d+1<u) {
		if (t[m] < val) { d=m;  } // drop the lower range
		else { u=m; } // drop the upper range
		m=(long)floor((double)(d+u)/2.0);
	}
	
	dd=fabs(t[d]-val);   du=fabs(t[u]-val);
	if (dd < du) id = d; else id=u;
	
	return id;
}
long cpeds_find_value(long val,long * t,long ts, long start, long num) {
	long d,m,u,id;
	long dd,du;
	
	d=start; u=d+num-1; m=(long)floor((double)(d+u)/2.0);
	
	while (d+1<u) {
		if (t[m] < val) { d=m;  } // drop the lower range
		else { u=m; } // drop the upper range
		m=(long)floor((double)(d+u)/2.0);
	}
	
	dd=labs(t[d]-val);   du=labs(t[u]-val);
	if (dd < du) id = d; else id=u;
	
	return id;
}


/****************************************************************************************************************/
bool cpeds_val_in_array(double val, double* t, long size) {
	long i;
	for (i=0;i<size;i++) { if (t[i] == val) return true; }
	return false;
}

/****************************************************************************************************************/
double* cpeds_point2array(long N, cpeds_point *Pxy, long coord) {
	double * t = new double[N];
	long i;
	if (coord==0)   for (i=0;i<N;i++) t[i]=Pxy[i].x;
	if (coord==1)   for (i=0;i<N;i++) t[i]=Pxy[i].y;
	return t;
}
/****************************************************************************************************************/
//void cpeds_MCMC_parameters_fit
/* void cpeds_matrix_io(cpeds_strarg whattodo, cpeds_strarg what) { */


/* } */
/****************************************************************************************************************/
void cpeds_save_matrix_lin(double *Dvec, long vec_size, long vec_num, string file, bool deleteDvec, bool append) {
	long k,l;
	FILE* f;
	char c;
	if (append) c='a';
	else c='w';
	//  f=fopen(file.c_str(),"w");
	f=fopen(file.c_str(),&c);
	fprintf(f,"matrix is %li %li\n",vec_size,vec_num);
	for (k=0;k<vec_size;k++) {
		for (l=0;l<vec_num;l++) {
			fprintf(f," %.12lE ",Dvec[k+l*vec_size]);
		}
	}
	fprintf(f,"\n");
	fclose(f);
	if (deleteDvec) delete [] Dvec;
}
/***************************************************************************************/
void cpeds_print_matrix(double *Dvec, long vec_size, long vec_num) {
	long k,l;
	printf("matrix is %li %li\n",vec_size,vec_num);
	for (k=0;k<vec_size;k++) {
		printf("|");
		for (l=0;l<vec_num;l++) {
			printf(" %.12lE ",Dvec[k+l*vec_size]);
		}
		printf("|\n");
	}
	printf("|\n");
	
}

/****************************************************************************************************************/
void cpeds_print_matrix(const matrix<double>& m) {
	long rn = m.RowNo();
	long cn = m.ColNo();
	long i,j;
	
	for (i=0;i<rn;i++) {
		for (j=0;j<cn;j++) {
			printf("%lf ",m(i,j));
		}
		printf("\n");
	}
	printf("\n");
}

/****************************************************************************************************************/
// saves the cpeds matrix in the rows-major ordering
void cpeds_save_matrix(double * M, long rows, long cols, string name, bool deleteM,bool append) {
	FILE * f;
	long k,i,l;
	if (append)   
		f = fopen(name.c_str(),"a");
	else
		f = fopen(name.c_str(),"w");
	l=0;
	for (k=0;k<rows;k++) {
		for (i=0;i<cols;i++) { fprintf(f,"%.12lE ",M[l]); l++; } fprintf(f,"\n"); }
	fclose(f);
	if (deleteM) delete [] M;
}
/***************************************************************************************/
void cpeds_save_matrix(float * M, long rows, long cols, string name, bool deleteM,bool append) {
	FILE * f;
	long k,i,l;
	if (append)   
		f = fopen(name.c_str(),"a");
	else
		f = fopen(name.c_str(),"w");
	l=0;
	for (k=0;k<rows;k++) {
		for (i=0;i<cols;i++) { fprintf(f,"%.12E ",M[l]); l++; } fprintf(f,"\n"); }
	fclose(f);
	if (deleteM) delete [] M;
}
/***************************************************************************************/
void cpeds_save(fftwComplex* M, long rows, string name, bool deleteM) {
	FILE * f;
	long k,i,l;
	f = fopen(name.c_str(),"w");
	l=0;
	for (k=0;k<rows;k++) {
		fprintf(f,"%.12lE %.12lE\n",M[l][0],M[l][1]); l++; } 
	fclose(f);
	if (deleteM) delete [] M;	
}

/****************************************************************************************************************/
matrix<double>& cpeds_shift_matrix(matrix<double>& m, long Nrows, long Ncols) {
	long rn = m.RowNo();
	long cn = m.ColNo();
	long i,j,ip,jp;
	matrix<double> tmp=m;
	
	for (i=0;i<rn;i++) {
		for (j=0;j<cn;j++) {
			ip=(i-Nrows)%rn; if (ip<0) ip=rn+ip;
			jp=(j-Ncols)%cn; if (jp<0) jp=cn+jp;
			tmp(i,j)=m(ip,jp);
		}
	}
	m=tmp;
	return m;
}
/****************************************************************************************************************/
matrix<double>& cpeds_pad_matrix(matrix<double>& m, long NrowPaddingAbove, long NrowPaddingBelow, long NcolPaddingBefore, long NcolPaddingAfter) {
	long rn = m.RowNo();
	long cn = m.ColNo();
	long i,j;
	if (NrowPaddingAbove<0) NrowPaddingAbove=0;
	if (NrowPaddingBelow<0) NrowPaddingBelow=0;
	if (NcolPaddingBefore<0) NcolPaddingBefore=0;
	if (NcolPaddingAfter<0) NcolPaddingAfter=0;
	
	matrix<double> tmp(rn+NrowPaddingAbove+NrowPaddingBelow,cn+NcolPaddingBefore+NcolPaddingAfter);
	tmp.Null();
	
	for (i=0;i<rn;i++) {
		for (j=0;j<cn;j++) {
			tmp(i,j)=m(i,j);
		}
	}
	m=cpeds_shift_matrix(tmp,NrowPaddingAbove,NcolPaddingBefore);
	return m;
}
/****************************************************************************************************************/
double* cpeds_matrix2array(const matrix<double>& m, long& size, bool rowMajor) {
	long rows=m.RowNo();
	long cols=m.ColNo();
	double *t = new double[rows*cols];
	long i,j,k=0;
	
	if (rowMajor) {
		for (i=0;i<rows;i++) {
			for (j=0;j<cols;j++) {
				t[k]=m(i,j);
				k++;
			}
		}
	}
	else {
		for (j=0;j<cols;j++) {
			for (i=0;i<rows;i++) {
				t[k]=m(i,j);
				k++;
			}
		}
	}
	
	size=k;
	return t;
}
/****************************************************************************************************************/
float* float_double_array(double* t, long size) {
	float* f = new float[size];
	for (long i=0;i<size;i++) { f[i]=(float)t[i]; }
	return f;
}

/***************************************************************************************/
void cpeds_shift_array(double*t, long size, long N, bool fwd) {
	// this is a generalized version of the implementation below for arbitrary shift size
	double tmp,tmpNext;
	long iNext,i,i0,j;
	if (N<0) { printf("cpeds_shift_array: N<0. use fwd=false instead\n"); exit(0); }
	N=N%size;
	if (N==0) return;
	if (fwd==false) N=size-N;
	
	long s2=size-2;
	i=0;
	j=0;
	i0=i;
	tmp=t[i];
	do {
		iNext=(i+N) % size;
		tmpNext=t[iNext];
		t[iNext]=tmp;
		tmp=tmpNext;
		i=iNext;
		j++;
		if (i==i0 and j<size) { i++; tmp=t[i]; i0=i;}
	} while (j<size);
	return;
	
	//#ifdef DEBUG
	//		printf("shifting array only shifts by one currently. FIX ME SOON");
	//#endif
	//		// this is a smarter version of the commented out implementation 
	//	double tmp;
	//	long s2=size-2;
	//	if (fwd) {
	//		tmp=t[1];
	//		t[1]=t[0];
	//		t[0]=t[size-1];
	//		for (long i = s2; i >= 2; i--) { t[i+1]=t[i]; }
	//		t[2]=tmp;
	//	}
	//	else {
	//		tmp=t[s2];
	//		t[s2]=t[size-1];
	//		t[size-1]=t[0];
	//		for (long i = 1; i < s2; i++) { t[i-1]=t[i]; }
	//		t[s2-1]=tmp;		
	//	}
	//	
	//	return t;
	
	//	double *t2=new double[size];
	//	long ip;
	//	if (fwd) {
	//		for (long i = 0; i < size; i++) {
	//			ip=(i+N) % size;
	//			t2[ip]=t[i];				
	//		}
	//	}
	//	else {
	//		for (long i = 0; i < size; i++) {
	//			ip=(i+N) % size;
	//			t2[i]=t[ip];				
	//		}
	//	}
	//	delete [] t;
	//	return t2;
}

/***************************************************************************************/
void cpeds_shift_array(fftwComplex *t, long size, long N, bool fwd) {
	// this is a generalized version of the implementation below for arbitrary shift size
	fftwComplex tmp,tmpNext;
	long iNext,i,i0,j;
	if (N<0) { printf("cpeds_shift_array: N<0. use fwd=false instead\n"); exit(0); }
	N=N%size;
	if (N==0) return;
	if (fwd==false) N=size-N;
	
	long s2=size-2;
	i=0;
	j=0;
	i0=i;
	tmp[0]=t[i][0];	tmp[1]=t[i][1];
	do {
		iNext=(i+N) % size;
		tmpNext[0]=t[iNext][0];	tmpNext[1]=t[iNext][1];
		t[iNext][0]=tmp[0];	t[iNext][1]=tmp[1];
		tmp[0]=tmpNext[0]; tmp[1]=tmpNext[1];
		i=iNext;
		j++;
		if (i==i0 and j<size) { i++; tmp[0]=t[i][0]; tmp[1]=t[i][1]; i0=i; }
	} while (j<size);
	return;
}

/****************************************************************************************************************/
const matrix<double> cpeds_array2matrix(const double* t, long size, long vecSize, bool rowMajor) {
	matrix<double> m;
	long rows;
	long cols;
	long i,j,k=0;
	
	if (rowMajor) {
		cols=vecSize;
		rows=size/vecSize;
		m.SetSize(rows,cols);
		for (i=0;i<rows;i++) {
			for (j=0;j<cols;j++) {
				m(i,j)=t[k];
				k++;
			}
		}
	}
	else {
		rows=vecSize;
		cols=size/vecSize;
		m.SetSize(rows,cols);
		for (j=0;j<cols;j++) {
			for (i=0;i<rows;i++) {
				m(i,j)=t[k];
				k++;
			}
		}
	}
	
	return m;
}


/****************************************************************************************************************/

void cpeds_change_matrix_ordering_from_rows_to_cols_major(double *M, long rows, long cols) {
	long size=rows*cols;
	double * tmp= new double[size];
	long i,j,k=0;
	
	// copy array
	for (i=0;i<size;i++) { tmp[i]=M[i]; }
	
	// change ordering
	for (i=0;i<rows;i++) {
		for (j=0;j<cols;j++) { M[k]=tmp[j*rows+i]; k++; } }
	
	delete [] tmp;
	
}

/****************************************************************************************************************/
/* Calculates the covariance matrix of the measured data in Dvec table, organized in vec_num vectors 1-row vectors, each of size vec_size and form (1___x___vec_size)  */
/* If i iterates vector index and j iterates variate in i'th vector then the ordering of the Dvec array is i-major: i.e. vector index major (i.e. variate-minor). */
/* Hence the Dvec is a set of vectors where each row-vector is a single measurement of all variates */
/* The resulting cov matrix is a square var_num___x___var_num size symmetric matrix */
/* given by pointer cov. */
/* THe NCov parameter is used to indicate from which measurement to start the calculation: can be 0..vec_num-1 -- NOT IMPLEMENTED YET*/

/* double * cpeds_calculate_covariance_matrix(double *Dvec, long vec_size, long vec_num, long NCovst) { */
double * cpeds_calculate_covariance_matrix(double *Dvec, long vec_size, long vec_num, bool diagonal) {
	long i,k,l;
	long idxC,idxD;
	double *cov;
	
	if (diagonal) cov= new double[vec_size];
	else cov= new double[vec_size*vec_size];
	
	double *av = new double[vec_size];
	double tmp,Nleo;
	/*   long vec_num_loc=vec_num-NCovst; */
	Nleo=double(vec_num-1);
	if (Nleo==0) { printf("too few (1) data vectors given for covariance matrix computations. Will stop now."); exit(0); }
	
	// calculate average vectors first
	for (k=0;k<vec_size;k++) {
		av[k] = 0;
		for (i=0;i<vec_num;i++) { av[k] += Dvec[i*vec_size+k]; }
		av[k]/=(double)vec_num;
	}
	
	if (diagonal) { // only the diagonal elements are stored
		for (k=0;k<vec_size;k++) {
			cov[k] = 0;
			for (i=0;i<vec_num;i++) { idxD=i*vec_size;  tmp=Dvec[idxD+k]-av[k]; cov[k] += tmp*tmp; }
			cov[k] /= Nleo;
		}
	}
	else {  
		for (k=0;k<vec_size;k++) {
			for (l=0;l<=k;l++) {
				idxC=k*vec_size+l;      cov[idxC] = 0;
				for (i=0;i<vec_num;i++) { idxD=i*vec_size;  cov[idxC] += (Dvec[idxD+k]-av[k])*(Dvec[idxD+l]-av[l]); }
				cov[idxC] /= Nleo;
				cov[l*vec_size+k] = cov[idxC]; //symetrize
			}
		}
	}
	delete [] av;
	return cov;
}
/****************************************************************************************************************/
//! Calculates the two-sided, upper-tail quantile probability of getting value x, based on values given in array t of size ts which probe the underlying PDF. */
/*! The probability is that of getting a measurment deviating in absolute value from second quartile (Q24) by more than the measured value x. */
/*! Too small and too big deviations will be unprobable (unlike in the Q probability) */
/*! NOTE: assumes that the data in t are sorted in rising order */

double cpeds_quantile_P(double x,double * t,long ts) {
	long I_Q24 = (long)round(ts/2);
	long xIDX = cpeds_find_value(x,t,ts,0,ts);
	double P;
	
	if (xIDX > I_Q24) P = (double)(2*(ts-xIDX))/(double)ts;
	else { if (xIDX==0) xIDX++; P = (double)(2*xIDX)/(double)ts; }
	return P;
}


/**************************************************************************************************************/
double cpeds_quantile_P_sort(double x,double * t,long ts) {
	gsl_sort(t,1,(size_t)ts);
	return cpeds_quantile_P(x,t,ts);
}
/**************************************************************************************************************/
//
// The sign indicates the measurment position with regard to 2nd quartile: - left side, + right side
// the probability is of getting a measurment deviating in abs value from Q24 by more than the measured value.
// CAUTION !! assumes that the data in t are sorted in rising order

double cpeds_quantile_P_signed(double x,double * t,long ts) {
	long I_Q24 = (long)round(ts/2);
	long xIDX = cpeds_find_value(x,t,ts,0,ts);
	double P;
	
	if (xIDX > I_Q24) P = (double)(2*(ts-xIDX))/(double)ts;
	else { if (xIDX==0) xIDX++; P = -(double)(2*xIDX)/(double)ts; }
	return P;
}
/***************************************************************************************/
double cpeds_quantile_P_sort_signed(double x,double * t,long ts) {
	gsl_sort(t,1,(size_t)ts);
	return cpeds_quantile_P_signed(x,t,ts);
}
/**************************************************************************************************************/
// the sign indicates the measurment position with regard to 2nd quartile: - left side, + right side
// the probability is of getting a measurment deviating in abs value from Q24 by more than the measured value.
// P=2 is a control value to detect in case something unpredicted happened.
// Routine can also extrapolate probability using gaussian exponential lower/upper tail extrapolation
// CAUTION !! assumes that the data in t are sorted in rising order

double cpeds_quantile_P_signed_interpolated(double x,double * t,long ts) {
	double I_Q24,Q24;
	double xIDX;
	double P=2,gap,f;
	double sqrt2=1.41421356237310;
	
	I_Q24 = round((double)ts/2)-1; Q24 = t[(long)I_Q24];
	
	if (x > t[ts-1]) { // do Gaussian extrapolation
		/*     f=-invCDFGauss(0.5/(double)ts)*fabs(x-Q24)/fabs(Q24-t[ts-1]); // old extrapolation, erroneous I believe */
		/*     f=-invCDFGauss(1.0/(double)ts)*fabs(x-Q24)/fabs(Q24-t[ts-1]); // this one does not account for the variance of the MC PDF */
		/*     P=1.0-gsl_sf_erf(f/sqrt2); */
		P=cpeds_extrapolate_gauss(x,t,ts);
	}
	else {
		if (x < t[0]) {  // do Gaussian extrapolation
			/*       f=-invCDFGauss(0.5/(double)ts)*fabs(x-Q24)/fabs(Q24-t[0]);  // old extrapolation, erroneous I believe */
			/*       f=-invCDFGauss(1.0/(double)ts)*fabs(x-Q24)/fabs(Q24-t[0]); //// this one does not account for the variance of the MC PDF */
			/*       P=gsl_sf_erf(f/sqrt2)-1.0; */
			P=cpeds_extrapolate_gauss(x,t,ts);
		}
		else { // do linear interpolation within the probed regime
			xIDX = (double)cpeds_find_value(x,t,ts,0,ts);
			gap=x-t[(long)xIDX];
			if (gap==0) { f=0; } else {
				if (gap > 0) { f=gap/(t[(long)xIDX+1]-t[(long)xIDX]); } else { f=gap/(t[(long)xIDX]-t[(long)xIDX-1]); }
			}
			/*       printf("*********f=%lE\n",f); */
			
			if (fabs(xIDX-I_Q24) <= 1) { // aviod P>1 case
				/* 	if ((ts % 2) == 0) { P = 2*((double)ts-(xIDX+f)-1)/(double)ts; } else { P = (2*((double)ts-(xIDX+f))-1)/(double)ts; } */
				P=((double)ts-1)/(double)ts; // I can change that some day to do better interpolation in this region
				if (x < t[(long)I_Q24]) P = -P;    }
			else { // don't mind the potential assymetry in placeing Q_24, since the undertainty of the true Q24 would probably have bigger impact on P than due to 1/ts miscalculation from even/odd problem in quantile approach.
				if (x >= t[(long)I_Q24]) {       P = (2*((double)ts-(xIDX+f)))/(double)ts; }    else {       P = -(2*(xIDX+1+f))/(double)ts; } }
		}
	}
	
	return P;
}

double cpeds_quantile_P_sort_signed_interpolated(double x,double * t,long ts) {
	gsl_sort(t,1,(size_t)ts);
	return cpeds_quantile_P_signed_interpolated(x,t,ts);
}
/***************************************************************************************/
double cpeds_quantile(double* t, long ts, double quantile, bool deleteAfter) {
	double ret;
	long idx;
	cpeds_sort_data(ts,t,12);
	idx=long(round(quantile*ts));
	if (idx<0) idx=0;
	if (idx>=ts) idx=ts-1;
	ret=t[idx];
	if (deleteAfter) delete [] t;
	return ret;
}
/**************************************************************************************************************/
// the signed routine for calculation gaussian extrapolation (interpolation is also possible) given set of data.
// the tail probabilities are caclulated separatelly depening on which tail the data sits. and
// to that tail the gauss curve is fitted. the sewing assumes that the CDF should reach 2/N value at the x[0]
// (and also at x[ts-1]). the sign of the probabilitiy is only indication on which tail the data sits.
// SO THIS ROUTUNE RETURNS THE TAIL PROBABILITY NOT UPPER OR LOWER TAIL PROB.
// IF YOU WANT A TAIL PROBABILITY -- IT'S JUST 0.5 OF THE RETURNED VALUE SINCE THE GAUSSIAN PDF IS SYMMETRICAL
// WARNING ! if this routine is used as fitting gauss curve within the MCPDF probed region, some |P|>1
// are possible around Q24. to be fixed !
// CAUTION !! assumes that the data in t are sorted in rising order
// the input data MUST BE SORTED IN ASCENDING ORDER
double cpeds_extrapolate_gauss(double x, double* t, long ts) {
	double s, nsig, nsigMC = fabs(invCDFGauss(1/(double)ts)); //printf("******* %lE\n",nsigMC);
	double P=-2; // err detection value in case something goes wrong
	long iQ24 = (long)round((double)ts/2);
	long i,j;
	double Q24 = t[iQ24];
	double sq2=1.41421356237310;
	double *tt = new double[ts];
	bool lower;
	
	// calculate tail variance of the MC PDF to find out n_sigma
	if (x <= Q24) { // will calculate according to lower tail
		/*     for (i=0;i<iQ24;i++) { j=2*i; tt[j]=t[i]-Q24; tt[j+1]=-tt[j]; } */
		/*     s=sqrt(cpeds_variance(tt,2*iQ24-1)); */
		s=fabs(t[0]-Q24)/nsigMC; //s/=2;
		nsig=nsigMC+(t[0]-x)/s;
		lower=true;
	}
	else { // will calculate according to upper tail
		/*     for (i=iQ24;i<ts;i++) { j=2*(i-iQ24); tt[j]=t[i]-Q24; tt[j+1]=-tt[j]; } */
		/*     s=sqrt(cpeds_variance(tt,2*(ts-iQ24)-1)); */
		s=fabs(t[ts-1]-Q24)/nsigMC; //s/=2;
		nsig=nsigMC+(x-t[ts-1])/s;
		lower=false;
	}
	
	
	/*   P=0.5*(1.0-gsl_sf_erf(nsig/sq2)); */
	P=1.0-gsl_sf_erf(nsig/sq2);
	
	if (lower) P=-P;
	/*   printf("tail sigma: %lE, nsig: %lE nsigMC: %lE, P: %lE \n",s,nsig,nsigMC,P); */
	delete [] tt;
	return P;
}

double cpeds_sort_extrapolate_gauss(double x, double* t, long ts) {
	gsl_sort(t,1,(size_t)ts);
	return cpeds_extrapolate_gauss(x,t,ts);
}

/**************************************************************************************************************/
double cpeds_extrapolate_linear(double x, double x1, double y1, double x2, double y2) {
	//y=(y2-y1)/(x2-x1)* (x-x1)+y1
	double f=(y2-y1)/(x2-x1);
	return f*(x-x1)+y1;
}
/**************************************************************************************************************/
/*  calculates the lower tail inverse gaussian CDF  */
/* invoke with P \in {0,p/2} for finding the numer of sigmas deviation */
/* e.g. invCDFGauss ( 0.0027/2) = -3 (sigma) (and remember about the signs)*/
double invCDFGauss(double p) {
	double q, r;
	
	errno = 0;
	
	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < _invCDFGauss_LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((_invCDFGauss_c[0]*q+_invCDFGauss_c[1])*q+_invCDFGauss_c[2])*q+_invCDFGauss_c[3])*q+_invCDFGauss_c[4])*q+_invCDFGauss_c[5]) /
				((((_invCDFGauss_d[0]*q+_invCDFGauss_d[1])*q+_invCDFGauss_d[2])*q+_invCDFGauss_d[3])*q+1);
	}
	else if (p > _invCDFGauss_HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((_invCDFGauss_c[0]*q+_invCDFGauss_c[1])*q+_invCDFGauss_c[2])*q+_invCDFGauss_c[3])*q+_invCDFGauss_c[4])*q+_invCDFGauss_c[5]) /
				((((_invCDFGauss_d[0]*q+_invCDFGauss_d[1])*q+_invCDFGauss_d[2])*q+_invCDFGauss_d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
		q = p - 0.5;
		r = q*q;
		return (((((_invCDFGauss_a[0]*r+_invCDFGauss_a[1])*r+_invCDFGauss_a[2])*r+_invCDFGauss_a[3])*r+_invCDFGauss_a[4])*r+_invCDFGauss_a[5])*q /
				(((((_invCDFGauss_b[0]*r+_invCDFGauss_b[1])*r+_invCDFGauss_b[2])*r+_invCDFGauss_b[3])*r+_invCDFGauss_b[4])*r+1);
	}
}
/**************************************************************************************************************/
/*
 * The standard normal CDF, for one random variable.
 *
 *   Author:  W. J. Cody
 *   URL:   http://www.netlib.org/specfun/erf
 *
 * This is the erfc() routine only, adapted by the
 * transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 */
/* double cpeds_CDFGauss(double u) { */
/*  const double a[5] = { */
/*   1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001, */
/*   1.887426188426510e+002,3.209377589138469e+003 */
/*  }; */
/*  const double b[5] = { */
/*   1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002, */
/*   1.813893686502485e+003,8.044716608901563e+003 */
/*  }; */
/*  const double c[9] = { */
/*   2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00, */
/*   6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02, */
/*   1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03 */
/*  }; */
/*  const double d[9] = { */
/*   1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02, */
/*   5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03, */
/*   4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03 */
/*  }; */
/*  const double p[6] = { */
/*   1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1, */
/*   1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4 */
/*  }; */
/*  const double q[6] = { */
/*   1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00, */
/*   5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3 */
/*  }; */
/*  register double y, z; */

/*  if (_isnan(u)) */
/*   return _Nan._D; */
/*  if (!_finite(u)) */
/*   return (u < 0 ? 0.0 : 1.0); */
/*  y = fabs(u); */
/*     if (y <= 0.46875*M_SQRT2) { */
/*   /\* evaluate erf() for |u| <= sqrt(2)*0.46875 *\/ */
/*   z = y*y; */
/*   y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4]) */
/*        /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]); */
/*   return 0.5+y; */
/*  } */
/*  z = exp(-y*y/2)/2; */
/*  if (y <= 4.0) { */
/*   /\* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 *\/ */
/*   y = y/M_SQRT2; */
/*   y = */
/* ((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8]) */


/* /((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]); */

/*   y = z*y; */
/*     } else { */
/*   /\* evaluate erfc() for |u| > sqrt(2)*4.0 *\/ */
/*   z = z*M_SQRT2/y; */
/*   y = 2/(y*y); */
/*         y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5]) */
/*     /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]); */
/*         y = z*(M_1_SQRTPI-y); */
/*     } */
/*  return (u < 0.0 ? y : 1-y); */
/* }; */

//***************************************************************************************************
// this function calculates the Wigner 3J symbol values for all ms=0
double cpeds_W3jm0(long l1, long l2, long l3) {
	long L=l1+l2+l3;
	long L12,L13,L23,d12,d13,d23,r,g;
	double tmpd, f1a1,f1a2,f1a3;
	long double A,frac1,frac2;
	long i,j,tmpl;
	
	// check the parity of all ls
	r=L%2;
	if (r == 1) { return 0; }
	else {
		// check the triangle relation
		L12=l1+l2; L13=l1+l3; L23=l2+l3;
		
		if (L12 < l3) return 0;
		if (L23 < l1) return 0;
		if (L13 < l2) return 0;
		
		d12=abs(l1-l2);    d13=abs(l1-l3);    d23=abs(l2-l3);
		if (d12 > l3) return 0;
		if (d23 > l1) return 0;
		if (d13 > l2) return 0;
		
		// ok, it's non-zero; derive it.
		
		//sort in rising order -- it doesn't change anything
		if (l3 < l2) { tmpl=l2; l2=l3; l3=tmpl; }
		if (l2 < l1) { tmpl=l1; l1=l2; l2=tmpl; }
		if (l3 < l2) { tmpl=l2; l2=l3; l3=tmpl; }
		if (l2 < l1) { tmpl=l1; l1=l2; l2=tmpl; }
		/*     printf("l1: %li, l2: %li, l3: %li \n",l1,l2,l3); */
		
		
		g=L/2;
		// sign factor
		A=(double)cpeds_pow_int(-1,g);
		
		// first fraction
		f1a1=L-2*l1;    f1a2=L-2*l2;    f1a3=L-2*l3;
		frac1=cpeds_factorial_ld(f1a2)*cpeds_factorial_ld(f1a3)/cpeds_factorial_frac_ld(L+1,f1a1+1);
		/*     m=cpeds_get_max(cpeds_get_max(f1a1,f1a2,&i),f1a3,&j); */
		/*     if (j==2) frac1=cpeds_factorial(f1a1)*cpeds_factorial(f1a2)/cpeds_factorial_frac(L+1,f1a3+1);  */
		/*     else { */
		/*       if (i==1) frac1=cpeds_factorial(f1a2)*cpeds_factorial(f1a3)/cpeds_factorial_frac(L+1,f1a1+1); */
		/*       if (i==2) frac1=cpeds_factorial(f1a1)*cpeds_factorial(f1a3)/cpeds_factorial_frac(L+1,f1a2+1); */
		/*     } */
		
		/*     frac1 = cpeds_factorial(L-2*l1)*cpeds_factorial(L-2*l2)*cpeds_factorial(L-2*l3)/cpeds_factorial(L+1); */
		
		// second fraction
		
		frac2 = cpeds_factorial_frac_ld(g,g-l1+1)/(cpeds_factorial_ld(g-l2)*cpeds_factorial_ld(g-l3));
		
		/*     frac2 = cpeds_factorial(g)/(cpeds_factorial(g-l1)*cpeds_factorial(g-l2)*cpeds_factorial(g-l3)); */
		tmpd = (double)(A*sqrt(frac1)*frac2);
		
	}
	return tmpd;
}


//***************************************************************************************************
// this function calculates the Wigner 3J symbol values for all ms=0
/* double cpeds_W3jm0_lidia(long l1, long l2, long l3) { */
/*   long L=l1+l2+l3; */
/*   long L12,L13,L23,d12,d13,d23,r,g,A; */
/*   double tmpd, f1a1,f1a2,f1a3; */
/*   bigfloat frac1,frac2; */
/*   bigint tmpi; */
/*   long i,j,tmpl; */

/*   // check the parity of all ls */
/*   r=L%2; */
/*   if (r == 1) { return 0; }  */
/*   else { */
/*     // check the triangle relation */
/*     L12=l1+l2; L13=l1+l3; L23=l2+l3; */

/*     if (L12 < l3) return 0; */
/*     if (L23 < l1) return 0; */
/*     if (L13 < l2) return 0; */

/*     d12=abs(l1-l2);    d13=abs(l1-l3);    d23=abs(l2-l3); */
/*     if (d12 > l3) return 0; */
/*     if (d23 > l1) return 0; */
/*     if (d13 > l2) return 0; */

/*     // ok, it's non-zero; derive it. */

/*     //sort in rising order -- it doesn't change anything */
/*     if (l3 < l2) { tmpl=l2; l2=l3; l3=tmpl; } */
/*     if (l2 < l1) { tmpl=l1; l1=l2; l2=tmpl; } */
/*     if (l3 < l2) { tmpl=l2; l2=l3; l3=tmpl; } */
/*     if (l2 < l1) { tmpl=l1; l1=l2; l2=tmpl; } */
/* /\*     printf("l1: %li, l2: %li, l3: %li \n",l1,l2,l3); *\/ */


/*     g=L/2;  */
/*     // sign factor */
/*     A=(double)cpeds_pow_int(-1,g); */
/* /\*     printf("********* A=%li ",A); *\/ */

/*     // first fraction */
/*     f1a1=L-2*l1;    f1a2=L-2*l2;    f1a3=L-2*l3; */
/* /\*     cout << "a1: "<< f1a1<< " a2: "<<f1a2 << " a3: "<<f1a3 <<"   ||| " << cpeds_factorial_lidia(f1a2)<<" "<<cpeds_factorial_lidia(f1a3)<<" "<<cpeds_factorial_frac_lidia(L+1,f1a1+1); *\/ */
/*     tmpi=cpeds_factorial_lidia(f1a2)*cpeds_factorial_lidia(f1a3); frac1=tmpi; frac1/=cpeds_factorial_frac_lidia(L+1,f1a1+1); */
/* /\*     cout << " frac1: " << frac1 << " tmpi "<< tmpi; *\/ */

/*     // second fraction */

/*     frac2 = cpeds_factorial_frac_lidia(g,g-l1+1); tmpi=cpeds_factorial_lidia(g-l2)*cpeds_factorial_lidia(g-l3); frac2/=tmpi; */
/* /\*     cout << " frac2: " << frac2; *\/ */
/*     frac1=sqrt(frac1)*frac2; frac1*=A; */
/*     frac1.doublify(tmpd); */
/* /\*     cout << "  tmpd " << tmpd << "\n"; *\/ */

/*   } */
/*   return tmpd; */
/* } */
//***************************************************************************************************
bool cpeds_W3jm0_non_zero(long l1, long l2, long l3) {
	long L=l1+l2+l3;
	long L12,L13,L23,d12,d13,d23,r,g,A;
	
	// check the parity of all ls
	r=L%2;
	if (r == 1) { return false; }
	else {
		// check the triangle relation
		L12=l1+l2; L13=l1+l3; L23=l2+l3;
		
		if (L12 < l3) return false;
		if (L23 < l1) return false;
		if (L13 < l2) return false;
		
		d12=abs(l1-l2);    d13=abs(l1-l3);    d23=abs(l2-l3);
		if (d12 > l3) return false;
		if (d23 > l1) return false;
		if (d13 > l2) return false;
	}
	return true;
}

bool cpeds_W3jm0_lcond(long l1, long l2, long l3, long i) {
	long L=l1+l2+l3;
	long L12,L13,L23,d12,d13,d23,r,g,A;
	
	// check the parity of all ls
	/*   r=L%2; */
	
	// check the triangle relation
	
	if (i == 0) { //checking lower boundary cond.
		d12=abs(l1-l2);    d13=abs(l1-l3);    d23=abs(l2-l3);
		if (d12 > l3) return false;
		if (d23 > l1) return false;
		if (d13 > l2) return false;
	}
	else {  //checking upper boundary cond.
		L12=l1+l2; L13=l1+l3; L23=l2+l3;
		if (L12 < l3) return false;
		if (L23 < l1) return false;
		if (L13 < l2) return false;
	}
	return true;
}

//***************************************************************************************************
// returns a^b; a and b must be natural numbers and a > 0 !!
long cpeds_pow_int(long a, long b) {
	long i,x=1;
	for (i=0;i<b;i++) { x*=a; }
	return x;
}


//***************************************************************************************************
// this functions measures the number of cols in the first line of the file fn
// it was tested and it works for the normal linux txt files
// the field delimiter is single space
// it works good whether or not the last line in file ends with \n (this doesn't matter here since it's only for the first line)
// if you run it on some strange stuff then you better check the results.

long cpeds_get_cols_num_first_ln(strarg fn,char * lastc) {
	char c;
	long r,cols=0;
	bool oneorzero=true;
	FILE *f;
	
	f=fopen(fn,"r"); if (f==NULL) { printf("error: no such file\n"); return 0; }
	
	r=fscanf(f,"%c",&c);
	if (c!=' ' && feof(f)==0) cols++;
	
	while (c!='\n' && feof(f)==0) {
		r=fscanf(f,"%c",&c);
		if (c==' ' && feof(f)==0) cols++;
		oneorzero=false;
	}
	
	if (r!=0) { if (lastc!=NULL) *lastc=c; } else cols=0;
	if (oneorzero && c=='\n') cols=0;
	
	fclose(f);
	
	return cols;
}
//***************************************************************************************************
long cpeds_get_cols_num_first_ln(strarg fn, bool commentedFile, long maxLineSize) {
	ifstream f;
	f.open(fn);
	
	if (f.fail()) {
		return -1;
	}
	
	//  filebuf fb;
	//  fb.open (fn,ios::in);
	//  istream in(&fb);
	char l[maxLineSize];
	bzero(l,maxLineSize);
	
	bool go=true;
	while (go) {   
		//	  in.getline(&l[0],maxLineSize); 
		f.getline(l,maxLineSize);
		if (l[0]!='#') go=false;   
		//	  printf("reading: %s\n",l); 
	} // printf("reading\n"); }
	f.close();
	//  fb.close();
	std::string line=l;
	//  printf("dupa: %s\n",l);
	//  cout << line << "\n";
	std::stringstream is(line);
	std::string s;
	int cols=0;
	while (is >> s) {
		cols++;
		//     std::cout << s << std::endl; 
	}
	/* std::cout << cols; */
	
	return cols;
	
}

//***************************************************************************************************
// this functions measures the number of cols starting from a position in the file indicated by the f pointer
// it works in the same way as the cpeds_get_cols_num_first_ln but does not colose the file.
// the position of the f pointer is increased by the line width since this routine reads the whole line untill the
// \n character

long cpeds_get_cols_num(FILE *f,char * lastc) {
	char c,cp;
	long r,cols=0;
	bool oneorzero=true;
	long i=0;
	
	/* OLD STUFF BELOW */
	/*   r=fscanf(f,"%c",&c);  */
	/*   if (c!=' ' && feof(f)==0) cols++;  */
	/*   cp=c; */
	
	/*   while (c!='\n' && feof(f)==0) { */
	/*     r=fscanf(f,"%c",&c);  */
	/*     if (c==' ' && feof(f)==0 && c!=cp) cols++;  */
	/*     if (cp==' ' && c=='\n' && cols>0) cols--; */
	/*     cp=c; */
	/*  /\*    printf("%c\n",c); *\/ */
	/*     oneorzero=false; */
	/*   } */
	
	/*   if (r!=0) *lastc=c; else cols=0; */
	/*   if (oneorzero && c=='\n') cols=0; */
	/* OLD STUFF ABOVE */
	
	r=fscanf(f,"%c",&c);
	while (c!='\n' && feof(f)==0) {
		// reach next column
		while ((c==' ' || c=='\t') && feof(f)==0) { r=fscanf(f,"%c",&c); i++; /* printf("moving to next col %li\n",i); */ }
		if (c!=' ' && c!='\t' && c!='\n') { cols++; /* printf("cols: %li\n",cols); */ }
		i=0;
		// going through the column
		while ((c!=' ' && c!='\t') && feof(f)==0 && c!='\n') { r=fscanf(f,"%c",&c); i++; /* printf("going threough the col %li\n",i); */}
		i=0;
		
	}
	
	if (r!=0) *lastc=c; else cols=0;
	
	return cols;
}
/***************************************************************************************/
long long cpeds_get_file_size(string fName) {
	struct stat filestatus;
	stat( fName.c_str(), &filestatus );
	long long fs=filestatus.st_size;
	return fs;
}
/***************************************************************************************/
long long cpeds_get_txt_file_lines_count(string fName) {
    unsigned int number_of_lines = 0;
    FILE *infile = fopen(fName.c_str(), "r");
    int ch;

    while (EOF != (ch=getc(infile)))
        if ('\n' == ch)
            ++number_of_lines;
    fclose(infile);
    return number_of_lines;
}
//***************************************************************************************************
// this routine checks the txt file fn and returns an cpeds_queue object that contains the
// number of columns in each row of the file and much other useful info

cpeds_queue<long>* cpeds_get_txt_file_cols_rows(strarg fn) {
	long cols=0,lines=0;
	char c;
	cpeds_queue <long> *q;
	FILE * f;
	f=fopen(fn,"r"); if (f==NULL) { printf("error: no such file\n"); return NULL; }
	q = new cpeds_queue <long>;
	
	if (cpeds_get_file_size(fn)==0) return q;
	
	while (feof(f)==0) {
		lines++;
		cols=cpeds_get_cols_num(f,&c);
		q->addq(cols);
	}
	
	if (cols==0) q->delq();
	fclose(f);
	
	return q;
}
/***************************************************************************************/

bool cpeds_fileExists(string fname) {
	if (fname=="") return false;
	FILE* f = fopen(fname.c_str(),"r");
	if (f==NULL) return false;
	fclose(f);
	return true;
	/*
	 * Comment:
	 * this block has been uncommented because checking if file exists by using stat only
	 * worked correctly in a single threaded environment, but it lead to std::bad_alloc
	 * when run in MPI jobs on cluster. Perhaps some sort of MPI_File_open could / should be 
	 * used in parallel runs. Anyways this solution works fine on parallel jobs, but is possibly
	 * slower.
	 * 
	 * 
	 * author: blew
	 * date: Nov 19, 2012 5:36:51 PM
	 *
	 */
	
	//	struct stat buffer;
	//	if ( stat( fname.c_str(), &buffer )==0 ) return true ;
	//	return false;
}
/************************************************************************************************/
const matrix<double> cpeds_matrix_load(string fileName, string how, long * result) {
	unsigned long i,j;
	FILE * f=NULL;
	long row,col;
	double tmp;
	char c;
	cpeds_queue<long> *finfo;
	bool header,binary,fourbyte;
	if (how.find("header",0) == string::npos) header=false; else header=true;
	if (how.find("binary",0) == string::npos) binary=false; else binary=true;
	if (how.find("float",0) == string::npos) fourbyte=false; else fourbyte=true;
	
	if (binary)  {f = fopen(fileName.c_str(),"rb");} else {  f = fopen(fileName.c_str(),"r"); }
	if (f == NULL) { printf("   -- ERROR: no file: %s\n",fileName.c_str()); matrix<double> m; if (result!=NULL) *result=-1; return m ; }
	
	if (header) {    fscanf(f,"%li %li",&row,&col); }
	else {
		if (binary) {
			printf("* binary read requested with no header information\n");
			struct stat64 s; stat64(fileName.c_str(),&s);
			col=1;
			if (fourbyte) {
				printf("* requested 4-byte float type read\n");
				row=(long)s.st_size/(long)sizeof(float);
			}
			else {
				printf("* assuming 8-byte float type read\n");
				row=(long)s.st_size/(long)sizeof(double);
			}
			printf("* assuming the matrix of size: rows %li cols %li\n",row,col);
		}
		else {
			col=0;
			row=0;
			finfo=cpeds_get_txt_file_cols_rows(fileName.c_str());
			if (finfo!=NULL) {
				if (finfo->get_size()>0) {
					col=(*finfo)(0);
					row=(*finfo).get_size();    	  
				}
			}
			printf("* assuming the matrix of size: rows %li cols %li\n",row,col);
		}
		
	}
	printf("  -- reading matrix of size: %li rows, %li cols from file %s\n",row,col,fileName.c_str()); //exit(0);
	matrix<double> M; 
	
	if (col==0 or row==0) { if (result!=NULL) *result=-1; return M; }
	M.SetSize(row,col);
	
	if (binary) {
		size_t ret,siz;
		float flt;
		if (header) fread(&c,sizeof(c),1,f); // read in the \n character-- 1 byte to get the right values !! this might be system dependent or library dependent !!
		if (fourbyte) {
			siz=sizeof(flt);
			for (i=0;i<M.RowNo();i++) {
				for (j=0;j<M.ColNo();j++)  { ret=fread(&flt,siz,1,f); M(i,j) = flt; }
			}
		}
		else {
			siz=sizeof(tmp);
			for (i=0;i<M.RowNo();i++) {
				for (j=0;j<M.ColNo();j++)  { ret=fread(&tmp,siz,1,f); M(i,j) = tmp; }
			}
		}
	}
	else {
		int ret;
		for (i=0;i<M.RowNo();i++) {
			for (j=0;j<M.ColNo();j++)  { ret=fscanf(f,"%lE ",&tmp); M(i,j) = tmp; }//printf("i=%li, j =%li\n",i,j); }
		}
	}
	
	
	fclose(f);
	if (result!=NULL) *result=0;
	return M;
}


/************************************************************************************************/
long cpeds_matrix_save(const matrix<long>& M, string fileName, string how, int precision) {
	matrix<double> md=matrix<double>(M);
	return cpeds_matrix_save(md, fileName,how);
}
long cpeds_matrix_save(const matrix<double>& M, string fileName, string how, int precision) {
	unsigned long i,j;
	FILE * f=NULL;
	double tmpd;
	bool header,binary,fourbyte,longtype;
	if (how.find("header",0) == string::npos) header=false; else header=true;
	if (how.find("binary",0) == string::npos) binary=false; else binary=true;
	if (how.find("float",0) == string::npos) fourbyte=false; else fourbyte=true;
	
	if (header) {
		f = fopen(fileName.c_str(),"w");
		if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",fileName.c_str());  return -1; }
		fprintf(f,"%li %li\n",(long)M.RowNo(),(long)M.ColNo());
		printf("  -- saving matrix of size: %li rows %li cols to file %s\n",(long)M.RowNo(),(long)M.ColNo(),fileName.c_str());
		fclose(f);
	}
	else {
		f = fopen(fileName.c_str(),"w");
		if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",fileName.c_str());  return -1; }
		fclose(f);
	}
	
	if (binary) {
		f = fopen(fileName.c_str(),"ab");
		if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",fileName.c_str());  return -1; }
		if (fourbyte) {
			printf("* requested 4-byte float type save\n");
			float flt;
			size_t siz=sizeof(float);
			for (i=0;i<M.RowNo();i++) {
				for (j=0;j<M.ColNo();j++)  { flt = (float)(M(i,j));  fwrite(&flt,siz,1,f); }
			}
		}
		else {
			printf("* assuming 8-byte float type save\n");
			size_t siz=sizeof(tmpd);
			for (i=0;i<M.RowNo();i++) {
				for (j=0;j<M.ColNo();j++)  { tmpd = (double)(M(i,j));  fwrite(&tmpd,siz,1,f); }
			}
		}
	}
	else {
		f = fopen(fileName.c_str(),"a");
		if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",fileName.c_str());  return -1; }
		stringstream fmt;
		fmt << "%." << precision << "lE ";
		for (i=0;i<M.RowNo();i++) {
			for (j=0;j<M.ColNo();j++)  { tmpd = M(i,j); fprintf(f,fmt.str().c_str(),tmpd); }   fprintf(f,"\n");
		}
	}
	
	fclose(f);
	return 0;
}
/***************************************************************************************/
double cpeds_get_memory_usage() {
	ifstream proc("/proc/self/status");
	string s;
	while(getline(proc, s), !proc.fail()) {
		if(s.substr(0, 6) == "VmSize") {
			break;
		}
	}
	//parse to number in kB
	int i = s.size();
	if (i< 7) return -1.0;
	double vmsize=strtod(&s.c_str()[7],NULL);
	return vmsize;
}



/************************************************************************************************/
// This class implements deriviations of the spherical bessel function

/********************************************************************/
/* ORIGINAL FORTRAN CODE FROM THE CMBFAST */
/********************************************************************/
/*       subroutine bjl(L,x,jl) */
/* c  Calculates the spherical bessel function j_l(x)  */

/* c  and optionally its derivative for real x and integer l>=0.  */

/* c  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from  */

/* c  G.N.Watson, A Treatise on the Theory of Bessel Functions, */
/* c  2nd Edition (Cambridge University Press, 1944). */
/* c  Higher terms in expansion for x near l given by */
/* c  Airey in Phil. Mag. 31, 520 (1916). */

/* c  This approximation is accurate to near 0.1% at the boundaries */
/* c  between the asymptotic regions; well away from the boundaries */
/* c  the accuracy is better than 10^{-5}. The derivative accuracy */
/* c  is somewhat worse than the function accuracy but still better */
/* c  than 1%. */

/* c  Point *jlp initially to a negative value to forego calculating */
/* c  the derivative; point it to a positive value to do the derivative */
/* c  also (Note: give it a definite value before the calculation */
/* c  so it's not pointing at junk.) The derivative calculation requires  */

/* c  only arithmetic operations, plus evaluation of one sin() for the */
/* c  x>>l region.   */


/* c  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu */
/* c  This fortran version only computes j_l(x) */
/* 	implicit double precision(a-h,o-z) */
/* 	integer L */
/* 	real*8 nu, nu2,ax,ax2,beta,beta2,beta4,beta6 */
/*   	real*8 sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1 */
/*  	real*8 cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b */
/*  	real*8 trigarg,trigcos,expterm,prefactor,llimit,ulimit,fl */
/* 	real*4 x,jl */

/* 	PI = 3.1415926536 */
/* 	ROOTPI = 1.772453851 */
/* 	GAMMA1 = 2.6789385347             !/\* Gamma function of 1/3 *\/ */
/* 	GAMMA2 = 1.3541179394             !/\* Gamma function of 2/3 *\/ */


/*   	ax = abs(dble(x)) */
/*  	fl = l */


/* 	beta = fl**0.325 */
/*  	llimit=1.31*beta   !/\* limits of asymptotic regions; fitted *\/ */
/*   	ulimit=1.48*beta */

/*  	 nu= fl + 0.5         */

/*  	 nu2=nu*nu */

/*   	if (l .lt. 0) then */
/* 		print*, 'Bessel function index < 0\n' */
/* 		stop */
/* 	end if */

/* c          /\************* Use closed form for l<6 **********\/ */

/*  	if (l .lt. 6) then                 */

/*     	sx=sin(ax) */
/* 	cx=cos(ax) */
/* 	ax2=ax*ax */

/*    	    if(l .eq. 0) then */
/*      		 if(ax .gt. 0.001) then */
/* 		 jl=real(sx/ax) */
/* 		 else  */

/* 		 jl=real(1.0d0-ax2/6.0d0) */
/* 		 end if				!   /\* small x trap *\/ */
/* 	    endif */


/*     	    if(l .eq. 1) then  */

/*      		 if(ax .gt. 0.001) then */
/* 		 jl=real((sx/ax -cx)/ax) */
/* 		 else  */

/* 		 jl=real(ax/3.0d0) */
/* 		 end if */
/* 	    endif */

/*     	    if(l .eq. 2) then */
/*     		  if(ax .gt. 0.001) then */
/* 		  jl=real((-3.0d0*cx/ax  */
/*      1                    -sx*(1.0d0-3.0d0/ax2))/ax) */
/*     		  else  */

/* 		  jl=real(ax2/15.0d0) */
/* 		  end if */
/* 	    endif */

/*             if(l .eq. 3) then */
/*       		if(ax .gt. 0.001) then */
/* 		jl=real((cx*(1.0d0-15.0d0/ax2) */
/*      1                    -sx*(6.0d0-15.0d0/ax2)/ax)/ax) */
/*      		else  */

/* 		jl=real(ax*ax2/105.0d0) */
/* 		endif */
/* 	    endif */

/*     	    if(l .eq. 4) then */
/*       		if(ax .gt. 0.001) then */
/* 	jl=real((sx*(1.0d0-45.0d0/(ax*ax)+105.0d0 */
/*      1       /(ax*ax*ax*ax)) +cx*(10.0d0-105.0d0/(ax*ax))/ax)/ax) */
/*       		else  */

/* 		jl=real(ax2*ax2/945.0d0) */
/* 		end if */
/* 	    endif */

/*    	     if(l .eq. 5) then  */

/*       		if(ax .gt. 0.001) then */
/* 	jl=real((sx*(15.0d0-420.0d0/(ax*ax)+945.0d0 */
/*      1    /(ax*ax*ax*ax))/ax -cx*(1.0-105.0d0/(ax*ax)+945.0d0 */
/*      1                                    /(ax*ax*ax*ax)))/ax) */
/*      		else  */

/* 		jl=real(ax2*ax2*ax/10395.0d0) */
/* 		endif */
/* 	     endif */


/* c          /\********************** x=0 **********************\/ */

/*   	else if (ax .lt. 1.d-30) then */
/*     	jl=0.0 */

/* c          /\*************** Region 1: x << l ****************\/ */

/* 	else if (ax .le. fl+0.5-llimit) then */


/* c       beta=acosh(nu/ax) */
/* 	if (nu/ax .lt. 1.d0) print*, 'trouble with acosh' */
/* 	beta = dlog(nu/ax + dsqrt((nu/ax)**2 - 1.d0) )  */
/* 		!(4.6.21) */
/* 	cotb=nu/sqrt(nu*nu-ax*ax)      ! /\* cotb=coth(beta) *\/ */
/*     	cot3b=cotb*cotb*cotb */
/*     	cot6b=cot3b*cot3b */
/*         secb=ax/nu */
/*     	sec2b=secb*secb */
/*     	sec4b=sec2b*sec2b */
/*     	sec6b=sec4b*sec2b */
/*     	sum1=2.0+3.0*sec2b */
/*     	expterm=sum1*cot3b/(24.0*nu) */
/*     	sum2=4.0+sec2b */
/*     	expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2) */
/*     	sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b */
/*    	expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2) */
/*     	sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b */
/*    	expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2) */
/*     	expterm=exp(-nu*beta+nu/cotb-expterm) */
/*     	prefactor=sqrt(cotb/secb)/(2.0*nu) */
/*     	jl=real(prefactor*expterm) */

/* c          /\**************** Region 2: x >> l ****************\/ */


/*   	else if (ax .ge. fl+0.5+ulimit) then          */


/*     	beta=acos(nu/ax) */
/*     	cotb=nu/sqrt(ax*ax-nu*nu)      !/\* cotb=cot(beta) *\/ */
/*     	cot3b=cotb*cotb*cotb */
/*     	cot6b=cot3b*cot3b */
/*     	secb=ax/nu */
/*     	sec2b=secb*secb */
/*     	sec4b=sec2b*sec2b */
/*     	sec6b=sec4b*sec2b */
/*     	trigarg=nu/cotb - nu*beta - PI/4.0 */
/*     	sum1=2.0+3.0*sec2b */
/*     	trigarg = trigarg - sum1*cot3b/(24.0*nu) */
/*     	sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b */
/*     	trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2) */
/*     	trigcos=cos(trigarg) */
/*     	sum2=4.0+sec2b */
/*     	expterm=sum2*sec2b*cot6b/(16.0*nu2) */
/*     	sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b */
/*     	expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2) */
/*     	expterm=exp(-expterm) */
/*     	prefactor=sqrt(cotb/secb)/nu */
/*     	jl=real(prefactor*expterm*trigcos) */

/* c          /\***************** Region 3: x near l ****************\/ */

/*   	else                        */



/*     	beta=ax-nu         */

/*     	beta2=beta*beta */
/*     	beta4=beta2*beta2 */
/*     	beta6=beta2*beta4 */
/*     	sx=6.0/ax */
/*     	sx2=sx*sx */
/*     	cx=sqrt(sx)                    */

/*     	secb=sx**0.333333333     */

/*     	sec2b=secb*secb */

/*     	deriv1=GAMMA1*secb */
/*     	deriv1= deriv1+ beta*GAMMA2*sec2b */
/*     	sum1=(beta2/6.0-1.0/15.0)*beta */
/*     	deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0 */
/*     	sum2=beta4/24.0-beta2/24.0+1.0/280.0 */
/*     	deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0 */
/*     	sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0 */
/*     	deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0 */
/*     	sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0- */
/*      2	13.0/31500.0)*beta */
/*     	deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0 */
/*     	sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0 */
/*      1               -121.0*beta2/907200.0 + 7939.0/232848000.0)*beta */
/*     	deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0     */

/*     	jl=real(deriv1*cx/(12.0*ROOTPI)) */

/* 	end if */

/* 	return */
/* 	end */

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */


// this is an implementation of the above fortran code into C.

double cpeds_spherical_bessel_fn(long l,double x) { // spherical bessel functions -- adopted from cmbfast code
	
	double GAMMA1 = 2.6789385347;             /* Gamma function of 1/3 */
	double GAMMA2 = 1.3541179394;             /* Gamma function of 2/3 */
	
	double nu, nu2,ax,ax2,beta,beta2,beta4,beta6;
	double sx,sx2,cx,sum1,sum2,sum3,sum4,sum5,deriv1;
	double cotb,cot3b,cot6b,secb,sec2b,sec4b,sec6b;
	double trigarg,trigcos,expterm,prefactor,llimit,ulimit,fl;
	double jl;
	
	ax = fabs(x);
	fl=(double)l;
	
	beta = pow(fl,0.325);
	
	llimit=1.31*beta; /* limits of asymptotic regions; fitted */
	ulimit=1.48*beta;
	
	nu= fl + 0.5;
	
	nu2=nu*nu;
	
	if (l < 0) {
		printf("Bessel function index < 0\n");
		exit(0);
	}
	
	
	/************* Use closed form for l<6 **********/
	
	if (l < 6) {
		
		sx=sin(ax);
		cx=cos(ax);
		ax2=ax*ax;
		
		if(l == 0) {
			if(ax > 0.001) jl=sx/ax;  else jl=1.0-ax2/6.0;
		}   /* small x trap */
		
		
		if(l == 1) {
			if(ax > 0.001) jl=(sx/ax-cx)/ax;      else jl=ax/3.0;
		}
		
		if(l == 2) {
			if(ax > 0.001) jl=(-3.0*cx/ax-sx*(1.0-3.0/ax2))/ax;
			else jl=ax2/15.0;
		}
		
		if(l == 3) {
			if(ax > 0.001) jl=(cx*(1.0-15.0/ax2)-sx*(6.0-15.0/ax2)/ax)/ax;
			else jl=ax*ax2/105.0;
		}
		
		if(l == 4) {
			if(ax > 0.001) jl=(sx*(1.0-45.0/(ax*ax)+105.0/(ax*ax*ax*ax)) +cx*(10.0-105.0/(ax*ax))/ax)/ax;
			else jl=ax2*ax2/945.0;
		}
		
		if(l == 5) {
			if(ax > 0.001) jl=(sx*(15.0-420.0/(ax*ax)+945.0/(ax*ax*ax*ax))/ax -cx*(1.0-105.0/(ax*ax)+945.0/(ax*ax*ax*ax)))/ax;
			else jl=ax2*ax2*ax/10395.0;
		}
		
		/********************** x=0 **********************/
	}
	else {
		
		if (ax < 1e-30) {      jl=0.0; }
		
		/*************** Region 1: x << l ****************/
		else {
			
			if (ax <= fl+0.5-llimit) {
				
				/* c	beta=acosh(nu/ax) */
				if (nu/ax < 1.0) printf("trouble with acosh\n");
				beta = log(nu/ax + sqrt((nu/ax)*(nu/ax) - 1.0) );
				/* 		!(4.6.21) */
				cotb=nu/sqrt(nu*nu-ax*ax);       /* cotb=coth(beta) */
				cot3b=cotb*cotb*cotb;
				cot6b=cot3b*cot3b;
				secb=ax/nu;
				sec2b=secb*secb;
				sec4b=sec2b*sec2b;
				sec6b=sec4b*sec2b;
				sum1=2.0+3.0*sec2b;
				expterm=sum1*cot3b/(24.0*nu);
				sum2=4.0+sec2b;
				expterm = expterm - sum2*sec2b*cot6b/(16.0*nu2);
				sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b;
				expterm = expterm - sum3*cot3b*cot6b/(5760.0*nu*nu2);
				sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b;
				expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2);
				expterm=exp(-nu*beta+nu/cotb-expterm);
				prefactor=sqrt(cotb/secb)/(2.0*nu);
				jl=prefactor*expterm;
			}
			/**************** Region 2: x >> l ****************/
			
			else {
				if (ax >= fl+0.5+ulimit) {
					
					beta=acos(nu/ax);
					cotb=nu/sqrt(ax*ax-nu*nu);      /* cotb=cot(beta) */
					cot3b=cotb*cotb*cotb;
					cot6b=cot3b*cot3b;
					secb=ax/nu;
					sec2b=secb*secb;
					sec4b=sec2b*sec2b;
					sec6b=sec4b*sec2b;
					trigarg=nu/cotb - nu*beta - PI/4.0;
					sum1=2.0+3.0*sec2b;
					trigarg = trigarg - sum1*cot3b/(24.0*nu);
					sum3=16.0-1512.0*sec2b-3654.0*sec4b-375.0*sec6b;
					trigarg = trigarg - sum3*cot3b*cot6b/(5760.0*nu*nu2);
					trigcos=cos(trigarg);
					sum2=4.0+sec2b;
					expterm=sum2*sec2b*cot6b/(16.0*nu2);
					sum4=32.0+288.0*sec2b+232.0*sec4b+13.0*sec6b;
					expterm = expterm - sum4*sec2b*cot6b*cot6b/(128.0*nu2*nu2);
					expterm=exp(-expterm);
					prefactor=sqrt(cotb/secb)/nu;
					jl=prefactor*expterm*trigcos;
				}
				
				/***************** Region 3: x near l ****************/
				
				else {
					
					beta=ax-nu;
					
					beta2=beta*beta;
					beta4=beta2*beta2;
					beta6=beta2*beta4;
					sx=6.0/ax;
					sx2=sx*sx;
					cx=sqrt(sx);
					
					secb=pow(sx,0.333333333);
					
					sec2b=secb*secb;
					
					deriv1=GAMMA1*secb;
					deriv1= deriv1+ beta*GAMMA2*sec2b;
					sum1=(beta2/6.0-1.0/15.0)*beta;
					deriv1 = deriv1 - sum1*sx*secb*GAMMA1/3.0;
					sum2=beta4/24.0-beta2/24.0+1.0/280.0;
					deriv1 = deriv1 - 2.0*sum2*sx*sec2b*GAMMA2/3.0;
					sum3=beta6/720.0-7.0*beta4/1440.0+beta2/288.0-1.0/3600.0;
					deriv1 = deriv1 + 4.0*sum3*sx2*secb*GAMMA1/9.0;
					sum4=(beta6/5040.0-beta4/900.0+19.0*beta2/12600.0-13.0/31500.0)*beta;
					deriv1 = deriv1 + 10.0*sum4*sx2*sec2b*GAMMA2/9.0;
					sum5=(beta4*beta4/362880.0-beta6/30240.0+71.0*beta4/604800.0-121.0*beta2/907200.0 + 7939.0/232848000.0)*beta;
					deriv1 = deriv1 - 28.0*sum5*sx2*sx*secb*GAMMA1/27.0;
					
					jl=deriv1*cx/(12.0*sqrtPI);
				}
			}
		}
	}
	
	return jl;
}

/*****************************************************************************************************************/
// this integrates the function given in tab fx sampled on grid in table x in size number of samples
// the x table and fx table must correspond to each other and the argumente of the functions must be given in increasing order
double cpeds_integrate_1d(long size, double *x, double *fx) {
	return cpeds_integrate_1d(size,x,fx,0,size-1);
}

// integrates the functions given in tab fx sampled on grid in table x in size number of samples
// in range from imin to imax.  size >= imax+1; imax>imin>=0
// the x table and fx table must correspond to each other and the argumente of the functions must be given in increasing order
double cpeds_integrate_1d(long size, double *x, double *fx, long imin, long imax) {
	long i;//,s=size-1;
	double r=0;
	
	if (imin>=imax) { printf("NOTE: cpeds_integrate_1d: minimal integration limit >= than upper limit: imin: %li imax: %li; will return 0\n",imin,imax); return 0; }
#pragma omp parallel for private(i) shared(imax,x,fx) reduction(+ : r)
	for (i=imin;i<imax;i++) {
		r=r+(x[i+1]-x[i]) * 0.5*(fx[i+1]+fx[i]);
	}
	return r;
}

/*****************************************************************************************************************/
double * cpeds_interpolate(double *X,double *Y, long N, double *Xint, long Nint, string type, bool check_points) {
	long i;
	double P;
	if (check_points && (type == "cspline_periodic" || type == "linear_periodic")) {
		P=X[N-1]-X[0];
		
		for (i=0; i<Nint;i++) {
			if (Xint[i]<X[0]) Xint[i]=P-Xint[i];
			if (Xint[i]>X[0]) Xint[i]=Xint[i]-P;
		}
	}
	return cpeds_interpolate(X,Y,N,Xint,Nint,type);
}
/***************************************************************************************/
double * cpeds_interpolate(double *X,double *Y, long N, double *Xint, long Nint, string type) {
	long i;
	//  long Nleo=N-1;
	double *Yint=NULL;
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline;
	
	if (type == "cspline") { spline = gsl_spline_alloc(gsl_interp_cspline, N); }
	if (type == "cspline_periodic") { spline = gsl_spline_alloc(gsl_interp_cspline_periodic, N); }
	if (type == "linear")  { spline = gsl_spline_alloc(gsl_interp_linear, N); }
	//  if (type == "linear_periodic")  { spline = gsl_spline_alloc(gsl_interp_linear_periodic, N+2); }
	if (type == "akima")  { spline = gsl_spline_alloc(gsl_interp_akima, N); }
	if (type == "akima_periodic")  { spline = gsl_spline_alloc(gsl_interp_akima_periodic, N); }
	if (type == "steffen")  { spline = gsl_spline_alloc(gsl_interp_steffen, N); }
	
	Yint = new double[Nint];
	
	
	if (type == "akima_periodic" or type == "cspline_periodic") {
		gsl_spline_init(spline, X, Y, N);
		for (i=0; i<Nint;i++) {
			Yint[i] = gsl_spline_eval(spline, Xint[i], acc);	  
		}
	}
	else {
		gsl_spline_init(spline, X, Y, N);
		for (i=0; i<Nint;i++) {
			if (Xint[i]<X[0]) { 
				Yint[i]=cpeds_extrapolate_linear(Xint[i],X[0],Y[0],X[1],Y[1]);
			}
			else
				if (Xint[i]>X[N-1]) {  
					Yint[i]=cpeds_extrapolate_linear(Xint[i],X[N-1],Y[N-1],X[N-2],Y[N-2]);
				}
				else
					Yint[i] = gsl_spline_eval(spline, Xint[i], acc);
		}	  
	}
	gsl_spline_free(spline);
	
	gsl_interp_accel_free(acc);
	
	return Yint;
}
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
void* cpeds_bicubic_interpolation_ccoef(double* y, double* y1, double* y2, double* y12, double d1, double d2, double (&c)[4][4]) {
	
	int l,k=0,j,i;
	double xx,d1d2=d1*d2;
	double cl[16];
	double x[16];
	//	static int wt[16][16]; for (i = 0; i < 16; i++) for (j = 0; j < 16; j++) { wt[i][j]=wt_d[k++]; }
	for (i=0;i<4;i++) {
		x[i]=y[i];
		x[i+4]=y1[i]*d1;
		x[i+8]=y2[i]*d2;
		x[i+12]=y12[i]*d1d2;
	}
	for (i=0;i<16;i++) {
		xx=0.0;
		for (k=0;k<16;k++) xx += cpeds_bicubic_interpolation_ccoef_wt_d[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) c[i][j]=cl[l++];
	
}


void* cpeds_bicubic_interpolation(double* y, double* y1, double* y2, double* y12, 
		const double x1l, const double x1u, const double x2l, const double x2u,
		const double x1, const double x2, double &ansy, double &ansy1, double &ansy2) {
	
	int i;
	double t,u,d1=x1u-x1l,d2=x2u-x2l;
	double c[4][4];
	cpeds_bicubic_interpolation_ccoef(y,y1,y2,y12,d1,d2,c);
	if (x1u == x1l || x2u == x2l)
		throw("Bad input in routine cpeds_bicubic_interpolation");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	ansy=ansy2=ansy1=0.0;
	for (i=3;i>=0;i--) {
		ansy=t*ansy+((c[i][3]*u+c[i][2])*u+c[i][1])*u+c[i][0];
		ansy2=t*ansy2+(3.0*c[i][3]*u+2.0*c[i][2])*u+c[i][1];
		ansy1=u*ansy1+(3.0*c[3][i]*t+2.0*c[2][i])*t+c[1][i];
	}
	ansy1 /= d1;
	ansy2 /= d2;
}


/***************************************************************************************/
/***************************************************************************************/
/*****************************************************************************************************************/
/* double* cpeds_interpolate_spline_on_sphere(long n, double* b, double* l, double* f,   long Nb, long Nl, double* B, double* L, int * err) { */

/*   double* F; */
/*   F=c_cssgridd((int)n, b,l,f,(int)Nb,(int)Nl, B,L,err); */
/*   return F; */
/* } */


/*****************************************************************************************************************/

/* yvals is the 1d array of size of the number of rows in the matrix that defines the yvalues for the corresponding  */
/* data in the matrix in each row. It doesn't have to be equally spaced, but it should be in increasing order. */


void cpeds_matrix_derivative(matrix<double>* m, double* yvals) {
	
	long i,j;
	long cols,rows,rowsleo,rowslet;
	
	cols = m->ColNo();
	rows = m->RowNo();
	rowsleo=rows-1;
	rowslet=rows-2;
	
	double * D = new double[cols];
	double dy1,dy2,Dj;
	
	// the first row
	dy1=yvals[1]-yvals[0];
	for (j=0;j<cols;j++) {
		D[j] = ( (*m)(1,j)-(*m)(0,j) )/ dy1;
	}
	
	// the i'th row
	
	for (i=1;i<rowsleo;i++) {
		dy1=yvals[i]-yvals[i-1];
		dy2=yvals[i+1]-yvals[i];
		for (j=0;j<cols;j++) {
			// store the derivative value from the previous row
			Dj=D[j];
			
			// derive the derivative value from the current row
			D[j] = 0.5*(  ( (*m)(i,j)-(*m)(i-1,j) )/ dy1 + ( (*m)(i+1,j)-(*m)(i,j) )/ dy2  );
			
			// overwrite the value from the row above back in the matrix with its derivative
			(*m)(i-1,j)=Dj;
			
		}
	}
	
	// the last row
	dy2=yvals[rowsleo]-yvals[rowslet];
	for (j=0;j<cols;j++) {
		// derive the derivative value from the current row
		(*m)(rowsleo,j) = ( (*m)(rowsleo,j)-(*m)(rowslet,j) )/ dy2;
		
		// overwrite the value from the row above back in the matrix with its derivative
		(*m)(rowslet,j)=D[j];
	}
	
	//  cpeds_print_matrix(*m);
	// clean up
	delete [] D;
	
	
}
/***************************************************************************************/
void cpeds_derivative(double* x, double* y, long size, double* periodX, double *periodY) {
	long i,j;
	long rows,rowsleo,rowslet;
	
	rows = size;
	rowsleo=rows-1;
	rowslet=rows-2;
	
	double D;
	double dx1,dx2,Dj;
	
	if (periodY==NULL) {
		if (periodX==NULL) {
			// the first value
			dx1=x[1]-x[0];
			D = ( y[1]-y[0] )/ dx1;
			// the i'th point
			for (i=1;i<rowsleo;i++) {
				dx1=x[i]-x[i-1];
				dx2=x[i+1]-x[i];
				// store the derivative value from the previous row
				Dj=D;
				// derive the derivative value from the current row
				D = 0.5*(  ( y[i]-y[i-1] )/ dx1 + ( y[i+1]-y[i] )/ dx2  );
				// overwrite the value from the row above back in the matrix with its derivative
				y[i-1]=Dj;
			}
			
			// the last row
			dx2=x[rowsleo]-x[rowslet];
			// derive the derivative value from the current row
			y[rowsleo] = ( y[rowsleo]-y[rowslet] )/ dx2;
			
			// overwrite the value from the row above back in the matrix with its derivative
			y[rowslet]=D;		 
		}
		else {
//			// the first value
//			dx1=x[1]-x[0];
			if (rows<2) return;
			long i,iprev,inext;
			double f0=y[0];
			double fprev = y[rowsleo];
			double f=y[0];
			double fnext=y[1];
//			double thres=1e-14;
			iprev=rowsleo;
			i=0;
			inext=1;
			dx1=x[0]-x[rowsleo]+*periodX;
//			if (dx1<thres) dx1=x[0]-x[rowsleo-1]+*periodX;
			dx2=x[1]-x[0];
			// derive the derivative value from the current row
			D = 0.5*(  ( f-fprev )/ dx1 + ( fnext-f )/ dx2  );
			i=0;
//			printf("%i x:%lf yprev:%lf y:%lf ynext:%lf dx1:%lf dx2:%lf D:%lf\n",i,x[i],fprev,f,fnext,dx1,dx2,D);
			// overwrite the value from the row above back in the matrix with its derivative
			y[i]=D;

			// the i'th point
			for (i=1;i<rowsleo;i++) {
				// store function values
				fprev=f;
				f=fnext;
				fnext=y[(i+1) % rows];

				dx1=x[i]-x[(i-1+rows) % rows];
				dx2=x[(i+1) % rows]-x[i];
				// derive the derivative value from the current row
				D = 0.5*(  ( f-fprev )/ dx1 + ( fnext-f )/ dx2  );
//				printf("%i x:%lf yprev:%lf y:%lf ynext:%lf dx1:%lf dx2:%lf D:%lf\n",i,x[i],fprev,f,fnext,dx1,dx2,D);
				// overwrite the value from the row above back in the matrix with its derivative
				y[i]=D;
			}

			// the last point
			i=rowsleo;
			fprev=f;
			f=fnext;
			fnext=f0;
			dx1=x[rowsleo]-x[rowslet];
			dx2=x[0]-x[rowsleo]+*periodX;
//			if (dx2<thres) dx2=x[1]-x[rowsleo]+*periodX;
			D = 0.5*(  ( f-fprev )/ dx1 + ( fnext-f )/ dx2  );
			i=rowsleo;
			// overwrite the value from the row above back in the matrix with its derivative
			y[i]=D;
//			printf("%i x:%lf yprev:%lf y:%lf ynext:%lf dx1:%lf dx2:%lf D:%lf\n",i,x[i],fprev,f,fnext,dx1,dx2,D);
						
		}
	}
	else {
		double D1,D2,d1,d2,d3;
		long i1,i2;
		// the first value
		dx1=x[1]-x[0];
		d1 = ( y[1]-y[0] )/ dx1;
		d2 = ( y[1]+*periodY-y[0] )/ dx1;
		d3 = ( y[1]-*periodY-y[0] )/ dx1;
		D=cpeds_getMinAbs3(d1,d2,d3);
#ifdef DEBUG
		printf("d1: %lE d2: %lE d3: %lE, D: %lE\n",d1,d2,d3,D);
#endif
		// the i'th point
		for (i=1;i<rowsleo;i++) {
			dx1=x[i]-x[i-1];
			dx2=x[i+1]-x[i];
			// store the derivative value from the previous row
			Dj=D;
			// derive the derivative value from the current row
			d1 = ( y[i]-y[i-1] )/ dx1;
			d2 = ( y[i]+*periodY-y[i-1] )/ dx1;
			d3 = ( y[i]-*periodY-y[i-1] )/ dx1;
			D1=cpeds_getMinAbs3(d1,d2,d3);
			d1 = ( y[i+1]-y[i] )/ dx2;
			d2 = ( y[i+1]+*periodY-y[i] )/ dx2;
			d3 = ( y[i+1]-*periodY-y[i] )/ dx2;
			D2=cpeds_getMinAbs3(d1,d2,d3);
			
			D = 0.5*(  D1 + D2  );
			// overwrite the value from the row above back in the matrix with its derivative
			y[i-1]=Dj;
		}
		
		// the last row
		dx2=x[rowsleo]-x[rowslet];
		// derive the derivative value from the current row
		d1 = ( y[rowsleo]-y[rowslet] )/ dx2;
		d2 = ( y[rowsleo]+*periodY-y[rowslet] )/ dx2;
		d3 = ( y[rowsleo]-*periodY-y[rowslet] )/ dx2;
		
		y[rowsleo] = cpeds_getMinAbs3(d1,d2,d3);
		
		// overwrite the value from the row above back in the matrix with its derivative
		y[rowslet]=D;
		
	}
	
	
}
/***************************************************************************************/
double cpeds_getMinAbs3(double v1, double v2, double v3) {
	double D=0;
	if (fabs(v1)<fabs(v2)) {
		if (fabs(v1)<fabs(v3)) D=v1; else D=v3;
	}
	else { 
		if (fabs(v2)<fabs(v3)) D=v2; else D=v3;
	}
	return D;
}
/*****************************************************************************************************************/
void cpeds_RT32_Model4(double AZ, double ZD, double *dAZ, double *dZD) {
	
	//const double p[13] = {
	//-0.032367, -0.049374, 0.051222, -0.054595, -0.000161, 0.002850,
	//0.001370, 0.020659, -0.012753, 0.008571, 0.004812, -0.006485, 0.028568 };
	
	const double p[13] = { -3.2367E-02,-4.9374E-02, 5.1222E-02,-5.4595E-02,
			-1.6113E-04, 2.8503E-03, 1.3701E-03, 2.0659E-02,-1.2753E-02,
			8.5710E-03, 4.8125E-03,-6.4854E-03, 2.8568E-02 };
	
	double pi180 = M_PI/180.;
	
	AZ = AZ * pi180;
	double Alt = (90.-ZD) * M_PI/180.;
	
	double xi = p[4]*pi180;
	double zeta = p[5]*pi180;
	double sigma = p[2]*pi180;
	double beta = p[3]*pi180;
	double sh = sin(Alt);
	double ch = cos(Alt);
	
	double sinxi = sin(xi);
	double sinzeta = sin(zeta);
	double se = sqrt(sinxi*sinxi + sinzeta*sinzeta);
	
	if (ch < 0.) { se = -se; }
	
	double ce = sqrt(1.-se*se);
	double alfa = atan2(sin(zeta), sin(xi));
	double AT = atan2(ch*sin(alfa-AZ), sh*se-ch*ce*cos(alfa-AZ)) -
			atan2(sin(alfa), -ce*cos(alfa));
	
	double hT = asin(ce*sh + se*ch*cos(alfa-AZ));
	double AT_Az = AT - AZ;
	if (ch < 0.) {
		AT_Az = M_PI - AT_Az;
		hT = M_PI - hT;
	}
	
	*dAZ = cpeds_RT4_arcsin((sin(sigma)*sin(hT) + sin(beta)) /
			(cos(hT)*cos(sigma)));
	
	double hb = atan2(sin(hT)*cos(sigma) + cos(hT)*sin(sigma)*sin(*dAZ),
			cos(hT)*cos(*dAZ));
	
	*dAZ = fmod((*dAZ + AT_Az)*180./M_PI, 360.) +
			p[0] + 
			p[8]*sin(2.*AZ) +
			p[9]*cos(2.*AZ) +
			p[11]*sin(3.*AZ)*sin(Alt) + 
			p[12]*cos(Alt)*cos(AZ/4.);
	
	*dZD = (hb-Alt) * 180./M_PI +
			p[1] +
			p[6]*ch +
			p[7]*sh +
			p[10]*sin(2.*AZ);
	*dZD = -(*dZD);
}
/***************************************************************************************/
void cpeds_RT32_Model4e(double AZ, double ZD, double *dAZ, double *dZD) {
	const double p[9] = { 
			3.4144712346E-03, -9.1488156466E-04, -2.2984989887E-03, 1.8860129779E-02, 
			-4.3813651976E-02, -1.0356952146E-02, 8.2631028153E-03, -7.4978342202E-03,
			5.0519549549E-03
	}; // deg
	const double q[7] = { 
			8.2681943479E-02, 4.2711366808E-05, 2.7049249587E-04, -8.7588060472E-04, 
			-3.4899753854E-02, -4.1424115143E-03, 3.6971967552E-03
	}; // deg
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	double MA=0,MZ=0;
	double A=AZ*PI180;
	double Z=ZD*PI180;
	double sinA=sin(A);
	double cosA=cos(A);
	double sin2A=sin(2.0*A);
	double cos2A=cos(2.0*A);
	double sinZ=sin(Z);
	double cosZ=cos(Z);
	
	// Model4e
	MA+=p[0]+( p[1]*sinA-p[2]*cosA+p[3] )/tan(Z) + p[4]/sinZ;
	MA+=p[5]*sin2A + p[6]*cos2A + p[7]*sin(3.0*A)*cosZ + p[8]*cos(A/4)*sinZ;
	
	MZ+=q[0]+q[1]*cosA+q[2]*sinA+q[3]*sinZ+q[4]*cosZ+q[5]*sin2A+q[6]*cos2A;
	
	*dAZ=MA;
	*dZD=MZ;
}
/***************************************************************************************/
/***************************************************************************************/
void cpeds_RT32_Model5(double AZ, double ZD, double *dAZ, double *dZD) {
	const double deltaA_A_phi_T[51][3] = {
			{9.1464679913E-06,0.0,0.0},
			{7.4812314678E-04,-1.0146440002E+00,3.5964783888E+02},
			{2.1099373994E-04,2.1395266021E+00,1.7982391944E+02},
			{7.5789593125E-04,3.7969801449E-01,1.1988261296E+02},
			{5.6486694523E-04,-2.5560328871E-01,8.9911959720E+01},
			{2.9594479076E-04,2.9213364814E+00,7.1929567776E+01},
			{1.5902554032E-03,2.0405774586E+00,5.9941306480E+01},
			{5.1419002198E-04,-2.6892085396E+00,5.1378262697E+01},
			{2.4555827774E-04,-2.6879404830E+00,4.4955979860E+01},
			{7.0900164890E-04,-2.4479912145E+00,3.9960870987E+01},
			{1.1665254860E-03,-1.8508451646E+00,3.5964783888E+01},
			{2.1747761860E-04,5.0597093916E-01,3.2695258080E+01},
			{3.4186615678E-04,9.6981714032E-01,2.9970653240E+01},
			{2.2048263759E-04,-1.0939918206E+00,2.7665218375E+01},
			{3.4358437983E-03,-2.6928190457E+00,2.5689131349E+01},
			{4.6753622473E-04,-2.8479383897E+00,2.3976522592E+01},
			{7.7632348788E-05,-1.4847739440E+00,2.2477989930E+01},
			{6.8024460573E-05,-2.2613095138E+00,2.1155755228E+01},
			{1.5780593336E-04,1.8409436599E+00,1.9980435493E+01},
			{1.5954792635E-04,-2.4045974145E+00,1.8928833625E+01},
			{1.5467645909E-04,-2.6642158778E+00,1.7982391944E+01},
			{1.8160142754E-04,2.3099465958E-01,1.7126087566E+01},
			{1.0035780153E-04,7.1399949804E-01,1.6347629040E+01},
			{1.0490219191E-04,5.5218331689E-01,1.5636862560E+01},
			{2.7213169735E-04,1.9217781072E+00,1.4985326620E+01},
			{2.8771359850E-05,2.2834929006E+00,1.4385913555E+01},
			{2.7055732417E-04,2.9936859119E+00,1.3832609188E+01},
			{1.1862895918E-04,-2.5581520809E+00,1.3320290329E+01},
			{7.9562328734E-05,-2.6941803972E+00,1.2844565674E+01},
			{1.8689729643E-04,-2.3185283681E+00,1.2401649617E+01},
			{8.0638670934E-05,-1.6240672909E+00,1.1988261296E+01},
			{1.8136793028E-04,-1.2375765810E+00,1.1601543190E+01},
			{1.6136708224E-04,-2.3978970641E+00,1.1238994965E+01},
			{1.2495378166E-04,6.0103896559E-01,1.0898419360E+01},
			{2.4145259184E-04,1.5398138470E+00,1.0577877614E+01},
			{2.2646069620E-04,2.0528858146E+00,1.0275652539E+01},
			{1.6499465983E-04,2.4674367246E+00,9.9902177466E+00},
			{1.2696608887E-04,-2.2342792246E+00,9.7202118616E+00},
			{1.8477887487E-04,-2.4500988264E+00,9.4644168126E+00},
			{2.1960431737E-04,-6.3960034023E-01,9.2217394584E+00},
			{1.0188446233E-04,-4.9474299296E-01,8.9911959720E+00},
			{2.1053922078E-04,-2.9518911738E+00,8.7718985092E+00},
			{1.4852096038E-03,2.2113244978E+00,8.5630437828E+00},
			{1.7162835158E-04,2.0271244979E+00,8.3639032297E+00},
			{1.1407557322E-04,1.5943235394E+00,8.1738145200E+00},
			{8.0313176267E-05,-2.9625652152E-01,7.9921741973E+00},
			{1.3482542279E-04,2.4574276126E+00,7.8184312800E+00},
			{1.3187051810E-04,-2.9775507141E+00,7.6520816783E+00},
			{1.3109424134E-04,-2.0719697657E+00,7.4926633100E+00},
			{1.3625561070E-04,-1.4991715456E+00,7.3397518139E+00},
			{6.3555475863E-05,1.5036934663E+00,7.1929567776E+00}};
	
	const double deltaZ_A_phi_T[51][3] = {
			{6.1031365843E-05,0.0,0.0},
			{5.8974717125E-04,1.7811473351E+00,3.5990974471E+02},
			{4.0854852247E-04,-2.6003551057E+00,1.7995487235E+02},                                                                                                 
			{3.4867701410E-03,1.5128086420E+00,1.1996991490E+02},                                                                                                  
			{7.8911821183E-04,1.9672588206E+00,8.9977436177E+01},                                                                                                  
			{6.9451672748E-04,-9.5592799286E-01,7.1981948942E+01},                                                                                                 
			{9.4651234414E-04,1.6532224066E+00,5.9984957452E+01},
			{4.0666916301E-04,2.6466373050E+00,5.1415677816E+01},
			{3.2267969308E-04,2.9303803378E+00,4.4988718089E+01},
			{1.1016303575E-03,-1.6497812365E+00,3.9989971634E+01},
			{1.0202613956E-03,-2.2047342894E+00,3.5990974471E+01},
			{2.3225751584E-04,-1.2424765126E+00,3.2719067701E+01},
			{2.5049337354E-04,2.2845075865E+00,2.9992478726E+01},
			{1.8276649393E-04,-2.6514459471E-01,2.7685364978E+01},
			{2.1991151848E-03,-2.6982071766E+00,2.5707838908E+01},
			{4.5762406162E-04,-3.9979919812E-01,2.3993982981E+01},
			{2.8537630491E-04,-4.3498146135E-01,2.2494359044E+01},
			{2.5601999128E-04,-9.5130402198E-01,2.1171161454E+01},
			{3.4131811378E-04,-2.6255431693E+00,1.9994985817E+01},
			{4.0785879696E-04,-2.9737380475E+00,1.8942618143E+01},
			{7.0419095424E-04,-2.9510473627E+00,1.7995487235E+01},
			{2.1776010725E-04,3.0084577024E+00,1.7138559272E+01},
			{3.6958604739E-05,-2.1961016080E+00,1.6359533850E+01},
			{3.6834571411E-04,4.9972426188E-02,1.5648249770E+01},
			{2.8645969469E-04,3.1083111914E-01,1.4996239363E+01},
			{3.0475897324E-04,9.4605913134E-02,1.4396389788E+01},
			{1.3821256068E-04,-1.2051938560E+00,1.3842682489E+01},
			{1.8788813573E-04,2.7052042092E+00,1.3329990545E+01},
			{1.5646459375E-04,-2.2659842917E+00,1.2853919454E+01},
			{4.7487733425E-04,-2.9503798986E+00,1.2410680852E+01},
			{1.7631710246E-04,-2.6475751136E+00,1.1996991490E+01},
			{8.5467498390E-05,-2.4013118071E+00,1.1609991765E+01},
			{3.3306771678E-04,-4.8102299618E-01,1.1247179522E+01},
			{3.2772469936E-04,-3.0997894318E-02,1.0906355900E+01},
			{3.2294873161E-04,2.0693484746E-01,1.0585580727E+01},
			{1.6459431559E-04,1.7563465081E+00,1.0283135563E+01},
			{2.0195800309E-04,-2.6410653936E+00,9.9974929086E+00},
			{3.8195320183E-04,-2.9463440179E+00,9.7272903976E+00},
			{3.3587917938E-04,-2.5961626323E+00,9.4713090713E+00},
			{1.7860306368E-04,-3.0561297896E+00,9.2284549926E+00},
			{3.1128186374E-04,2.3303732166E-01,8.9977436177E+00},
			{9.5328259266E-05,-3.1221250565E+00,8.7782864563E+00},
			{7.4697876787E-04,1.8867612798E+00,8.5692796359E+00},
			{2.0636941761E-04,-2.3087059333E+00,8.3699940630E+00},
			{1.4737174254E-04,-6.7120776279E-01,8.1797669252E+00},
			{2.8705226183E-05,-1.0659823844E+00,7.9979943269E+00},
			{9.9876318620E-05,9.8926442465E-01,7.8241248850E+00},
			{2.0844522051E-04,-2.9464022958E+00,7.6576541428E+00},
			{2.3603680020E-04,2.4934551141E+00,7.4981196815E+00},
			{4.7268543178E-05,2.2828282702E+00,7.3450968308E+00},
			{1.6503420594E-04,6.5779151054E-01,7.1981948942E+00}};
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	double MA=0,MZ=0;
	
	// Model4e
	cpeds_RT32_Model4e(AZ,ZD,&MA,&MZ);
	
	// Model5
	MA+=deltaA_A_phi_T[0][0];
	MZ+=deltaZ_A_phi_T[0][0];
	for (unsigned long i = 1; i < 51; i++) {
		MA+=deltaA_A_phi_T[i][0]*sin(-twoPI/deltaA_A_phi_T[i][2]*(AZ-180)+deltaA_A_phi_T[i][1]);
		MZ+=deltaZ_A_phi_T[i][0]*sin(-twoPI/deltaZ_A_phi_T[i][2]*(AZ-180)+deltaZ_A_phi_T[i][1]);
	}
	
	*dAZ=MA;
	*dZD=MZ;
}

/***************************************************************************************/
/***************************************************************************************/
void cpeds_RT32_Model6(double AZ, double ZD, double *dAZ, double *dZD) {
	static const double p[9] = { 
			2.0926054359E-03, -6.5575253165E-04, -2.3636332835E-03, 2.2555977599E-02, 
			-4.7203651029E-02, -1.0162365385E-02, 8.3838275454E-03, -7.2000506137E-03, 
			9.1692472305E-03
	}; // deg
	static const double q[7] = { 
			7.9005272694E-02, -5.3791653813E-04, 5.0630943699E-04, -8.5484780147E-04, 
			-2.9886032652E-02, -3.1614987640E-03, 3.0010137421E-03
	}; // deg
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	double PI180=M_PI/180l;

	double MA=0,MZ=0;
	double A=AZ*PI180;
	double Z=ZD*PI180;
	double sinA=sin(A);
	double cosA=cos(A);
	double sin2A=sin(2.0*A);
	double cos2A=cos(2.0*A);
	double sinZ=sin(Z);
	double cosZ=cos(Z);
	
	// Model6 (same as 4e)
	MA+=p[0]+( p[1]*sinA-p[2]*cosA+p[3] )/tan(Z) + p[4]/sinZ;
	MA+=p[5]*sin2A + p[6]*cos2A + p[7]*sin(3.0*A)*cosZ + p[8]*cos(A/4)*sinZ;
	
	MZ+=q[0]+q[1]*cosA+q[2]*sinA+q[3]*sinZ+q[4]*cosZ+q[5]*sin2A+q[6]*cos2A;
	
	*dAZ=MA;
	*dZD=MZ;
}


/***************************************************************************************/
void cpeds_RT32_Model6r(double AZ, double ZD, double *dAZ, double *dZD) {
	static const double model6r_deltaA_A_phi_T[51][3] = {
	{-5.4052091778E-06,0.0,0.0},
	{1.0894708728E-04,4.1341562864E-01,3.5992281292E+02},
	{2.5215619290E-04,5.6822121996E-01,1.7996140646E+02},
	{4.6095613516E-04,2.9441955486E-01,1.1997427097E+02},
	{4.1367769878E-04,-4.8621443899E-01,8.9980703231E+01},
	{3.0088509691E-04,3.0252058169E+00,7.1984562584E+01},
	{1.0620559418E-03,2.1745347939E+00,5.9987135487E+01},
	{3.8878223208E-04,-3.1313025514E+00,5.1417544703E+01},
	{2.7951425081E-04,-2.5798792344E+00,4.4990351615E+01},
	{5.3999757599E-04,-2.1778876462E+00,3.9991423658E+01},
	{7.2679600946E-04,-1.6936917634E+00,3.5992281292E+01},
	{2.5529784004E-04,2.8282366612E-01,3.2720255720E+01},
	{2.9106472682E-04,7.6485700787E-01,2.9993567744E+01},
	{2.7525716401E-04,-2.0080083260E+00,2.7686370225E+01},
	{2.1933041120E-03,-2.6542419960E+00,2.5708772352E+01},
	{5.3872506774E-04,-2.8539218640E+00,2.3994854195E+01},
	{1.0181699361E-05,-9.6905888687E-01,2.2495175808E+01},
	{7.1469068319E-05,2.8442997769E-01,2.1171930172E+01},
	{8.9047483028E-05,1.1656911441E+00,1.9995711829E+01},
	{1.0655760963E-04,-2.1767394894E+00,1.8943305943E+01},
	{9.9651970291E-05,-2.5115746428E+00,1.7996140646E+01},
	{1.2821948351E-04,2.1994039634E-01,1.7139181568E+01},
	{6.5351241669E-05,5.4780523997E-01,1.6360127860E+01},
	{6.1162895620E-05,8.6691756528E-01,1.5648817953E+01},
	{1.8855743883E-04,1.9798472634E+00,1.4996783872E+01},
	{6.7719475044E-05,1.9927536319E+00,1.4396912517E+01},
	{1.8931207812E-04,2.9568950411E+00,1.3843185112E+01},
	{1.0078220794E-04,-2.8920413310E+00,1.3330474553E+01},
	{4.0767554300E-05,-2.5376680403E+00,1.2854386176E+01},
	{1.0928749466E-04,-2.0401718229E+00,1.2411131480E+01},
	{9.0257528881E-05,-1.0519689560E+00,1.1997427097E+01},
	{1.5212451485E-04,-1.1138011133E+00,1.1610413320E+01},
	{9.2411836450E-05,-2.3740420852E+00,1.1247587904E+01},
	{7.6153459087E-05,1.5302000410E+00,1.0906751907E+01},
	{1.9626085038E-04,1.8300021595E+00,1.0585965086E+01},
	{1.8727756508E-04,2.1621727848E+00,1.0283508941E+01},
	{1.4232890629E-04,2.4931938421E+00,9.9978559145E+00},
	{9.4138802171E-05,-2.2083334463E+00,9.7276435925E+00},
	{1.4084332031E-04,-1.9481298436E+00,9.4716529716E+00},
	{1.9536033543E-04,-6.8838359495E-01,9.2287900749E+00},
	{8.2446588601E-05,-9.6259838995E-01,8.9980703231E+00},
	{2.4400134467E-04,2.9526307151E+00,8.7786051932E+00},
	{9.8072023132E-04,2.3307532858E+00,8.5695907839E+00},
	{2.1390897307E-04,2.1074334733E+00,8.3702979749E+00},
	{6.6522008905E-05,1.0801414425E+00,8.1800639300E+00},
	{9.9529905510E-05,-5.5180598001E-01,7.9982847316E+00},
	{3.5209705412E-05,-2.7263349507E+00,7.8244089766E+00},
	{9.3834835918E-05,-2.6139615990E+00,7.6579321898E+00},
	{8.9500640427E-05,-1.9104066915E+00,7.4983919359E+00},
	{9.8419481445E-05,-1.3570099126E+00,7.3453635290E+00},
	{3.8416705956E-05,1.9380291669E+00,7.1984562584E+00}};	

	static const double model6r_deltaZ_A_phi_T[51][3] = {
	{-8.1365002104E-05,0.0,0.0},
	{6.8153143625E-04,1.8566471516E+00,3.5993330163E+02},
	{8.0276553853E-04,6.2346768146E-01,1.7996665081E+02},
	{3.1353249053E-03,1.5060743322E+00,1.1997776721E+02},
	{5.7128028648E-04,2.2463673038E+00,8.9983325407E+01},
	{7.3862779417E-04,-1.0521020455E+00,7.1986660326E+01},
	{1.0171517085E-03,1.6184551146E+00,5.9988883605E+01},
	{4.8066628699E-04,2.4825157208E+00,5.1419043090E+01},
	{4.3295070971E-04,2.6460429499E+00,4.4991662704E+01},
	{9.9303708243E-04,-1.6620022346E+00,3.9992589070E+01},
	{9.9181240767E-04,-2.1467225385E+00,3.5993330163E+01},
	{3.1478679972E-04,-1.2304888317E+00,3.2721209239E+01},
	{2.0794967254E-04,2.5071517278E+00,2.9994441802E+01},
	{1.8700613324E-04,-4.3572740395E-01,2.7687177048E+01},
	{2.2634330276E-03,-2.6627200568E+00,2.5709521545E+01},
	{4.5512255245E-04,-4.8007446007E-01,2.3995553442E+01},
	{2.9639338714E-04,-4.0903472673E-01,2.2495831352E+01},
	{2.2743062840E-04,-7.5328898543E-01,2.1172547155E+01},
	{2.9219228125E-04,-2.6014615213E+00,1.9996294535E+01},
	{3.7548801499E-04,-3.0930441145E+00,1.8943857980E+01},
	{6.9911991950E-04,-2.9613359633E+00,1.7996665081E+01},
	{2.3590137335E-04,-3.0445841103E+00,1.7139681030E+01},
	{1.2020535175E-04,-1.8541619814E+00,1.6360604619E+01},
	{3.5662490523E-04,6.0632793772E-02,1.5649273984E+01},
	{2.7604153874E-04,3.1457129769E-01,1.4997220901E+01},
	{3.1017334266E-04,2.5615070172E-01,1.4397332065E+01},
	{9.7878427064E-05,-9.6559042031E-01,1.3843588524E+01},
	{1.6903054283E-04,2.5837365446E+00,1.3330863023E+01},
	{1.0566498470E-04,-2.4298101658E+00,1.2854760772E+01},
	{4.4697467042E-04,-3.0075360508E+00,1.2411493160E+01},
	{1.8852485924E-04,-2.5695966992E+00,1.1997776721E+01},
	{1.1373460908E-04,-2.2696712066E+00,1.1610751665E+01},
	{3.3902886963E-04,-4.4427765759E-01,1.1247915676E+01},
	{3.1193327869E-04,-1.1247843592E-01,1.0907069746E+01},
	{3.1384620483E-04,2.5157359262E-01,1.0586273577E+01},
	{1.3616519105E-04,1.9158236019E+00,1.0283808618E+01},
	{2.0379036342E-04,-2.6192127362E+00,9.9981472675E+00},
	{3.5173645627E-04,-2.8788011663E+00,9.7279270710E+00},
	{3.0866427479E-04,-2.6351402381E+00,9.4719289902E+00},
	{1.4946135630E-04,3.0691900119E+00,9.2290590161E+00},
	{3.1927234141E-04,3.0976592517E-01,8.9983325407E+00},
	{9.5800191188E-05,-2.7916637092E+00,8.7788610153E+00},
	{7.3610353401E-04,1.9615116112E+00,8.5698405150E+00},
	{2.3859212363E-04,-2.1458609443E+00,8.3705418983E+00},
	{1.7218889751E-04,-7.0497940155E-01,8.1803023097E+00},
	{3.8685203394E-05,-5.1909926658E-01,7.9985178140E+00},
	{9.2599695927E-05,1.2572882738E+00,7.8246369919E+00},
	{1.9018201078E-04,-2.9726178182E+00,7.6581553538E+00},
	{2.3462698474E-04,2.4735954629E+00,7.4986104506E+00},
	{3.2524418222E-05,1.7518975356E+00,7.3455775843E+00},
	{1.6871184968E-04,5.9613845640E-01,7.1986660326E+00}};
	
	double twoPI=M_PI*2.0;

	double MA=0,MZ=0;

	int i;
	
	//
	// check the input
	//
	if (AZ>180) AZ-=360;
	if (AZ<-180) AZ+=360;
	if (ZD>90) ZD=90;
	if (ZD<0 ) ZD=0.1;
	
	// Model6
	cpeds_RT32_Model6(AZ,ZD,&MA,&MZ);
	
	// Model6r
	MA+=model6r_deltaA_A_phi_T[0][0];
	MZ+=model6r_deltaZ_A_phi_T[0][0];
	for (i = 1; i < 51; i++) {
		MA+=model6r_deltaA_A_phi_T[i][0]*sin(-twoPI/model6r_deltaA_A_phi_T[i][2]*(AZ-180)+model6r_deltaA_A_phi_T[i][1]);
		MZ+=model6r_deltaZ_A_phi_T[i][0]*sin(-twoPI/model6r_deltaZ_A_phi_T[i][2]*(AZ-180)+model6r_deltaZ_A_phi_T[i][1]);
	}
	
	*dAZ=MA;
	*dZD=MZ;
}

/***************************************************************************************/
/*!
	\brief this is a function returning the corrections in A,h for requested A,h which are used in the control system meteofix
	\details 
	@param
	@return

	\date Jun 5, 2013, 6:57:13 PM
	\author Bartosz Lew
 */
double cpeds_RT4_arcsin(double x) {
	if (x > 1.) x = 1;
	if (x < -1.) x = -1.;
	return asin(x);
}
