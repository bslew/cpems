#include "cpeds-math.h"
#include "cpeds-direction.h"
#include <libnova/transform.h>
#include <libnova/sidereal_time.h>
#include <libnova/aberration.h>
#include <libnova/nutation.h>
#include "novas.h"

/* ************************************************************************************************************************************************************************************ */
DirectionAh::DirectionAh(const DirectionRaDec& d, const cpedsDirection& observer, double JD, bool localTime, bool refract, bool checkDirection) {
  *this=d.toAh(observer,JD,localTime,refract);
  if (checkDirection) check();
}

/* ************************************************************************************************************************************************************************************ */
cpedsDirection DirectionRaDec::toGal() const {
//#ifdef WITH_NOVAS
	double glon,glat;
	double iradec[2];
	double iradecICRS[2];
	iradec[0]=ra()*PI180inv/15.0; // convert to hours
	iradec[1]=dec()*PI180inv; // convert to degrees
	frame_tie(iradec,-1,iradecICRS); // convert to ICRS system
	equ2gal(iradecICRS[0], iradecICRS[1],&glon, &glat);
//	equ2gal(iradec[0], iradec[1],&glon, &glat);
	return cpedsDirection(glon,glat)*PI180;  
//#else
//  ln_equ_posn equ;
//  equ.ra=ra()*PI180inv;
//  equ.dec=dec()*PI180inv;
//
//  ln_gal_posn gal;
//
//  ln_get_gal_from_equ2000(&equ,&gal);
//  return cpedsDirection(gal.l,gal.b,val())*PI180;  
//#endif
}
/* ************************************************************************************************************************************************************************************ */
const DirectionRaDec& DirectionRaDec::fromGal(const cpedsDirection& galactic) {

  ln_equ_posn equ;
  ln_gal_posn gal;
  gal.l=galactic.lon()*PI180inv;
  gal.b=galactic.lat()*PI180inv;

  ln_get_equ2000_from_gal(&gal,&equ);
  // ln_get_equ_from_gal(&gal,&equ);
  setRa(equ.ra*PI180);
  setDec(equ.dec*PI180);
  return *this;
}

/* ************************************************************************************************************************************************************************************ */
DirectionRaDec& DirectionRaDec::aberrate(double JD) {
  ln_equ_posn mean;
  mean.ra=ra()*PI180inv;
  mean.dec=dec()*PI180inv;
  ln_equ_posn tmp;
  
  ln_get_equ_aber(&mean,JD, &tmp);
  setRa(tmp.ra*PI180);
  setDec(tmp.dec*PI180);
  return *this;
}


/* ************************************************************************************************************************************************************************************ */
DirectionRaDec& DirectionRaDec::unaberrate(double JD) {
	  ln_equ_posn mean;
	  mean.ra=ra()*PI180inv;
	  mean.dec=dec()*PI180inv;
	  ln_equ_posn tmp,delta;
	  
	  ln_get_equ_aber(&mean,JD, &tmp);
	  delta.ra=tmp.ra*PI180-ra();
	  delta.dec=tmp.dec*PI180-dec();
	  setRa(ra()-delta.ra);
	  setDec(dec()-delta.dec);
	  return *this;
}

/* ************************************************************************************************************************************************************************************ */
DirectionRaDec& DirectionRaDec::toJD(double jd, bool nutate) {
	  ln_equ_posn mean;
//	  double ratmp,dectmp;
	  mean.ra= ra()*PI180inv;
	  mean.dec=dec()*PI180inv;
	  ln_equ_posn pos;
	  
	  if (nutate) {
//		  ln_get_equ_prec2(&mean, epoch(), 2433283.0, &pos);
		  ln_get_equ_prec2(&mean, epoch(), jd, &pos);
		  ln_lnlat_posn ecl;
		  equ2ecl(jd,0,0,pos.ra/15.0,pos.dec,&ecl.lng,&ecl.lat);
//		  ln_get_ecl_from_equ(&pos,jd,&ecl);

		  ln_nutation nutation;
//		  ln_get_nutation(jd,&nutation);
//		  ecl.lng+=nutation.longitude;
//		  ecl.lat+=nutation.obliquity;		  
//		  printf("libnova nutation [\"]: %lf %lf\n",nutation.longitude*3600, nutation.obliquity*3600);
		  iau2000a(jd,0,&nutation.longitude,&nutation.obliquity);
		  ecl.lng+=nutation.longitude*PI180inv;
		  ecl.lat+=nutation.obliquity*PI180inv;		  
//		  printf("novas nutation [\"]: %lf %lf\n",nutation.longitude*PI180inv*3600, nutation.obliquity*PI180inv*3600);
//		  printf("mean obliquity [deg]: %lf\n",nutation.ecliptic);
//		  ln_get_nutation(jd+3000.0,&nutation);
//		  printf("mean obliquity [deg]: %lf\n",nutation.ecliptic);
//		  ln_get_equ_from_ecl(&ecl,jd,&pos);
//		  mean=pos;
//		  ln_get_equ_prec2(&mean, 2433283.0, jd, &pos);
//		  ecl2equ_vec()
	  }
	  else {
		  ln_get_equ_prec2(&mean, epoch(), jd, &pos);		  
	  }
	  	
	  set(pos.ra*PI180,pos.dec*PI180);
	  setEpoch(jd);
	  return *this; 
}
/***************************************************************************************/
DirectionRaDec& DirectionRaDec::toJD(double jd) {
	  double pos1[3];
	  double pos2[3];
	  radec2vector(ra()*PI180inv/15.0,dec()*PI180inv,1.0,pos1);
	  precession(epoch(),pos1,CPEDS_JD2000,pos2);
	  precession(CPEDS_JD2000,pos2,jd,pos1);
	  double ra,dec;
	  vector2radec(pos1,&ra,&dec);
	  setRa(ra*15.0*PI180);
	  setDec(dec*PI180);
	  setEpoch(jd);
	  return *this;
}
/***************************************************************************************/
DirectionRaDec& DirectionRaDec::nutate(int direction) {
	  double pos1[3]; // geocentric equatorial rectangular coordinates
	  double pos2[3]; // geocentric equatorial rectangular coordinates
	  ln_lnlat_posn ecl;
	  ln_nutation n;
	  double RA,DEC;
	  double Elon,Elat;
	  double jd=epoch();

//      Convert equatorial spherical coordinates to a vector (equatorial
//      rectangular coordinates).

	  radec2vector(ra()*PI180inv/15.0,dec()*PI180inv,1.0,pos1);
//      Nutate equatorial rectangular coordinates from mean equator and
//      equinox of epoch to true equator and equinox of epoch. Inverse
//      transformation may be applied by setting flag 'direction'.
	  
	  nutation(jd,direction,0,pos1,pos2); 
//		  ln_get_nutation(jd,&nutation);
//		  ecl.lng+=nutation.longitude;
//		  ecl.lat+=nutation.obliquity;		  
//		  printf("libnova nutation [\"]: %lf %lf\n",nutation.longitude*3600, nutation.obliquity*3600);
//	  iau2000b(jd,0,&n.longitude,&n.obliquity);
//	  nu2000k(jd,0,&n.longitude,&n.obliquity);
//	  printf("JD: %.15lf\n",jd);
//	  printf("JD2000: %.15lf\n",CPEDS_JD2000);
	  vector2radec(pos2,&RA,&DEC);
//	  gcrs2equ(jd,1,0,RA,DEC,&RA,&DEC);
	  
/*
	  equ2ecl(jd,0,0,ra()*PI180inv/15.0,dec()*PI180inv,&Elon,&Elat);
//	  radec2vector(ra()*PI180inv/15.0,dec()*PI180inv,1.0,pos1);
//	  equ2ecl_vec(jd,0,0,pos1,pos2);

	  nutation_angles((jd-CPEDS_JD2000)/(36525.0),0,&n.longitude,&n.obliquity); // ["]
	  double mObliq=mean_obliq(jd);
	  printf("novas nutation [\"]: %lf %lf, obliquity [\"]: %lE\n",n.longitude, n.obliquity,mObliq);
//	  ln_nutation lnnutation;
//	  ln_get_nutation(jd,&lnnutation);
//	  printf("libnova obliquity: %lE\n",lnnutation.ecliptic*3600);

	  
	  Elon+=n.longitude/3600.0;
	  Elat+=n.obliquity/3600.0;
	  radec2vector(Elon/15.0,Elat,1.0,pos1);
	  ecl2equ_vec(jd,0,0,pos1,pos2);
	  vector2radec(pos2,&RA,&DEC);
*/
	  
	  
	  
	  
	  setRa(RA*15.0*PI180);
	  setDec(DEC*PI180);
	  return *this;

}


/* ************************************************************************************************************************************************************************************ */
DirectionAh DirectionRaDec::toAh(const cpedsDirection& observer, double JDutc, double ut1_utc, double polar_x, double polar_y, double P, double T, bool refract) const {
//	double tt = JDutc + double(CPEDS_TAI_UTC + CPEDS_TT_TAI_OFFSET)/86400.0;
	// ut1 = utc + dAT +32.184 - dt
//		printf("tt: %15.8lf\n",tt);
	double jd_ut1 = JDutc + ut1_utc/86400.0;
	double dt=CPEDS_TAI_UTC + CPEDS_TT_TAI- ut1_utc;
//	printf("jd_ut1: %15.8lf\n",jd_ut1);
//	printf("CPEDS_TAI_UTC: %lf\n",CPEDS_TAI_UTC);
	on_surface obs;
	obs.longitude=observer.lon()*PI180inv;
	obs.latitude=observer.lat()*PI180inv;
	obs.height=observer.val();
	obs.pressure=P;
	obs.temperature=T;
//	printf("%lf %lf\n",obs.longitude, obs.latitude);
//	printf("ra: %lf, dec: %lf\n",ra()*PI180inv, dec()*PI180inv);
	short int refr=0;
	if (refract) refr=2; else refr=0;
	double zd;
	double az;
	double rar,decr;
	equ2hor(jd_ut1,dt,0,polar_x,polar_y,&obs,ra()*PI180inv/15.0, dec()*PI180inv,refr,&zd,&az,&rar,&decr);	
//	printf("az %lf el: %lf\n",az, 90.0-zd);
	DirectionAh ah(az*PI180, PIsnd-zd*PI180);
	return ah;
}


/***************************************************************************************/
DirectionAh DirectionRaDec::toAh(const cpedsDirection& observer, double JD, bool localTime, bool refract) const {
  if (localTime) JD-=observer.lon()/twoPI; // convert to UT time

  
  ln_equ_posn object;
  object.ra = ra()*PI180inv;
  object.dec=dec()*PI180inv;

  ln_lnlat_posn obs;
  obs.lng=observer.lon()*PI180inv;
  obs.lat=observer.lat()*PI180inv;
  
  ln_hrz_posn position;
  
  double sidereal = ln_get_mean_sidereal_time(JD); // BLcomment (Jan 15, 2017, 4:57:56 PM): changed to mean LST

  ln_get_hrz_from_equ_sidereal_time(&object,&obs,sidereal,&position);
  if (refract) DirectionAh(position.az*PI180-PI, position.alt*PI180).refract();
  return DirectionAh(position.az-180, position.alt)*PI180;
}
/***************************************************************************************/
/***************************************************************************************/
double DirectionRaDec::VradialLSR(double JDutc) {
	DirectionRaDec radec(*this); 
	double rah,ram,ras,decd,decm,decs;
//	double JD=cpeds_julian_time(2014,1,5,14.0);
	double JD0=radec.epochYr();

//	radec*=PI180;
	radec.toJD(JDutc,true);
	
	double aLMST, Vsun,Vobs,Vdop, ra,dec;
	ra=radec.ra();
	dec=radec.dec();
	cpeds_vdrt32_(&JDutc,&ra,&dec,&aLMST,&Vsun,&Vobs,&Vdop);
	
//	printf("JD: %.10lf\n",JD);
//	radec.print_direction("requested direction",true);
//	printf("aLMST [deg]: %.10lf\n",aLMST*PI180inv);
//	printf("Vsun [km/s]: %.10lf\n",Vsun);
//	printf("Vobs [km/s]: %.10lf\n",Vobs);
//	printf("Vdop [km/s]: %.10lf\n",Vdop);
	
	return Vdop;
}

/* ************************************************************************************************************************************************************************************ */
DirectionRaDec DirectionAh::toRaDec(const cpedsDirection& observer, double JD, bool localTime, bool aberrate) const {
	  if (localTime) JD-=observer.lon()/twoPI; // convert to UT time

	  ln_hrz_posn object;
	  object.az=(A()-PI)*PI180inv;
	  object.alt=h()*PI180inv;

	  ln_lnlat_posn obs;
	  obs.lng=observer.lon()*PI180inv;
	  obs.lat=observer.lat()*PI180inv;
	  
	  ln_equ_posn position;

	  ln_get_equ_from_hrz(&object,&obs,JD,&position);
	  if (aberrate) return DirectionRaDec(position.ra*PI180, position.dec*PI180).aberrate(JD);
	  return DirectionRaDec(position.ra, position.dec,JD)*PI180;
}
/* ************************************************************************************************************************************************************************************ */
DirectionRaDec DirectionAh::toRaDec(const cpedsDirection& observer, double JD, double ut1_utc, double DeltaAT, bool localTime, double polar_x, double polar_y, double P, double T, bool refract, int LSTtype) const {
//	double tt_ut1=DeltaAT+CPEDS_TT_TAI_OFFSET-ut1_utc;

	if (localTime) JD-=observer.lon()/twoPI; // convert to JD utc

//	  double REC,DEC,ST;

	  cpedsDirection hadec = cpeds_horizontalToEquatorialFirst(A(),h(),observer.lat());
	  double lon_sec = observer.lon()*PI180inv/15*3600;
	  double lst=cpeds_local_sidereal_time_novas(JD+ut1_utc/86400.0,ut1_utc,DeltaAT,lon_sec,LSTtype)/86400*twoPI; // [rad]
//	  printf("lst novas: %.15lf\n",lst*PI180inv/15*3600);
//	  lst=cpeds_local_sidereal_time(JD+lon_sec/86400+ut1_utc/86400.0,observer.lon()*PI180inv); // this is consistent with novas LST
//	  printf("lst cpeds: %.15lf\n",lst*PI180inv/15*3600);
//	  lst=ln_get_apparent_sidereal_time(JD+lon_sec/86400+ut1_utc/86400.0)*3600/86400*twoPI;
//	  printf("lst nova app: %.15lf\n",lst*PI180inv/15*3600);
//	  lst=ln_get_mean_sidereal_time(JD+lon_sec/86400+ut1_utc/86400.0)*3600/86400*twoPI;
//	  printf("lst nova mean: %.15lf\n",lst*PI180inv/15*3600);
/*
	  cpeds lst is highly (<0.01 s) compatible with novas mean sidereal time with equonix method.
	  if localJulianTime is given with dUT1 correction.
	  The two are strongly (~2s) inconsistent with the ln_get_mean_sidereal_time of libnova.
*/
	  
	  cpedsDirection radec =cpeds_equatorialFirstToEquatorialSecond(hadec.get_direction(),lst);
	  
	  return DirectionRaDec(radec.lon(), radec.lat(),JD);
	  
//	  on_surface obs;
//	  obs.longitude=observer.lon()*PI180inv;
//	  obs.latitude=observer.lat()*PI180inv;
//	  obs.height=observer.val();
//	  obs.pressure=P;
//	  obs.temperature=T;
//	  
//	  double gst;
//	  sidereal_time(JD,0,delta_t,1,0,0,&gst);
//	  double pos[3],vel[3],pos2[3];
//	  terra(&obs,gst,pos,vel); // convert geodetic to geocentric location
//	  ter2cel(JD,0,delta_t,1,0,1,polar_x,polar_y,pos,pos2); // rotate from earth-fixed system to space fixed system applying corrections for polar motion, earth rotation, nutation, precession etc.
//	  double ra,dec;
//	  vector2radec(pos2,&ra,&dec);
//	  return DirectionRaDec(ra*15.0, dec,JD)*PI180;
}
/***************************************************************************************/
cpedsDirection& cpedsDirection::toGeographic(double a, double b) {
	lat()=atan(a*a/b/b*tan(lat()));
	return *this;
}
/* ************************************************************************************************************************************************************************************ */
string cpedsDirection::print_direction(string comment, bool show, int format) const {
  double h1,m1,s1;
  double h2,m2,s2;
  if (format==0 or format==1) {
	  if (format==0) {
		  cpeds_angToDMS(lon()*PI180inv,&h1,&m1,&s1);
		  cpeds_angToDMS(lat()*PI180inv,&h2,&m2,&s2);
	  }
	  if (format==1) {
		  cpeds_angToHMS(lon()*PI180inv,&h1,&m1,&s1);
		  cpeds_angToDMS(lat()*PI180inv,&h2,&m2,&s2);
	  }
  }
  if (format==2 or format==3) {
	  if (format==2) {
		  cpeds_angToDMS(lon(),&h1,&m1,&s1);
		  cpeds_angToDMS(lat(),&h2,&m2,&s2);
	  }
	  if (format==3) {
		  cpeds_angToHMS(lon(),&h1,&m1,&s1);
		  cpeds_angToDMS(lat(),&h2,&m2,&s2);		  
	  }
  }
  string s;
  char tmpch[1000];
  if (format==0 or format==2)
	  if (format==0)
		  sprintf(tmpch,"%s> lon [deg]: %lf lat [deg]: %lf val: %lE    ||| lon %li [d] %li [m] %lE [s], lat  %li [d] %li [m] %lE [s] \n", comment.c_str(), cpeds_check_phi(lon())*PI180inv, cpeds_check_b(lat())*PI180inv, val(),long(h1),long(m1),s1,long(h2),long(m2),s2);  
	  else
		  sprintf(tmpch,"%s> lon [deg]: %lf lat [deg]: %lf val: %lE    ||| lon %li [d] %li [m] %lE [s], lat  %li [d] %li [m] %lE [s] \n", comment.c_str(), cpeds_check_phi(lon()*PI180)*PI180inv, cpeds_check_b(lat()*PI180)*PI180inv, val(),long(h1),long(m1),s1,long(h2),long(m2),s2);  
  if (format==1 or format==3)
	  if (format==1)
		  sprintf(tmpch,"%s> lon [deg]: %lf lat [deg]: %lf val: %lE    ||| lon %li [h] %li [m] %lE [s], lat  %li [d] %li ['] %lE [\"] \n", comment.c_str(), cpeds_check_phi(lon())*PI180inv, cpeds_check_b(lat())*PI180inv, val(),long(h1),long(m1),s1,long(h2),long(m2),s2);  
	  else
		  sprintf(tmpch,"%s> lon [deg]: %lf lat [deg]: %lf val: %lE    ||| lon %li [h] %li [m] %lE [s], lat  %li [d] %li ['] %lE [\"] \n", comment.c_str(), cpeds_check_phi(lon()*PI180inv)*PI180, cpeds_check_b(lat()*PI180)*PI180inv, val(),long(h1),long(m1),s1,long(h2),long(m2),s2);  
  if (show) printf("%s",tmpch);
  s=tmpch;
  return s;
}
cpedsDirection& cpedsDirection::subtract(const cpedsDirection& rhs) {
	lon()-=rhs.lon();
	lat()-=rhs.lat();
	return *this;
}
cpedsDirection& cpedsDirection::add(const cpedsDirection& rhs) {
	lon()+=rhs.lon();
	lat()+=rhs.lat();
	return *this;
}
/* ************************************************************************************************************************************************************************************ */

string DirectionAh::print_direction(string comment, bool show) const { 
  double d1,m1,s1;
  double d2,m2,s2;
  cpeds_angToDMS(A()*PI180inv,&d1,&m1,&s1);
  cpeds_angToDMS(h()*PI180inv,&d2,&m2,&s2);
  string s;
  char tmpch[1000];
  sprintf(tmpch,"%s> A [deg]: %lf h [deg]: %lf   ||| A %li [deg] %li ['] %lE [\"], h %li [deg] %li ['] %lE [\"]\n", comment.c_str(), A()*PI180inv, h()*PI180inv,long(d1),long(m1),s1,long(d2),long(m2),s2);  
  if (show) printf("%s",tmpch);
  s=tmpch;
  return s;
}
/* ************************************************************************************************************************************************************************************ */
string DirectionRaDec::print_direction(string comment, bool show) const {
  double h1,m1,s1;
  double h2,m2,s2;
  cpeds_angToHMS(ra()*PI180inv,&h1,&m1,&s1);
  cpeds_angToDMS(dec()*PI180inv,&h2,&m2,&s2);
  string s;
  char tmpch[1000];
  sprintf(tmpch,"%s> Ra [deg]: %lf Dec [deg]: %lf Epoch: %lf  val: %lE ||| Ra %li [h] %li [m] %lE [s], Dec %li [deg] %li ['] %lE [\"]\n", comment.c_str(), ra()*PI180inv, dec()*PI180inv, epoch(), val(), long(h1),long(m1),s1,long(h2),long(m2),s2);  
  if (show) printf("%s",tmpch);
  s=tmpch;
  return s;
}
/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
