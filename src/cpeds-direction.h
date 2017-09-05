/*!
  \file implements the class cpeds_direction
*/
/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <libnova/refraction.h>
#include <libnova/precession.h>
#include <QtCore/QPointF>
#include "cpeds-common.h"
#include "cpeds-consts.h"
#include "cpeds-math.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */
class  DirectionRaDec;

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
#ifndef __directionRADEC__
#define __directionRADEC__
/*!
  \struct directionRADEC
  \brief stores a direction in equatorial 2nd system
  \details 
  
  \date 2009/05/26 21:21:42 
  \author Bartosz Lew
*/
typedef struct { 
  double ra; //!< right ascension in equatorial system
  double dec; //!< declination in equatorial system
} directionRADEC;

#endif
/* **************************************************************************************************** */
#ifndef __directionAh__
#define __directionAh__
/*!
  \struct directionAh
  \brief stores a direction in horizontal system
  \details 
  
  \date 2009/05/26 21:22:31 
  \author Bartosz Lew
*/
typedef struct { 
  double A; // azimuth in horizontal system
  double h; // elevation in horizontal system
} directionAh;
#endif

/* **************************************************************************************************** */
#ifndef __directionLONLAT__
#define __directionLONLAT__
/*!
  \struct directionLONLAT
  \brief geographical location container, or a general direction container
  \details 
  
  \date 2009/05/26 21:23:20 
  \author Bartosz Lew
*/
typedef struct { 
  double lon; //!< longitude
  double lat; //!< latitude
} directionLONLAT;

//extern directionLONLAT cpeds_directionLONLAT_default;//={0,0};

#endif


/* **************************************************************************************************** */
#ifndef __directionRaDecAh__
#define __directionRaDecAh__
/*!
  \struct directionRaDecAh
  \brief a general directions container for two coordinate systems
  \details 
  The azimuth is assumed to be zero in the north.
  h is assumed to be zero at the horizon.
  The convention is that the angles inside the library are all in radians
  Degrees are sometimes used with communication with the user.

  \date 2009/05/26 21:24:08 
  \author Bartosz Lew
*/
typedef struct { 
  double ra; //!< right ascension in equatorial system
  double dec; //!< declination in equatorial system
  double A; //!< azimuth in horizontal system
  double h; //!< elevation in horizontal system
} directionRaDecAh;
#endif



#ifndef CPEDS_DIRECTION_CLASS
#define CPEDS_DIRECTION_CLASS



/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsDirection
  \brief Encapsulates  cpeds_direction container %class

  \details 

  \date 2009/06/30 20:01:54 
  \author Bartosz Lew
*/
class cpedsDirection {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:
  
  typedef cpeds_direction cpeds_direction_type;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  cpedsDirection() { lon()=0.0; lat()=0.0; setVal(0.0); }
  cpedsDirection(const cpedsDirection& d) { _d=d._d; setVal(d.val()); }
  cpedsDirection(cpeds_direction d) { _d=d; setVal(0.0); }
  cpedsDirection(double lon, double lat) { _d.l=lon; _d.b=lat; setVal(0.0); }
  cpedsDirection(double lon, double lat,  double v) { _d.l=lon; _d.b=lat; setVal(v); }
  cpedsDirection(directionLONLAT n) { _d.l=n.lon; _d.b=n.lat; setVal(0.0);}
  virtual ~cpedsDirection() {};
  

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  double lon() const { return _d.l; }
  double l() const { return _d.l; }
  double& lon() { return _d.l; }
  double setLon(const double lon) { return _d.l=lon; }
  double lat() const { return _d.b; }
  double b() const { return _d.b; }
  double& lat() { return _d.b; }
  double setLat(const double lat) { return _d.b=lat; }
  void setVal(const double v) { _v=v; }
  double& value() { return _v; }
  double val() const { return _v; }
//  double getVal() const { return _v; }
  virtual cpedsDirection& set(double lon, double lat, double val=0) { _d.l=lon; _d.b=lat; _v=val; return *this; }
  QPointF toQPointF() const { return QPointF(l(),b()); }

  virtual cpeds_direction get_direction() const { return _d; }
  virtual cpeds_direction set_direction(cpeds_direction d) {  return _d=d; }
  directionLONLAT getLONLAT() const { directionLONLAT n; n.lon=lon(); n.lat=lat(); return n; }


  /*!
    \brief This and the following routines rotate the direction p about the X axis by Ax [rad]
    \details 
    @param p - direction to be rotated
    @param Ax - angle in rad by which to rotate
    @return rotated direction
    
    This routine does not change the object direction.

    \date 2009/08/06 00:06:49 
    \author Bartosz Lew
  */
  cpeds_direction Rx(cpeds_direction p,double Ax) const {
    return cpeds_Rx(p,Ax);
  }

  cpeds_direction Ry(cpeds_direction p,double Ay) const {
    return cpeds_Ry(p,Ay);
  }
  
  cpeds_direction Rz(cpeds_direction p,double Az) const {
    return cpeds_Rz(p,Az);
  }

  /*!
    \brief This and the following routines rotate the direction p about the X axis by Ax [rad]
    \details 
    The X axis is defined as the one that is in direction l=90 deg and b=0 deg i.e. it points to the East at the horizon in the Ah CS
    @param p - a reference to an object to be rotated; The p object will not be changed
    @param Ax - angle in rad by which to rotate
    @return The rotated result reference and the result itself is stored in the object
    
    \date 2009/08/06 00:09:13 
    \author Bartosz Lew
  */
  cpedsDirection& Rx(cpedsDirection& p,double Ax) { set_direction( Rx(p.get_direction(),Ax) ); return *this; }
  cpedsDirection& Ry(cpedsDirection& p,double Ay) { set_direction( Ry(p.get_direction(),Ay) ); return *this; }
  cpedsDirection& Rz(cpedsDirection& p,double Az) { set_direction( Rz(p.get_direction(),Az) ); return *this; }

  /*!
    \brief This and the following routines rotate the direction about the X axis by Ax [rad]
    \details 
    @param Ax - angle in rad by which to rotate
    @return The rotated result reference and the result itself is stored in the object. The information
    about the initial direction is lost.
    
    \date 2009/08/06 00:09:13 
    \author Bartosz Lew
  */
  cpedsDirection& Rx(double Ax) { set_direction( Rx(get_direction(),Ax)); return *this; }
  cpedsDirection& Ry(double Ay) { set_direction( Ry(get_direction(),Ay)); return *this; }
  cpedsDirection& Rz(double Az) { set_direction( Rz(get_direction(),Az)); return *this; }


  /*!
  \brief calculates the cross product of two directions in the sky.
  \details 
  @param n1, n2 - directions in the sky
  @return The result is the direction perpendicular to n1 and n2 directions; If n1==n2 then (l,b)=(0,0) is returned.
  
  \date 2009/10/23 16:53:55 
  \author Bartosz Lew
  */
  const cpedsDirection crossProduct(const cpedsDirection& n1, const cpedsDirection& n2) const {
    if (n1==n2) return cpedsDirection();
    double x1=cpeds_sph2cart(0, PIsnd-n1.lat(), n1.lon());
    double y1=cpeds_sph2cart(1, PIsnd-n1.lat(), n1.lon());
    double z1=cpeds_sph2cart(2, PIsnd-n1.lat(), n1.lon());
    double x2=cpeds_sph2cart(0, PIsnd-n2.lat(), n2.lon());
    double y2=cpeds_sph2cart(1, PIsnd-n2.lat(), n2.lon());
    double z2=cpeds_sph2cart(2, PIsnd-n2.lat(), n2.lon());
    double x3,y3,z3;
    /* printf("xyz: %lE, %lE, %lE\n",x1,y1,z1); */
    /* printf("xyz: %lE, %lE, %lE\n",x2,y2,z2); */
    cpeds_cross_product_3d(x1,y1,z1,x2,y2,z2,&x3,&y3,&z3);
    /* printf("xyz: %lE, %lE, %lE\n",x3,y3,z3); */
    cpedsDirection(cpeds_cart2sph(1,x3,y3,z3),PIsnd-cpeds_cart2sph(0,x3,y3,z3)).print_direction();
    return cpedsDirection(cpeds_cart2sph(1,x3,y3,z3),PIsnd-cpeds_cart2sph(0,x3,y3,z3));
  }
  
  cpedsDirection& subtract(const cpedsDirection& rhs);
  cpedsDirection& add(const cpedsDirection& rhs);
  
  /*!
  \brief calculates the dot product of two directions in the sky.
  \details 
  @param n1, n2 - directions in the sky
  @return The result is the direction perpendicular to n1 and n2 directions
  
  */
  double dotProduct(const cpedsDirection& n1, const cpedsDirection& n2) const {
    return cpeds_dot_product_3d(n1.get_direction(),n2.get_direction());
  }


/*   /\*! */
/*     \brief calculates the intersection between the two great circles defined by n1 and n2 */
/*     \details  */
/*     @param n1,n2 - directions perpendicular to the great circles in the sky */
/*     @return */
    
/*     \date 2009/10/26 14:00:16  */
/*     \author Bartosz Lew */
/*   *\/ */
/*   cpedsDirection great_circles_intersection(const cpedsDirection n1, const cpedsDirection n2); */

  /*!
    \brief derives the angle between spherical line elements
    \details 
    @param
    @return returns the smaller of the two angles defined as angles between the intersecting lines, tangent to the sphere in the point of intersection

    All angles are in radians

    \date 2009/10/23 17:22:34 
    \author Bartosz Lew
  */
  double angle_between_great_circles(const cpedsDirection n11, const cpedsDirection n12, const cpedsDirection n21, const cpedsDirection n22) {
    cpedsDirection p1(crossProduct(n11,n12));
    cpedsDirection p2(crossProduct(n21,n22));
    
    p1.print_direction();
    p2.print_direction();
/*     // generate great circle that include n11 and n12 */
/*     skyMotion circ1(); */
/*     // generate great circle that include n21 and n22 */
/*     skyMotion circ2(); */

/*     //find the crossing point  -- no this is done smarter now, look below*/
    if (p1==p2) return 0;
    return cpeds_ang_n1n2(p1.get_direction(),p2.get_direction());
  }

  
  /*!
	\brief convert the direction (assumed to be geocentric longitude,latitude) to geographic longitude,latidute on the surface of ellipsoid
	\details 
	@param a - semi-major axis of the Earth ellipsoid
	@param b - semi-minor axis of the Earth ellipsoid
	@return

	The conversion does not account for possible altitude above the ellipsoid.
	The longitude in both coordinate systems is the same.
	See also:
	http://www.usno.navy.mil/USNO/earth-orientation/eo-info/general/conv-1996
	
	
	\date Dec 6, 2012, 1:39:24 PM
	\author Bartosz Lew
   */
  cpedsDirection& toGeographic(double a=6378137, double b=6356752.314);


  /*!
    \brief returns the angle between this and n direction
    \details 
    @param reference direction for angle calculation [rad]
    @return the angular separation in radians
    
    \date 2009/12/11 22:30:27 
  */
  virtual double angle(const cpedsDirection& n) const { return cpeds_ang_n1n2(PIsnd-lat(),lon(),PIsnd-n.lat(),n.lon()); }
  /*!
	\brief 
	\details 
	@param commnet- a string to print alongside the direction
	@param show - can be set false to suppress printout if only a returned string is needed
	@param format - \n
		0 - DMSDMS\n
		1 - HMSDMS\n
		2 - DMSDMS input in degrees\n
		3 - HMSDMS input in degrees\n
		
	@param 
	@return

	\date Dec 13, 2012, 1:49:32 PM
	\author Bartosz Lew
   */
  virtual string print_direction(string comment="", bool show=true, int format=0) const;

  virtual bool operator==(const cpedsDirection &rhs) const { if (lon()!=rhs.lon()) return false;  if (lat()!=rhs.lat()) return false; return true; } 
  virtual bool operator!=(const cpedsDirection &rhs) const { return !(*this==rhs); } 
  virtual const cpedsDirection& operator=(double rhs) { setLon(rhs);   setLat(rhs);  return *this; }
  virtual const cpedsDirection& operator=(const cpeds_direction rhs) { setLon(rhs.l);   setLat(rhs.b);  return *this; } 
  virtual const cpedsDirection& operator=(const cpedsDirection& rhs) { setLon(rhs.lon()); setLat(rhs.lat()); setVal(rhs.val());  return *this; } 
  virtual const cpedsDirection& operator*= (const double rhs) { setLon(lon()*rhs); setLat(lat()*rhs); return *this; } 
  const cpedsDirection operator* (const double rhs) const { return cpedsDirection(lon(), lat(),val())*=rhs; }
  virtual double operator* (const cpedsDirection& rhs) const { return dotProduct(*this,rhs); }
  virtual const cpedsDirection operator^ (const cpedsDirection& rhs) const { return crossProduct(*this,rhs); } //!< This cross product operator. The direction vecrot is cross multiplied with rhs and the result is returned. 
  virtual cpedsDirection operator- (const cpedsDirection& rhs) const { cpedsDirection tmp(*this); tmp.subtract(rhs); return tmp; }
  virtual cpedsDirection operator+ (const cpedsDirection& rhs) const { cpedsDirection tmp(*this); tmp.add(rhs); return tmp; }
  virtual const cpedsDirection& operator^= (const cpedsDirection& rhs)  { return *this=crossProduct(*this,rhs); } //!< This cross product operator. The direction vecrot is cross multiplied with rhs and the result is returned. The object itself gets cross-multiplied
  virtual const cpedsDirection& operator/= (const double rhs) { setLon(lon()/rhs); setLat(lat()/rhs); return *this; } 
  const cpedsDirection operator/ (const double rhs) { return cpedsDirection(lon(), lat(),val())/=rhs; }

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:
  double _v;
  cpeds_direction_type _d;

/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:

  double x,y,z;
/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };
#endif /* CPEDS_DIRECTION_CLASS */ 
























#ifndef DIRECTIONS
#define DIRECTIONS

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
//#include "scan-common.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */




































/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class DirectionAh
  \brief Encapsulates directionAh structure
  \details 
  
  \date 2009/07/01 17:46:20 
  \author Bartosz Lew
*/
class DirectionAh : public cpedsDirection {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
 DirectionAh() : cpedsDirection() { }
 DirectionAh(double A, double h) : cpedsDirection(A,h) { }
 DirectionAh(cpeds_direction n) : cpedsDirection(n) { }
 DirectionAh(const cpedsDirection& n) : cpedsDirection(n) { }
 DirectionAh(directionAh n) : cpedsDirection(n.A, n.h) { }
  //! performs an implicit conversion from RaDec to Ah for a given JD and observer's location
  DirectionAh(const DirectionRaDec& d, const cpedsDirection& observer, double JD, bool localTime=false, bool refract=false, bool checkDirection=true);
 DirectionAh(const DirectionAh& n) : cpedsDirection(n.A(),n.h(),n.val()) { }
  virtual ~DirectionAh() {};

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  directionAh getAh() const { directionAh n; n.A=lon(); n.h=lat(); return n; }
  double A() const { return lon(); } 
  double h() const { return lat(); } 
  double setA(double A) { return setLon(A); } 
  double seth(double h) { return setLat(h); } 

  /*!
    \brief returns the angle between this and n direction
    \details 
    @param reference direction for angle calculation
    @return the angular separation in radians
    
    \date 2009/12/11 22:30:27 
  */
  virtual double angle(const DirectionAh& n) const { return cpeds_ang_n1n2(PIsnd-h(),A(),PIsnd-n.h(),n.A()); }


  /*!
    \brief convert the direction to RaDec from A,h
    \details 
    @param observer - a location of the observer in [rad]
    @param JD - julian date (corresponding to lon=0 - UT time I think)
    @param localTime - specifies whether the supplied JD time is local (true) or UT (false)
    @param aberrate - if true, then the result will be corrected for aberration effect
    @return cpedsDirection object that stores the Ra and Dec of the stored A,h in the current object for the indicated location and JD [rad]

    The conversion takes into account the precession and nutation.\n
    It does not account for the atmospheric refraction nor aberration etc.
    
    \date 2009/12/10 14:48:12 
    \author Bartosz Lew
  */
  DirectionRaDec toRaDec(const cpedsDirection& observer, double JD, bool localTime=false, bool aberrate=false) const;
  /*!
	\brief convert the direction to RaDec from A,h using topo_star function from novas library
	\details 
    @param observer - a location of the observer in [rad]
    @param JD - JD_utc time scale (corresponding to observers location or UT, depending on localTime option)
    @param ut1_utc - UT1-UTC in seconds of time
	@param DeltaAT - TAI-UTC [integer s] - number of leap seconds since the beginning of TAI
    @param localTime - specifies whether the supplied JD time is local (true) or UT (false)
	@param polar_x - polar motion - currently not used
	@param polar_y - polar motion - currently not used
    @param P - pressure [mbar]
    @param T - temperature [degC]
    @param refract - if true then the direction should be first unrefracted usisng the P,T data but this is not 
    implemented now. So the default option is refract false and it is assumed that the A,h direction is in-space
	@return

	\date Jan 10, 2013, 1:55:06 PM
	\author Bartosz Lew
   */
  DirectionRaDec toRaDec(const cpedsDirection& observer, double JD, double ut1_utc=0, double DeltaAT=37, bool localTime=false, double polar_x=0, double polar_y=0, double P=1012, double T=0, bool refract=false) const;

  
  //! pressure is given in millibars=100 Pa=100N/m^2 and temperature in Celsius degrees
  double getRefraction(double pressure=1013, double temperature=20) const { return PI180*ln_get_refraction_adj(h()*PI180inv,pressure,temperature); }
  //! include the effect of the atmospheric refraction to the Ah direction.
  DirectionAh& refract(double pressure=1013, double temperature=20) { seth(h()-getRefraction(pressure,temperature)); return *this; }
  //! remove the effect of the atmospheric refraction to the Ah direction.
  DirectionAh& unrefract(double pressure=1013, double temperature=20) { seth(h()+getRefraction(pressure,temperature)); return *this; }

  //! checks the direction for the right values: l in <0,twoPI>, b in <-PIsnd,PIsnd>
  virtual DirectionAh& check() { double atmp,htmp; atmp=A(); htmp=h(); cpeds_check_bl(&htmp,&atmp); setA(atmp); seth(htmp); return *this; }

  virtual string print_direction(string comment="", bool show=true) const;

  virtual const DirectionAh& operator=(const directionAh rhs) { setA(rhs.A); seth(rhs.h); return *this; } 
  virtual const DirectionAh& operator=(const DirectionAh& rhs) { setA(rhs.A()); seth(rhs.h()); return *this; } 
  virtual const DirectionAh& operator=(const cpedsDirection& rhs) { setA(rhs.lon()); seth(rhs.lat()); return *this; } 
  virtual const DirectionAh& operator*= (const double rhs) { setLon(lon()*rhs); setLat(lat()*rhs); return *this; } 
  const DirectionAh operator* (const double rhs) { return DirectionAh(lon(), lat())*=rhs; }
  virtual const DirectionAh& operator/= (const double rhs) { setLon(lon()/rhs); setLat(lat()/rhs); return *this; } 
  const DirectionAh operator/ (const double rhs) { return DirectionAh(lon(), lat())/=rhs; }


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
//  directionAh _Ah;

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };
























/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class  DirectionRaDec
  \brief Encapsulates  directionRaDec structure
  \details 
  
  \date 2009/07/01 17:46:20 
  \author Bartosz Lew
*/
class  DirectionRaDec : public cpedsDirection {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
 DirectionRaDec() : cpedsDirection() { }
 DirectionRaDec(directionRADEC d) : cpedsDirection(d.ra,d.dec) { }
 DirectionRaDec(const cpedsDirection& n) : cpedsDirection(n) { }
 DirectionRaDec(double ra, double dec) : cpedsDirection(ra,dec) { }
 /*!
	\brief constructor for initializing ra dec and epoch in one line.
	\details 
	@param ra - ra in rad
	@param dec - dec in rad
	@param ep - epoch in JD

	\date Sep 5, 2012, 2:19:36 PM
	\author Bartosz Lew
  */
 DirectionRaDec(double ra, double dec, double ep) : cpedsDirection(ra,dec) { setEpoch(ep); }
 DirectionRaDec(double ra, double dec, double ep, double val) : cpedsDirection(ra,dec,val) { setEpoch(ep); }
 DirectionRaDec(cpeds_direction d) : cpedsDirection(d) { }
  //! performs implicit conversion from Ah to RaDec given the observer's location for a given UTC expressed in JD. All angles in radians
 DirectionRaDec(const DirectionAh& d, const cpedsDirection& observer, double JD, bool localTime=false, bool aberrate=false) { *this=d.toRaDec(observer,JD,localTime,aberrate); check();  }
 DirectionRaDec(const DirectionRaDec& n) : cpedsDirection(n.ra(),n.dec(),n.val()) { setEpoch(n.epoch()); }
  virtual ~  DirectionRaDec() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  double ra() const { return lon(); } 
  double dec() const { return lat(); } 
  double setRa(double ra) { return setLon(ra); } 
  double setDec(double dec) { return setLat(dec); }
  virtual DirectionRaDec& set(double ra, double dec) { cpedsDirection::set(ra,dec,val()); return *this; }


  //! These are provided for convenience so that this object could also serve as lon lat container
  directionRADEC getRADEC() const { directionRADEC n; n.ra=lon(); n.dec=lat(); return n; }
  cpeds_direction getlb() const { return _d; }

  //! set the epoch [year] for which the direction coordinates are given 
  void setEpochYr(long yy) { _epoch=cpeds_julian_time(yy,1,1,12.0); }
  //! set the epoch [JD] for which the direction coordinates are given 
  void setEpoch(double e) { _epoch=e; }
  //! get the epoch [JD] for which the direction coordinates are given 
  double epoch() const { return _epoch; }
  double epochYr() const { return cpeds_JDToYear(_epoch); }

  
  //! converts the direction from the current epoch to the requested epoch (given in JD) and position in [rad]
  // nutation is not implemented correctly
  DirectionRaDec& toJD(double jd, bool nutate); 
  /*!
	\brief novas implementation of precession
	\details 
	@param jd - target JD
	@return

	\date Jun 12, 2016, 9:42:09 PM
	\author Bartosz Lew
   */
  DirectionRaDec& toJD(double jd); 

  /*!
	\brief calculate equatorial position corrected for nutation
	\details 
	@param direction - 0 - converts from mean to true coordinates;
	otherwise, converts from true to mean coordinates
	@return returns nutation corrected direction (rad)

	This function was also tested for compatibility with KB's simplified formulas
	for nutation (test program cosmocalc)
	
	BTW: KB's formulas are sufficiently accurate for RT32 but are much faster.
	
	\date Jun 12, 2016, 6:47:00 PM
	\author Bartosz Lew
   */
  DirectionRaDec& nutate(int direction); 
  
  
  /*!
    \brief convert the direction to Ah
    \details 
    @param observer - a location of the observer in [rad] (geodetic coordinates)
    @param JD - julian date (corresponding to lon=0 - UT time I think)
    @param localTime - specifies whether the supplied JD time is local (true) or UT (false)
    @param refract - if true, then the result will be corrected for refraction effect (using default values)
    @return cpedsDirection object that stores the A and h of the stored Ra,Dec in the current object for the indicated location and JD [rad]

    The conversion takes into account the precession and nutation.\n
    It does not account for the atmospheric refraction nor annual aberration etc.
    
    \date 2009/12/10 14:48:12 
    \author Bartosz Lew
  */
  DirectionAh toAh(const cpedsDirection& observer, double JD, bool localTime=false, bool refract=false) const;

  /*!
	\brief convert the  ra,dec [rad] to geodetic A,h
	\details 
    @param observer - a location of the observer in [rad] (geodetic coordinates). Value should hold the altitude above ellipsoid [m]
    @param JD - julian date (UTC time in JD)
	@param ut1_utc - difference between ut1 and utc time. the default value is assumed for 2012/12/06 [s]
	@param polar_x - polar motion in x (see: http://maia.usno.navy.mil/ser7/ser7.dat and http://maia.usno.navy.mil/ser7/mark3.out)
	@param polar_y - polar motion in y
	@param P - pressure [mbar]
	@param T - temperature [Celsius]
	@return

	\date Dec 10, 2012, 11:04:30 AM
	\author Bartosz Lew
   */
  DirectionAh toAh(const cpedsDirection& observer, double JDutc, double ut1_utc=0.300522, double polar_x=0, double polar_y=0, double P=1012, double T=0, bool refract=false) const;
  /*!
    \brief returns the galactic coordinates from RADec given for epoch J2000 (in J2000 system)
    \details 
    @param
    @return
    
    All angles of the directions are given and returned in radians
    The value is also passed to the returned object.
    If compiled with WITH_NOVAS compilation option, then first conversion form J2000 to ICRS system is done
    and then to galactic coordinates.
    
    Dec 12, 2012, 3:36:29 PM
    	Novas is now normally compiled in. 

    \date 2009/12/10 16:54:57 
    \author Bartosz Lew
  */
  cpedsDirection toGal() const;

  /*!
    \brief sets the RaDec coordinates for epoch J2000 from galactic coordinates 
    \details 
    @param
    @return thes method returns reference to this
    
    All angles of the directions are given and returned in radians
    The value is NOT passed to the returned object

    \date 2009/12/10 16:54:57 
    \author Bartosz Lew
  */
  const DirectionRaDec& fromGal(const cpedsDirection& galactic);

  /*!
    \brief recalculate the position accounting for the annual aberration effect
    \details 
    @param JD - Julian Date of the coordinates (this should be the same as the value returned by epoch())
    @return returns the direction corrected for the aberration effects.
    This is still libnova implementation and does not account for relativistic effects.
        
    This mehtod was tested against KB's aberration routines and they agree at 0.05" level 
    (for a single tested direction and time).
    
    \date 2009/12/11 21:01:38 
    
    updated on: Jun 14, 2016, 1:52:00 PM -- the test against consistency with KB's aberration
    was done after correction the RT32 control system  implementation of that routine 
    for a number of mistakes due to rad/deg conventions.
     
  */
  DirectionRaDec& aberrate(double JD);

  /*!
    \brief subtract the effects of annual aberration form the current direction
    \details 
    @param JD - Julian Date of the coordinates (this should be the same as the value returned by epoch())
    @return
    
    It is assumed that the current direction is the apparent direction. 
    The application of aberrate() followed by unaberrate() will generally not return
    the input direction but the difference in practical applications are probably negligible,
    at the level of few mas.
    
    \date 2009/12/11 21:01:38 
  */
  DirectionRaDec& unaberrate(double JD);

  /*!
    \brief returns the angle between this and n direction
    \details 
    @param reference direction for angle calculation
    @return the angular separation in radians
    
    \date 2009/12/11 22:30:27 
  */
  virtual double angle(const DirectionRaDec& n) const { return cpeds_ang_n1n2(PIsnd-dec(),ra(),PIsnd-n.dec(),n.ra()); }

  //! checks the direction for the right values: RA in <0,twoPI>, Dec in <-PIsnd,PIsnd>
  virtual DirectionRaDec& check() { double ratmp,dectmp; ratmp=ra(); dectmp=dec(); cpeds_check_bl(&dectmp,&ratmp); setRa(ratmp); setDec(dectmp); return *this; }


  /*!
  	\brief convenience function for cpeds_vdrt32
  	\details 
  	@param JDutc - time of interest
  	@return velocity component of the observer towards requested radec direction wrt LSR [km/s]

	It is assumed that the direction is setup in radians (the default standard)
  	\date Jan 2, 2015, 4:20:10 PM
  	\author Bartosz Lew
  */
  double VradialLSR(double JDutc);

  
  
  virtual string print_direction(string comment="", bool show=true) const;

  virtual const DirectionRaDec& operator=(const cpedsDirection& rhs)  { cpedsDirection::operator=(rhs); setEpoch(rhs.val()); return *this; } 
  virtual const DirectionRaDec& operator=(const directionRADEC rhs)  { setRa(rhs.ra);   setDec(rhs.dec);   return *this; } 
  virtual const DirectionRaDec& operator=(const DirectionRaDec& rhs) { cpedsDirection::operator=(rhs); setEpoch(rhs.epoch()); return *this; } 
  
//  virtual const DirectionRaDec& operator*= (const double rhs) { setLon(lon()*rhs); setLat(lat()*rhs); return *this; }
  virtual const DirectionRaDec& operator*= (const double rhs) { cpedsDirection::operator *=(rhs); return *this; }
  const DirectionRaDec operator* (const double rhs) const { return DirectionRaDec(lon(), lat(),epoch(),val())*=rhs; }
//  virtual const DirectionRaDec& operator/= (const double rhs) { setLon(lon()/rhs); setLat(lat()/rhs); return *this; } 
  virtual const DirectionRaDec& operator/= (const double rhs) { cpedsDirection::operator /=(rhs); return *this; }
  const DirectionRaDec operator/ (const double rhs) { return DirectionRaDec(lon(), lat(),epoch(),val())/=rhs; }
  
  
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
//   directionRADEC _raDec;

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */
	 double _epoch; //!< given in JD

 };











/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class DirectionRaDecAh
  \brief Encapsulates directionRaDecAh structure
  \details 
  
  \date 2009/07/01 17:46:20 
  \author Bartosz Lew
*/
class DirectionRaDecAh : public DirectionAh, public DirectionRaDec {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  DirectionRaDecAh() { setA(0.0); seth(0.0); setRa(0.0); setDec(0.0); }
  DirectionRaDecAh(double ra, double dec, double A, double h) { setRa(ra); setDec(dec); setA(A); seth(h); }
  DirectionRaDecAh(directionRADEC n1, directionAh n2) { setRa(n1.ra); setDec(n1.dec); setA(n2.A); seth(n2.h); }
  DirectionRaDecAh(const DirectionRaDec& n1, const DirectionAh& n2) { setRa(n1.ra()); setDec(n1.dec()); setA(n2.A()); seth(n2.h()); }
/*   DirectionRaDecAh(double ra, double dec) { _d.ra=ra; _d.dec=dec; _d.A=0.0; _d.h=0.0; } */
/*   DirectionRaDecAh(double A, double h) { _d.ra=0.0; _d.dec=0.0; _d.A=A; _d.h=h; } */
  virtual ~ DirectionRaDecAh() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  directionRaDecAh getRaDecAh() const { directionRaDecAh d; d.A=A(); d.h=h(); d.ra=ra(); d.dec=dec(); return d; }

  const DirectionRaDecAh& operator=(const DirectionRaDecAh& rhs) { setA(rhs.A()); seth(rhs.h()); setRa(rhs.ra()); setDec(rhs.dec()); return *this; } 

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
/*     directionRaDecAh _d; */

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };


#endif /* Directions */ 
