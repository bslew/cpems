/*!
  \file implements the projections of the given set of directions onto a tangent plane
*/

#ifndef CPEDSPROJECT
#define CPEDSPROJECT

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "cpeds-math.h"
#include "cpeds-direction.h"
#include "cpeds-point_set.h"
#include "cpeds-direction_set.h"
//#include <stl/deque.h>

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsProject
  \brief Encapsulates projections of a set of directions onto a tangent plane
  \details 

  The input and output directions are assumed to be given in radians in a coordinate system
  consistent with (l,b) - galactic CS - with l in (0,2PI) b in (-PI,PI)
  
  
  \date 2009/11/05 23:16:51 
  \author Bartosz Lew
*/
class cpedsProject {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:

	typedef struct {
		bool invertLongitudeOnProject;
	} cpeds_project_t;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  /* cpedsProject(QList<cpedsDirection> &dirs, string projection); */
  /* cpedsProject(QList<cpedsPoint2D> &points, string projection); */
  /*!
    \brief class constructors
    \details 
    @param dirs - set of directions [rad]
    @return
    Use the names of the projections as in proj package
    
    \date 2010/03/10 18:14:54 
  */
  cpedsProject(const cpedsDirectionSet &dirs, string projection="stere");
  cpedsProject(const cpedsPointSet3D &points, string projection="stere");

  ~cpedsProject();


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  void setInvertLon(bool tf) { _cpeds_project.invertLongitudeOnProject=tf; }
  string getProjection() const { return _projection; }
  /*!
    \brief returns the direction of the tangent plane onto/from which projection takes place
  */
  cpedsDirection& projectionDirection() { return _plane; }
  const cpedsDirectionSet getDirections() const { return cpedsDirectionSet(_ds); }
  cpedsDirectionSet& directions() { return _ds; }
  const cpedsPointSet3D getPoints() const { return cpedsPointSet3D(_ps); }
  cpedsPointSet3D& points() { return _ps; }
  long pointsCount() { return _ps.count(); }
  long directionsCount() { return _ds.count(); }
  /*!
	\brief calculates a length scale in for the current projection corresponding to given angle
	\details 
	@param ang - angular scale to be converted to linear scale [rad]
	@return linear scale corresponding to ang on the projection plane

	\date Jan 4, 2013, 3:16:43 PM
	\author Bartosz Lew
   */
  double getLengthScaleFromAng(double ang);

  /*!
    \brief projects the set of directions onto a plane
    \details 
    @param n - direction of the tangent plane [deg] where longitude changes within [0,360) and latitude within [-90, 90]
    @return - returns a set of projected points

    The plane is assumed to be tangent to the sphere in point towards direction n
    The object must be initialized with the set of directions of interest

	The z coordinate of each of the output point will be set to the corresponding input direction's value.

    \date 2009/11/06 14:45:47 
    \author Bartosz Lew
  */
  cpedsPointSet3D projectOnPlane(const cpedsDirection& n, string projFormat="");

//  mscsFunction3dregc projectOnRectangularGrid(const cpedsDirection& n, string projFormat="");

  /*!
    \brief projects the set of points onto a sphere
    \details 
    @param n - direction of the tangent plane

    The plane is assumed to be tangent to the sphere in point towards direction n
    @return - returns a set of projected points

    The plane is assumed to be tangent to the sphere in point towards direction n
    The object must be initialized with the set of points of interest
    
    \date 2009/11/06 14:45:47 
    \author Bartosz Lew
  */
  cpedsDirectionSet projectOnSphere(const cpedsDirection& n, string projFormat="",bool usePointVals=false);

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
  cpeds_project_t _cpeds_project;

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

  cpedsDirectionSet _ds; //!< the directions to be projected - lon,lat [rad]
  cpedsPointSet3D _ps; //!< the projected points on the tangent plane

  cpedsDirection _plane; // tangent plane direction [rad] - lon, lat
  string _projection;


};
#endif /* CPEDSPROJECT */ 

