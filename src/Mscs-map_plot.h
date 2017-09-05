/*!
  \file implements plotting maps
*/

#ifndef MSCSMAPPLOT
#define MSCSMAPPLOT

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "Mscs-map.h"
#include <proj_api.h>
#include "cpeds-list.h"
#include "cpeds-point_set.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */
class mscsMapProj;

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
  \class mscsMapPlot
  \brief Encapsulates plotting of the spherical maps in various projections
  \details 
  
  \date 2010/03/10 15:33:08 
  \author Bartosz Lew
*/
class mscsMapPlot : public mscsMap {


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
 mscsMapPlot(string name="") : mscsMap(name) {}
 mscsMapPlot(const mscsMapPlot& parent) : mscsMap(parent) {}
 mscsMapPlot(const mscsMap& parent) : mscsMap(parent) {}
 ~mscsMapPlot() {}
/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

 long plot(long sizeX=0, long sizeY=0, string projection="moll", double l=0, double b=0, double Tmin=0, double Tmax=0, long color_num=50, string color_scheme="spectralRGB", string title="Temperature map", string bgcolor="white", double window_size=8, int output=1);

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


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */
 float* replot(mscsMapProj& flat, double& Tmin, double& Tmax, matrix<double>& m, matrix<long>& mask);
 /*!
   \brief return the matrix coordinate of a pixel with world coodrinate v
   \details 
   @param min - minimal world coordinate
   @param max - maximal world coordinate
   @param size - number of pixels in that direction
   @return position of element in matrix (should be integer)

   World coordinates are in the same sense as in the PGplot library - 
   these are coordinates consistent with those that are plotted in the axes in the plot.
   \date 2010/03/12 22:30:07 
 */
 double getMatrixCoord(double v,double min,double max, long size);

 void generateMapColors(string color_scheme,long color_num, double Tmin, double Tmax);

/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };
#endif /* MSCSMAPPLOT */ 

