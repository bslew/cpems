#include <proj_api.h>
#include "cpeds-math.h"
#include "cpeds-consts.h"
#include "cpeds-direction.h"
#include "Mscs-object.h"
#include "Mscs-map.h"

#ifndef MSCS_MAP_PROJ
#define MSCS_MAP_PROJ



/* #ifndef MSCS_MAP */
/* #include "Mscs-map.h" */
/* #endif */


/*!
  \class mscsMapProj
  \brief Encapsulates projections of the mscsMap maps onto a plane and management of the projected flat maps
  \details 
  
  This class accepts an Mscs class map object as a source data and performs requested projections of this map onto a flat  plane
  It can also read/save the flat image from/to disk.
  It also supports visualizations of the projected maps using imlib2 library, and basic interactive operations of this map.
  This class works in two coordinate systems. One is related to the
  matrix coordinates that stores the projected map and each pixel
  has integer coordinates.
  The allowed coordinates values range from 0 to size[XY]-1
  The other coordinate system is connected to the proj library
  in which the projections are performed.
  The allowed coordinates values depend on particular projection.
  The projected plane from proj library can have various sizes
  depending on the projection type.
  Eg. for the Mollweide projection the the matrix 0,0 coordinates
  correspond to the -1,0.5 coordinates in the proj library coodrinate system

  \date 2009/05/28 15:23:15 
  \author Bartosz Lew
*/
class mscsMapProj : public mscsObject {

 public:

  //! Constructor for the projected map class. Not really in use now
  mscsMapProj();
  //! Constructor for the projected map class. Not really in use now
  mscsMapProj(long x, long y);
  //! Constructor for the projected map class. 
  /*! Parameters define the map to be projected, X,Y sizes of he projected map in pixels, */
  /*! projection type, and the longitudal and lattitudal offsets of the central point in the projected map */
  mscsMapProj(const mscsMap& sphT, long X, long Y, string proj, double lonOff, double latOff);
  ~mscsMapProj();
  //! initiates the on-plane coordinates ranges for a given projection
  void initiate_range_variables();


  
  mscsMap& mapSph() { return sphMap; }
  /* projection methods */

  //! Central projections engine. Projects from the spherical map coordinates to flat projected map 
  void project_map();

  void setProjection(string p) { _projection = p; initiate_range_variables(); }
  string projection() { return _projection; }
/* IO handlers */
  void setMap(const mscsMap& sphT) {   sphMap=sphT; }
  void resetZoom() { setZoomWindow(0,0,getSizeX()-1,getSizeY()-1); zoom(1.0); }
  //! Returns the projected map value, given the matrix coordinates of the projected pixel.
  double get_T(long x, long y);
  double zoom() { return _zoomLevel; }
  void zoom(double z) { _zoomLevel=z; }

  /*!
    \brief Returns the direction in galactic coordinates (l,b) in radians extracted from the projected map.
    \details 
    @param x - x coordinate in the matrix
    @param y - y coordinate in the matrix
    @return direction corresdpoding to the pixel in the matrix field
    The correct values will come only from the region where the plot actually 
    exists. From outside of the plot region...well, I don't care. 
    The parametes x,y are the coordinates of a pixel in the map matrix
  
    \date 2010/03/12 16:22:42 
  */
  cpedsDirection getCoordLB(long x,long y);

  void reset() { set_LonLatOffset(0,0); resetZoom(); }
  void setLonLat1(const cpedsDirection& n) { setLon1(n.l()); setLat1(n.b()); }
  void setLonLat2(const cpedsDirection& n) { setLon2(n.l()); setLat2(n.b()); }

  void setLonLat0(const cpedsDirection& n) { setLon0(n.l()); setLat0(n.b()); }

  //! Defines the offset in longitude and latitude to be a central point in the projected map
  void set_LonLatOffset(double lon, double lat);
  //! set central meridian [rad] - changes from -PI, PI
  void setLon0(double l) { lonShift=l; } //-fmod(l, twoPI);  }
  //! set central parallel [rad] - changes from -PI/2, PI/2
  void setLat0(double b) { latShift=b; }//fmod(b, PIsnd); }

  void setLon1(double l) { _lon1=l; }//fmod(b, PIsnd); }
  void setLat1(double b) { _lat1=b; }//fmod(b, PIsnd); }
  void setLon2(double l) { _lon2=l; }//fmod(b, PIsnd); }
  void setLat2(double b) { _lat2=b; }//fmod(b, PIsnd); }

  double lon0() { return lonShift; }
  double lat0() { return latShift; }
  double lon1() { return _lon1; }
  double lat1() { return _lat1; }
  double lon2() { return _lon2; }
  double lat2() { return _lat2; }

  //! Handler returning the matrix size X
  long getSizeX() { return sizeX; }
  //! Handler returning the matrix size Y
  long getSizeY() { return sizeY; }
  //! Sets the pixel value in the flat map matrix at the given matrix coordinates 
  void set_T(long x, long y, double T);
  //! Sets a flag in the flat map.
  /*! The flags are: for either: ok- good projected pixel, bad - pixel without a counterpart in the  */
  /*! spherical coodrinates, masked-masked pixel */
  void set_ctl(long x, long y, long v);

  //! Defines the new window in matrix coordinates for the region to be magnified. 
  /*! The parameters must be given in the matrix coordinates */
  //! The "from" values must be smaller than "to" values
  void setZoomWindow(long fromX,long fromY,long toX,long toY);
  //! Answers whether the pixel in the flat map is masked or not.
  bool isMasked(long x, long y);
  //! Answers true is the pixel is not "bad"
  bool isSet(long x, long y);
  //! Zeros the control map (with flags) and the flat projected map
  void clear();
  matrix<double> getFlatMap(matrix<long>& mask) { mask=*ctl; return *map; }
  double* getFlatArray(matrix<long>& mask, long& n) { mask=*ctl; matrix<double> m=getFlatMap(mask); double* dtab = cpeds_matrix2array(m, n,true); return dtab; }


  //! returns the minimal X coordinate value in the projected plane for a given projection
  double getMinX() { return minX; }
  double getMinY() { return minY; }
  double getMaxX() { return maxX; }
  double getMaxY() { return maxY; }

/*   long save_as(string ext, string name); */

/* flat map operations */

/*   void show_equator(long width, MscsColor color); */
/*   void show_ecliptic(long width, MscsColor color); */
/*   void show_ecliptic_poles(long radius, MscsColor color); */
/*   void show_equator_poles(long radius, MscsColor color); */
/*   void show_meridians(long num, long width, MscsColor color); */
/*   void show_parallels(long num, long width, MscsColor color); */
/*   void show_meridians_labels(long num, MscsColor color, string font, long fontSize); */
/*   void show_parallels_labels(long num, MscsColor color, string font, long fontSize); */

/*   void show_colorbar(long width, long height, long xoff, long yoff, long colornum, long fontSize, string label); */
 

  /*!
    \brief zooms in the map around direction indicated by x,y element of the matrix
    \details 
    @param
    @return
    Ths is done by first calculating the direction corresponding to x,y cell in the matrix,
    setting the zoom box, and reprojecting map centered at that direction.
    The level of zoom in is controllable via zoomLevel.
    It defines by how much the current zoom window is shrunk in each direction.

    \date 2010/03/12 22:33:09 
    \author Bartosz Lew
  */
  void zoomIn(long x,long y, double zoomLevel);

 protected:
  string getProjInitString();


/* visualization */


 private:
  matrix<double> *map;
  matrix<long> *ctl;
  long sizeX,sizeY,pix_num; // size of the plot in pixels
  double maxX,maxY; // maximal values of flat coordinates to be given to the proj program for inverse projection onto sphere
  double minX,minY; // minimal values of flat coordinates to be given to the proj program for inverse projection onto sphere
  double lonShift,latShift; // offset in  longitude and lattitude of the projected map (in deg)
  double _lon1,_lat1, _lon2, _lat2; //!<< define the range of conic projections for the secant plane intersecting with sphere
  string _projection, _projInit;
  bool initiated;
  double acc;
  long zoomFromX, zoomFromY, zoomToX, zoomToY; // defines the state of the current zoom window
  double _zoomLevel;
  mscsMap sphMap;

  
  enum flags { masked=0, ok=1, bad=-1 };

};

#endif
