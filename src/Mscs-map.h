/*!
  \file A definition of the mscsMap class
*/

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


#ifndef MSCS_MAP
#define MSCS_MAP

//#define _GLIBCXX_EXPORT_TEMPLATE

//#ifndef MSCS_MAP_PROJ
//#endif

// standalone headers
//#include <fitsio.h>
#include <fitsio.h>
#include "cpeds-math.h"
#include "cpeds-templates.h"
#include "cpeds-cosmo.h"
#include "cpeds-list.h"
#include "cpeds-direction_set.h"
#include "Mscs-common.h"
#include "Mscs-map-window_function.h"
//#include "Mscs-minkowski-f.h"
/* #include "Mscs-power_spectrum.h" */
#include "Mscs-correlation_function.h"
//#include "Mscs-map_proj.h"

// interdependent headers
/* #include "Mscs-colormap.h" */
#include "Mscs-alms.h"
#include "ccSHT3.h"
#include <vector>



typedef struct {
  double n,m,s,S,K; // which stands for number of points, mean, variance, skewness, kurtosis
} stat_info;


/* **************************************************************************************************** */
/* USING NAMESPACES */
/* **************************************************************************************************** */

using namespace std; // for inside header file implementations

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
  \class mscsMap
  \brief the main class of the Mscs package: encapsulates two-dimentional spherical map
  \details
  Implements all I/O operations, maps handling, basic map operations, map generation,
  gaussian simulations, spherical harmonic transformations, masks management, all necessary handlers...

  \note This is a very big class, and most definitelly should be split into a number of smaller pieces
  eg. almsClass etc

  \note This class is rather badly written - public structures and variables,...
  too big and too clumsy.

  \date 2009/05/27 14:30:08
  \author Bartosz Lew
*/
class mscsMap : public mscsObject {

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:

  /*!
    \struct map_structure
    \brief the main constituent of the map
    \details

    \date 2009/05/28 15:06:28
    \author Bartosz Lew
  */
  typedef enum { ring=0, nested=1, unknownMapOrdering=-1  } mapOrderings;

  typedef struct {
    cpedsList<double> T; //!< temerature array
    cpedsList<double> Q; //!< stokes Q parameter - linear polarization
    cpedsList<double> U; //!< stokes U parameter - linear polarization
    cpedsList<double> V; //!< stokes U parameter - circular polarization
    cpedsList<double> N; //!< number of observations array
    cpedsList<double> m; //!< mask array
    cpedsDirectionSet n; //!< direction array
  } map_structure;

  typedef struct {
    bool T,Q,U,V,N,m,n; // map structures
    bool rotation; //
    long n2r,r2n; //!< conversion arrays between nested and ring orderings
  } mapAllocationFlags;

  typedef struct {
    bool merged;
    long masked_pix_num;
    double f_sky;
    long type;
    long multi_mask_lreg, multi_mask_breg, multi_mask_nside, multi_mask_reg_num;
  } maskInfo;

  typedef struct {
    long int pix_num, rows_num; //map pixel number
    long int coord_num;
    long int nside; //!< map nside
    double meanT,minT,maxT,varianceT, skewnessT, kurtosisT;
    long int iminT, imaxT;
    long n2r_conv_tab_loaded, r2n_conv_tab_loaded; //!< define if the conversion tables for ring and nest were loaded; the flag value indicated the nside for which they were loaded. 0 valus indicates that they weren't loaded.
    mapOrderings ordering; //!< the ordering of the map written in  a given coordinate system (1 - nested 2 - ring )(default is :1) \note this is to be replaced with enum types
    /* int pix_system; //! pixeliztion system of the map: 1 - helapix (default); 2 - tetrahedron; \note this is to be replaced with enum types */
  } mapInfoStructure;
/* ------------- */
/* CLASS FRIENDS */
/* ------------- */

  // ***************************************************************************************************************


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  //! generates an empty map object
  mscsMap(long ns=0);
  //! generates an empty map object with a name
  mscsMap(string _object_name, long ns=0);
  //! a cloning object
  mscsMap(const mscsMap &orig);
  //! a cloning object with a given name of the clone
  mscsMap(const string _object_name, const mscsMap &orig);


  //! Destructor
  virtual ~mscsMap();


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  //! a help function to initiate the member variables
  void set_object_initial_variables();

  /*!
    \brief associates a file name with the map object
    \details
    @param name - preferred file name

    \date 2009/05/28 15:27:04
    \author Bartosz Lew
  */
  void setLoadedMapFileName(string name);
  //! returns the file name associated with this map
  string getLoadedMapFileName() const;
  //! associates a file name with the loaded mask
  void setLoadedMaskFileName(string name);
  //! returns the associated mask file name
  string getLoadedMaskFileName() const;

  //! cloning function
  void clone(const mscsMap &orig);


  //--------------------------------------------------------
  // MAP INFORMATION AND I/O METHODS
  //--------------------------------------------------------



  // FILES METHODS
  //TEMPERATURE
  long savebinT (string map_file); //!< saves temperatures  to a file
  long loadbinT (string map_file, mapOrderings ordering=nested); //!< read in  temperatures from a file
  long savetxtT (string map_file); //!<file ending "-T"
  long loadtxtT (string map_file, mapOrderings ordering=nested);

  /*!
    \brief loads columns from fits file to structures in the map object
    \details
    @param fineName - name of the fits file
    @param colNames - list of column names to be read
    @param dstStructures - list of names of structures into which the corresponding columns from colNames should be loaded into.
    @return status code: 0 - success, negative values would ...

    colNames and dstStructures must be vectors of the same sizes.\n
    colNames must correspond to the column names in the fits file.

    dstStructures can contain one of the following:\n
    "T" - temperature map\n
    "Q" - stokes Q map\n
    "V" - stokes V map\n
    "U" - stokes U map\n
    "m" - mask map\n
    "N" - number of observations map\n
    "n" - coordinates map - NOT IMPLEMENTED YET -- NEED TO BE DONE SEPARATELY\n

    \date 2010/03/12 11:21:23
    \author Bartosz Lew
  */
  long loadfits(const string& fileName, const QList<string>& colNames, const QList<string>& dstStructures);

  /*!
    \brief load a single column of a given name into a single structure
    \details
    details as above
    \date 2010/03/12 11:36:17
  */
  long loadfits(const string& fileName, const string colName, const string dstStructure);
  long loadPLANCK_temp_fits(const string& fileName, string hduName);
  
  long savefits(string fileName);
  void printtxtT(); //!< prints a map to the screen

  long loadbinm (string fileName); //!< read in  mask from a file
  long savebinm (string fileName); //!< read in  mask from a file

  //COORDINATES
  /*!
    \brief saves coordinates from the coordinates structure of the map to file
    \details
    @param coord_file - file name
    @param how - introduced from backward compatibility reasons. Can define extra parameters to the load such as:\n
    "float" - to indicate that the coordinates are stored in 4-byte float numbers rather than default 8-byte numbers.
    The default value is "" - which will read 8-byte number.
    @return
  */
  long savebinC(string coord_file, string how=""); //saves coordinates to a file
  /*!
    \brief loads coordinates of the map from file into the coordinates structure
    \details
    @param coord_file - file name
    @param how - introduced from backward compatibility reasons. Can define extra parameters to the load such as:\n
    "float" - to indicate that the coordinates are stored in 4-byte float numbers rather than default 8-byte numbers.
    The default value is "" - which will read 8-byte number.
    @return
  */
  long loadbinC(string coord_file, string how=""); // read in coordinates from a file
  long savetxtC(string coord_file); // ( file ending "-Cg" )
  long loadtxtC(string coord_file); // how - 1 -- galactic coordinate system: l,b
                                                    //       2 -- flat coordinate system from proj program: x,y
                                                    //       3 -- galactic coordinates system l[0,360] b[+90,-90] - for the proj program
  void printtxtC(string what); // prints coordinates to the screen
                                     // prints the coordinates from the database as requested
 // what = g - 'galactic' coordinates eg. from the healpix map (from the map - data)
      // what = f - flat coordinates from the flatcoord data - this is eg. after making a projection
      // what = g proj - prints the galactic coordinates but in a format acceptable by 'proj' program:
      //                i.e. the lattidutes are calculated from equator at 0 to + and - 90 deg at the poles

  // PIXEL NUMBER OF OBSERVATIONS
  long savebinN(string  map_file);
  long loadbinN(string  map_file);
  long savetxtN(string  map_file);
  long loadtxtN(string  map_file);



  /*!
    \brief the main space allocation manager for the map structure
    \details
    @param whattodo - defines whether a structure should be allocated or deleted (assumes make (or load) or kill)
    @param for_what - defines the target structure:<br>
    "T" - temperature map<br>
    "Q" - stokes Q map<br>
    "U" - stokes U map<br>
    "V" - stokes V map<br>
    "m" - mask<br>
    "N" - number of observations of pixels in the map<br>
    "C" - coordinates for the map
    "rot" - allocates space for a previously saved rotation
    \date 2009/06/12 15:36:51
    \author Bartosz Lew
  */
  /* void  makekill_space_manager(string whattodo, string for_what, mapOrderings ord=nested, long how=0); */
  void  makekill_space_manager(string whattodo, string for_what, long how=0);

           // how parameter is usefull for loading coordinates since there are two
           // dynamical structures connected with them.  -- for flatcoord and the one from the map. so the manager
           // must know which one should be initiated

  /*!
    \brief the main map loader
    \details
    @param
    @return

    \date 2009/06/12 15:51:58
    \author Bartosz Lew
  */
  void loadsave_manager(string whattodo, string fileformat, string what, int how, string where, long whatmultipole);
  //void loadsave_manager(string whattodo, string fileformat, string what, int how, string where, int whatmultipole);
  //void loadsave_manager(char *whattodo, char *fileformat, char *what, int how, string where, int whatmultipole);

  /*!
    \brief loads the pixel trnasfer function from file
    \details
    @param fn - name of the transfer function. If "" given (default) then
    the standard name will be created accoding to naming definitions of the Mscs package
    defined in Mscs-global-defs.
    @param sc - load status code pointer. Will indicate whether the operation succeeded
    @return pixel transfer function - loaded from file, or unitary when loading failed

    \date 2010/02/26 14:00:35
    \author Bartosz Lew
  */
  mscsWindowFunction readPixTf(string fn, cpedsStatusCodes *sc);

  /* cpedsStatusCodes readPixTf(mscsWindowFunction& pixtf, string fn=""); */



  /**************************************************************************************************/
  /* //  MAP HANDLERS - these set or get temperatures from map in an ***ACTIVE*** coordinate system */
  /**************************************************************************************************/

  bool isMasked(long i);
  bool maskAllocated();
  void resetMaskInfo();
  bool maskLoaded() const { return loaded.m; }
  bool mapLoaded() const { return loaded.T; }
  bool coordLoaded() const { return loaded.n; }
  bool isRing();
  bool isNested();

  //! returns the state of the nest to ring ordering conversion array flags
  long n2rConvTabLoaded() const { return loaded.n2r; }
  //! returns the state of the ring to nested ordering conversion array flags
  long r2nConvTabLoaded() const { return loaded.r2n; }
  //! sets the state of the nested to ring ordering conversion array flags
  void n2rConvTabLoaded(long n) { loaded.n2r=n; }
  //! sets the state of the ring to nested ordering conversion array flags
  void r2nConvTabLoaded(long n) { loaded.r2n=n; }

  //! returns the nside value
  long int nside() const { return mapInfo.nside; } //!< returns the actual nside of the map
  long int pixNum() const { return mapInfo.pix_num; } //!< returns the actual pixels number in the map
  long maskedPixNum() const { return mask.masked_pix_num; }
  long multi_mask_reg_num() const { return mask.multi_mask_reg_num; }
  long int coordNum() const { return mapInfo.coord_num; } //!< returns the actual pixels number in the map
  mapOrderings ordering() const { return mapInfo.ordering; }

  /* long get_rows_num(); */
  long ringNum() const { return cpeds_get_ring_num_healpix(nside()); }

  /*!
    \brief  Sets the map nside pixel number and coordinates number variables, but it doesn't allocate the map space
    \details
    @param n - nside to be set

    This routine does not allocate memory space for maps since it doesn't know which maps you want to process.
    (temperature, mask, Nobs, etc). To make an empty map  (to allocate space for a given map) use
    make_kill_space_manager member function.

    \date 2010/02/24 10:32:11
    \author Bartosz Lew
  */
  void set_nside(long n);

  /*!
    \brief sets map pixel number variable and updates corresponding nside
    \details
    @param n - number of pixels to be set
    If this donesn't cleanly indicate some integer nside value then warning message is issued

    \date 2010/03/11 09:37:39
    \author Bartosz Lew
  */
  void set_pixNum(long n);

  /*!
    \brief sets ordering flag
    \details
    @param ord - ordering

    This function only sets the ordering flag but doesn't really perform any ordering conversions so
    one must use this function wisely- eg in order to read ring ordered map as nested etc.
  */
  void set_ordering(mapOrderings ord);

  // the direction is passed in the system consistent with THIS class  - i.e. in l,b
  const cpedsList<double>& get_T() const { return map.T; }
  const cpedsList<double>& get_Q() const { return map.Q; }
  const cpedsList<double>& get_U() const { return map.U; }
  const cpedsList<double>& get_V() const { return map.V; }
  const cpedsList<double>& get_N() const { return map.N; }
  const cpedsList<double>& get_m() const { return map.m; }
  const cpedsDirectionSet& get_n() const { return map.n; }

  double get_T(long i) const { return map.T.value(i); } // returns the temperature i'th pixel
  void   set_T(long i, double t) { map.T[i]=t; } // sets the temeperature of i'th pixe
  void set_Treg(long reg, double t);
  void set_Treg(double l, double b, double t);
  //! sets temperature of pixels with temperature T to val
  void set_Ttoval(double t, double val);
  void set_Treg_if_gr_abs(long reg, double t);
  void set_Treg_if_le_abs(long reg, double t);
  double get_T(const cpedsDirection& nn) const; //!< returns the temperature from requested direction
  double get_T(const cpedsDirection& nn, long *i) const; //!< returns the temperature from requested direction and the corresponding pixel number
  double get_T(double l, double b) const; //!< returns the temperature from requested direction
  void set_T(double l, double b, double val); //!< sets the temperature from requested direction
  void   set_T(const cpedsDirection& nn, double t); //!< sets the temeperature at requested direction
  cpedsDirection get_TminmaxCoord(bool max);

  /*!
    \brief get the modifiable reference to the pixel with number num
    \details
    @param num - pixel number
    @return a modifiable reference to the temperature in the pixel
  */
  double& T(long num) { return map.T[num]; }
  double& Q(long num) { return map.Q[num]; }
  double& U(long num) { return map.U[num]; }
  double& V(long num) { return map.V[num]; }
  double& N(long num) { return map.N[num]; }
  double& m(long num) { return map.m[num]; }
  cpedsDirection& n(long num) { return map.n[num]; }

  /*!
    \brief get the modifiable reference to the temperature object
    \details
    @return a modifiable reference to the temperature object

    Note that if you use this to assign new values to this vector, then remember that the assignment operator does
    not copy any statistical information on the temperature map (i.e mean variance ,etc). To be up to date, you need to
    either copy that info youself or run calculate_map_stats()
  */
  cpedsList<double>& T() { return map.T; } // sets the temeperature of i'th pixe
  //! same as for T() but for Q()
  cpedsList<double>& Q() { return map.Q; } // sets the temeperature of i'th pixe
  //! same as for T() but for U()
  cpedsList<double>& U() { return map.U; } // sets the temeperature of i'th pixe
  //! same as for T() but for V()
  cpedsList<double>& V() { return map.V; } // sets the temeperature of i'th pixe
  //! same as for T() but for N()
  cpedsList<double>& N() { return map.N; } // sets the temeperature of i'th pixe
  //! same as for T() but for m()
  cpedsList<double>& m() { return map.m; } // sets the temeperature of i'th pixe
  /*!
    \brief same as for T() but for n()
    \details
    Note that the naming here is somewhat inconsistent since usually letter C is used to deal with coordinates not n.
    This is possibly to be changed in future.
  */
  cpedsDirectionSet& n() { return map.n; }

  double get_Q(long num) const { return map.Q.value(num); } //!< returns the Q stoks parameter at i'th pixel
  double get_U(long num) const { return map.U.value(num); } //!< returns the Q stoks parameter at i'th pixel
  double get_V(long num) const { return map.V.value(num); } //!< returns the Q stoks parameter at i'th pixel

  //! returns a shortened list of points in from the map from the requested structure defined by "what"
  cpedsList<double> get_nonMasked(string what);
  cpedsDirectionSet get_nonMaskedDirs();
  cpedsDirectionSet get_maskedDirs();


  void set_m(long num,double m) { map.m[num]=m; } //!< sets mask at pixel num to value m
  double get_m(double l, double b) const; //!< returns the mask value at coordinates (l,b) [rad]
  double get_m(long i) const { return map.m.value(i); } //!< returns the mask value in the pixel with number num

  cpedsDirection get_C(long int num) const { return map.n.value(num); } //!< returns the coordinates of i'th pixel
  void get_C(long int num, double * l, double * b)  const; //!< returns the coordinates of i'th pixel
  void set_C(long int num, const cpedsDirection& nn) { map.n[num]=nn; } //!< sets the coordinates of i'th  pixel
  long  get_Ci(const cpedsDirection&  nn) const; //!< returns the coordinates of i'th pixel
  cpedsDirection  get_C(const cpedsDirection&  nn) const; //!< returns the coordinates of i'th pixel

  double get_N(long i) const { return map.N.value(i); } //!< returns the number of i'th pixel observations
  void set_N(long i, double n) { map.N[i]=n; } //!< sets the number of i'th pixel observations to n
  void clear_map() { T()=0.0; } //!< clears whole temperature map

  /*!
    \brief populates the map structure with a value x
    \details
    @param - x value to be inserted into the map
    @param - what:  defines which part of the map_pixel structure should be populated
    T - stands for temperature map
    m - stands for mask
    N - stands for pixel numbe observations

    \date 2009/05/28 22:50:47
    \author Bartosz Lew
  */
  void setValue(double x,string what);
  /* void setValue(const cpedsDirection& n);  */ // this routine doesn't make much sense

  /* void clear_map(map_structure * maptoclear); //!< clears whole temperature map */


  /*!
    \brief kills the map structure and all it's components
    \details
    @param pointer to the map structure to be deleted. Formally it only makes sense to use own maps to delete.

    This routime frees allocated memory int the map structures, but does NOT change the flags status.<br>
    This routine is to be used by some map transformations; it is only for special applications, be careful
    of how you use it.

    \date 2009/06/12 15:05:12
    \author Bartosz Lew
  */
  void kill_map_space(map_structure * map); //!<
  /*!
    \brief removes all maps and resets the corresponding variables
  */
  void killAll(); // formerly known as:  kill_map(); //!< clears the whole map_structure structure
  void kill_window_functions();
  /*!
    \brief Allocates space for map structure and all it's components.
    \details
    @param
    @return

    \note THIS IS NOT IMPLEMENTED. CHECK THE CODE
    \date 2009/06/12 15:08:07
    \author Bartosz Lew
  */
  void make_map_space(map_structure * map);

  /* /\*! */
  /*   \brief Clones the map structure space (only) and returns the cloned map structure. */
  /*   \details  */
  /*   @param map - pointer to the source map to be cloned */
  /*   @param pix_num - number of pixels to be allocated in the clone structure */
  /*   @return allocated clone map_structure structure  */

  /*   Clones the map structure space in the new area of memory with a given size */
  /*   by pix_num parameter. This routine allowes to allocate space for the same structures */
  /*   but in different resolution; it does NOT copy the data */

  /*   \date 2009/06/12 15:10:01  */
  /*   \author Bartosz Lew */
  /* *\/ */
  /* map_structure clone_map_space(map_structure * map,long pix_num); */
  /* void copy_map(map_structure *map_from, map_structure *map_to, cpedsList<long> *conv_tab); */
  void copy_map(mscsMap *map_from, mscsMap *map_to, cpedsList<long> *conv_tab);
  void copy_map(cpedsList<double> *from, cpedsList<double> *to, cpedsList<long> *conv_tab);
  void import_map_data(mscsMap &from_here, string what, int how);
  void import_map_data_to(mscsMap &from_here, string what, string where, int how);
  /* mscsPowerSpectrum* export_C_l(mscsPowerSpectrum * Cl); */
  int load_conv_tabs(long nside, string type); // loades the ordering conversion arrays from file for a given nside
  void set_conv_tabs(long nside); // sets the ordering conversion arrays for a given nside
  void kill_conv_tabs();




  /***********************/
  /* // PLOTTING METHODS */
  /***********************/

  //void plot_Tmap(); // plots the pgplot map of simulated or physical map
  //void plot_Tcircle(cpeds_direction n, double r); //plots the tempetature on a circle in cpeds_direction n of angular radius r

  //  long plot_map(int projection, double dl, double db, int output,  long color_num, string color_scheme,string title,double window_size); // plots the flat map in a requested projection






  /***************************/
  /* MAP STATISTICS HANDLERS */
  /***************************/
  double get_meanT() const { return mapInfo.meanT; }
  double get_varianceT() const { return mapInfo.varianceT; }
  double get_skewnessT() const { return mapInfo.skewnessT; }
  double get_kurtosisT() const { return mapInfo.kurtosisT; }




  //-------------------------------------------------------------
  // MAP DRAWING/GENERATION  METHODS
  //-------------------------------------------------------------
  void make_multi_maskLB(long l_reg_num, long b_reg_num, string mask_info_filename,  double Ath, double Aphi,double Al, bool rot_short);
  void make_multi_maskHP(long nside_loc, string mask_info_filename="", double Ath=0, double Aphi=0,double Al=0, bool rot_short=false);
  void mask_region(long reg);
  int make_circle_dot(double l, double b, double r, double v, string what, string region_type, long region_type_dot_points, string operation);
  void mk_ring(double l, double b, double s, double v, string what);
  void draw_circles(const matrix<double>& lbrv, string dstMap, string overplot_region_type, long overplot_region_type_dot_points, string operation);
  void make_equatorial_mask(double b);

  void markEcliptic(double val, double width);
  void markEquator(double val, double width);
  void markEclipticPoles(double val, double size);
  void markNSPoles(double val, double size);
  void markGalacticMeridians(double val, long N, double width);
  void markGalacticParallels(double val, long N, double width);
  void markEquatorialMeridians(double val, long N, double width);
  void markEquatorialParallels(double val, long N, double width);

  /*!
	\brief create dipole component in the full sky with axis along direction m 
	\details 
	@param m - dipole axis direction in l,b [deg]
	@return
	
	assumes that the map's ns is already set.

	\date Mar 19, 2013, 11:51:07 AM
	\author Bartosz Lew
   */
  void makeDipole(cpedsDirection m);



  //-------------------------------------------------------------
  // MAP OPERATION  METHODS
  //-------------------------------------------------------------
  //! adds value to the map without changing values under the mask
  void TaddVal(double val);
  //! subtracts value from the map without changing values under the mask
  void TsubVal(double val);
  //! multiplies  the  map by the value without changing values under the mask
  void TmulVal(double val);
  //! divides  the  map by the value without changing values under the mask
  bool TdivVal(double val);
  void cutAbove(double val);
  void cutBelow(double val);
  void replaceIfGreater(double val, double val2);
  void replaceIfLess(double val, double val2);
  void replaceIfAbsGr(double val, double val2);
  void replaceIfAbsLe(double val, double val2);
/*   void xorVal(long val); */
/*   void orVal(long val); */
  void TaddMap(const mscsMap& map);
  void TsubMap(const mscsMap& map);
  void TmulMap(const mscsMap& map);
  void TdivMap(const mscsMap& map);

  //! adds t to i'th pixel
  void Tadd(long i, double t);
  //! subtracts t from i'th pixel
  void Tsub(long i, double t);
  //! multiplies  i'th pixel by t
  void Tmul(long i, double t);
  //! divides i'th pixel by t
  void Tdiv(long i, double t);


  void change_map_resolution(long int new_nside);
  void invert_mask();
  cpeds_queue<long int>* get_circle(double l, double b, double r,double step, string region_type, long region_type_dot_points);
  //void make_multi_maskLB(long l_reg_num, long b_reg_num, string mask_info_filename);
  void clear_multimask();
  void mask_map_merge();
  void multimask_map_merge();
  void mask2map();
  void check_mask();
  void clean_mask();
  void scramble_over_nsigma(double n);
  void scramble_over_minmax(double min,double max);
  void save_reference_orientation();



  /*!
    \brief performs a gaussian smoothing of the temperature map using scalar spherical harmonic transformation
    \details
    @param beamFWHM - FWHM of the beam to be used for smoothing the map [rad]
    If debeamFWHM>0 then the beam kernel will be formed and SHT coefficients will be multiplied by this kernel before SH synthesis. [rad]
    @param debeamFWHM - FWHM of the beam to be used for desmoothing the map - i.e. for the first SHT. Default -1 which means that no desmoothing is done.
    If debeamFWHM>0 then the beam kernel will be formed and SHT coefficients will be divided by this kernel after SH analysis
    @param lmax - maximal multipole number up to which SHT is done. Default: -1 which indicates that the SHT is done up to lmax=2*nside()
    @param pixtf - if true then the pixel transfer function will be used.
    @param method - indicates the SHT method to be used.

    The FWHM is given in radians.
    \date 2010/02/23 14:25:24
    \author Bartosz Lew
  */

  void gaussianSmoothT(double beamFWHM, double debeamFWHM, long lmax=-1, bool pixtf=false, long method=0);

  /*!
	\brief interpolate on sphere using tangent plane projections and radial basis functions interpolation scheme
	\details 
	@param toNs - target Healpix resolution parameter
	@param rbfName - radial basis function to be used for interpolation. 
	Valid names are: \n
	"gaussian" - phi(r) = exp(-0.5*(r/r0)^2)
	
	@param r0 - interpolation scale parameter [deg]
	@param dirsMask - a pointer to an mscsMap object of resolution toNs with mask structure initialized to indicate which directions should have
	values interpolated (0 - don't interpolate, 1 - interpolate); If not supplied (or null given) then all directions will be interpolated.
	@param R - maximal range (angular separation) [deg] within which neighboring points are taken for interpolation process. The selection process currently costs N^2 operations where N is the number 
	of directions in the data to be interpolated. Larger r1 will select more points for interpolation (resulting possibly in a better interpolation) but it will cost M^3 operations (where
	M is the number of selected directions) to perform interpolation using these M points.
	@return returns a new mscsMap object with the interpolated field

	The data directions and function values are taken from this object's temperature vector. All directions are considered unless they are masked.


	\date Nov 7, 2012, 3:13:48 PM
	\author Bartosz Lew
   */
  void interpolate(mscsMap& outMap, double R, double r0, string rbfName="gaussian");

  /*!
    \brief White noise, gaussian map generator
    \details
    @param m - mean of the gaussian distribution
    @param s - variance of the gaussian distribution
    @param method - cpeds method to be used for the RN generation<br>
    method: 1 - the invers CDF from cpeds is used for generagion of random nubers<br>
    method: 2 - gsl library is used for generagion of random nubers
    method: 3 - central limit theorem with 100 numbers
    @param rns - random numbers object pointer; if null given (default) than a new one will be allocated and destroyed when it is done.\n
    if non-null is given then it will be used and will not be destroyed and can be re-used
    for further calls. In this case after all is done you need to delete it from outside of this method.
    This is useful for very frequent calls to this method (faster than 1Hz) to avoid generation of random numbers from
    the same seed.


    \note FOR THE MOMENT THE method ARGUMENT IS IGNORED AND CENTRAL LIMIT THEOREM IS USED FOR RN GENERATION
    \date 2009/06/05 15:05:24
    \author Bartosz Lew
  */
  void make_gaussian_map(double m=0, double s=1, int method=1, cpedsRNG *rns=NULL);

  /*!
    \brief White noise, map generator (unifirm distribution)
    \details
    @param min - minimal value of the noise
    @param max - maximal value of the noise

    THE rest of the details as in generate_gaussian_vector

    \date 2009/06/05 15:05:24
    \author Bartosz Lew
  */
  void generate_uniform_vector(double min, double max, int method, cpedsRNG *rns);
  /*!
    \brief returns an argument sorted function with map temperatures and corresponding pixel numbers as function values
    \details
    @param direction defines how the temperatures should be sorted: 12 - ascending, 21 - descending
    @return mscsFunction
    Quick sort is used for sorting

    \date 2010/02/16 09:42:20
    \author Bartosz Lew
  */
  mscsFunction get_value_sorted_pixel_numbers(long direction);
                                              // whatisgaussian: "D" - generates gaussian distributed moduls and uniform random phis of alms
                                              //                 "RI" - generates gaussian real and unreal parts of alms
// whatkind: 1 - dont use the information about the power spectrum in C_l table
// whatkind: 2 - DO use that information and generate the gaussian alms but the assumption is that the C_l are not normalized by factor l(l+1)/2pi
// whatkind: 3 - DO use that information and generate the gaussian alms but the assumption is that the C_l are normalized by factor l(l+1)/2pi
// method: 1 -  the invers CDF from cpeds is used for generagion of random nubers (iterative)
// method: 2 -  the invers CDF from cpeds is used for generagion of random nubers (by half division)
// method: 3 -  gsl library is used   for generagion of random nubers
// method: 4 -  the invers CDF from cpeds is used for generagion of random nubers (by half division) but make just one call to the cpeds function for all alms_num at once.
// in this method the variance throughout all multipoles is the same. nevertheless the underlying power spectrum is respected

  //void gaussian_smoothing(); // performs the gaussian smoothing on a map
  void average_map_in_rings(cpeds_queue<double>** qp);
  /*   void set_cosmic_variance(); */

  /*!
	\brief returns i'th ring from the map as a mscsFunction
	\details
	@param i - number of the map ring to return (starts from 0)
	@return mscsFunction containing the longitudes as arguments and temperature values as function values

	THe map must be in ring ordering and the coordinates must be loaded.

	\date May 19, 2010, 10:59:56 AM
	\author Bartosz Lew
  */
  mscsFunction getRing(long i);
  /*!
    \brief calculates the correlation function of the loaded map
    \details
    @param theta_min - the minimal angle in the correlation function [deg]
    @param theta_max - the maximal angle in the correlation function [deg]
    @param resolution - step in resolving the correlation function [deg]
    @return the correlation function pointer is returned and also C_th protected variable set; old contents will be deleted
    The correlation function arguments are stored in radians.

    \date 2009/06/03 16:51:29
    \author Bartosz Lew
  */
  mscsCorrelationFunction calculate_C_th(double theta_min, double theta_max, double resolution); // calculates the correlation function on a map // if resolution is -1 then it's taken from values set by read_binC_th_parameters (or txt of course too :)

  /*!
	\brief calculate S correlation statistic on the map
	\details 
    @param theta_min - the minimal angle in the correlation function [deg]
    @param theta_max - the maximal angle in the correlation function [deg]
    @param resolution - step in resolving the correlation function [deg]
    @return calculates S correlation statistic
    
    S = 2 * <Ti mi Tj mj> / ( <Ti^2 mi mj> + <Tj^2 mi mj> )

	\date Jun 11, 2018, 3:28:45 PM
  */
  mscsCorrelationFunction calculate_Sth(double theta_min, double theta_max, double resolution); // calculates the correlation function on a map // if resolution is -1 then it's taken from values set by read_binC_th_parameters (or txt of course too :)

  /* /\*! */
  /*   \brief extracts the power spectrum from the loaded map using MASTER method */
  /*   \details  */
  /*   @param */
  /*   @return */

  /*   \date 2009/06/03 16:55:28  */
  /*   \author Bartosz Lew */
  /* *\/ */
  /* mscsPowerSpectrum * extract_C_l(const mscsMap long lmax_loc, mscsPowerSpectrum * Blsq, mscsPowerSpectrum * pseudoN, mscsPowerSpectrum * pseudoW, matrix <double>* M, matrix <double>* C, mscsPowerSpectrum * Nfs, string Mllinfo); */
  double calculate_circ_statistics(cpeds_direction C, double r); // calculates the circle correlation function around a given dircection -- dedicated to TSZ
  //void calculate_3P_corr_function();
  //void calculate_bispectrum();
/*   void initiate_dodecahedron(); */
/*   double * calculate_dodecahedron_circle_stistic(double l, double b, double g, double a, double s); */

  /* void test_gaussianity_chisq(string what); */
  double calculate_meanT();
  double calculate_varianceT();
  double calculate_rmsT();
  double calculate_minT();
  double calculate_maxT();
  double get_integralT();
  double getMinT() const { return mapInfo.minT; }
  double getMaxT() const { return mapInfo.maxT; }
  int rotate_map(double Ath, double Aphi,double Al, bool rot_short, string what);
  int prepare_rotation(double Ath, double Aphi,double Al);
  int remove_rotation();
  int rotate_quick(string what, bool unrotate=false);
  void calculate_map_stats(int output=0);
  stat_info * calculate_map_stats_with_multi_mask(int output);
  stat_info * calculate_map_variance_with_multi_mask(int output=0);
  long * count_multi_mask_pix_num();
  double calculate_map_momentum(); // calculates the "angular momentum" of the map in real space: i.e. returns M = sum_p T^2*cos(|b|)^2


  /*!
    \brief find maximal/minimal map momentum
    \details
    @param nsst - the nside from which the initial set of directions is taken
    based on the healpix pixelization (for b>0)
    @param nsen -- the nside at which the nested indepth
    search is supposed to end the maximization, thus limiting the accuracy.
    @param acc - parameter returns the accuracy of the result in radians
    @param Mmax - pointer to allocated double where max/minimal momentum is returned
    @param maximize - assumes either true or false for  maximalization or minimalization respectively
    @return direction of the lmax bmax in gal coords. in rad. of the axis orientation

    routine finds the orientation of the axis that maximizes/minimizes the angular momentum of the map
    using calculate_map_momentum() method, and nested healpix indepth search.
    the input parameters are:

    \date 2010/03/12 15:01:02
    \author Bartosz Lew
  */
  cpedsDirection maximize_map_momentum(long nsst, long nsen, double *Mmax, double *acc, bool maximize);

  /*!
    \brief fits the dipole in the map
    \details
    @param range - defines the maximal(minimal) amplitude of the dipole
    @param range - defines the maximal step in the dipole amplitude search
    @param AMP - is filled with the amplitude of the fitted dipole
    @param remove_dipole - if true the fitted dipole is removed from the map
    @return returns the direction of the fitted dipole

    \note THIS ROUTINE SHOULD BE RETESTED BEFORE USING
    \date 2010/02/23 15:20:09
    \author Bartosz Lew
  */
  cpedsDirection fit_map_dipole(double range, double acc, double *AMP, bool remove_dipole=false);
  double calculate_map_to_map_chisq(mscsMap* map);













  /*!
    \brief calculate the spherical haromonic transformation: derives a map scalar from alms
    \details
    @param method- defines which SH method is to be used
    how = 1 - ccSHTx
    how = 2 - spherepack ( need to include -DINCLUDE_SPHEREPACK at compilation time )
    how = 7,8,9 - old and slow routines, if you want to use these you better check the code for all the requirements

    @param a - reference to alms object with alms to be used for the map generation. nside of the map must be preset; if NULL is given here, then by default the alm pointer will be tried
    @param lmax - maximal multipole number upto which to perform synthesis
    @param beamtf - is a beam transfer function file name to be used in the SH space; the alms will be multiplied by this function before synthesis
    @param pixtfype - pixel transfer function type; if "" is given then no pixel transfer function is used; Currently also "healpix" pix tf is supported
    @param smoothGauss - is a gaussian smoothing function kernel full width at half maximum in given in degrees; if =0 than no smoothing is performed

    @return a synthesized map

    A map is generated in this object. If the map already exists then it will be lost.

    The alms must be antysymetrical for this to work.

    Remember to set nside to the requested value before calling this function


    the smoothing is possible while composing the map since the gaussian smoothing is a convolution of a function with a smoothing kernel and the convolution of two
    functions is a product of their fourier transforms. Since we create a map from it's fourier coefficients it's natural to do the smoothing in the fourier space.
    there are two ways of smoothing the map:
    1) the gaussian smoothing with the gaussian kernel : thFWHM parameter: if negative no smoothing is done
    2) with a beam window function of an experiment: used to be gaussian but can be replaced with WMAP DA depentant window functions:
    smoothing = "" - no smoothing is done
    smoothing = "filename" - smoothing according to wmap window functions read from file filename with a window function

    \date 2009/06/04 14:18:24
    \author Bartosz Lew
  */
  mscsMap& SH_synthesis(const mscsAlms& a, long lmax, const mscsWindowFunction& beamtf, string pixtfType, double smoothGauss=0, long method=1);
  mscsMap& SH_synthesis(const mscsAlms& a, long lmax, const mscsWindowFunction& beamtf, const mscsWindowFunction&  pixtf, double smoothGauss=0, long method=1);

  /*!
    \brief shorter version of the routine above
    \details

    The beamtf and pixtf is assumed to be unitary.
  */
  mscsMap& SH_synthesis(const mscsAlms& a, long lmax) {   mscsWindowFunction bl("beam tf",lmax);  bl.make_unit_kernel();  return SH_synthesis(a,lmax,bl,""); }





  /*!
    \brief calculates the spherical harmonic transformation using the current map (returns alms)
    \details
    @param how - defines which SH method is to be used
    how = 1 - ccSHTx
    how = 2 - spherepack ( need to include -DINCLUDE_SPHEREPACK at compilation time )
    @param a - if this is NULL then a pointer to the alms object is returned, and the same pointer is assigned to the protected variable alm.
    If the alm is not NULL before the transformation then a warning message is issued. The contendts of the object that
    lives under that pointer will be lost ! A warning message will be issued though. If alm==NULL then a new object
    will be generated to contain the map's alms.
    @param lmax - the maximal number of multipole upto which to calculate the alms
    @param pixtfFile - pixel transfer function file name; if "" is given then no pixel transfer function is used
    @param smoothing - is a beam transfer function file name to be used in the SH space; the alms will be divided by this function
    @param thFWHM - is a gaussian smoothing function kernel full width at half maximum in given in degrees; if <0 than no smoothing is performed
    @param save_partial: 1 - save the partial files with temperature maps according to names given by partial_files_prefix
    save_partial = 0 - do not save partial files
    @param partial_files_prefix - prefix for saving the partial files for maps with increasing number of multipoles included and individual multipoles maps: THIS IS DEPRECIATED, YOU CAN IGNORE THIS AND GIVE ""
    @param plms_file - a file with Legendre polynomials: THIS IS DEPRECIATED, YOU CAN IGNORE THIS and put ""
    @param show_process - setting graphical visualization of the progress - THIS IS DEPRECIATED AND DOESN'T WORK WITH THE ccSHT library; you can ignere this and use ""

    @return - a pointer to an mscsAlms object with the corresponding alms
    pix_system - 1 - HEALPIX pixelization system - the transform is performed and then the map is stored in a given pixelization system
    pix_system - 2 - TETRAHEDRON pixelization system (not implemented yet)

    Most of the parameters are not very much important and are kept from backward compatibility reasons but
    in fact only one method for SH calculations at the moment is in use (ccSHT) because it's the fastest method implemented so far.

    \date 2009/06/04 13:53:14
    \author Bartosz Lew
  */
  const mscsAlms SH_analysis(long lmax, const mscsWindowFunction& beamtf, string pixtfType, double desmoothGauss=0, long method=1);
  const mscsAlms SH_analysis(long lmax, const mscsWindowFunction& beamtf, const mscsWindowFunction& pixtf, double desmoothGauss=0, long method=1);

  /*!
    \brief shorter version of the routine above
    \details

    The beamtf and pixtf is assumed to be unitary.
  */
  const mscsAlms SH_analysis(long lmax) {   mscsWindowFunction bl("beam tf",lmax,getVerbosityLevel());  bl.make_unit_kernel();  return SH_analysis(lmax,bl,bl); }

  /*!
    \brief create an inverse noise weighed co-added map
    \details
    @param daNum - total number of DAs to be co-added
    @param sigma0 - a pointer to an array of size daNum with standard deviations per pixel
    @param da - a pointer to an array of maps to be co-added (They may have Nobs initiated)
    @param weight - a string that defines what INC should be done:<br>
    "Nobs" - uses the number of pixel observations in the da maps to create the INC<br>
    "simple" - assumes that the number of pixel observations is always 1 for all pixels in the map


    \date 2009/06/05 16:15:23
    \author Bartosz Lew
  */
  //void mk_INC_DA(long daNum, double* sigma0, mscsMap* da, string weight);
/*   void mk_INC_DA(long DAst,long DAen,mscsMap* DAs,string weight); */  //commented out during transition into version-1.0


  void shift_mean_to(double val, bool calc_stats=true );


/*   void calculate_minkowski_area(double sigmaDA, long thres_num, double th_min, double th_max); */
/*   void calculate_minkowski_circ(double sigmaDA, long thres_num, double th_min, double th_max); */
/*   void calculate_minkowski_genus(double sigmaDA, long thres_num, double th_min, double th_max); */
  //void calculate_minkowski_v0(long thres_num, double nsigma);
  mscsFunction calculate_minkowski_v0(long thres_num, double min, double max);
  mscsFunction calculate_minkowski_v1(long thres_num, double min, double max);
  mscsFunction calculate_minkowski_v2(long thres_num, double min, double max);
/*   minkowski_fs* calculate_minkowski_v0v1v2(long thres_num, double min, double max); */
  vector<mscsFunction> calculate_minkowski_v0v1v2(long thres_num, double min, double max);

/*   void test_minkowski_circ(); */


  /* void conv_ring2nest(map_structure * map_loc=NULL); */
  void conv_ring2nest();
  void conv_nest2ring();



  // fills the map coordinates within the actuall coordinates pixelizations scheme and sets the coor_loaded status
  // it's done with regard to the current settings of nside and pix_system
  void set_map_coord(double db=0.0,double dl=0.0);
  void norm_by_maxT();
  void norm_by_mean();
  void norm_by_variance(); // divides the map by it's variance
  void norm_by_stddev(); // divides the map by it's standard deviation
  void norm_by_rms(); // divides the map by it's variance
  /* bool isHealpix(); */
  void mk_abs();
  void mk_sqrtAbs();
  void invertT();
  void logarithmT(double base);
  void powerT(double exponent);





  /*********************************************/
  /* // MINKOWSKI FUNCTIONAL RELATED FUNCTIONS */
  /*********************************************/
  /*!
	\brief count pixels in temperature map above threshold mu
	\details 
	@param mu - threshold
	@param nsig - not used
	@param sigma0 - not used
	@return pixel count

	\date Jun 7, 2018, 7:34:14 PM
   */
  long Ngrmu(double mu, double nsig, double sigma0);
  long Nlemu(double mu);
  void Ngrmu(double mu, long* counts);
  void Ngrmu(double *t, long thres_num, long** counts);
  void Ninthresmu(double *mapl, double *t, long thres_num, double** counts);


  /**************************************************/
  /* // METHODS FOR DIFFERENTIAL CALCUS ON THE MAPS */
  /**************************************************/

  mscsMap cov_derivative_phi();
  mscsMap cov_derivative_th();
  mscsMap partial_derivative_th();
  mscsMap partial_derivative_phi();



  /*************/
  /* OPERATORS */
  /*************/
  const mscsMap& operator=(const mscsMap& rhs) {   
	  if (this != &rhs) { 
		  map=rhs.map; 
		  mask=rhs.mask; 
		  mapInfo=rhs.mapInfo; 
		  loaded=rhs.loaded; 
		  r2n_conv=rhs.r2n_conv;
		  n2r_conv=rhs.n2r_conv;
		  setLoadedMaskFileName(rhs.getLoadedMaskFileName());
		  setLoadedMapFileName(rhs.getLoadedMapFileName());

		  this->mscsObject::operator=(rhs); 
	  }
	  return *this; 
  }

 protected:
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */

/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  void read_binmap_parameters(string map_file); //!< reads the map information from a binary map file
  void read_txtmap_parameters(string map_file); //!< reads the map information from a text map file
  void read_bincoord_parameters(string coord_file); //!< reads the coordinates information from a binary coordinates file
  void read_txtcoord_parameters(string coord_file); //!< reads the coordinates information from a text coordinates file



  // METHODS FOR PIXEL ROTATIONS
  cpedsDirection RzRy(const cpedsDirection& p,double Az, double Ay);
  cpedsDirection Rx(const cpedsDirection& p,double Ax);
  cpedsDirection Rx(const cpedsDirection& p,double sinAx, double cosAx);
  cpedsDirection Ry(const cpedsDirection& p,double Ay);
  cpedsDirection Ry(const cpedsDirection& p,double sinAy, double cosAy);
  cpedsDirection Rz(const cpedsDirection& p,double Az);
  cpedsDirection Rz(const cpedsDirection& p,double sinAz, double cosAz);

  int rotate_map_directions(cpedsDirectionSet& X, int axis, double ang);

  // optimized math routines for SKregstat project
  void get_moments_of_distribution(const cpedsList<double>& tab, double* mean, double * variance, double * skewness, double * kurtosis);

/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
/* CLASS VARIABLES DECLARATION */
  string loadedMapFileName;
  string loadedMaskFileName;


  // MAP VARIABLES
  map_structure map, *map_ptr;//, *mask; // a mapw
  cpedsList<long> n2r_conv, r2n_conv; //!< ring2nest, nest2ring conversion tables pointers
  cpedsList<double> rotation_map;

  // MAP RELATED VARIABLES AND FLAGS
  maskInfo mask;
  mapInfoStructure mapInfo;
  mapAllocationFlags loaded;
  /* long int masked_pix_num, multi_mask_lreg, multi_mask_breg, multi_mask_nside, multi_mask_reg_num; */
  /* double f_sky; */
  /* int mask_type; */

  //cpedsDirection reference_dir; //!< a direction of the first pixel in the map, kept to save the original information about the orientation of the map if some rotations are performed
  int default_fourier_method;
  /* long seed_offset; //!< offset to be applied to the seed taked from the current time (no of secs from 1980) in case of parallel runs \note - this is to be removed once RNG is placed in objects */







 private:
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */

/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */
  //! read map nside from fits file
  long get_fits_map_nside(fitsfile **fptr);
  //! read map ordering from fits file
  mapOrderings get_fits_map_ordering(fitsfile **fptr);
  /*!
    \brief read a column by index
    \details
    @param colno - indicates the column number to read
    @param fitsfile - must be a valid opern fits file pointer pointer
    @param result - is set to zero if operation succeded, and to non-zero if it failed
    @return
  */
  cpedsList<double> get_fits_map_data_column(fitsfile **fptr, long colno, int* result);
  /*!
    \brief read a column by name
    \details
    @param colName - indicates the column name to be read
    @param fitsfile - must be a valid opern fits file pointer pointer
    @param result - is set to zero if operation succeded, and to non-zero if it failed
    @return
  */
  cpedsList<double> get_fits_map_data_column(fitsfile **fptr, string colName, int* result);
  //! list fits file headers
  void list_fits_map_header(string fileName);


  void precalculate_fourier_stuff(long8 ring_num, const mscsMap& mapRING, flt8 *ThVals,long8 *ThBreaks, flt8* phi0, flt8 * dphi);

//  void precalculate_fourier_stuff(long ring_num, const mscsMap& mapRING, double *ThVals,long *ThBreaks, double* phi0, double * dphi);


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */
};


#endif

/* #define NSIDE 1024 //these are the definitions of the maximal possible values */
/* #define HEALPIX_PIXNUM 12*NSIDE*NSIDE */
/* #define MAX_L 2048 // maximal multipole number */
/* #define MAX_M 2*MAX_L+1 */
/* #define MAX_C_TH_POINTS 10000  // maximal number of points that tabulate  correlation function */
/* #define MINKOWSKI_LEVEL_NUMER 200 */
/* #define NUMBER_OF_WINDOW_FUNCTIONS 10 */

/* typedef struct { // moved to Mscs-common.h since 2009-02-03 */
/*     double R; // real part */
/*     double I; // unreal part */
/*     double D; // module */
/*     double P; // phase */
/* } a_lm; */
  //typedef a_lm struct znumber;
//typedef a_lm a_lms[MAX_L][MAX_M];



/* typedef struct { */
/*   int ok; */
/*   double A,s,m; */
/*   double Astart,sstart,mstart; */
/*   double dA,ds,dm; */
/*   double xmin,xmax; */
/*   long int num; */
/*   double chisq; */
/*   double mean,skew,kurt; */
/* } fit_params; */

/* typedef struct { */
/*   fit_params circ,area,genus; */
/* } fit_info; */



  /* //! generates an empty map object and passes the info on various directories */
  /* /\*! \note This is depreciated. It is not needed to be done that way since this information is global *\/ */
  /* mscsMap(string _object_name, string data_path_str); */

  /* //! generates an empty map object and passes the info on various directories */
  /* /\*! \note This is depreciated. It is not needed to be done that way since this information is global *\/ */
  /* mscsMap(string _object_name, package_dirs dirs); */
/*   mscsMap(string data_path_str); */

/*   //! gives a name to the object */
/*   void set_object_name(string name); */
/*   //! returns the object's name */
/*   string get_object_name(); */

/*   void read_binC_th_parameters(string C_th_file); */
/*   void read_txtC_th_parameters(string C_th_file); */
/*   void read_binC_l_parameters(string C_l_file); */
/*   void read_txtC_l_parameters(string C_l_file, int colnum); */
/*   void read_txtM_parameters(string M_file); */
/*   void read_binM_parameters(string M_file); */
/*   void read_txtbl_parameters(string _blfile); */ //commented out during transition into version-1.0








/*   double * get_map_data_addr(string towhat); */
  //pixel get_Pi(long int num);
  //pixel get_PC(direction n);



  //void loadtxtCflat(string  coord_file); // reads the coordinates from an output file from the "proj" program and puts it into a *flatmap viariable
/*   void killflatcoord(); // kills the themparature map and coordinates */
  /* void killflatmap(); //!< kills the flat, projected map */


  /* int plot_temperature_histogram(string filename); */

/*   int plot_minkowski_area(int output);    // plots the minkowski functional area */
/*   int plot_minkowski_circ(int output);    // plots the minkowski functional circumference */
/*   int plot_minkowski_genus();    // plots the minkowski functional circumference */
/*   int plot_minkowski_v(int whichone, int output); */



  /* void plot_mm_reg_number(long reg); */
  // projection: 1 - moll - Mollweide
  // xsize - size of the map in pixels in x cpeds_direction

  //void load_map_colors(string color_scheme, int color_num);  //commented out during transition into version-1.0
  //void show_progress(string what, string filename);  //commented out during transition into version-1.0

  // CONVERTING METHIDS
  //void convert_lb2xy(); // calls externam mapping program
  //void convert_xy2lb();



  /* bool almsLoaded();   */


  /* void gaussianSmooth(double beam_size, string plms_file); */




/*   double calculate_minkowski_area_min(); */
/*   double calculate_minkowski_area_max(); */
/*   double calculate_minkowski_circ_min(); */
/*   double calculate_minkowski_circ_max(); */
/*   double calculate_minkowski_genus_min(); */
/*   double calculate_minkowski_genus_max(); */
/*   double calculate_minC_l(); */
/*   double calculate_maxC_l(); */
/*   double calculate_minC_th(); */
/*   double calculate_maxC_th(); */

/*   void calculate_flat_coord(int projection_type, double dl, double db); // convert the corrdinates data into some XY projection for plotting purposes */

  //void calculate_minkowski_genus();
/*   void fit_function_3parameters_by_MCMC(int f, long int num, double *data, double Xmin, double Xmax, double ndelta, long int chainnum2, double* params); */


/*   void set_seed_offset(long offset); // commented out on 2009-01-22 -- moved to cpeds and individual applications  */
/*   long get_seed_offset(); // commented out on 2009-01-22 -- moved to cpeds and individual applications  */


/*   double* GSL_random_uniform_numbers(double min, double max, long int num);   */   //commented out during transition into version-1.0
/*   double* GSL_random_gauss_numbers(double fm, double fs, long int num);   */  //commented out during transition into version-1.0
  //  void calculate_flat_map_parameters(int printinfo); // calculates the extremal values of the map coordinates
/*   float* create_resized_flat_map(int printinfo); // fits the converted to desires screen size */
  // returns the pointer to the dynamical table of size x pixels and y corresponding to the map y/x ratio
  // also sets the map size variables xpix_num_flat,ypix_num_flat1



  //commented out during transition into version-1.0
/*   int  makekill_space_manager(string whattodo, string fileformat, string for_what, string what_file, int how); */





// ***************************************************************************************************************
/* CLASS DATA CONTROL FLAGS */


  /* int map_loaded,coord_loaded; */
  /* int mask_loaded, mask_merged; */
  /* int Nobs_loaded; */
  /* int rotation_loaded; */

  //int alms_loaded; //////// ***************************** temporary //commented out during transition into version-1.0

  /* int flat_coord_loaded, flat_map_loaded; */ //commented out during transition into version-1.0
  /* int flat_coord_ordering; // the ordering of the map written a given coordinate system (1 - nested 2 - ring )(default is :1) */ //commented out during transition into version-1.0

// ***************************************************************************************************************
  // FLAT MAP RELATED VARIABLES
/*   directionxy *flatcoord; //!< projected coordinates \note this is to be removed */
/*   long int flat_coord_num; //!< number of coordinates (pixels) in this pixelisation */
  /* mscsMapProj* flatmap; */ //commented out during transition into version-1.0
/*   float *flatmap; //!< this is the map used for plotting purposes: \note this should be changed */
/*   double minx_flat,maxx_flat,miny_flat,maxy_flat,x2yratio_flat; */
/*   int xpix_num_flat,ypix_num_flat; */
/*   long int pix_num_flat; //!< this is the rectangular number of points in the flat map (including the points from outside of the sphere projection region */
/*   filenamestr flatmap_file_name; */
/*   double flatmap_minT,flatmap_maxT; */
/*   bool flatmap_change_zrange; */

// ***************************************************************************************************************
  //SH SPACE RELATED VARIABLES
  /* mscsAlms *alm;//,aSIM,aTH; */ //commented out during transition into version-1.0
  /* long lmax; */
/*   long int alms_num; // number (amount) of alms together, not the element number starting from 0. ////////\***************************** temporary */

// ***************************************************************************************************************
  //! window functions
  //mscsWindowFunction * bl_tab[NUMBER_OF_WINDOW_FUNCTIONS]; //static array of pointers to objects keeping all information about window functions in all DAs of the WMAP experiment   //commented out during transition into version-1.0
  // first 8 are assumed to represent the WMAP window functions. the remaining 2 are free for other use.
  /* mscsWindowFunction * bl; */
  // long bl_lnum;   //commented out during transition into version-1.0
/*   window_function *wf1, *wf2; */
// ***************************************************************************************************************
// topology
//  dodecahedron * dodec;
// ***************************************************************************************************************
  // spectra
  /* mscsCorrelationFunction* C_th; //!< correlation function pointer */
  //mscsPowerSpectrum * C_l; //!< power spectrum pointer //commented out during transition into version-1.0
/*   double C_l[MAX_L][3]; // power spectrum of the map; 0 - l , 1 - C_l, 2 - cosmic variance stuff  */
/*   double C_th[MAX_C_TH_POINTS][2]; // correlation function of the map */

/*   long lmax_C_l,lmin_C_l,l_num_C_l; // actuall: lmax and lmin in the fourier data, - these are the limits for which the C_l will be stored calculated or loaded -  */
/*                                    //not the multipoles of the minimal and maximal value of the C_l. */
/* /\*   long int alm_max_num; // number (amount) of alms together, not the element number starting from 0. *\/ */
/*   double thetamin_C_th, thetamax_C_th, dtheta_C_th; //  minimal and maxiamal theta angle in corralation function and the angle step */
/*   double C_l_min,C_l_max, C_th_min, C_th_max; */
/*   int point_num_C_th, C_l_mini, C_l_maxi, C_th_maxi, C_th_mini; // number of points that tabulate the correlation function */

// ***************************************************************************************************************
  // minkowski functionals
/*   int minkowski_level_num_circ, minkowski_level_num_area,minkowski_level_num_genus, minkowski_level_num; */
/*   double minkowski_area[MINKOWSKI_LEVEL_NUMER][6]; // 0 - nu, 1 - C(nu) 2,3 - lower/upper x error 4,5 - lower/upper y error */
/*   double minkowski_circ[MINKOWSKI_LEVEL_NUMER][6]; */
/*   double minkowski_genus[MINKOWSKI_LEVEL_NUMER][6]; */
/*   double area_min, area_max, circ_min, circ_max, genus_min, genus_max; */
/*   fit_info mink_fit_info; */
/*   long int area_mini, area_maxi, circ_mini, circ_maxi, genus_mini, genus_maxi; */
  //! minkowski functionals pointers
  //minkowski_f *mink_v0,*mink_v1,*mink_v2; // this is probably a bit too much to have it all in one object but whatever //commented out during transition into version-1.0

// ***************************************************************************************************************
  //! WMAP maps related information
  //removed during transition into version-1.0
// ***************************************************************************************************************
  // HELP VARIABLES
  /* double Klm_prev; // Klm  =  (l-m)!/(l+m)! */
  /* long l_prev,m_prev, ok_prev,rek_depthl; */
  /* filenamestr data_path; *///commented out during transition into version-1.0
/*   double * exclude_values_list; */

// experimentall
  /* bool cff; */
// experimentall

// ***************************************************************************************************************
  // VISUALIZATION VARIABLES
/*   int vis_sizeofxy; */
/*   int vis_iter; */
/*   int plot_line_color; */
/* /\*   int minkowski_plot_from // flag valriable indicating whether the plotting data come from minkowski arraych of from the objects - this is temporary I think !!! *\/ */
/* /\*     // by default 0 - plot from minkowski_area etc, if 1 - plot using  *\/ */
/*   float *vis_xtab, *vis_ytab; */
/*   float vis_x,vis_y; */
/*   float vis_x_prev,vis_y_prev; */
  /* bool first_time; */ //commented out during transition into version-1.0



  /* a_lm a(long l,long m); */

  /* FILE* f; */

