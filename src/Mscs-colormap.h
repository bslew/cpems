#include <string.h>
//#include <stdlib.h>

#ifndef MSCS_COLORMAP
#define MSCS_COLORMAP

#include "Mscs-common.h"
#include "cpeds-msgs.h"


//using namespace std;


/*!
  \class mscsColormap
  \brief generaes colormaps for 3-D plots.
  \details 
  This class handles the generation and assignment of colors to map values
  for 3-D plots of projected maps.
  Some predefined color schemes can be used.
  
  \date 2009/05/27 11:55:35 
  \author Bartosz Lew
*/
class mscsColormap {

 public:

  /*!
    \brief Generates color pallete with num colors in it, given the name of the color
    scheme in "colors" argument.
    \details 
    @param colors - a requested color scheme; can assume one of the following values:
    spectralBGR, spectralRGB, spectral, color, grayscale
    @param num - total valsCount of colors to be used
    @param minV - the minimal value to be plotted
    @param maxV - the maximal value to be plotted
    
    The colors are linearly assigned and distributed beteween the minV and
    maxV values
    
    \date 2009/05/27 11:58:33 
    \author Bartosz Lew
  */
  mscsColormap(string colors="spectralRGB", long num=50, double minV=0, double maxV=1);
  ~mscsColormap();
  
  //! Zeros thresholds and color tables.
  void clear_variables();

  long colorNum() const { return valsCount; }
  //! Generates thresholds for the colors assignments
  void mk_thres();
  //! Generates the color scheme
  void generateColors();

  //! Returns an array of thresholds as double
  double* export_thres();
  //! Returns an array of thresholds as float
  float* export_thresf();
  //! Returns a normalized array of thresholds as float - from 0 for lowest threshold  to 1 - for highest threshold
  float* export_thresf_norm();


  float* export_rf();
  float* export_gf();
  float* export_bf();

  //! Returns an MscsColor given a value
  MscsColor getC(double v);

  //! Sets the gamma for a colorscheme to g
  /*! \notice currently not used */
  void setGamma(double g);

  //! Defines the standard deviaiton for the color schemes that base on the RGB Gaussian shaped fileters
  /*! The red filter peaks on the lowes value, Blue on the highest value */
  /*! and green on the value in the middle */
  void setSigma(double s);

  //! Defines the standard deviaiton for the color schemes that base on the RGB Gaussian shaped fileters
  /*! Each sigma value is controlled seperatelly */
  /*! The red filter peaks on the lowes value, Blue on the highest value */
  /*! and green on the value in the middle */
  void setSigma(double sr, double sg, double sb);

  //! Defines the standard deviaiton for the color schemes that base on the RGB Gaussian shaped fileters
  /*! The red filter peaks on the lowes value, Blue on the highest value */
  /*! and green on the value in the middle */
  void setSigma(MscsColorD c);

  //! Defines the predefined RGB sigma values
  MscsColorD getSigmaRGB0();

  //! Returns the current sigma value, common for all three colors.
  double getSigma0();

  //! Prints information on the defined colors to the standard output.
  void print_colors();

  //! Sets the current color scheme for colors generation
  void setColorScheme(string c);

  //! returns the next color scheme name given the current colors name
  static long getColorSchemesNumber() { return 5; }
  static string cycleColorSchemes(string colors);
  static string getColorSchemeName(long id);
  static long getSchemeID(string colors);

 protected:
  cpedsMsgs *msgs; //!< standard CPEDS user communication and logging object
  double min,max; //!< minimal/maximal projected value
  long valsCount; //!< total number of colors in the color scheme
  string colorScheme; //!< generated color scheme

  double *thres; //!< thresholds array of size valsCount 
  long *r; //!< array with intensities of the red color for each threshold value
  long *g; //!< array with intensities of the green color for each threshold value
  long *b; //!< array with intensities of the blue color for each threshold value

  //! standard deviations that help defining the gaussian color filters
  double gamma,sigma,sigma0,sigmaR,sigmaG,sigmaB,sigmaR0,sigmaG0,sigmaB0;
  

  //  enum colorScheme { rainbow, grayscale }

};

#endif
