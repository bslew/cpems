#ifndef __cpedsVerbosityLevel__
#define __cpedsVerbosityLevel__
/*!
  \enum cpeds_messageImportance
  \brief list of available verbosity levels
  \details
  Zero by default does not produce any output, while Top prints out all messages.

  \date 2009/05/27 00:25:05
  \author Bartosz Lew
*/
typedef enum { Zero=0, Low, Medium, High, Top, MustSee } cpeds_messageImportance;
typedef cpeds_messageImportance cpeds_VerbosityLevel;
#define CPEDS_defaultVerbosityLevel High
#endif

#ifndef __cpedsCoordinateSystemType__
#define __cpedsCoordinateSystemType__
/*!
  \enum cpedsCStype
  \brief List of coordinate systems used in CPEDS
  \details

  \date 2009/05/27 00:26:35
  \author Bartosz Lew
*/
typedef enum { cpedsHorizontalCS, cpedsEquatorialFirstCS, cpedsEquatorialSecondCS, cpedsEclipticalCS, cpedsGalacticCS, cpedsPolarAzimuthalCS } cpedsCStype;
#endif

#ifndef filenamestr
typedef char filenamestr[1000];
#endif

#ifndef strarg
typedef const char strarg[1000];
#endif

#ifndef __cpeds_point__
#define __cpeds_point__
/*!
  \struct cpeds_point
  \brief A point container
*/
typedef struct {
    double x;
    double y;
} cpeds_point;


typedef struct {
    long x;
    double y;
} cpeds_pointL;
#endif

#ifndef __cpeds_direction__
#define __cpeds_direction__
/*!
  \struct cpeds_direction
  \brief A direction container
  \details
*/
typedef struct {
    double l; //!< galactic longitude
    double b; //!< galactic lattitude
} cpeds_direction;
#endif

#ifndef __cpeds_pixel__
#define __cpeds_pixel__
/*!
  \struct cpeds_pixel
  \brief A map pixel container.
  \details
  This structure is here defined only for some old compatibility reasons with
  Mscs and should be removed at some point.

  \note This is rather to be removed (depreciated code)
  \date 2009/05/27 00:28:12
  \author Bartosz Lew
*/
typedef struct {
    double *T; // temerature
    double *N; // number of observations
    double *m; // mask
    cpeds_direction *n;
} cpeds_pixel;


/*!
  \enum
  \brief  cpeds status codes list
  \details

  \date 2009/05/29 00:24:00
  \author Bartosz Lew
*/
typedef enum { cpedsSuccess=0, cpedsNoSuchFile, cannotOpenFile, wrongFormat, cpedsCannotWriteToFile, cannotConnect, cpedsError, cpedsNoSuchHDF5dataset } cpedsStatusCodes;


#endif


#ifndef __cpeds_common__
#define __cpeds_common__
#include <vector>
#include <string>

std::vector<std::string> cpeds_strToVec(std::string s,char delim=',');

#endif
