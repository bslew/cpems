#ifndef CPEDS_CLASS_DEFS
#define CPEDS_CLASS_DEFS
#include <string.h>
#include "cpeds-math.h"
#include "cpeds-msgs.h"
#include "cpeds-templates.h"

using namespace std;

/*!
  \class cpeds_PDF
  \brief A probability distribution function (PDF) container
  \details 

  \date 2009/05/27 00:30:02 
  \author Bartosz Lew
*/
class cpeds_PDF { 
 public:
  // PDF info variables
  double EX; //!< expectation value
  double MODE; //!< mode value
  double norm;
  bool EXdone,MODEdone,normdone;
  cpeds_queue_double *MODES; //!< list of mode values in the pdf function
  cpeds_queue_long *MODESidx; //!< corresponding list of the arguments in the pdf array (x arguments)

  cpeds_point* pdf; //!< PDF to work on: x - argument, y - PDF
  bool pdf_fields_swapped; //!< true means that x is swapped with y with the switch_pdf_fields() method;
  long N; //!< size of PDF
  double DX; //!< size of step in argumets inferred from two first args 0 and 1

  string object_name;



  // CLs info variables
  double *CL2side_down; //!< CL lower range
  double *CL2side_up; //!< CL upper range
  double *CL; //!< CLs
  long NCL; //!< number of CLs
  cpeds_queue_double * CLRpts; //!< array of cpeds queues, each defining the set of points within a requested CL range

  


  /*!
    \brief 
    \details 
    The input arrays of the pdf function should be resonably resonable eg. sorter in rising abscissa order etc :P
    but itdoesn't need to be normalized I guess
    p1 p2 p3 CLs in rising order to the corresponding returned CL ranges
    
    the confidence limits can be derived in a couple of different ways. Here, I derive the CL ranges as the ranges
    integrating from the ML point (the mode) saidways untill a requested CL is reached.

    @param x,y - list of point coordinates
    @param n - size of x,y arrays
    @param cl - WHAT IS THIS ? please update 
    \date 2009/05/27 00:33:54 
    \author Bartosz Lew
  */
  cpeds_PDF(double *x, double *y, long n, cpeds_queue_double *cl);
  cpeds_PDF(double *x, double *y, long n);
/*   cpeds_PDF(long ncl); */
  ~cpeds_PDF();
  
  void initiate_CLs_info(cpeds_queue_double *cl);
  void initiate_CLs_info(long NCLs);
  void destroy_CLs_info();

  void set_name(string name);
  string get_name();

  void get_info();


  // some handlers

  double get_EX();
  double get_MODE();
  void get_PDF_modes(cpeds_queue_double** PDF_modes, cpeds_queue_long** PDF_modes_idx);
  double get_CL(long i); // returns the i'th CL from internal CL array
  double get_CL_range(bool lower, long i, bool* cont);  // cont flag indicates whether there are gaps in the returned CLrange due to possible multi-modality of the PDF. The continuity is checked by looking whether the points in the CLr are separated by the same, constant separation d, where d=pdf[1].x-pdf[0].x. If the separation is larger then it can indicate a gap in the PDF at requested i'th CL

  void derive_CLs();  // derives the CLranges for requested values of CL given via the constructor
  void derive_CLs(cpeds_queue_double *cl);  // derives the CLranges for requested values of CL given via parameter object

  void integrate_to_CLr(double CLr, bool gaps, long nthmode);

  void save_CLr(string f,long n);

 protected:
  cpedsMsgs *msgs;
  long i,j,k;
  void derive_EX();
  void derive_MODE();
  void find_normalization();

  void switch_pdf_fileds();
  bool CLr_gaps(long j);
  
};

#endif

