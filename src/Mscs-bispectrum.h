#include "Mscs-map.h"

#ifndef MSCS_BISPECTRUM
#define MSCS_BISPECTRUM

using namespace std;

//!This class derives the angle averaged bispectrum:
/*! This is it calculates: 

 \f$B_{l_1 l_2 l_3} = \sum_{m_1 m_2 m_3}\Biggl(  
 \begin{array}{ccc}l_1& l_2& l_3\\
 m_1& m_2& m_3\end{array}\Biggr) B_{l_1 m_1, l_2 m_2, l_3 m_3}\f$ 

 where the \f$B_{l_1 m_1, l_2 m_2, l_3 m_3}\f$ is the 
full bispectrum defined as :

 \f$B_{l_1 m_1, l_2 m_2, l_3 m_3}= \langle a_{l_1 m_1}a_{l_2 m_2}a_{l_3 m_3}\rangle\f$ */

/*! This class includes parallel implementation of the bispectrum calculation. */

template <class T>
class Mscs_bispectrum {
  // CLASS PUBLIC METHODS
 public:

  typedef struct { 
    T R; // real part
    T I; // unreal part
  } a_lm_Bispectrum;


  //! constructor
  Mscs_bispectrum(string name, long L1min,long L1max,long L2min,long L2max,long L3min,long L3max, map_class* src_alms);
  ~Mscs_bispectrum();

  //! class handler
  T get(long l1,long l2,long l3);
  //! class handler
  void set(long l1,long l2,long l3,T v);

/*   T operator (long l1, long l2, long l3) { return get(l1,l2,l3); } */
/*   Mscs_bispectrum & operator = (T v) { set(l1,l2,l3 */
  

  void zeroBispectrum();
  void deriveBispectrum();

 protected:
  // CLASS VARIABLES AND STRUCTURES
  string objectName;
  long l1min,l1max, l2min,l2max, l3min,l3max;
  long sizel1, sizel2,sizel3;
  map_class *alm;
  a_lm_Bispectrum ***b;

/*   T ***bR; */
/*   T ***bI; */

  // CLASS PROTECTED METHODS
  bool rangesOK();
  void allocateBispectrum(bool allocate);

  bool isEven(long l1,long l2,long l3);


 private:

  



};

#endif
