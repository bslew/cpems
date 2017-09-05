#include "Mscs-common.h"

using namespace std;

class dodecahedron {

 public:
  dodecahedron(double a,double s, map_class * T); // creates a dodec. in "0" reference position
  //dodecahedron(double a,double s, long res, ); // creates a dodec. in "0" reference position
  dodecahedron(double l, double b, double g, double a,double s,  map_class * T); // creates a dodec with requested parameters
  ~dodecahedron();

/*   int Sstat(map_class T); // calculates S statistics on map T with dodec. orientation defined by Di (D) */
/*   int Sstat(double l, map_class T); // calculates S statistics on map T with dodec. in requested orientation defined by Di (D) */
  double * Sstat(double l, double b, double g, double a,double s); // calculates S statistics on map T with dodec. in requested orientation defined by Di (D)
/*   int D2Di(); */

/*   int save_D0(); */
/*   int save_Di0(); */



 private:
  long circ_pix_num,nside,pix_num,ring_num;
  direction Dface0[12]; // definiton of dodecahedron's faces pointing vectors in "0" reference position.
  double Dl,Db,Dg,Da,Ds; // position (l,b), orientation (g) and circle size (a) and shift parameters (s).
  double Dth;

  //  direction *c; // pointer to array of directions to points on a circle
  direction **D0, **D; //
  long **Di;
  map_class *map;

  double S[7]; // array of S statistics on the dodecahedron for a given map (0-5 for individual pairs of circles and 6 - all circles
  double K[7]; // array of K statistics on the dodecahedron for a given map (0-5 for individual pairs of circles and 6 - all circles)

  double *btab; // array of all b values in the map
  double *ltab; // array of all l values in the map pix. by pix.
  long * pirtab; // pix. in ring tab; keeps number of pixels in each ring
  long * ringtab; // array keeps the pix. numers which start a new ring

  void initiate_dodecahedron_structures(string what);
  void define_Dfaces0();
  void derive_D0(double a, double s);
  void derive_D(double l, double b, double g, double a,double s);
  void derive_Di();

  direction Rx(direction p,double Ax);
  direction Rz(direction p,double Az);
  direction RzRx(direction p,double Az, double Ax);
  direction RzRxRz(direction p, double Azz,double Ax, double Az);
  void RzD(direction **D, double Az);
  void RzRxRzD(direction **D, double Azz, double Ax, double Az);
  void copyD(direction **D1, direction **D2);
  //  RotC(long i,
  //help variables
  long i,j;


};

