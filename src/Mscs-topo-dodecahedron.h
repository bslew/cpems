#include "Mscs-common.h"
#include "Mscs-map.h"

using namespace std;

class dodecahedron {

 public:
  double **S; // array of S statistics on the dodecahedron for a given map (0-5 for individual pairs of circles and 6 - all circles
  //  double K[7]; // array of K statistics on the dodecahedron for a given map (0-5 for individual pairs of circles and 6 - all circles)


  //dodecahedron(double a,double s, map_class * T); // creates a dodec. in "0" reference position
  //dodecahedron(double a,double s, long res, ); // creates a dodec. in "0" reference position
  dodecahedron(double l, double b, double g, double a,double s,  map_class * T, double _amax); // creates a dodec with requested parameters
  ~dodecahedron();

/*   int Sstat(map_class T); // calculates S statistics on map T with dodec. orientation defined by Di (D) */
/*   int Sstat(double l, map_class T); // calculates S statistics on map T with dodec. in requested orientation defined by Di (D) */
  double ** Sstat(double l, double b, double g, double a,double s); // calculates S statistics on map T with dodec. in requested orientation defined by Di (D)
  double ** Sstat_cur();
  void shift_circlesi(double p);
  void save_dodecahedron_circles(string fname);
  void derive_D(double l, double b, double g, double a,double s);
  void derive_D(double l1, double b1, double l2, double b2, double a,double s);
  void set_rho(long rho_loc);
  double get_Dl();
  double get_Db();
  double get_Dg();
  double get_Dg(double l1, double b1, double l2, double b2);
  double get_Da();
  double get_Ds();

/*   int D2Di(); */

/*   int save_D0(); */
/*   int save_Di0(); */



 private:
  long circ_pix_num,nside,pix_num,ring_num, rho, Drho;
  direction Dface0[12]; // definiton of dodecahedron's faces pointing vectors in "0" reference position.
  double Dl,Db,Dg,Da,Ds,Damax; // position (l,b), orientation (g) and circle size (a) and shift parameters (s).
  bool Damax_b;
  double Dth;
  double Clifford_rot,Clifford_rot2;

  //  direction *c; // pointer to array of directions to points on a circle
  direction **D0, **D, **Drotg; // D0 keeps the initial dodecahedron 0;90;0;a;0, D keeps the rotated dodecahedron by g,b,l angles and could be named Drotgbl; Drotg keeps the dodecahedron rotated from initial position by g
  long **Di;
  double *map,*mask,*msq;
  double *w;

  double *btab, *thtab; // array of all b values in the map
  double *ltab; // array of all l values in the map pix. by pix.
  long * pirtab; // pix. in ring tab; keeps number of pixels in each ring
  long * ringtab; // array keeps the pix. numers which start a new ring

  void initiate_dodecahedron_structures(string what,map_class *T);
  void kill_dodecahedron_structures();
  void init_dodecahedron_structures(double a);
  void init_dodecahedron_Dstructures();
  long initiate_rho();
  long initiate_rho(long rho_loc);
  void calculate_circ_pix_num(long rho, double a);
  void define_Dfaces0();
  void match_circles(direction **D, double p);
  void derive_D0(double a, double s);
  void derive_Di();

  direction Rx(direction p,double Ax);
  direction Ry(direction p,double Ay);
  direction Rz(direction p,double Az);
  direction RzRx(direction p,double Az, double Ax);
  direction RzRy(direction p,double Az, double Ay);
  direction RyRz(direction p,double Ay, double Az);
  direction RzRxRz(direction p, double Azz,double Ax, double Az);
  direction RzRyRz(direction p, double Azz,double Ay, double Az);
  void RzD(direction **D, double Az);
  void RyRzD(direction **D, double Ay, double Az);
  void RzRyD(direction **D, double Az, double Ay);
  void RzRxRzD(direction **D, double Azz, double Ax, double Az);
  void RzRyRzD(direction **D, double Azz, double Ay, double Az);
  void copyD(direction **D1, direction **D2);
  //  RotC(long i,
  //help variables
  long i,j;


};

