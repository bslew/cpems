#include "cpeds-point_set.h"


int main() {

  cpedsPointSet3D ps;

  ps.append(cpedsPoint3D(1,2,3));
  ps.append(cpedsPoint3D(4,5,6));
  matrix<double>m=ps.exportAsField(10,10);
  cpeds_matrix_save(m,"field.map");

  ps.print();

  cpedsPoint2D p1(1,2);
  p1.print_point();

  return 0;

}
