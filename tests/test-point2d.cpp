#include "cpeds-math.h"
#include "cpeds-point2d.h"
#include "cpeds-project.h"


int main() {

  double *lon,*lat;
  long N=10;

  cpeds_random_uniform_directions_on_sphere(N, &lat, &lon, 0, 5*PI180, 0, 5*PI180);
  cpedsDirectionSet n(N,lon,lat);
  n.print();
  n.save("dirs.txt");

  cpedsProject P(n);
  cpedsPointSet2D p=P.projectOnPlane(cpedsDirection(2.5,2.5)*PI180);
  printf("\n");
  p.print();
  p.save("points.txt");
}
