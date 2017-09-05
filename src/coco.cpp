#include "cpeds-math.h"
#include "cpeds-direction.h"


int main() {

  DirectionRaDec n;

  n.fromGal(cpedsDirection(209,-54)*PI180).print_direction("cold spot ra-dec");

  n=cpedsDirection(cpeds_HMSToAng(15,22,11.47)*PI180,cpeds_DMSToAng(28,54,6.2)*PI180);
  n.print_direction("ra-dec");
  n.toGal().print_direction("lb");

}
