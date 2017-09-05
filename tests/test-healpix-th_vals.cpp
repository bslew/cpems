#include <stdlib.h>
#include <stdio.h>
//#include <cpgplot.h>
#include <math.h>
//#include <fitsio.h>
#include <string.h>
//#include <time.h>
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"
#include "Mscs-function.h"

int main(int argc, char** argv) {
  Mscs_initiate_global_variables();
  // Mscs_print_global_variables();

  long ns=strtol(argv[1],NULL,10);

  mscsMap test("test map");
  test.setVerbosityLevel(Zero);

  test.set_nside(ns);
  test.makekill_space_manager("make","T");
  test.set_map_coord();
  
  mscsFunction f("th vals");
  for (long i=0;i<test.pixNum();i++) {
    f.newPoint(test.get_C(i).b()*PI180inv,0);
  }
  f.unique();
  f.sortFunctionArgAscending();
  f.print();

  // .savebinT("test-map-lmax2-Tn-bin");
}
