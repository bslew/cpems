#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
//#include "chealpix.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"

void print_program_usage();


int main(int ARGC, char *ARGV[]) {
  Mscs_initiate_global_variables();
  mscsMap map("map");

  long int nside;
  filenamestr filename;

  if (ARGC == 1) { print_program_usage(); exit(1);}
  //if (ARGC != 4)  { print_program_usage(); exit(1);}

  nside = strtol(ARGV[1],NULL,10);
  //strcpy(ordering,ARGV[2]);

  map.set_nside(nside);
  map.makekill_space_manager("make","T",1);
  map.set_map_coord();
  sprintf(filename,"%s%s%li%s",MSCS_DATA_DIR.c_str(),MSCS_GLOBAL__NESTED_COORDINATES_PREF.c_str(),nside,MSCS_GLOBAL__NESTED_COORDINATES_SUFF.c_str());
  map.savebinC(filename);

}

void print_program_usage() {
/*   printf(" incorrect syntax !!! usage: precalculate_n2rr2n_conv nside ordering [n2r/r2n] \n\n");  */
  printf(" incorrect syntax !!! usage: precalculate_coords nside \n\n"); 
}
