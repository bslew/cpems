#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"

void print_program_usage();


int main(int ARGC, char *ARGV[]) {
  map_class ilc,mask;
  filenamestr infile,mask_file,mask_prefix;
  long i;
  long th_num = 15;
  double nsigma;

  Mscs_initiate_global_variables();

  if (ARGC == 1) { print_program_usage(); exit(1);}
  strcpy(infile,ARGV[1]);
  strcpy(mask_file,ARGV[2]);
  th_num = strtol(ARGV[3],NULL,10);

  Mscs_initiate_global_variables();
//--------------------------------------------------------------------------------------------------------------
//  CALCULATE THE MINKOWSKI FUNCTIONALS FROM THE SAMPLE OF GAUSSIAN RANDOM SIMULATIONS
//--------------------------------------------------------------------------------------------------------------  
  // load mask
  if (strcmp(mask_file,"nomask") != 0) { sprintf(mask_prefix,"%s%s",program_dir,mask_file);      mask.loadfitsT(mask_prefix,2);   mask.change_map_resolution(256); }
  // load file
  ilc.loadbinT(infile,1);  ilc.change_map_resolution(256);    

  if (strcmp(mask_file,"nomask") != 0) { 
    for (i=0;i<mask.pix_num;i++) { ilc.map[i].m = mask.map[i].m; } // copy mask
    ilc.mask_loaded = 1; // !!!!!!!!!!!!!!!!!!!!!! THIS IS TERRIBLE WHAT I DO HERE !!!!!!!!!!!!!!!!!!! THIS IS NOT THE OBJECTIVE ORIENTED PROGRAMMING !!!!!!!!!!!!!!!
    ilc.check_mask();
    ilc.mask_map_merge();
  }
  ilc.calculate_map_stats(1); nsigma=4.0;//}
  ilc.calculate_minkowski_circ(-1,th_num,-nsigma*sqrt(ilc.varianceT),nsigma*sqrt(ilc.varianceT)); ilc.plot_minkowski_circ(3);
  ilc.calculate_minkowski_area(-1,th_num,-nsigma*sqrt(ilc.varianceT),nsigma*sqrt(ilc.varianceT)); ilc.plot_minkowski_area(3);
/*   make_detailed_file_name(&outfile,WMAPcalc_dir,256,512,"wmap","",0,"",0,"kp0","smooth",8); */
  ilc.savetxtM(infile,1);  // save area    
  ilc.savetxtM(infile,2);  // save circ.

}

void print_program_usage() {
  printf(" incorrect syntax !!! usage: calculate_minkowski_functionals map_file (prefix) mask_file (prefix) thres_num\n\n"); 
}
