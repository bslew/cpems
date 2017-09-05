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

  long int nside,i,j;
  filenamestr filename,cmd;
  FILE *f;

  if (ARGC == 1) { print_program_usage(); exit(1);}
  //if (ARGC != 4)  { print_program_usage(); exit(1);}

  nside = strtol(ARGV[1],NULL,10);
  //strcpy(ordering,ARGV[2]);

  map.set_nside(nside);
  map.makekill_space_manager("make","T",1);

/*   sprintf(cmd,"if [ -f %sn2r-convtab-%li ]; then rm %sn2r-convtab-%li; fi\n",map.data_path,nside,map.data_path,nside); */
/*   system(cmd); */
  for (i=0;i<map.pixNum();i++) {  map.set_T(i,(double)i); }

  map.conv_nest2ring();
//  sprintf(filename,"%s%s%li%s",MSCS_DATA_DIR.c_str(),MSCS_GLOBAL__N2R_CONV_TAB_PREF.c_str(),nside,MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN.c_str());
//  printf("  -- writting to file: %s\n",filename);
//  f = fopen(filename,"w");
//  for (i=0;i<map.pixNum();i++) { 
//    j= (long)(map.get_T(i));
//    fprintf(f,"%li\n",j);
//  }
//  fclose(f);

//  sprintf(cmd,"if [ -f %sr2n-convtab-%li ]; then rm %sr2n-convtab-%li; fi\n",MSCS_DATA_DIR.c_str(),nside,MSCS_DATA_DIR.c_str(),nside);
//  system(cmd);
  for (i=0;i<map.pixNum();i++) { map.set_T(i,(double)i); }
  map.conv_ring2nest(); 
//  sprintf(filename,"%s%s%li%s",MSCS_DATA_DIR.c_str(),MSCS_GLOBAL__R2N_CONV_TAB_PREF.c_str(),nside,MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN.c_str());
//  printf("  -- writting to file: %s\n",filename);
//  f = fopen(filename,"w");
//  for (i=0;i<map.pixNum();i++) { 
//    j= (long)(map.get_T(i));
//    fprintf(f,"%li\n",j);
//  }
//  fclose(f);

}

void print_program_usage() {
/*   printf(" incorrect syntax !!! usage: precalculate_n2rr2n_conv nside ordering [n2r/r2n] \n\n");  */
  printf(" incorrect syntax !!! usage: precalculate_n2rr2n_conv nside \n\n"); 
}
