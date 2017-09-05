#include <stdlib.h>
#include <stdio.h>
//#include <cpgplot.h>
#include <math.h>
//#include <iostream>
//#include <fstream>
#include <string>
//#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"

using namespace std;

int main(int argc, char **argv) {
  long i,num;
  string infile, cmdstring;
  double l,b,ang,a1,a2,neg;
  double refl1,refb1,refl2,refb2;
  long n,type;
  FILE *f=NULL;

  cmdstring="statSKregstat SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100 -d SKregstat--wmap3_ILC___x___mm-LB_64_16-1_100 --stats --with_data --C  --alpha 0.1 --force_mm_num 77 --force_map_num 9000 --mm_list_file list-mm-LB_64_16-1_100 --chisq_file SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100___CL90.chisq -o SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100___CL90 --Cov --fast --calculate_MC_MOD_PDF --MC_MOD_PDF_sim_num 1000 --MC_MOD_PDF_st_sim 9001   --data_dir /home/blew/external300/SKregstats/ --mm_start 77"; // for removing -r list

  cmdstring="statSKregstat SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100 -d SKregstat--wmap3_ILC___x___mm-LB_64_16-1_100 --stats --with_data --C  --alpha 0.1 --force_mm_num 77 --force_map_num 9000 --mm_list_file list-mm-LB_64_16-1_100 --chisq_file SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100___CL90.chisq -o SKregstat--sim-wmap3_ILC-1_10000___x___mm-LB_64_16-1_100___CL90 --Cov --fast --calculate_MC_MOD_PDF --MC_MOD_PDF_sim_num 1000 --MC_MOD_PDF_st_sim 9001   --data_dir /home/blew/external300/SKregstats/ --mm_start 77 --negate_r"; // for removing all but -r list

  if (argc == 1) { printf("usage: calc_ang_dist "
		  "[infile eg mmID77info] "
		  "[number of directions 1 - for l1,bl1, and l2,b2 not used, "
		  "2 - for l1,b1,l2,b2 , "
		  "3 - for l1,b1 with cut in b only] "
		  "l1 b1 [l2 b2] "
		  "ang "
		  "negate[1/0]\n\n\n "
		  "infile - file with region number l b; \n"
		  "if the filename STDIN is given then the l1 b1 l2 b2 coordinates from standard input are taken; \n"
		  "the negate option also applies; \n"
		  "all other command line parameters then are not used \n\n\n "
		  "negate - in case of 2nd mode (without infile) this option defines:\n "
		  "0 - return the cos(n1,n2),\n "
		  "1 - return the angle n1n2 in deg,\n "
		  "2 - angle n1 n2 in rad,\n "
		  "3 - sin (n1 n2),\n "
		  "4 - cos(n1,n2) mod 180  ie. the cos(n1,n2) for ang(n1,n2) <= 90 and cos(PI-ang(n1,n2)) for ang(n1,n2) > 90 (for odd_symmetries project)\n "
		  "5 - sin(n1,n2) mod 180  ie. the sin(n1,n2) for ang(n1,n2) <= 90 and sin(PI-ang(n1,n2)) for ang(n1,n2) > 90 (for odd_symmetries project)\n "
		  "6 - ang(n1,n2) mod 180  ie. the ang(n1,n2) for ang(n1,n2) <= 90 and PI-ang(n1,n2) for ang(n1,n2) > 90 (for odd_symmetries project) \n\n\n ang - unused in case of 2nd more; it is to be used with 1st mode only and with a special use for SKregstat project only\n\n The data passed to STDIN are assumed to be in deg in lon lat format"); exit(0); }
  infile=argv[1];

  
/*   refl1=226; refb1=-32; ang = 10; */
/*   refl2=226; refb2=-32;  */

  refl1 = strtod(argv[3],NULL); refb1 = strtod(argv[4],NULL);
  refl2 = strtod(argv[5],NULL); refb2 = strtod(argv[6],NULL);
  ang =  strtod(argv[7],NULL);
  neg =  strtod(argv[8],NULL);


  refl1 *= PI180;  refb1 *= PI180;   refb1 = PIsnd-refb1; ang *=PI180;
  refl2 *= PI180;  refb2 *= PI180;   refb2 = PIsnd-refb2;

/*   a2=cpeds_ang_n1n2(refb1,refl1, refb2,refl2); */
/*   printf("%lf \n",180/PI*a2); */
/*   exit(0); */

  if (infile!="STDIN") {

    f = fopen(infile.c_str(),"r"); 
    if (f !=NULL) {
      printf("reading infine: %s\n", infile.c_str());
      i=0;
      type=strtol(argv[2],NULL,10); printf("type: %li\n",type);
      printf("%s ",cmdstring.c_str());
      
      while (fscanf(f,"%li %lE %lE",&n,&l,&b) != EOF) {
	if (type == 1 || type == 2) a1=cpeds_ang_n1n2(refb1,refl1, PI180*(90-b), PI180*l);
	if (type == 2) a2=cpeds_ang_n1n2(refb2,refl2, PI180*(90-b), PI180*l);
	//printf("%li %lE %lE %lf %lf \n",n,l,b,180/PI*a1,180/PI*a2);
	
	if (neg==0) {      
	  if (type == 1 || type == 2) { if ( a1 <= ang) { printf(" -r %li ",n); } } // for n1 direction
	  if (type == 2 && a2 <= ang) { printf(" -r %li ",n); }  // for n2 direction
	  if (type == 3 && fabs(b) < 180/PI*(PIsnd-refb1)) { printf(" -r %li ",n); } // for galactic cut only in in the first numers given are important
	}
	else {
	  if (type == 1) { if ( a1 > ang) { printf(" -r %li ",n); } } // for n1 direction
	  if (type == 2 && a1 > ang && a2 > ang) { printf(" -r %li ",n); }  // for n2 direction
	  if (type == 3 && fabs(b) > 180/PI*(PIsnd-refb1)) { printf(" -r %li ",n); } // for galactic cut only in in the first numers given are important
      }
      }
      fclose(f);
      printf("\n");
    }
    else  {
      if (neg == 0 || neg==100) printf("%lE ",cpeds_cosang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 1 || neg==100) printf("%lE ",PI180inv*cpeds_ang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 2 || neg==100) printf("%lE ",cpeds_ang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 3 || neg==100) printf("%lE ",sin(cpeds_ang_n1n2(refb1,refl1, refb2, refl2)));
      if (neg == 4 || neg==100) { 
	if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",cpeds_cosang_n1n2(refb1,refl1, refb2, refl2)); } 
	else { printf("%lE ",cos(PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); }
      }
      if (neg == 5 || neg==100) { 
	if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",sin(cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); } 
	else { printf("%lE ",sin(PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); }
      }
      if (neg == 6 || neg==100) { 
	if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",cpeds_ang_n1n2(refb1,refl1, refb2, refl2)); } 
	else { printf("%lE ",PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2)); }
      }
      printf("\n");  

    }
  }
  else { // reading from STDIN
    while (scanf("%lE %lE %lE %lE",&refl1,&refb1, &refl2, &refb2)!=EOF) {
      refl1 *= PI180;  refb1 *= PI180;   refb1 = PIsnd-refb1; ang *=PI180;
      refl2 *= PI180;  refb2 *= PI180;   refb2 = PIsnd-refb2;

      if (neg == 0 || neg==100) printf("%lE ",cpeds_cosang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 1 || neg==100) printf("%lE ",PI180inv*cpeds_ang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 2 || neg==100) printf("%lE ",cpeds_ang_n1n2(refb1,refl1, refb2, refl2));
      if (neg == 3 || neg==100) printf("%lE ",sin(cpeds_ang_n1n2(refb1,refl1, refb2, refl2)));
      if (neg == 4 || neg==100) { 
      if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",cpeds_cosang_n1n2(refb1,refl1, refb2, refl2)); } 
      else { printf("%lE ",cos(PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); }
      }
      if (neg == 5 || neg==100) { 
	if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",sin(cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); } 
	else { printf("%lE ",sin(PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2))); }
      }
      if (neg == 6 || neg==100) { 
	if (cpeds_cosang_n1n2(refb1,refl1, refb2, refl2) >=0) { printf("%lE ",cpeds_ang_n1n2(refb1,refl1, refb2, refl2)); } 
	else { printf("%lE ",PI-cpeds_ang_n1n2(refb1,refl1, refb2, refl2)); }
      }
      printf("\n");
    };
  }
}

