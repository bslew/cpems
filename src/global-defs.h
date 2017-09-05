#include <string.h>

// definitions of directories where usefull stuff is kept

filenamestr WMAPdata_dir = "/home/blew/cosmo-data/wmap/";
filenamestr tmpWMAPdata_dir = "/home/blew/programy/map-scs/";
filenamestr program_dir = "/home/blew/programy/map-scs/";
filenamestr WMAPcalc_dir = "/home/blew/programy/map-scs/WMAPstuff/";
filenamestr WMAPcalc_dir1 = "/home/blew/programy/map-scs/WMAPstuff/cut-sky-test/";
filenamestr data_dir = "/home/blew/programy/map-scs/data/";
filenamestr WMAPsim_dir = "/scratch1/gaussian_simulations/wmap/";


// definitions of files to be used in programs

filenamestr GLOBAL_path;
long GLOBAL_nside;
long GLOBAL_lmax;
filenamestr GLOBAL_file_prefix;
filenamestr GLOBAL_DA1;
int GLOBAL_DA1i;
filenamestr GLOBAL_DA2;
int GLOBAL_DA2i;
filenamestr GLOBAL_mask_file;
filenamestr GLOBAL_smoothing_method;
int GLOBAL_fourier_method_num;

filenamestr GLOBAL_Plmsfile_forward;
filenamestr GLOBAL_Plmsfile_inverse;


// definitions of common routines for programs of map-scs package


void  make_detailed_file_name(filenamestr *name,
				    strarg pathP,
				    long nsideP,
				    long lmaxP,
				    strarg file_prefixP,
				    strarg DA1P,
    				    int DA1iP,
				    strarg DA2P,
				    int DA2iP,
				    strarg mask_fileP,
				    strarg smoothing_methodP,
				    int fourier_method_numP) {

  strcpy(GLOBAL_path,pathP);
  GLOBAL_nside = nsideP;
  GLOBAL_lmax = lmaxP;
  strcpy(GLOBAL_file_prefix,file_prefixP);
  strcpy(GLOBAL_DA1,DA1P); GLOBAL_DA1i = DA1iP;
  strcpy(GLOBAL_DA2,DA2P); GLOBAL_DA2i = DA2iP;
  strcpy(GLOBAL_mask_file,mask_fileP);
  strcpy(GLOBAL_smoothing_method,smoothing_methodP);
  GLOBAL_fourier_method_num = fourier_method_numP;

  sprintf(GLOBAL_Plmsfile_forward,"%sPlms-%li-%li",data_dir,GLOBAL_nside,2*GLOBAL_nside);
  if ( GLOBAL_fourier_method_num == 7 ) { sprintf(GLOBAL_Plmsfile_inverse,"%sPlms-%li-%li",data_dir,GLOBAL_nside,2*GLOBAL_nside); }
  if ( GLOBAL_fourier_method_num == 8 ) { sprintf(GLOBAL_Plmsfile_inverse,"%sPmthls-%li-%li",data_dir,GLOBAL_nside,2*GLOBAL_nside); }

  sprintf(*name,"%s%li-%li-%s-%s%i-%s%i-%s-sm_%s-Fmet_%i",pathP,nsideP,lmaxP,file_prefixP,DA1P,DA1iP,DA2P,DA2iP,mask_fileP,smoothing_methodP,fourier_method_numP);

}

