#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include "Mscs-common.h"

/* using namespace std; */
/* using namespace math; */
#ifndef MSCS_PLMS_TAB
#define MSCS_PLMS_TAB

class Plms_tab {
  // this object must also handle the partial files for the case for 512-1024 -- the file exceeds 2Gb and so it must be divided into two or more smaller files
  // and the object must handle them to continously read from them.
 public:

/*   Plms_tab(long nside, long lmax, string plms_file); */
  Plms_tab(string _object_name, long nside, long lmax, strarg plms_file);
  Plms_tab(long nside, long lmax, strarg plms_file);
  ~Plms_tab();
  void load_new(); 
  double Pnext(); // make this smarter
  void Pskip(long skip);
  long get_lmax();
/*   double P(int l,int m,int x) { return 0; // this isn't used yet */
/*   } */

 private:
  string object_name;
  struct stat Pfile_info;
  double * Ptab, tmp; 
  unsigned long int size; // size ~160 Mb
  long int lmin,lmax,mmin,mmax,l,m,x;
  unsigned long int loaded_num,how_much,from_i,to_i,i,j,call_num;
  long int nside_map,nside_file;
  long int lmax_map,lmax_file;
  FILE * fplm;

  void initiate_data(long nside, long lmax, strarg plms_file);

};

#endif

