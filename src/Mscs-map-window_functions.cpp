// THIS IS NOT USED AT THE MOMENT
// Map-scs USES THE window_function CLASS INSTEAD AND AN ARRAY OF THESE OBJECTS
// ALTHOUGH THIS OBJECT IS COOL :)
// for now it's just wasted worl :(

#include "Mscs-map-window_functions.h"

// CONSTRUCTOR
window_functions::window_functions(long num, long lsize) { 
  long i;
  // initiate the window function space

  wf = new w_f[num*lsize]; // initiate the array of window functions 
  wf_name = new string[num];

  min_l = new double[num];  max_l = new double[num];
  min_wl = new double[num];  max_wl = new double[num];

  min_li = new long[num];  max_li = new long[num];
  min_wli = new long[num];  max_wli = new long[num];

  
  // initiate w_f parameters
  w_f_num = num;
  size = lsize;
  min_i = 0;
  max_i = lsize-1;

  // zero other info
  for (i=0;i<lsize;i++) { min_li[i] = max_li[i] = min_wli[i] = max_wli[i] = 0; }
}

// DESTRUCTOR
window_functions::~window_functions() {
  delete wf;
}


// HANDLERS
double window_functions::get_l(long j, long i) {
  double ret;
  if (check_irange(j,i))  ret = wf[j*size+i].l; else ret = 0;
  return ret;
}

double window_functions::get_wl(long j, long i) {
  double ret;
  if (check_irange(j,i)) ret = wf[j*size+i].wl; else ret = 0;
  return ret;
}

void window_functions::set_l(long j, long i, double l) {
  wf[j*size+i].l = l;
}

void window_functions::set_wl(long j, long i, double wl) {
  wf[j*size+i].wl = wl;
}

string window_functions::get_name(long j) {
  return wf_name[j];
}

void window_functions::set_name(long j, string wfname) {
  wf_name[j]  = wfname;
}

// this returns the range of table indexes that are allowed
void window_functions::get_irange(long *_min_i, long *_max_i) {
  *_min_i = min_i; *_max_i = max_i;
}

long window_functions::get_size() {
  return size;
}

void window_functions::clear(long j) {
/*  to be implemented */
}

void window_functions::get_wfrange(long j, double * _min_l, double * _max_l, double * _min_wl, double * _max_wl) {
  calculate_w_f_minmax_values(j);
  *_min_l = min_l[j];   *_max_l = max_l[j];
  *_min_wl = min_wl[j];   *_max_wl = max_wl[j];
}


//--------------------------------------------------------------------------------
// prrivate methods
//--------------------------------------------------------------------------------

// checking index range passed to handlers
bool window_functions::check_irange(long j, long i) {
  bool ret;
  if ((j < w_f_num) && (j >= 0) && (i <= max_i) && (i >= min_i)) ret = true;   else ret = false;
  return ret;
}

void window_functions::calculate_w_f_minmax_values(long j) {
  long i;

  min_l[j] = wf[j*size].l;  max_l[j] = wf[j*size].l;
  min_wl[j] = wf[j*size].wl;  max_wl[j] = wf[j*size].wl;

  //for (j=0;j<w_f_num;j++) { 
    for (i=min_i;i<=max_i;i++) { 
      if (wf[j*size+i].l < min_l[j]) { min_l[j] = wf[j*size+i].l; min_li[j] = i; }
      if (wf[j*size+i].l > max_l[j]) { max_l[j] = wf[j*size+i].l; max_li[j] = i; }
      if (wf[j*size+i].wl < min_wl[j]) { min_wl[j] = wf[j*size+i].wl; min_wli[j] = i; }
      if (wf[j*size+i].wl > max_wl[j]) { max_wl[j] = wf[j*size+i].wl; max_wli[j] = i; }
    }
  //}

}
