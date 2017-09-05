// THIS IS NOT USED AT THE MOMENT
// Map-scs USES THE window_function CLASS INSTEAD AND AN ARRAY OF THESE OBJECTS
// ALTHOUGH THIS OBJECT IS COOL :)
// for now it's just wasted worl :(

#include <stdlib.h>
#include <stdio.h>
#include <string>

using namespace std;
class window_functions {

 public:
  window_functions(long num, long lsize); // number of windows functions to initialize and its common size
  ~window_functions();

  double get_l(long j, long i);
  double get_wl(long j, long i);
  void set_l(long j, long i, double nu);
  void set_wl(long j, long i, double vnu);
  string get_name(long j);
  void set_name(long j, string wfname);

  void get_irange(long *_min_i, long *_max_i);
  long get_size(); // they all have the same size
  void clear(long j);
  void get_wfrange(long j, double * _min_l, double * _max_l, double * _min_wl, double * _max_wl);

 private:

  typedef struct {
    double l;
    double wl;
  } w_f;

  long w_f_num;
  w_f * wf; // pointer to array holding the functional
  string * wf_name;
  double *min_l, *max_l, *min_wl, *max_wl; // minimal and maximal values of the functional values in threshold and the value respectively
  long min_i, max_i, size; // minimal and maximal indexes of the functional array and its size
  long *min_li,*max_li,*min_wli,*max_wli; // the indexes of the minimal and maximal nu and vnu values respectively in the object


  bool check_irange(long j, long i);
  void calculate_w_f_minmax_values(long j);

};
