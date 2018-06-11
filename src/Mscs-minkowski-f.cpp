#include "Mscs-minkowski-f.h"
#include "Mscs-map.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF THE CLASS class minkowski_f
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// CONSTRUCTOR
minkowski_f::minkowski_f(string ftype, long th_num) {
	//  long i;
	// initiate the functional space
	//  f = new v[th_num];
	// initiate functional parameters
	size = th_num;
	min_i = 0;
	max_i = th_num-1;
	func_type = ftype;
	
	// zero other info
	min_nu = max_nu = min_vnu = max_vnu = 0;
	//  for (i=min_i;i<=max_i;i++) {     f[i].nu=0.0; f[i].vnu=0.0;  }
	
	if (th_num>0) setPointsNum(th_num);
	//  printf("|%s> * initiating minkowski functional object. size: %li\n",func_type.c_str(),th_num);
}

// DESTRUCTOR
minkowski_f::~minkowski_f() {
	//  delete [] f;
}


// HANDLERS
//double minkowski_f::get_nu(long i) {
//  double ret;
//  if (check_irange(i))  ret = f[i].nu; else ret = 0;
//  return ret;
//}
//
//double minkowski_f::get_vnu(long i) {
//  double ret;
//  if (check_irange(i)) ret = f[i].vnu; else ret = 0;
//  return ret;
//}
//
//void minkowski_f::set_nu(long i, double nu) {
//  f[i].nu = nu;
//}
//
//void minkowski_f::set_vnu(long i, double vnu) {
//  f[i].vnu = vnu;
//}

string minkowski_f::get_type() {
	return func_type;
}

void minkowski_f::set_type(string ftype) {
	func_type  = ftype;
}

// this returns the range of table indexes that are allowed
//void minkowski_f::get_irange(long *_min_i, long *_max_i) {
//  *_min_i = min_i; *_max_i = max_i;
//}
// returns the number of thresholds
long minkowski_f::get_size() {
	return pointsCount();
}

void minkowski_f::clear() {
	/*  to be implemented */
	clearFunction();
}

//void minkowski_f::get_frange(double * _min_nu, double * _max_nu, double * _min_vnu, double * _max_vnu) {
//  calculate_functional_minmax_values();
//  *_min_nu = min_nu;   *_max_nu = max_nu;
//  *_min_vnu = min_vnu;   *_max_vnu = max_vnu;
//}

// normalize the minkowski functional
void minkowski_f::normalize_V0() {
	long i;
	long sizele1=get_size()-1;
	for (i=0;i<get_size();i++) setf(i,f(i)/f(sizele1));
	
}

//void minkowski_f::save(string name) {
//  FILE *f;
//  long i;
//  f=fopen(name.c_str(),"w");
//  for (i=0;i<size;i++) {    fprintf(f,"%lE %lE\n",get_nu(i),get_vnu(i)); }
//  fclose(f);
//}

//--------------------------------------------------------------------------------
// prrivate methods
//--------------------------------------------------------------------------------

// checking index range passed to handlers
/*
bool minkowski_f::check_irange(long i) {
  bool ret;
  if ((i <= max_i) && (i >= min_i)) ret = true;   else ret = false;
  return ret;
}

void minkowski_f::calculate_functional_minmax_values() {
  long i;

  min_nu = f[0].nu;  max_nu = f[0].nu;
  min_vnu = f[0].vnu;  max_vnu = f[0].vnu;

  for (i=min_i;i<=max_i;i++) { 
    if (f[i].nu < min_nu) { min_nu = f[i].nu; min_nui = i; }
    if (f[i].nu > max_nu) { max_nu = f[i].nu; max_nui = i; }
    if (f[i].vnu < min_vnu) { min_vnu = f[i].vnu; min_vnui = i; }
    if (f[i].vnu > max_vnu) { max_vnu = f[i].vnu; max_vnui = i; }
  }

}
 */






















////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF THE CLASS class minkowski_fs
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// general constructor
minkowski_fs::minkowski_fs() {
	fnum=0;
	mfs=NULL;
}

// some other general constructor
minkowski_fs::minkowski_fs(long num) {
	long i;
	fnum = num;
	mfs=new mfs_type[fnum];
	
	for (i=0;i<fnum;i++) {
		mfs[i].v0=NULL;
		mfs[i].v1=NULL;
		mfs[i].v2=NULL;
		mfs[i].regid=0;
		mfs[i].size=0;
	}
	
}

// constructor for a given common number of thresholds and number of set of functionals
minkowski_fs::minkowski_fs(long num, long thnum) {
	long i;
	string name;
	filenamestr tmpch;
	
	fnum = num;
	mfs=new mfs_type[fnum];
	
	for (i=0;i<fnum;i++) {
		sprintf(tmpch,"F%li_%s",i,"v0"); name=tmpch; mfs[i].v0 = new minkowski_f(name,thnum);
		sprintf(tmpch,"F%li_%s",i,"v1"); name=tmpch; mfs[i].v1 = new minkowski_f(name,thnum);
		sprintf(tmpch,"F%li_%s",i,"v2"); name=tmpch; mfs[i].v2 = new minkowski_f(name,thnum);
		mfs[i].regid=i;
		mfs[i].size=0;
	}
	
	object_name="minkowski_functionals";
}

/****************************************************************************************************/


minkowski_fs::~minkowski_fs() {
	fnum=0;
	delete [] mfs;
	mfs=NULL;
}

// returns the number of sets of minkowski functionals stored in the object
long minkowski_fs::get_size() { return fnum; }
void minkowski_fs::set_name(string name) { object_name=name; }
string minkowski_fs::get_name() { return object_name; }

// saves the functionals to a file
void minkowski_fs::save(string name, long how) {
	FILE * f;
	long thres_num,i,j;
	
	if (how == 1) {
		name+=".mink";
		printf("|%s> * saving the minkowski functionals to file %s\n",object_name.c_str(),name.c_str());
		thres_num=mfs[0].v0->get_size(); // assume the same number of threshods for all functional types and all regions
		f = fopen(name.c_str(),"w");
		fprintf(f,"# %li %li\n",fnum, thres_num);
		
		for (j=0;j<fnum;j++) { 
			for (i=0;i<thres_num;i++) {
				fprintf(f,"%li %li %lE %lE %lE %lE\n",mfs[j].regid,mfs[j].size, mfs[j].v0->get_nu(i),mfs[j].v0->get_vnu(i),mfs[j].v1->get_vnu(i),mfs[j].v2->get_vnu(i));
			}
		} 
		
		fclose(f);
	}
	
}

void minkowski_fs::normalize_V0s() {
	long i;
	for (i=0;i<fnum;i++) mfs[i].v0->normalize_V0();
	
}

double minkowski_fs::get_mf(long i,long imf,long inu) {
	double ret=0;
	if (imf==0) { ret=mfs[i].v0->get_vnu(inu); }
	if (imf==1) { ret=mfs[i].v1->get_vnu(inu); }
	if (imf==2) { ret=mfs[i].v2->get_vnu(inu); }
	return ret;
}
























