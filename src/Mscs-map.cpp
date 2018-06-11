#define _ISOC99_SOURCE
#include <features.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
//#include <string.h>
#include <cpgplot.h>
//#include <fitsio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_coupling.h>
#include <proj_api.h>
//#include "mpspecfun.h" // BLcomment (Apr 14, 2011, 6:33:46 PM): the usage of W3J symbols will be done in SH space or in individual programs; this functionality is not yet restored in this version of Mscs
#include <complex>
#include "matrix.h"

#include "Mscs-map.h"
#include "Mscs-global-defs.h"
//#include "Mscs-power_spectrum.h"
#include "cpeds_angular_correlation_fn.h"
//#include "cpeds-consts.h"

//#include "matpack.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
//using namespace MATPACK; // BLcomment (Apr 14, 2011, 6:35:43 PM): see above comment
/* using namespace LiDIA; */
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

// Implementation Of The Map Class
// most general routines

// CONSTRUCTOR
//************************************************************************
mscsMap::mscsMap(long ns) :
	mscsObject("mscsMap",Zero) {
	set_object_initial_variables();
	set_nside(ns);
}
//************************************************************************
mscsMap::mscsMap(string _object_name, long ns) :
	mscsObject(_object_name,Zero) {
	set_object_initial_variables();
	set_nside(ns);
}

//************************************************************************
// CONSTRUCTOR FOR CLONING OBJECT
mscsMap::mscsMap(const mscsMap &orig) 
:	mscsObject(orig) 
{
	set_object_initial_variables();
	clone(orig);
}

//************************************************************************
mscsMap::mscsMap(const string _object_name, const mscsMap &orig) :
	mscsObject(orig) {
	setName(_object_name);
	set_object_initial_variables();
	clone(orig);
}
//************************************************************************
// CONSTRUCTOR FOR INITIALIZING THE DATA DIRECTORIES INFO
/* mscsMap::mscsMap(strarg data_path_str) { */
/*   set_object_initial_variables(""); */
/*   strcpy(data_path,data_path_str);  */
/* } */
/* mscsMap::mscsMap(string _object_name, strarg data_path_str) { */
/*   set_object_initial_variables(_object_name); */
/*   strcpy(data_path,data_path_str);  */
/* } */

/* mscsMap::mscsMap(string _object_name, package_dirs dirs) { */
/*   set_object_initial_variables(_object_name); */
/*   strcpy(data_path,dirs[7]);  */
/* } */

void mscsMap::resetMaskInfo() {
	mask.merged = false;
	mask.masked_pix_num = 0;
	mask.f_sky = 1.0;
	//    long type;
	mask.multi_mask_lreg = mask.multi_mask_breg = mask.multi_mask_nside
			= mask.multi_mask_reg_num = 0;
}

//************************************************************************
void mscsMap::set_object_initial_variables() {
	// int i;

	// initialization of various variables
	//************************************************************************

	loadedMapFileName = "";
	loadedMaskFileName = "";
	
	map_ptr = &map;
	
	
	/*   MAP INFORMATION BLOCK INITIALIZATION */
	mapInfo.nside = 0;
	mapInfo.pix_num = 0;
	mapInfo.rows_num = 0;
	mapInfo.ordering = nested;
	/* alms_num=0; */
	/* lmax=0; */

	//map = new pixel;  //map_ptr = new pixel; *map_ptr = *map;
	/* map->T = NULL;  map->m = NULL;   map->N = NULL;   map->n = NULL;  // REAL SPACE DATA INITIALIZATION */
	/* map->Q = NULL; map->U = NULL; map->V = NULL; */
	//alm = NULL; // SH SPACE DATA INITIALIZATION

	//flatmap = NULL; flatcoord = NULL;
	/* n2r_conv_tab = NULL; r2n_conv_tab = NULL;  // CONVERSION OF ORDERING TABLES INITIALIZATION */
	/* rotation_map = NULL; // INITIALIZATION OF THE ROTATION MAP */
	/*   exclude_values_list = NULL; */

	//************************************************************************
	// DATA ALLOCATION FLAGS INITIALIZATION
	loaded.T = false;
	/* alms_loaded = 0; */
	/* flat_coord_loaded = 0; */
	loaded.m = false;
	loaded.N = false;
	loaded.n = false;
	loaded.rotation = false;
	//seed_offset = 0;

	n2rConvTabLoaded(0);
	r2nConvTabLoaded(0);
	
	resetMaskInfo();
	/* ok_prev = 0; rek_depthl=1; // this is not used */
	default_fourier_method = 8;
	
	
	//************************************************************************
	// program directories stuff
#pragma omp critical
	Mscs_initiate_global_variables();
	/* strcpy(data_path,PROGRAM_DIRS[7]); // hope this works */

	//************************************************************************
	// window functions stuff
	/* for (i=0;i<NUMBER_OF_WINDOW_FUNCTIONS;i++) bl_tab[i] = NULL;  */

	//************************************************************************
	// spectra stuff
	/* lmax_C_l = 0; lmin_C_l = 0; l_num_C_l = 0; */
	/* for (i=0;i<MAX_L;i++) { C_l[i][0]=C_l[i][1]=C_l[i][2] = 0; } */

	//************************************************************************
	// visualization stuff -- not sure if this is used
	/* vis_sizeofxy = 1000; vis_iter=0; */
	/* vis_xtab = new float[vis_sizeofxy];   vis_ytab = new float[vis_sizeofxy]; */
	/* plot_line_color = 4; */
	/* for (i=0;i<vis_sizeofxy;i++) { vis_xtab[i]=vis_ytab[i] = 0; } */
	/* first_time = true; */
	/* flatmap_minT=flatmap_maxT=0; */
	/* flatmap_change_zrange = false; */
	/* reference_dir.l=0; reference_dir.b=0; */

	//************************************************************************
	// gsl and other external stuff
	// srand(time(0)-seed_offset); // commented out on 2009-01-22 -- moved to cpeds and individual applications
#pragma omp critical
	gsl_set_error_handler_off();
	
	
	//************************************************************************
	//************************************************************************
	// physical stuff
	/* C_2WMAP = 0.12844E-09; // normalization factor;  */

	// definition of temperature uncertainties in different frequency channels
	// WMAP I yr data.
	/* sigma0_K1 = 1.439E-3; sigma0_K2 = 1.464E-3; // K1=K1 K2=Ka1 these values are actually from wmap3 TO BE VERIFIED */
	/* sigma0_Q1 = 2.26677E-3; sigma0_Q2 = 2.15567E-3; */
	/* sigma0_V1 = 3.28789E-3; sigma0_V2 = 2.93683E-3; */
	/* sigma0_W1 = 5.85196E-3; sigma0_W2 = 6.53276E-3; sigma0_W3 = 6.88032E-3; sigma0_W4 = 6.72537E-3; */

	/* sigma0_WMAP[0] = sigma0_K1;   sigma0_WMAP[1] = sigma0_K2; */
	/* sigma0_WMAP[2] = sigma0_Q1;   sigma0_WMAP[3] = sigma0_Q2; */
	/* sigma0_WMAP[4] = sigma0_V1;   sigma0_WMAP[5] = sigma0_V2; */
	/* sigma0_WMAP[6] = sigma0_W1;   sigma0_WMAP[7] = sigma0_W2; */
	/* sigma0_WMAP[8] = sigma0_W3;   sigma0_WMAP[9] = sigma0_W4; */

	// sigma0 in [K] for various DAs of the WMAP
	// WMAP III yr data. from astro-ph/0603451 and astro-ph/0603452
	/* sigma0_K1 = 1.439E-3; sigma0_K2 = 1.464E-3; // K1=K1 K2=Ka1 */
	/* sigma0_Q1 = 2.245E-3; sigma0_Q2 = 2.135E-3; */
	/* sigma0_V1 = 3.304E-3; sigma0_V2 = 2.946E-3; */
	/* sigma0_W1 = 5.883E-3; sigma0_W2 = 6.532E-3; sigma0_W3 = 6.885E-3; sigma0_W4 = 6.744E-3; */

	/* sigma0_WMAP[10] = sigma0_K1;   sigma0_WMAP[11] = sigma0_K2; */
	/* sigma0_WMAP[12] = sigma0_Q1;   sigma0_WMAP[13] = sigma0_Q2; */
	/* sigma0_WMAP[14] = sigma0_V1;   sigma0_WMAP[15] = sigma0_V2; */
	/* sigma0_WMAP[16] = sigma0_W1;   sigma0_WMAP[17] = sigma0_W2; */
	/* sigma0_WMAP[18] = sigma0_W3;   sigma0_WMAP[19] = sigma0_W4; */

	// sigma0 in [K] for various DAs of the WMAP
	// WMAP 5 yr data. from http://lambda.gsfc.nasa.gov/product/map/dr3/pub_papers/fiveyear/basic_results/wmap5basic.pdf
	/* sigma0_K1 = 1.436E-3; sigma0_K2 = 1.470E-3; // K1=K1 K2=Ka1 */
	/* sigma0_Q1 = 2.254E-3; sigma0_Q2 = 2.141E-3; */
	/* sigma0_V1 = 3.314E-3; sigma0_V2 = 2.953E-3; */
	/* sigma0_W1 = 5.899E-3; sigma0_W2 = 6.565E-3; sigma0_W3 = 6.926E-3; sigma0_W4 = 6.761E-3; */

	/* sigma0_WMAP[20] = sigma0_K1;   sigma0_WMAP[21] = sigma0_K2; */
	/* sigma0_WMAP[22] = sigma0_Q1;   sigma0_WMAP[23] = sigma0_Q2; */
	/* sigma0_WMAP[24] = sigma0_V1;   sigma0_WMAP[25] = sigma0_V2; */
	/* sigma0_WMAP[26] = sigma0_W1;   sigma0_WMAP[27] = sigma0_W2; */
	/* sigma0_WMAP[28] = sigma0_W3;   sigma0_WMAP[29] = sigma0_W4; */

	// WMAP Beam sizes - the Full Width at Half Maximum in degrees for different DA in degrees
	/* thFWHM_Q = 0.49;  thFWHM_V = 0.33; thFWHM_W = 0.21;  */
	//strcpy(data_path,"/home/blew/programy/Mscs/data/"); // THIS IS WRONG --- THIS INFORMATION SHOULD BE GIVEN THROUGH PROGRAM_DIRS[] VARIABLE TO THE CONSTRUCTOR
	// THIS IS TO BE REMOVED SOON - when all programs are updated to use the new constructor
	//************************************************************************
}

/* void mscsMap::set_object_name(string name) { */
/*   printf("|%s> * changing name to \"%s\":\n",object_name.c_str(),name.c_str()); */
/*   object_name=name; */
/* } */
/* string mscsMap::get_object_name() { */
/*   return object_name; */
/* } */

//************************************************************************
string mscsMap::getLoadedMaskFileName() const {
	return loadedMaskFileName;
}

void mscsMap::setLoadedMaskFileName(string name) {
	loadedMaskFileName = name;
}

string mscsMap::getLoadedMapFileName() const {
	return loadedMapFileName;
}

void mscsMap::setLoadedMapFileName(string name) {
	loadedMapFileName = name;
}

//************************************************************************
// THIS METHOD IS TO CLONE OBJECT'S DATA INTO A NEW OBJECT. IS EVERYTHING IS COPIED ?
void mscsMap::clone(const mscsMap &orig) {
	// long i;

	msgs->say("cloning from object: " + orig.getName(), Medium);
//
//	mapInfo = orig.mapInfo;
//	mask = orig.mask;
//	map = orig.map;
	
	(*this)=orig;

}

//************************************************************************

// DESTRUCTOR

//************************************************************************
mscsMap::~mscsMap() {
	/* long i; */

	/* kill_conv_tabs(); */

	/* kill_map();  //delete map_ptr;    */
	/* if (rotation_map != NULL) { delete [] rotation_map; loaded.rotation = false; } */

}
//************************************************************************


//-------------------------------------------------------------
// MAP OPERATION  METHODS
//-------------------------------------------------------------
//************************************************************************
// loads the ordering conversion arrays from file for a given nside
// if it's not possible to open the array files for read then -1 value is returned
int mscsMap::load_conv_tabs(long ns, string type) {
	string fileName;
	// cpedsList<long> q;
	long ret = 0;
	
	msgs->say("loading arrays for N2R and R2N conversions", Medium);
	msgs->say(
			"current status: n2r_conv_tab_loaded = " + msgs->toStr(
					n2rConvTabLoaded()) + " r2n_conv_tab_loaded = "
					+ msgs->toStr(r2nConvTabLoaded()), Low);
	
	
	/* printf("|%s>  -- current status: n2r_conv_tab_loaded = %li, r2n_conv_tab_loaded = % li\n",object_name.c_str(),n2r_conv_tab_loaded,r2n_conv_tab_loaded); */
	// kill old tabs if they exist.
	if (type == "n2r") if (n2rConvTabLoaded() != 0) {
		n2r_conv.clear();
	}
	if (type == "r2n") if (r2nConvTabLoaded() != 0) {
		r2n_conv.clear();
	}

	// load new tabs
	// reading for nested2ring conversion
	if (type == "n2r") {
		fileName = MSCS_DATA_DIR + MSCS_GLOBAL__N2R_CONV_TAB_PREF
				+ msgs->toStr(ns) + MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN;
		msgs->say("trying to load conversion arrays from file: " + fileName,
				Low);
		ret = n2r_conv.load(fileName, true, "long", false);
		if (ret != 0) {
			msgs->error("There was problem in loading the file:" + fileName,
					Low);
			return ret;
		}
		msgs->say("success", Low);
		if (n2r_conv.size() != cpeds_get_healpix_pix_num(nside())) {
			msgs->error(
					"The conversion table size is invalid (" + msgs->toStr(
							n2r_conv.size()) + ", should be: " + msgs->toStr(
							cpeds_get_healpix_pix_num(nside()))
							+ "). Will not load", Low);
			n2r_conv.clear();
			return -1;
		}
		n2rConvTabLoaded(ns);
	}

	// reading for ring2nested conversion
	if (type == "r2n") {
		fileName = MSCS_DATA_DIR + MSCS_GLOBAL__R2N_CONV_TAB_PREF
				+ msgs->toStr(ns) + MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN;
		msgs->say("trying to load conversion arrays from file: " + fileName,
				Low);
		ret = r2n_conv.load(fileName, true, "long", false);
		if (ret != 0) {
			msgs->error("There was problem in loading the file:" + fileName,
					Low);
			return ret;
		}
		msgs->say("success", Low);
		if (r2n_conv.size() != cpeds_get_healpix_pix_num(nside())) {
			msgs->error(
					"The conversion table size is invalid (" + msgs->toStr(
							r2n_conv.size()) + ", should be: " + msgs->toStr(
							cpeds_get_healpix_pix_num(nside()))
							+ "). Will not load", Low);
			r2n_conv.clear();
			return -1;
		}
		r2nConvTabLoaded(ns);
	}
	msgs->say(
			"current status: n2r_conv_tab_loaded = " + msgs->toStr(
					n2rConvTabLoaded()) + " r2n_conv_tab_loaded = "
					+ msgs->toStr(r2nConvTabLoaded()), Low);
	return ret;
}

//************************************************************************
// sets the ordering conversion arrays for a given nside
void mscsMap::set_conv_tabs(long nside) {
	// for now not used: to be implemented
}

//************************************************************************
// kills the conversion tables
void mscsMap::kill_conv_tabs() {
	msgs->say("deleting conversion arrays", Medium);
	/* if (n2r_conv_tab != NULL ) { delete n2r_conv_tab; n2r_conv_tab = NULL; }  */
	n2r_conv.clear();
	n2rConvTabLoaded(0);
	/* if (r2n_conv_tab != NULL ) { delete r2n_conv_tab; r2n_conv_tab = NULL; }  */
	r2n_conv.clear();
	r2nConvTabLoaded(0);
	
}

//************************************************************************
// converts the map from the ring ordering to nested ordering using CPEDS functions

// void mscsMap::conv_ring2nest(mscsMap * map_loc) {
void mscsMap::conv_ring2nest() {
	
	if (ordering()==mscsMap::nested) return;
	long int i;
	mscsMap maptmp(*this);
	// if (map_loc==NULL) map_loc=this;

	msgs->say("converting to nested ordering", Medium);
	// maptmp = *map_loc;

	if (r2nConvTabLoaded() == nside()) {
		msgs->say("R2N tabs are loaded - converting...", Low);
	} else {
		if (r2nConvTabLoaded() == 0) {
			msgs->say("R2N tabs are NOT loaded - trying to load from file...",
					Low);
		} else {
			msgs->say(
					"wrong R2N tabs loaded - trying to load new ones from file",
					Low);
			kill_conv_tabs();
		}
		if (load_conv_tabs(nside(), "r2n") == 0) {
			msgs->say("R2N conversion tabs loaded", Low);
		} else {
			msgs->say(
					"The is no file for R2N coversion. Generating conversion table (this is SLOW).",
					Low);
			// generating the conversion table
			double *mapRING, *r2n_conv_tabd;
			long pix_num = pixNum();
			r2n_conv_tabd = new double[pix_num];
			mapRING = new double[pix_num];
			for (i = 0; i < pix_num; i++) {
				mapRING[i] = (double) i;
			}
			cpeds_conv_ring2nest_healpix(nside(), mapRING, r2n_conv_tabd);
			r2n_conv.clear();
			for (i = 0; i < pix_num; i++) {
				r2n_conv.append((long) r2n_conv_tabd[i]);
			}
			r2nConvTabLoaded(nside());
			delete r2n_conv_tabd;
			delete mapRING;
			r2n_conv.save(
					MSCS_DATA_DIR + MSCS_GLOBAL__R2N_CONV_TAB_PREF
							+ msgs->toStr(nside())
							+ MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN, true, "long");
			r2n_conv.save(
					MSCS_DATA_DIR + MSCS_GLOBAL__R2N_CONV_TAB_PREF
							+ msgs->toStr(nside())
							+ MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT, false, "long");
		}
	}
	
	
	// converting map structures
	copy_map(&maptmp, this, &r2n_conv);
	/* kill_map_space(map_loc); // kill the old map space; map will now point to new momory areas */
	/* if (map_loc == map) { map_ordering = 1; *map = maptmp; } // if the change of ordering was asked for the object's map then set the flat accordingly */
	/* *map_loc = maptmp;  */

	// if (map_loc==this) *this=maptmp;
	// *map_loc=maptmp;
	set_ordering(mscsMap::nested);
}
//************************************************************************
// converts the map from the nest ordering to ring ordering using CPEDS functions or precalculated conversion arrays
// void mscsMap::conv_nest2ring(map_structure * map_loc) {
void mscsMap::conv_nest2ring() {
	if (ordering()==mscsMap::ring) return;

	long int i;
	mscsMap maptmp(*this);
	// if (map_loc==NULL) map_loc=&map;

	msgs->say("converting to ring ordering", Medium);
	// maptmp = *map_loc;

	// if we have conversion tables ready
	if (n2rConvTabLoaded() == nside()) {
		msgs->say("N2R tabs are loaded - converting..", Low);
	} else { // if not then try to load
		if (n2rConvTabLoaded() == 0) {
			msgs->say("N2R tabs are NOT loaded - trying to load from file", Low);
		} else {
			msgs->say(
					"wrong N2R tabs loaded - trying to load new ones from file",
					Low);
			kill_conv_tabs();
		}
		if (load_conv_tabs(nside(), "n2r") == 0) {
			msgs->say("N2R conversion tabs loaded", Low);
		} // if there is no file then we have to do all the job
		else {
			msgs->say(
					"The is no file for N2R coversion. Converting in the usual way (this is SLOW).",
					Low);
			// generating the conversion table
			double * mapNEST, *n2r_conv_tabd;
			long pix_num = pixNum();
			n2r_conv_tabd = new double[pix_num];
			mapNEST = new double[pix_num];
			for (i = 0; i < pix_num; i++) {
				mapNEST[i] = (double) i;
			}
			cpeds_conv_nest2ring_healpix(nside(), mapNEST, n2r_conv_tabd);
			n2r_conv.clear();
			for (i = 0; i < pix_num; i++) {
				n2r_conv.append((long) n2r_conv_tabd[i]);
			}
			n2rConvTabLoaded(nside());
			delete n2r_conv_tabd;
			delete mapNEST;
			n2r_conv.save(
					MSCS_DATA_DIR + MSCS_GLOBAL__N2R_CONV_TAB_PREF
							+ msgs->toStr(nside())
							+ MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN, true, "long");
			n2r_conv.save(
					MSCS_DATA_DIR + MSCS_GLOBAL__N2R_CONV_TAB_PREF
							+ msgs->toStr(nside())
							+ MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT, false, "long");
		}
	}
	
	
	// converting map structures
	copy_map(&maptmp, this, &n2r_conv);
	/* kill_map_space(map_loc); // kill the old map space; map will now point to new momory areas */
	/* if (map_loc == map) { map_ordering = 0; *map = maptmp; } // if the change of ordering was asked for the object's map then set the flat accordingly */

	// if (map_loc==map_ptr) map=maptmp;
	// *map_loc = maptmp;
	set_ordering(mscsMap::ring);
}

//************************************************************************
// This routine does not check whtether the mask is even allocated/loaded etc. so
// be careful; check first with maskAllocated();
bool mscsMap::isMasked(long i) {
	if (get_m(i) >= 1.0) return false;
	return true;
}
//************************************************************************
bool mscsMap::maskAllocated() {
	if (get_m().count() > 0) return false;
	return true;
}

/* **************************************************************************************************** */
void mscsMap::TaddVal(double val) {
	long i, pix_num = pixNum();
	if (val != 0.0) {
		if (maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) += val;
			}
		} else {
			T() += val;
			/* for (i=0;i<pix_num;i++) { map->T[i]+=val; } */
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::TsubVal(double val) {
	long i, pix_num = pixNum();
	if (val != 0.0) {
		if (maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) -= val;
			}
		} else {
			T() -= val;
			/* for (i=0;i<pix_num;i++) { map->T[i]-=val; } */
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::TmulVal(double val) {
	long i, pix_num = pixNum();
	if (val != 1.0) {
		if (maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) *= val;
			}
		} else {
			T() *= val;
			/* for (i=0;i<pix_num;i++) { map->T[i]*=val; } */
		}
	}
}

/* **************************************************************************************************** */
// this returns false if you try to divide by 0.
bool mscsMap::TdivVal(double val) {
	long i, pix_num = pixNum();
	bool r = true;
	
	if (val != 0.0) {
		if (maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) /= val;
			}
		} else {
			T() /= val;
			/* for (i=0;i<pix_num;i++) { map->T[i]/=val; } */
		}
	} else {
		r = false;
	}
	return r;
}

/* void mscsMap::xorVal(long val) { */
/*   long i; */
/*   if (val!=0.0) { */
/*     if (mask_loaded == 1) { */
/*       for (i=0;i<pix_num;i++) { if (!isMasked(i)) map->T[i]+=val; } */
/*     } */
/*     else {        */
/*       for (i=0;i<pix_num;i++) { map->T[i]+=val; } */
/*     } */
/*   } */

/* } */

/* **************************************************************************************************** */
void mscsMap::cutAbove(double val) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && get_T(i) > val) {
				set_T(i, val);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) > val) {
				set_T(i, val);
			}
		}
	}
}

/* **************************************************************************************************** */

void mscsMap::cutBelow(double val) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && get_T(i) < val) {
				set_T(i, val);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) < val) {
				set_T(i, val);
			}
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::replaceIfGreater(double val, double val2) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && get_T(i) > val) {
				set_T(i, val2);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) > val) {
				set_T(i, val2);
			}
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::replaceIfLess(double val, double val2) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && get_T(i) < val) {
				set_T(i, val2);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) < val) {
				set_T(i, val2);
			}
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::replaceIfAbsGr(double val, double val2) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && fabs(get_T(i)) > val) {
				set_T(i, val2);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (fabs(get_T(i)) > val) {
				set_T(i, val2);
			}
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::replaceIfAbsLe(double val, double val2) {
	long i, pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && fabs(get_T(i)) < val) {
				set_T(i, val2);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (fabs(get_T(i)) < val) {
				set_T(i, val2);
			}
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::TaddMap(const mscsMap& mapP) {
	long i, pix_num = pixNum();
	mscsMap Clone("tmpclone", mapP);
	Clone.change_map_resolution(nside());
	
	if (maskLoaded()) {
		if (Clone.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i) && !Clone.isMasked(i)) T(i) += Clone.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) += Clone.get_T(i);
			}
		}
	} else {
		if (Clone.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!Clone.isMasked(i)) T(i) += Clone.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				T(i) += Clone.get_T(i);
			}
		}
	}
}
/* **************************************************************************************************** */
void mscsMap::TsubMap(const mscsMap& mapP) {
	long i, pix_num = pixNum();
	mscsMap copy("tmpcopy", mapP);
	copy.change_map_resolution(nside());
	
	if (maskLoaded()) {
		if (copy.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i) && !copy.isMasked(i)) T(i) -= copy.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) -= copy.get_T(i);
			}
		}
	} else {
		if (copy.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!copy.isMasked(i)) T(i) -= copy.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				T(i) -= copy.get_T(i);
			}
		}
	}
	
}

/* **************************************************************************************************** */
void mscsMap::TmulMap(const mscsMap& mapP) {
	long i, pix_num = pixNum();
	mscsMap copy("tmpcopy", mapP);
	copy.change_map_resolution(nside());
	
	if (maskLoaded()) {
		if (copy.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i) && !copy.isMasked(i)) T(i) *= copy.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				if (!isMasked(i)) T(i) *= copy.get_T(i);
			}
		}
	} else {
		if (copy.maskLoaded()) {
			for (i = 0; i < pix_num; i++) {
				if (!copy.isMasked(i)) T(i) *= copy.get_T(i);
			}
		} else {
			for (i = 0; i < pix_num; i++) {
				T(i) *= copy.get_T(i);
			}
		}
	}
	
}

/* **************************************************************************************************** */
// the division will not be performed on pixels that have 0 value.
void mscsMap::TdivMap(const mscsMap& mapP) {
	long i, pix_num = pixNum();
	mscsMap copy("tmpcopy", mapP);
	copy.change_map_resolution(nside());
	// double T;

	// prepare the divisor;
	if (copy.maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (copy.T(i) == 0) copy.m(i) = 0.0;
		}
	} else {
		copy.makekill_space_manager("make", "m");
		copy.clean_mask();
		for (i = 0; i < pix_num; i++) {
			if (copy.T(i) == 0) copy.m(i) = 0.0;
		}
	}
	//  printf(" **** divisor=%lE mask=%lE\n",copy.map->T[0],copy.map->m[0]);
	

	// divide maps
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i) && !copy.isMasked(i)) T(i) /= copy.get_T(i);
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (!copy.isMasked(i)) {
				T(i) /= copy.get_T(i);
			}
		} //printf(" result=%lE\n",map->T[i]); } }
	}
}

//************************************************************************
// new_nside must be number 2^N where N=1,2,3,4,5,6..../
// this changes the map ordering to nested if it was ring originally !!!
void mscsMap::change_map_resolution(long int new_nside) {
	long int i, j, l, average_pix, factor;
	bool wasring = false;
	// long int new_pix_num = 12*new_nside*new_nside;
	map_structure new_map;
	double Tl, Nl, Ql, Vl, Ul, ml, mmask_num; // l - indicates local values
	long pix_num = pixNum();
	
	
	//  if (nside() != new_nside) { new_map = map;  }

	if (nside() > new_nside) { // doing degrading resloution
		msgs->say(
				"degrading map resolution from " + msgs->toStr(nside())
						+ " to: " + msgs->toStr(new_nside), Medium);
		if (ordering() == ring) {
			conv_ring2nest();
			wasring = true;
		}
		factor = (long int) (round((double) nside() / (double) new_nside));
		average_pix = factor * factor;
		for (i = 0; i < pix_num; i++) {
			Tl = Ql = Ul = Vl = ml = Nl = mmask_num = 0;
			l = i + average_pix;
			if (loaded.T) for (j = i; j < l; j++) {
				Tl += get_T(j);
			}
			if (loaded.Q) for (j = i; j < l; j++) {
				Ql += get_Q(j);
			}
			if (loaded.U) for (j = i; j < l; j++) {
				Ul += get_U(j);
			}
			if (loaded.V) for (j = i; j < l; j++) {
				Vl += get_V(j);
			}
			if (loaded.N) for (j = i; j < l; j++) {
				Nl += get_N(j);
			}
//			if (loaded.m) for (j = i; j < l; j++) {
//				if (get_m(i) == 0) mmask_num++;
//				else mmask_num = get_m(i);
//			}
			if (loaded.m) for (j = i; j < l; j++) {
				mmask_num += get_m(j);
			}

			/* pix /= (double)average_pix; // calculate the average temperature in the future pixel and new number of observations */
			if (loaded.T) {
				new_map.T.append(Tl / (double) average_pix);
			}
			if (loaded.Q) {
				new_map.Q.append(Ql / (double) average_pix);
			}
			if (loaded.U) {
				new_map.U.append(Ul / (double) average_pix);
			}
			if (loaded.V) {
				new_map.V.append(Vl / (double) average_pix);
			}
			if (loaded.N) {
				new_map.N.append(Nl / (double) average_pix);
			} // change w.r.t to versions <1.0
			if (loaded.m) {
				if ((long) (mmask_num) == average_pix) new_map.m.append(1.0);
				else new_map.m.append(0.0);
			}
			i += average_pix - 1;
			//      k++;
			//printf("nside=%li, newnside=%li i=%li, j=%li, k=%li\n",nside, new_nside,i,j,k);
		}
	} else
		
		if (nside() < new_nside) { // doing progradig resloution
			msgs->say(
					"prograding map resolution from " + msgs->toStr(nside())
							+ " to: " + msgs->toStr(new_nside), Medium);
			if (ordering() == ring) {
				conv_ring2nest();
				wasring = true;
			}
			factor = new_nside / nside();
			average_pix = factor * factor;
			for (i = 0; i < pix_num; i++) {
				/* Tl=Ql=Ul=Vl=ml=Nl=mmask_num = 0; */
				for (j = 0; j < average_pix; j++) {
					if (loaded.T) {
						new_map.T.append(get_T(i));
					}
					if (loaded.Q) {
						new_map.Q.append(get_Q(i));
					}
					if (loaded.V) {
						new_map.V.append(get_V(i));
					}
					if (loaded.U) {
						new_map.U.append(get_U(i));
					}
					if (loaded.N) {
						new_map.N.append(get_N(i) / (double) average_pix);
					}
					if (loaded.m) {
						new_map.m.append(get_m(i));
					}
					//k++;
				} // copy temeratures to the new pixels from old pixel
			}
		}
	
	if (nside() != new_nside) {
		set_nside(new_nside);
		/* pix_num = new_pix_num; */
		/* rows_num = cpeds_get_ring_num_healpix(nside); */
		/* coord_num = pix_num; */
		/* kill_map_space(map); // but leave the flags unchanged */
		map = new_map; //*map_ptr = *map;
		if (loaded.n) {
			set_map_coord();
		} // if the coords were set, then we must set them again for a new resolution
		if (wasring) {
			conv_nest2ring();
		}
		check_mask();
		msgs->say(
				"new map nside: " + msgs->toStr(nside()) + " pix_num: "
						+ msgs->toStr(pixNum()), Low);
		// delete new_map;
	}
}
//************************************************************************
// rotates the map by requested angles in th,phi,l coordinates given in degr. in gal. CS
// the rotation to th,phi direction is performed along the geodesic by default (if rot_short is true)
// it is done by the following rotation:
// Rz(PI/2-phi)Rx(th)Rz(-PI/2+phi) * Rz(l)
// first part is the transformation of the CS to orientation siutable for a single rotation by th around x axis
// and followed by rotation to restore the old CS.
// the first Rz(l)  is the initial rotation
// If the rot_short=false then the first rotation in phi is not performed and so the whole rotation is not
// equivalent to "draging" the North Pole to the requested orientation
// In order to avoid the holes in the output map all rotations are done in the inversed order by -angles
// given to the procedure. This is only a technical trick.

int mscsMap::rotate_map(double Ath, double Aphi, double Al, bool rot_short,
		string what) {
	long i, j;
	/* direction n,*Xn; */
	cpedsDirectionSet n;
	cpedsList<double> Xmap, *XmapP;
	long pix_num = pixNum();
	
	msgs->say(
			"rotating map: th: " + msgs->toStr(Ath) + " phi: " + msgs->toStr(
					Aphi) + " l: " + msgs->toStr(Al), Medium);
	if (Ath != 0 || Aphi != 0 || Al != 0) {
		if (!coordLoaded()) {
			set_map_coord(0, 0);
		}

		n = get_n(); // copy the map directions; these will be rotated

		// define what is to be rotated
		if (what == "m") {
			XmapP = &map.m;
			Xmap = map.m;
		}
		if (what == "T") {
			XmapP = &map.T;
			Xmap = map.T;
		} // set the pointer of data to be rotated
		if (what == "rot") {
			XmapP = &rotation_map;
			Xmap = rotation_map;
		}
		if (what == "all") {
			prepare_rotation(Ath, Aphi, Al);
		}
		
		if (what != "all") {
			Al *= PI180;
			Ath *= PI180;
			Aphi = (90 - Aphi) * PI180; // convert to radians

			// rotate directions
			if (Ath != 0) { // rotate in b (around y axis defined as l=90, b=0, x axix: l=0 b=0, z axis: b=90)
				if (Aphi != 0) {
					rotate_map_directions(n, 3, Aphi);
				} // Rz(-phi) ; everything is reversed to avoid the holes in the output map
				if (Ath != 0) {
					rotate_map_directions(n, 1, Ath);
				} // Rx(-th)
				if (Aphi != 0 && rot_short) {
					rotate_map_directions(n, 3, -Aphi);
				} // Rz(phi)
			}
			if (Al != 0 && rot_short) {
				rotate_map_directions(n, 3, -Al);
			} // Rz(l)
			

			for (i = 0; i < pix_num; i++) { // make rotated data
				cpeds_ang2pix_healpix(nside(), &j, PIsnd - n[i].b(), n[i].l(),
						long(1)); // check the range and find new pix number
				(*XmapP)[i] = Xmap[j];
			}
		} else {
			rotate_quick("T");
			rotate_quick("Q");
			rotate_quick("V");
			rotate_quick("U");
			rotate_quick("N");
			rotate_quick("m");
		}
	}
	
	
	//if  (what == "all") map->m = Xn; // to be done
	return 0;
}
//************************************************************************
// prepares the convertion map for pixels for a requested rotation. the array is stored in the rotation_map structure for further use with rotate_quick
int mscsMap::prepare_rotation(double Ath, double Aphi, double Al) {
	long i;
	long pix_num = pixNum();
	msgs->say("PREPARING ROTATION.", High);
	if (!loaded.rotation) makekill_space_manager("make", "rot", 1);
	for (i = 0; i < pix_num; i++) {
		rotation_map[i] = (double) i;
	}
	rotate_map(Ath, Aphi, Al, true, "rot");
	loaded.rotation = true;
	msgs->say("ROTATION LOADED.", Low);
	return 1;
}
//************************************************************************
// performs a quick rotation according to the information stored in the rotation_map
// unrotate=false  - normal rotation
// unrotate=true  - back rotation
// what defines the date to be rotated: currently T and m
int mscsMap::rotate_quick(string what, bool unrotate) {
	long i;
	cpedsList<double> *Xmap, tmp;
	long pix_num = pixNum();
	
	if (loaded.rotation) {
		msgs->say("Quick rotating the map of " + what, Low);
		// set the pointer of data to be rotated
		if (what == "T") Xmap = &map.T;
		if (what == "Q") Xmap = &map.Q;
		if (what == "V") Xmap = &map.V;
		if (what == "U") Xmap = &map.U;
		if (what == "N") Xmap = &map.N;
		if (what == "m") Xmap = &map.m;
		
		tmp.makeLength(pixNum());
		if (unrotate) for (i = 0; i < pix_num; i++) {
			tmp[(long) rotation_map[i]] = (*Xmap)[i];
		}
		else for (i = 0; i < pix_num; i++) {
			tmp[i] = (*Xmap)[(long) rotation_map[i]];
		}
		
		(*Xmap) = tmp;
	} else {
		msgs->error("rotation was not loaded nor prepared, rotation not done.",
				High);
		return -1;
	}
	
	return 1;
}
//************************************************************************
int mscsMap::remove_rotation() {
	makekill_space_manager("kill", "rot", 1);
	return 1;
}
//************************************************************************
// sets the equatorial mask on the map. this cannot be pulled back unless you reload the map.
// b is given in degrees from equator to cut boundary
void mscsMap::make_equatorial_mask(double b) {
	long int i;
	long pix_num = pixNum();
	
	msgs->say(
			"setting map mask: symmetrical equatorial cut below |b| < "
					+ msgs->toStr(b) + " degr.", High);
	makekill_space_manager("make", "m", 1);
	if (coordLoaded()==false) set_map_coord();
	
	b *= PI180;
	for (i = 0; i < pix_num; i++) {
		if (fabs(get_C(i).b()) <= b) {
			m(i) = 0.0;
			mask.masked_pix_num++;
		}
		else {
			m(i) = 1.0;
		}
	}
	loaded.m = true;
	check_mask();
}
//************************************************************************
void mscsMap::invert_mask() {
	long i;
	long pix_num = pixNum();
	for (i = 0; i < pix_num; i++) {
		if (get_m(i) != 0) m(i) = 0.0;
		else m(i) = 1.0;
	}
	mask.masked_pix_num = pixNum() - mask.masked_pix_num;
	mask.f_sky = (double) mask.masked_pix_num / (double) pixNum();
	
}

//************************************************************************
//! Returns a queue object containing the numbers of pixels (in the current regionalization scheme) lyining in on the circle with the center in (l,b) and with radius r.
/*! @param l,b,r - galactic longitude, latitude and radius of the circle in degrees.*/
/*! @param step - parameter defines the accuracy in search for the mathing points. Expressed in fraction of the pixel size. eg. 0.5 gives step of 0.5 pixel size and r>1 will make the circle more sparse.*/
/*! @param strgion_type - the type of the region to be processed: "dot" gives a solid circle and "emptydot" gives just an annuli of thickness region_type_dot_points points */
/*! @param region_type_dot_points - is unused for the "dot" type */

cpeds_queue<long int>* mscsMap::get_circle(double l, double b, double r,
		double step, string region_type, long region_type_dot_points) {
	double pa = fourPI / (double) pixNum(); // pixel area in [rad^2]
	double ma = r * r * PI180 * PI180; // area to be masked [rad^2]
	double ps = sqrt(pa);
	double ms = sqrt(ma);
	double th, phi, dth, dphi, thE, thS, phiE;
	cpedsDirection n;
	// int retval;
	long num;
	
	cpeds_queue<long>* cq = new cpeds_queue<long> ;
	
	l *= PI180;
	b *= PI180;
	
	if (ms <= ps) {
		cpeds_ang2pix_healpix(nside(), &num, PIsnd - b, l, 1);
		cq->addq(num);
	} else {
		dth = dphi = step * ps; // let's make the step step factor smaller than the pixel size of the map
		phiE = twoPI;
		phi = 0;
		thE = PIsnd - ms;
		if (region_type == "dot") thS = 0;
		if (region_type == "emptydot") thS = ms
				- (double) region_type_dot_points * ps;
		th = PIsnd - thS;
		
		while (th >= thE) {
			phi = 0;
			while (phi < phiE) {
				n.lon() = phi;
				n.lat() = th;
				n = RzRy(n, l, PIsnd - b);
				cpeds_ang2pix_healpix(nside(), &num, PIsnd - n.b(), n.l(), 1);
				cq->addq(num);
				
				phi += dphi;
			}
			th -= dth;
		}
		
	}
	return cq;
}
//************************************************************************
// puts a round region on mask or map (depending on what) towards direction l,b (given in galactic coordinates in deg.) of size r [deg]
// with value v
// if the pixel size (4PI/pix_num) is twice smaller than the circle size (r^2), then the region is not masked, otherwise the
// whole pixel is masked.
int mscsMap::make_circle_dot(double l, double b, double r, double v,
		string what, string region_type, long region_type_dot_points,
		string operation) {
	// double pa = fourPI/(double)pixNum(); // pixel area in [rad^2]
	// double ma = r*r*PI180*PI180; // area to be masked [rad^2]
	// double ps = sqrt(pa);
	// double ms = sqrt(ma);
	// double th,phi,dth,dphi,thE,thS,phiE;
	cpedsDirection n;
	int retval;
	long num, i;
	cpedsList<double>* Dptr;
	cpeds_queue<long>* cq;
	
	
	/*   long i,is; // experimentall */
	/*   FILE* f; // experimental */

//	if (what == "m") {
//		makekill_space_manager("make", "m", 1);
//		if (loaded.m == false) clean_mask();
//	}
//	if (what == "T") {
//		makekill_space_manager("make", "T", 1);
//		if (loaded.T == false) clear_map();
//	}

	if (what == "m") {
		Dptr = &map.m;
		msgs->say(
				"making circular mask dot around l [deg] = " + msgs->toStr(l)
						+ " b [deg] = " + msgs->toStr(b) + " size [deg] = "
						+ msgs->toStr(r), Low);
	}
	if (what == "T") {
		Dptr = &map.T;
		msgs->say(
				"T making circular dot around l [deg] = " + msgs->toStr(l)
						+ " b [deg] = " + msgs->toStr(b) + " size [deg] = "
						+ msgs->toStr(r) + " T = " + msgs->toStr(v), Low);
	}
	/* if (what == "Q") { Dptr = &map.Q; msgs->say("Q making circular dot around l [deg] = "+msgs->toStr(l)+" b [deg] = "+msgs->toStr(b)+" size [deg] = "+msgs->toStr(r)+" T = "+msgs->toStr(v),Low); } */
	/* if (what == "U") { Dptr = &map.U; msgs->say("U making circular dot around l [deg] = "+msgs->toStr(l)+" b [deg] = "+msgs->toStr(b)+" size [deg] = "+msgs->toStr(r)+" T = "+msgs->toStr(v),Low); } */
	/* if (what == "V") { Dptr = &map.V; msgs->say("V making circular dot around l [deg] = "+msgs->toStr(l)+" b [deg] = "+msgs->toStr(b)+" size [deg] = "+msgs->toStr(r)+" T = "+msgs->toStr(v),Low); } */

	cq = get_circle(l, b, r, 0.5, region_type, region_type_dot_points);
	num = cq->get_size();
	
	if (operation == "add") {
		for (i = 0; i < num; i++) {
			Dptr[cq->getq(i)] += v;
		}
	} 
	else {
		for (i = 0; i < num; i++) {
			(*Dptr)[cq->getq(i)] = v;
		}
	}

	if (what == "m") {
		loaded.m = true;
	}
	if (what == "T") {
		loaded.T = true;
	}

	delete cq;
	retval = 1;
	
	return retval;
}

//************************************************************************
// makes a ring belt with poles at direction l,b with value v and size s about the
// equator on requested structure: now T or m
// at least nside must be initiated

void mscsMap::mk_ring(double l, double b, double s, double v, string what) {
	mscsMap tmpmap("mk_ring");
	cpedsDirection p;
	double s2;
	long i;
	long pix_num = pixNum();
	
	
	// initiate
	s2 = s * PI180 / 2;
	
	tmpmap.set_nside(nside());
	tmpmap.makekill_space_manager("make", "T", 1);
	tmpmap.set_map_coord(0, 0);
	tmpmap.clear_map();
	
	
	// make ring
	for (i = 0; i < pix_num; i++) {
		p = tmpmap.get_C(i);
		if (fabs(p.b()) <= s2) tmpmap.set_T(i, v);
	}

	// rotate it
	tmpmap.rotate_map(90.0 - b, l, 0.0, true, "T");
	
	
	// copy it onto the final structure
	if (what == "T") {
		if (!mapLoaded()) {
			makekill_space_manager("make", "T", 1);
			clear_map();
		}
		for (i = 0; i < pix_num; i++) {
			if (tmpmap.get_T(i) == v) {
				set_T(i, v);
			}
		}
		loaded.T = true;
	}
	
	if (what == "m") {
		if (!maskLoaded()) {
			makekill_space_manager("make", "m", 1);
			clean_mask();
		}
		for (i = 0; i < pix_num; i++) {
			if (tmpmap.get_T(i) == v) set_m(i, v);
		}
		loaded.m = true;
	}
	
}

//************************************************************************
//produces a mask which divides the sky into l_reg_num*b_reg_num adjacent regions
//with mask value at the pixels in the map numbering a region.
// the regions are defined by dividing latitudes and longitudes of values
// given by equations:
// b_i = +_ i * PI/2 / b_reg_num
// l_j = +_ j * PI   / l_reg_num
// where i = 0,..., b_reg_num, j = 0,.. l_reg_num
// regions ordering is consistent with the ring ordering in healpix pixelization scheme - i.e.
// the first region has number 1 and the numbering increases with increasing longitude
// and increasing latitude starting from the north pole and ending in the south pole.
// The number of regions is constant in all rings.
// the mask values are all grater then one (1) since the range [0,1] is reserved for
// calculations and information about the transparency of the mask the numberig of
// separate regions is the natural number it starts with 1 and ends with l_reg_num*b_reg_num
// with step by one. this is the convention which should be respected for all calculations
// involving the multi-mask stuff.
// the rotation parameters are given in degrees
// rot_short tells how the perform the rotation: if true along the geodesics if not along the axes of CS

void mscsMap::make_multi_maskLB(long l_reg_num, long b_reg_num,
		string mask_info_filename, double Ath, double Aphi, double Al,
		bool rot_short) {
	//void mscsMap::make_multi_maskLB(long l_reg_num, long b_reg_num, string mask_info_filename) {
	long int i, j;
	filenamestr tmpch;
	long total_reg_num = l_reg_num * b_reg_num;
	double th, b, l, ls, bs, reg_num;
	FILE *ff;
	cpedsList<double> M;
	cpedsDirection p_c, p_lmin, p_lmax, p_bmin, p_bmax;
	long pix_num = pixNum();
	
	msgs->say(
			"Setting regional mask: with " + msgs->toStr(l_reg_num)
					+ " longitudal divisions and " + msgs->toStr(b_reg_num)
					+ " lattitudal divisions and " + msgs->toStr(total_reg_num)
					+ " total number of regions", High);
	mask.multi_mask_lreg = l_reg_num;
	mask.multi_mask_breg = b_reg_num;
	mask.multi_mask_reg_num = total_reg_num; // set object multi-mask information variables

	makekill_space_manager("make", "m", 1);
	set_map_coord(0, 0);
	//conv_nest2ring(map);
	ls = 2 * PI / ((double) l_reg_num);
	bs = PI / ((double) b_reg_num);
	
	if (maskLoaded()) {
		M = get_m();
	} // save current mask

	if (l_reg_num == 1 && b_reg_num == 2) { // bypass (reduce) the problem of unequal pixel distribution between the hemispheres
		conv_nest2ring();
		j = pix_num / 2;
		for (i = 0; i < j; i++) {
			m(i) = 1;
		}
		for (i = j; i < pix_num; i++) {
			m(i) = 2;
		}
		conv_ring2nest();
	} else { // usual division of shpere for other cases
		for (i = 0; i < pix_num; i++) {
			l = get_C(i).l();
			th = PIsnd - get_C(i).b();
			reg_num = round(
					floor(th / bs) * (double) l_reg_num + floor(l / ls) + 1.0);
			m(i) = reg_num; // dont overwrite the mask that is already loaded (eg. kp0)
			//printf("i=%li bs=%lf bl=%lf, l=%lf th=%lf ||| floor(th/bs) = %lf floor(l/ls)=%lf, reg_num=%lf\n",i,180/PI*bs,180/PI*ls,180/PI*l,180/PI*th,floor(th/bs),floor(l/ls),reg_num);
		}
	}
	rotate_map(Ath, Aphi, Al, rot_short, "m"); //rotate multi mask if requested

	if (M.count() == pix_num) {
		for (i = 0; i < pix_num; i++) {
			if (m(i) < 1) m(i) = M[i];
		}
	} // apply previous mask (but not multi-mask)
	loaded.m = true;
	
	
	// print/save multi-mask information
	if (mask_info_filename.size() != 0) { // save the mask information to file
		sprintf(tmpch,
				"multi_mask_info-%s-LB%li_%li-Ath%.3lf-Aphi%.3lf-Al%.3lf",
				mask_info_filename.c_str(), l_reg_num, b_reg_num, Ath, Aphi, Al);
		ff = fopen(tmpch, "w");
		//fprintf(f,"%li %li\n",l_reg_num,b_reg_num);
		reg_num = 0;
		b = PIsnd;
		Al *= PI180;
		Ath *= PI180;
		Aphi = (90 - Aphi) * PI180; // convert to radians
		while (fabs(b + PIsnd) > 1e-3) {
			l = 0;
			while (fabs(l - twoPI) > 1e-3) {
				reg_num++;
				p_c.set((2 * l + ls) / 2, (2 * b - bs) / 2); // this formula assumes that there must be as least 2 regions in the longituad divisions or else it won't work correctly; NOTE that: if there's only one division in l then for regions other than polar caps there's no center point.
				if (l_reg_num == 1) {
					if (reg_num == 1) {
						p_c.setLat(PIsnd);
					}
					if (reg_num == total_reg_num) {
						p_c.setLat(-PIsnd);
					}
				}
				p_lmin.set(l, p_c.b());
				p_lmax.set(l + ls, p_c.b());
				p_bmin.set(p_c.l(), b - bs);
				p_bmax.set(p_c.l(), b);
				p_c = Rz(Rx(Rz(Rz(p_c, Al), Aphi), -Ath), -Aphi);
				p_lmin = Rz(Rx(Rz(Rz(p_lmin, Al), Aphi), -Ath), -Aphi);
				p_lmax = Rz(Rx(Rz(Rz(p_lmax, Al), Aphi), -Ath), -Aphi);
				p_bmin = Rz(Rx(Rz(Rz(p_bmin, Al), Aphi), -Ath), -Aphi);
				p_bmax = Rz(Rx(Rz(Rz(p_bmax, Al), Aphi), -Ath), -Aphi);
				fprintf(ff, "%.0lf %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE\n",
						reg_num, 180 / PI * p_c.l(), 180 / PI * p_c.b(),
						180 / PI * p_lmin.l(), 180 / PI * p_lmin.b(),
						180 / PI * p_lmax.l(), 180 / PI * p_lmax.b(),
						180 / PI * p_bmin.l(), 180 / PI * p_bmin.b(),
						180 / PI * p_bmax.l(), 180 / PI * p_bmax.b());
				l += ls;
			}
			b -= bs;
		}
		fclose(ff);
	}
	
}
//************************************************************************
void mscsMap::make_multi_maskHP(long nside_loc, string mask_info_filename,
		double Ath, double Aphi, double Al, bool rot_short) {
	//void mscsMap::make_multi_maskHP(long nside_loc, string mask_info_filename) {
	long i;
	mscsMap newmask("multi_mask_healpix");
	FILE* f;
	filenamestr tmpch;
	cpedsList<double> M;
	cpedsDirection p;
	long pix_num = pixNum();
	
	msgs->say(
			"Setting regional healpix based mask: with nside " + msgs->toStr(
					nside_loc) + " and " + msgs->toStr(
					12 * nside_loc * nside_loc) + " regions", High);
	newmask.set_nside(nside_loc);
	newmask.makekill_space_manager("make", "m", 1);
	newmask.clean_mask();
	newmask.set_map_coord(0, 0); //newmask.conv_nest2ring(newmask.map);
	newmask.conv_nest2ring();
	for (i = 0; i < newmask.pixNum(); i++) {
		newmask.set_m(i, (double) i + 1);
	}

	if (mask_info_filename.size() != 0) { // save the mask information to file
		sprintf(tmpch, "multi_mask_info-%s-HP%li-Ath%.3lf-Aphi%.3lf-Al%.3lf",
				mask_info_filename.c_str(), nside_loc, Ath, Aphi, Al);
		Al *= PI180;
		Ath *= PI180;
		Aphi = (90 - Aphi) * PI180; // convert to radians
		f = fopen(tmpch, "w");
		for (i = 0; i < newmask.pixNum(); i++) {
			p = newmask.get_C(i);
			p = Rz(Rx(Rz(Rz(p, Al), Aphi), -Ath), -Aphi);
			fprintf(f, "%li %lE %lE\n", i + 1, 180 / PI * p.l(),
					180 / PI * p.b());
		}
		fclose(f);
		Al *= PI180inv;
		Ath *= PI180inv;
		Aphi = (PIsnd - Aphi) * PI180inv; // convert to degr.
	}
	newmask.conv_ring2nest();
	newmask.change_map_resolution(nside());
	if (maskLoaded()) {
		M = get_m();
	} // save current mask // carefull here !!! if the multimask is loaded then things become bad.
	import_map_data(newmask, "m", 1); // copy multi mask
	mask.multi_mask_nside = nside_loc;
	mask.multi_mask_reg_num = 12 * nside_loc * nside_loc; // remember that import_map imports also the information on multi_mask_reg_num; if not set will result in error; hence set here
	rotate_map(Ath, Aphi, Al, rot_short, "m"); //rotate multi mask if requested
	if (M.size() > 0) {
		for (i = 0; i < pix_num; i++) {
			if (M[i] < 1) m(i) = M[i];
		}
	} // apply previous mask but not multimask
	loaded.m = true;
	/* newmask.kill_map(); */
}
//************************************************************************
void mscsMap::mask_region(long reg) {
	long i;
	long pix_num = pixNum();
	for (i = 0; i < pix_num; i++) {
		if (get_m(i) == (double) reg) {
			m(i) = 0;
		}
	}
}
//************************************************************************
// cleans the mask from the map - this is the default state - m value = 1.0 - no mask
void mscsMap::clean_mask() {
	// long int i;
	msgs->say("clearing the map mask", Medium);
	makekill_space_manager("make", "m", 1);
	m() = double(1.0);
	// setValue(1.0,"m");
	loaded.m = true;
	mask.masked_pix_num = 0;
	mask.f_sky = 0.0;
	mask.multi_mask_reg_num = 1;
}
//************************************************************************
// this routing scrambles all temperatures that deviate from the mean by more than n sigma
void mscsMap::scramble_over_nsigma(double n) {
	// long i;
	double zmin = mapInfo.meanT - n * sqrt(mapInfo.varianceT);
	double zmax = mapInfo.meanT + n * sqrt(mapInfo.varianceT);
	//calculate_map_stats(0);
	msgs->say(
			"scrambling the temperature map to fit within +_ " + msgs->toStr(n)
					+ " sigma of the map.\n", High);
	
	scramble_over_minmax(zmin, zmax);
	
}
//************************************************************************
// this routing scrambles all temperatures that exceed above and belov max, min respectively
void mscsMap::scramble_over_minmax(double min, double max) {
	long i;
	long pix_num = pixNum();
	
	
	//calculate_map_stats(0);
	msgs->say(
			"scrambling the temperature map below " + msgs->toStr(min)
					+ " and over " + msgs->toStr(max), High);
	
	for (i = 0; i < pix_num; i++) {
		if (get_T(i) < min) {
			T(i) = min;
		}
		if (get_T(i) > max) {
			T(i) = max;
		}
	}
	calculate_map_stats(1);
	
}
//************************************************************************
mscsFunction mscsMap::get_value_sorted_pixel_numbers(long direction) {
	long i;
	mscsFunction f("get_value_sorted_pixel_numbers");
	/* cpeds_point * t = new cpeds_point[pix_num]; */
	/* long *p = new long[pix_num]; */
	double x, y;
	long pix_num = pixNum();
	
	for (i = 0; i < pix_num; i++) {
		if (get_m(i) != 0) x = get_T(i);
		else x = 0;
		y = (double) i;
		f.newPoint(x, y);
	}
	if (direction == 12) f.sortFunctionArgAscending();
	if (direction == 21) f.sortFunctionArgDescending();
	
	return f;
}

//************************************************************************
// checks the status of the map mask
// valid states of the mask_loaded flag:
// 0 - there is no mask - all values on the mask are 1
// 1 - some values of the mask in the map are not 1
// !!!!!!!!!!!! development -- add checking for correct number of regions if multimask is defined

void mscsMap::check_mask() {
	long int i, mini, maxi;
	double min, max;
	long pix_num = pixNum();
	
	msgs->say(
			"mask information: loaded: " + msgs->toStr(loaded.m) + " merged: "
					+ msgs->toStr(mask.merged), High);
	if (get_m().size() > 0) {
		mask.masked_pix_num = 0;
		for (i = 0; i < pix_num; i++) {
			if (m(i) == 0) {
				mask.masked_pix_num++;
			}
		}
		if (mask.masked_pix_num == 0) {
			loaded.m = false;
		}
		mask.f_sky = (double) mask.masked_pix_num / (double) pix_num;
		msgs->say(
				"masked_pix_num: " + msgs->toStr(mask.masked_pix_num)
						+ " unmasked_pix_num: " + msgs->toStr(
						pix_num - mask.masked_pix_num)
						+ " masked sky fraction: " + msgs->toStr(mask.f_sky),
				Medium);
		/* cpeds_find_minmax_value(map->m,pix_num,&min,&max,&mini,&maxi); */
		get_m().getMinMaxValues(&min, &max, &mini, &maxi);
		msgs->say(
				"mask min. value: " + msgs->toStr(min) + " mask max. value: "
						+ msgs->toStr(max), Medium);
		if (max == 1.0) {
			mask.multi_mask_reg_num = 1;
		}
		if (max > 1) {
			mask.multi_mask_reg_num = (long) max;
			msgs->say(
					"NOTE: the maximal mask value > 1. setting the multi_mask_reg_num to (long)max value in the mask:  "
							+ msgs->toStr(mask.multi_mask_reg_num)
							+ ". Hope this is what you want. Note that the number of regions in the multimask does not include the region 0.",
					Medium);
		}
	}
}

//************************************************************************
void mscsMap::clear_multimask() {
	long int i;
	long pix_num = pixNum();
	msgs->say("clearing the map multimask", Medium);
	makekill_space_manager("make", "m", 1);
	for (i = 0; i < pix_num; i++)
		if (get_m(i) > 1) {
			m(i) = 1.0;
		}
	//mask_loaded = 1;
}

//************************************************************************
void mscsMap::mask_map_merge() {
	if (mapLoaded()==false) {
		mask2map();
	}
	else {
		if (maskLoaded()) {
			msgs->say("applying mask to the temperature map", High);
//			printf("%li %li\n",T().size(),m().size());
			T() *= m();
			mask.merged = true;
			check_mask();
			calculate_map_stats(0);
		} else {
			msgs->error("mask is not loaded", Medium);
		}		
	}
}
//************************************************************************
// This merges only the pixels that are not to be considered in analysis
// and skips all other pixels. Otherwise in case of the multimask the
// multiplication might make no sense.
void mscsMap::multimask_map_merge() {
	long int i;
	long pix_num = pixNum();
	if (maskLoaded()) {
		msgs->say("applying multi-mask to the temperature map", Medium);
		for (i = 0; i < pix_num; i++) {
			if (get_m(i) < 1.0) T(i) *= get_m(i);
		}
		mask.merged = true;
		check_mask();
		calculate_map_stats(0);
	} else {
		msgs->say("mask is not loaded and merge requested", Low);
	}
}

//************************************************************************
void mscsMap::mask2map() {
	// long i;
	msgs->say("copying mask to the temperature map", Medium);
	makekill_space_manager("make", "T", 1);
	//makekill_space_manager("make","m",1);
	if (maskLoaded()) {
		T() = get_m();
		loaded.T = true;
	} else {
		msgs->say("mask is not loaded", Low);
		loaded.T = false;
	}
}

//************************************************************************
// this function saves the orientation of the first pixel in the map in case if you want to do some rotations with the map
void mscsMap::save_reference_orientation() {
	/*   reference_dir.l = map[0].n.l; */
	/*   reference_dir.b = map[0].n.b; */
}
//************************************************************************
// rotates the map by a given angle; l and b are given in deg.
//void mscsMap:rotate_map(double l, double b) {

//}
//************************************************************************
void mscsMap::gaussianSmoothT(double beamFWHM, double debeamFWHM, long lmax,
		bool pixtf, long method) {
	if (lmax == -1) lmax = 2 * nside();
	mscsWindowFunction beamTf("beamtf", lmax);
	beamTf.make_gaussian_kernel(beamFWHM);
	mscsWindowFunction debeamTf("beamtf", lmax);
	if (debeamFWHM>0) debeamTf.make_gaussian_kernel(debeamFWHM);
	
	mscsWindowFunction pixTf("pixtf", lmax);
	cpedsStatusCodes sc;
	if (pixtf) {
		pixTf = readPixTf("", &sc);
		if (sc != cpedsSuccess) {
			msgs->warning(
					"Count not load the pixel window function. Will use a flat transfer function.",
					High);
			pixTf.make_unit_kernel();
		}
	} else {
		pixTf.make_unit_kernel();
	}

	mscsAlms a;
	if (debeamFWHM>0) a=SH_analysis(lmax, debeamTf, pixTf, method);
	else a=SH_analysis(lmax);
	
	SH_synthesis(a, lmax, beamTf, pixTf, method);
}

//************************************************************************
void mscsMap::average_map_in_rings(cpeds_queue<double>** qp) {
	long r, Nr, j, i, k, kp;
	cpeds_queue<double> q, *qm = NULL;
	double m;
	if (qp != NULL) {
		qm = new cpeds_queue<double> ;
	}

	msgs->say("averaging the map in rings", High);
	
	conv_nest2ring();
	
	Nr = cpeds_get_ring_num_healpix(nside());
	msgs->say("ring num " + msgs->toStr(Nr), Medium);
	
	k = 0;
	for (r = 0; r < Nr; r++) {
		j = cpeds_get_pixnum_in_ring_healpix(nside(), r);
		
		kp = k;
		if (maskLoaded()) {
			for (i = 0; i < j; i++) {
				if (get_m(k) == 1) {
					q.addq(get_T(k));
				}
				k++;
			}
		}
		else {
			for (i = 0; i < j; i++) {					q.addq(get_T(k));				k++;			}
		}
		
		m = q.mean();
		if (qp != NULL) qm->addq(m);
		k = kp;

		if (maskLoaded()) {
			for (i = 0; i < j; i++) {
				if (get_m(k) == 1) {
					set_T(k, m);
				}
				k++;
			}			
		}
		else {
			for (i = 0; i < j; i++) {	set_T(k, m);				k++;			}
		}
		q.delete_all_queue();
	}

	conv_ring2nest();
	if (qp != NULL) {
		(*qp) = qm;
	}
	
}
//************************************************************************
void mscsMap::calculate_map_stats(int output) {
	long int i, j = 0;
	cpedsList<double> tmp;
	cpedsList<long> idx;
	
	double th1, phi1, th2, phi2;
	// double var2;
	long pix_num = pixNum();
	
	if (maskLoaded()) { // we calculate the statistics ONLY on the unmasked regions
	//	  printf("!!!!!!!!!DEBUG:   the mask is still loaded \n");
		for (i = 0; i < pix_num; i++) {
			if (get_m(i) != 0) {
				tmp.append(get_T(i));
				idx.append(i);
			}
		} //exit(0);
		tmp.getMinMaxValues(&mapInfo.minT, &mapInfo.maxT, &mapInfo.iminT,
				&mapInfo.imaxT);
		mapInfo.iminT = idx[mapInfo.iminT];
		mapInfo.imaxT = idx[mapInfo.imaxT];
		get_moments_of_distribution(tmp, &mapInfo.meanT, &mapInfo.varianceT,
				&mapInfo.skewnessT, &mapInfo.kurtosisT);
		mask.masked_pix_num = pix_num - idx.size();
		mask.f_sky = (double) mask.masked_pix_num / (double) pix_num;
		
		
		//meanT = cpeds_mean_value(tmp2,j); varianceT = cpeds_variance(tmp2,j);   skewnessT = cpeds_skewness(tmp2,j); kurtosisT = cpeds_kurtosis(tmp2,j);
		// delete tmp2; delete idx2; delete idx;
	} else {
		get_T().getMinMaxValues(&mapInfo.minT, &mapInfo.maxT, &mapInfo.iminT,
				&mapInfo.imaxT);
		
		get_moments_of_distribution(get_T(), &mapInfo.meanT,
				&mapInfo.varianceT, &mapInfo.skewnessT, &mapInfo.kurtosisT);
		
	}
	if (output != 0) {
		msgs->say(
				"map stats: nside: " + msgs->toStr(nside())
						+ " masked_pix_num: "
						+ msgs->toStr(mask.masked_pix_num)
						+ " mask_loaded is: " + msgs->toStr(loaded.m)
						+ " merged: " + msgs->toStr(mask.merged), High);
		msgs->say(
				"appriximated pixel size [deg]: " + msgs->toStr(
						PI180inv * cpeds_pix_size_healpix(nside())), Medium);
		msgs->say(
				"minT: " + msgs->toStr(mapInfo.minT) + " maxT: " + msgs->toStr(
						mapInfo.maxT) + " imin: " + msgs->toStr(mapInfo.iminT)
						+ " imax: " + msgs->toStr(mapInfo.imaxT), Medium);
		msgs->say(
				"mean: " + msgs->toStr(mapInfo.meanT)
						+ " sqrt(variance) (sigma): " + msgs->toStr(
						sqrt(mapInfo.varianceT)) + " skewness: " + msgs->toStr(
						mapInfo.skewnessT) + " kurtosis: " + msgs->toStr(
						mapInfo.kurtosisT), Medium);
		cpeds_pix2ang_healpix(nside(), mapInfo.iminT, &th1, &phi1, 1); // CAUTION !!! this is for nested and NO check is performed about the map's ordering
		msgs->say(
				"Tmin: l = " + msgs->toStr(180 / PI * phi1) + ", b = "
						+ msgs->toStr(180 / PI * (PIsnd - th1)) + " iminT = "
						+ msgs->toStr(mapInfo.iminT), Medium);
		cpeds_pix2ang_healpix(nside(), mapInfo.imaxT, &th2, &phi2, 1);
		msgs->say(
				"Tmax: l = " + msgs->toStr(180 / PI * phi2) + ", b = "
						+ msgs->toStr(180 / PI * (PIsnd - th2)) + " iminT = "
						+ msgs->toStr(mapInfo.imaxT), Medium);
		msgs->say(
				"Tmin - Tmax relative angle: " + msgs->toStr(
						180 / PI * cpeds_ang_n1n2(th1, phi1, th2, phi2)),
				Medium);
	}
	
	if ((coordLoaded()) && (mapInfo.ordering == nested)) {
		if (output != 0) {
		}
	}
	/* delete [] tmp;   */
	//exit(0);
}

//************************************************************************
// this routine calculates the m,s,S,K statistics on a map but in different
// regions separatelly defined by multimask (created eg. by make_multi_mask routine)
// the address of array with statistics is returned

stat_info * mscsMap::calculate_map_stats_with_multi_mask(int output) {
	long i, reg, reg_tot = mask.multi_mask_reg_num;
	stat_info *reg_stat = new stat_info[reg_tot];
	double **t = new double*[reg_tot];
	long * mm_pix_num = new long[reg_tot];
	long * itab = new long[reg_tot];
	double mean_loc, variance_loc, skewness_loc, kurtosis_loc;
	long pix_num = pixNum();
	cpedsList<double> cl;
	
	
	/* temporary commented out */
	/*   if (output > 0) { printf("|%s> * calculating the map statistis in separate regions on sky: reg_num = %li \n",object_name.c_str(),reg_tot); } */
	/* temporary commented out */

	mm_pix_num = count_multi_mask_pix_num();
	/*   for (i=0;i<reg_tot;i++) { mm_pix_num[i] = 0; } // initialize mm_pix_num array */
	/*   for (i=0;i<pix_num;i++) { if (map->m[i] >= 1.0) { mm_pix_num[(long)(map->m[i])-1] += 1; }} // count points in different regions but outside the mask */

	/* for (j=0;j<64;j++) { printf("%li ",mm_pix_num[j]); } printf(" mask: %lE\n",map->m[i]); */
	/*     if (map->m[i] >= 64) printf("dupa\n"); */
	/*     if (map->m[i] < 0) printf("DUPA\n"); */

	for (i = 0; i < reg_tot; i++) {
		itab[i] = 0;
		if (mm_pix_num[i] > 0) {
			t[i] = new double[mm_pix_num[i]];
		} else t[i] = NULL;
	} // initiate space for individual maps in regions

	for (i = 0; i < pix_num; i++) { // copy the map information from regions into separate tables
		if (get_m(i) >= 1.0) {
			reg = (long) (get_m(i)) - 1;
			(t[reg])[itab[reg]] = get_T(i);
			itab[reg]++;
		}
	}

	/*   for (i=0;i<reg_tot;i++) {  // calculate statistics in regions */
	/*     reg_stat[i].n = mm_pix_num[i]; */
	/*     if (reg_stat[i].n > 0) reg_stat[i].m = cpeds_mean_value(t[i],mm_pix_num[i]); else reg_stat[i].m = 0; */
	/*     if (reg_stat[i].n > 0) reg_stat[i].s = sqrt(cpeds_variance(t[i],mm_pix_num[i])); else reg_stat[i].s = 0; */
	/*     if (reg_stat[i].n > 1) reg_stat[i].S = cpeds_skewness(t[i],mm_pix_num[i]); else reg_stat[i].S = 0; */
	/*     if (reg_stat[i].n > 1) reg_stat[i].K = cpeds_kurtosis(t[i],mm_pix_num[i]); else reg_stat[i].K = 0; */
	/*   } */

	/* FROM OPTIMALIZATION REASOSNS THE DISTRIBUTIONS MOMENTS ROUTINES ARE MOVED TO MSCS AND OPTIMIZED FOR THE PURPOSE OF THE SKREGSTAT PROJECT */
	for (i = 0; i < reg_tot; i++) { // calculate statistics in regions
		cl.clear();
		cl.fromCarray(t[i], mm_pix_num[i]);
		get_moments_of_distribution(cl, &mean_loc, &variance_loc,
				&skewness_loc, &kurtosis_loc);
		/*     mean_loc = 1; variance_loc=1; skewness_loc=1;kurtosis_loc=1; */
		reg_stat[i].n = mm_pix_num[i];
		reg_stat[i].m = mean_loc;
		reg_stat[i].s = sqrt(variance_loc);
		reg_stat[i].S = skewness_loc;
		reg_stat[i].K = kurtosis_loc;
	}
	/* FROM OPTIMALIZATION REASOSNS THE DISTRIBUTIONS MOMENTS ROUTINES ARE MOVED TO MSCS AND OPTIMIZED FOR THE PURPOSE OF THE SKREGSTAT PROJECT */

	/* FROM OPTIMALIZATION REASOSNS */
	/* temporary commented out */
	/*   if (output > 0) { */
	/*     for (i=0;i<reg_tot;i++) { */
	/*       printf("|%s>  -- region: %li n = %.0lf, m = %lE, s = %lE, S = %lE, K = %lE\n",object_name.c_str(),i+1,reg_stat[i].n,reg_stat[i].m,reg_stat[i].s,reg_stat[i].S,reg_stat[i].K); */
	/*     } */
	/*   } */
	/* temporary commented out */
	/* FROM OPTIMALIZATION REASOSNS */

	//free memory
	for (i = 0; i < reg_tot; i++) {
		if (t[i] != NULL) delete[] t[i];
	}
	delete[] t;
	delete[] mm_pix_num;
	delete[] itab;
	return reg_stat;
}

//************************************************************************
// this routine calculates the variance statistics on a map but in different
// regions separatelly defined by multimask (created eg. by make_multi_mask routine)
// the address of array with statistics is returned
// The fields of Skewness and Kurtosis are used to store the total variance of the whole unmasked map
// and the ratio of the local in-region variance to the total variance in the unmasked map.
stat_info * mscsMap::calculate_map_variance_with_multi_mask(int output) {
	long i, j, reg, reg_tot = mask.multi_mask_reg_num;
	stat_info *reg_stat = new stat_info[reg_tot];
	double **t = new double*[reg_tot];
	long * mm_pix_num = new long[reg_tot];
	long * itab = new long[reg_tot];
	long unmasked_pix_num_loc = 0;
	double *map_loc = NULL;
	double variance_loc;
	long pix_num = pixNum();
	
	mm_pix_num = count_multi_mask_pix_num();
	
	for (i = 0; i < reg_tot; i++) {
		itab[i] = 0;
		if (mm_pix_num[i] > 0) {
			t[i] = new double[mm_pix_num[i]];
			unmasked_pix_num_loc += mm_pix_num[i];
		} else t[i] = NULL;
	} // initiate space for individual maps in regions

	// copy the map information from regions into separate tables
	j = 0;
	map_loc = new double[unmasked_pix_num_loc];
	for (i = 0; i < pix_num; i++) {
		if (get_m(i) >= 1.0) {
			reg = (long) (get_m(i)) - 1;
			(t[reg])[itab[reg]] = get_T(i);
			itab[reg]++;
			map_loc[j] = get_T(i);
			j++;
		}
	}

	variance_loc = sqrt(cpeds_variance(map_loc, unmasked_pix_num_loc));
	
	for (i = 0; i < reg_tot; i++) { // calculate statistics in regions
		reg_stat[i].n = mm_pix_num[i];
		reg_stat[i].m = 1;
		reg_stat[i].s = sqrt(cpeds_variance(t[i], mm_pix_num[i]));
		reg_stat[i].S = variance_loc;
		reg_stat[i].K = reg_stat[i].s / variance_loc;
	}

	//free memory
	for (i = 0; i < reg_tot; i++) {
		if (t[i] != NULL) delete[] t[i];
	}
	delete[] t;
	delete[] mm_pix_num;
	delete[] itab;
	delete[] map_loc;
	return reg_stat;
}

//************************************************************************
// counts number of pixels in regions defined by multi-mask
// mask must be loaded
// allocates the space and returns the array of size multi_mask_reg_num with number of pixels in regions
// the i'th index in the array corresponds to the region number i+1
// this does not calculate the pixels with values < 1 in the mask
long * mscsMap::count_multi_mask_pix_num() {
	long * mm_pix_num = new long[mask.multi_mask_reg_num];
	long i;
	long pix_num = pixNum();
	
	for (i = 0; i < mask.multi_mask_reg_num; i++) {
		mm_pix_num[i] = 0;
	} // initialize mm_pix_num array
	for (i = 0; i < pix_num; i++) {
		if (get_m(i) >= 1.0) {
			mm_pix_num[(long) (get_m(i)) - 1]++;
		}
	} // count points in different regions but outside the mask
	return mm_pix_num;
}

//************************************************************************
/* method for calculation the angular momentum of the map: ie. M = sum_i T^2_i * cos(|b_i|)^2 */
/* in the current orientation of the map */

double mscsMap::calculate_map_momentum() {
	long i;
	double M = 0, c, Nl = 0, Tl;
	long pix_num = pixNum();
	
	if (!coordLoaded()) set_map_coord(0, 0);
	if (!maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			c = cos(get_C(i).b());
			c *= c;
			Tl = get_T(i);
			Tl *= Tl;
			M += Tl * c;
			Nl += Tl;
		}
	} else {
		for (i = 0; i < pix_num; i++)
			if (map.m[i] != 0) {
				c = cos(get_C(i).b());
				c *= c;
				Tl = get_T(i);
				Tl *= Tl;
				M += Tl * c;
				Nl += Tl;
			}
	}
	/*   for (i=0;i<pix_num;i++) { c=cos(abs(map->n[i].b)); c*=c; T=map->T[i]; T=abs(T); M+=T*c; N+=T;} */

	return M / Nl;
	
	
	/*   conv_nest2ring(map.T); */

}
//************************************************************************
cpedsDirection mscsMap::maximize_map_momentum(long nsst, long nsen,
		double *MMAX, double *acc, bool maximize) {
	mscsMap *copy = new mscsMap("maximize map momentum");
	long dir_num, ns, nsmap;
	long i, j, imax, ist, ien;
	double l, b, th, lmax, bmax, thmax, M = 0, Mmax;
	cpedsList<double> backup;
	long pix_num = pixNum();
	bool rot_short = true; // define how the rotations are to be performed using rotate_map method
	// this must be true since otherwise he may not hit the right rotation to get the right orientation etc.

	// clone the object
	copy->clone(*this);
	
	
	/*   nsen = nside; */
	ns = nsst;
	if (maximize) Mmax = 0;
	else Mmax = 100;
	
	msgs->say("maximizing the angular momentum of the map", High);
	
	do {
		ns *= 4; //if (ns > nsen) ns=nsen;

		//set the map resolution for search in dir_num directions in northern sky
		if (ns * 4 >= copy->nside()) nsmap = copy->nside();
		else nsmap = 4 * ns;
		change_map_resolution(copy->nside());
		import_map_data(*copy, "T", 1);
		change_map_resolution(nsmap);
		backup = get_T();
		
		
		// define the directions in which we look for the maximal angular momentum axis
		if (ns == 4 * nsst) {
			dir_num = 12 * ns * ns;
			ist = 0;
			ien = dir_num;
		} else {
			ist = imax << 4;
			ien = ist + 16;
			// update the Mmax in the map in the new resolution
			rotate_map(PI180inv * thmax, PI180inv * (PI + lmax), 0, rot_short,
					"T");
			Mmax = calculate_map_momentum();
			for (j = 0; j < pix_num; j++) {
				T(j) = backup[j];
			}
		}
		
		
		// define a set of directions for first search over the northern sky
		for (i = ist; i < ien; i++) {
			cpeds_pix2ang_healpix(ns, i, &th, &l, 1);
			b = PIsnd - th; //dir[i].l=l; dir[i].b=PIsnd-th;
			/*       if (b>0) { */
			rotate_map(180 / PI * th, 180 / PI * (PI + l), 0, rot_short, "T");
			
			
			// perform the search and find the maximal value
			M = calculate_map_momentum();
			if (maximize) {
				if (M > Mmax) {
					Mmax = M;
					lmax = l;
					bmax = b;
					thmax = th;
					imax = i;
				}
			} else { // minimize
				if (M < Mmax) {
					Mmax = M;
					lmax = l;
					bmax = b;
					thmax = th;
					imax = i;
				}
			}

			// restore the old map orientation
			T() = backup;
			char ctmp[1000];
			sprintf(
					ctmp,
					"ns=%li nsmap=%li, i=%li, l=%lf b=%lf, M=%lE lmax =%lf, bmax= %lf, Mmax=%lE\n",
					ns, nsmap, i, l * 180 / PI, b * 180 / PI, M,
					lmax * 180 / PI, bmax * 180 / PI, Mmax);
			msgs->say(ctmp, Low);
			/*       } // b>0 */
		}
		/*     printf("koniec petli: ns=%li nsen=%li\n",ns,nsen); */
	} while (4 * ns < nsen);

	delete copy;
	
	*acc = cpeds_pix_size_healpix(nsen);
	*MMAX = Mmax;
	return cpedsDirection(lmax, bmax);
}

//************************************************************************
/* routine finds and returns the orientation of the HOT TIP of the dipole compenent in the map in the real space*/
/*  acc - parameter that defines the accuracy of the fit given in the raito of variances of the map- it defines the step in which the amplitude of the fitted dipole
 components are changed in each of the coordinates */
/* range - parameter that defines the upper range of the search over amplitudes in variances, the search will be done from -range to range */
/* the routine returns also the amplitude of the fitted dipole in the temperature units of the input map */

cpedsDirection mscsMap::fit_map_dipole(double range, double acc, double *AMP,
		bool remove_dipole) {
	mscsMap *d = new mscsMap("dipole");
	mscsAlms a("dipole alms");
	cpedsList<double> backup;
	cpedsList<double> backupD;
	double A, dxyz[3], Xmin, X, DX;
	double stdevM, stdevMsqrt3;
	long j;
	double l, b;
	cpedsDirection n;
	
	msgs->say("fitting dipole component of the map", High);
	
	
	// prepare the variance calibrated temperature map
	calculate_map_stats(0);
	stdevM = sqrt(mapInfo.varianceT);
	stdevMsqrt3 = stdevM * sqrt(3);
	backup = T();
	TdivVal(stdevM);
	check_mask();
	
	d->set_nside(nside());
	a.lmax(1);
	
	for (j = 0; j < 3; j++) {
		// make the dipole map for X axsis (l,b)=(0,0) (for the cold tip)
		if (j == 0) {
			a.set(0, 0, 0, 0);
			a.set(1, 0, 0, 0);
			a.set(1, 1, 0, 1);
		}
		// make the dipole map for Y axis (l,b) = (90,0) (for the cold tip)
		if (j == 1) {
			a.set(0, 0, 0, 0);
			a.set(1, 0, 0, 0);
			a.set(1, 1, -1, 0);
		}
		// make the dipole map for vertical orientation
		if (j == 2) {
			a.set(0, 0, 0, 0);
			a.set(1, 0, 1, 0);
			a.set(1, 1, 0, 0);
		}

		a.antysymmetrize();
		/*     printf("C_1 = %lE\n",d->calculate_C_l(1)); */
		d->SH_synthesis(a, 1);
		
		d->calculate_map_stats(1);
		d->norm_by_stddev(); // should here be norm_by_variance() ? - TODO: please check this bit
		// stdev=sqrt(d->get_varianceT);
		// d->TdivVal(stdev);
		backupD = d->get_T();
		
		A = -range;
		/* for (i=0;i<=pix_num;i++) { d->map->T[i]=backupD[i]*A; } */
		d->T() = backupD * A;
		X = calculate_map_to_map_chisq(d);
		Xmin = X;
		DX = 0.0;
		while (DX <= 0) {
			// set the dipole amplitude
			A += acc;
			/* for (i=0;i<=pix_num;i++) { d->map->T[i]=backupD[i]*A; } */
			d->T() = backupD * A;
			
			
			// calculate chisq between the map and dipole
			X = calculate_map_to_map_chisq(d);
			printf("A: %lE, range: %lE DX=%lE\n", A, range, DX);
			
			DX = X - Xmin;
			if (X <= Xmin) {
				Xmin = X;
				dxyz[j] = A;
			}
			
			
			/*       printf("A=%lE X= %lE Xmin= %lE DX=%lE dx=%lE\n",A,X,Xmin,DX,dx); */
		}
		
		
		/*     if (remove_dipole) { */
		/*       for (i=0;i<=pix_num;i++) { d->map->T[i]=backupD[i]*dxyz[j]; } */
		/*       for (i=0;i<=pix_num;i++) { backup[i]-=d->map->T[i]*stdevMsqrt3; } */
		/*     } */

	}

	printf("dx %lE, dy %lE, dz %lE\n", dxyz[0], dxyz[1], dxyz[2]);
	(*AMP) = sqrt(dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2])
			* stdevM * sqrt(3);
	
	b = cpeds_cart2sph(0, dxyz[0], dxyz[1], dxyz[2]); // this is theta
	l = PIsnd - cpeds_cart2sph(1, dxyz[0], dxyz[1], dxyz[2]);
	n.set(l, b);
	n.print_direction();
	cpeds_check_thphi(&b, &l);
	n.lat() = PIsnd - n.lat(); //TODO check the sense of these operations
	n.set(l, b); //TODO check the sense of these operations
	n.print_direction();
	
	if (remove_dipole) {
		a.set(0, 0, 0, 0);
		a.set(1, 0, 1, 0);
		a.set(1, 1, 0, 0);
		a.antysymmetrize();
		d->SH_synthesis(a, 1);
		d->rotate_map(PI180inv * (PIsnd - n.b()), PI180inv * n.l(), 0.0, true,
				"T");
		d->norm_by_maxT();
		/* for (i=0;i<=pix_num;i++) { backup[i]-=d->map->T[i]*(*AMP); } */
		backup -= d->get_T() * (*AMP);
	}

	// set back the original temperature map
	T() = backup;
	
	delete d;
	
	return n;
}

//************************************************************************
// this calculates X = 1/Npix Sum_pix (T1_i-T2_i)^2
// no check for the size consistency between the maps is done
// only the mask in this object is used
double mscsMap::calculate_map_to_map_chisq(mscsMap* map) {
	long i;
	double tmp1, tmp2, X = 0;
	long pix_num = pixNum();
	
	for (i = 0; i <= pix_num; i++) {
		if (get_m(i) != 0) {
			tmp1 = get_T(i);
			tmp2 = map->get_T(i);
			X += (tmp1 - tmp2) * (tmp1 - tmp2);
		}
	}
	return X / (double) (pix_num - mask.masked_pix_num);
}
//************************************************************************
double mscsMap::calculate_minT() {
	double min, max;
	long imin, imax;
	get_nonMasked("T").getMinMaxValues(&min, &max, &imin, &imax);
	return min;
}

//************************************************************************
double mscsMap::calculate_maxT() {
	double min, max;
	long imin, imax;
	get_nonMasked("T").getMinMaxValues(&min, &max, &imin, &imax);
	return max;
}
//************************************************************************
double mscsMap::calculate_meanT() {
	return get_nonMasked("T").mean();
}
//************************************************************************
double mscsMap::calculate_varianceT() {
	return get_nonMasked("T").variance();
}
//************************************************************************
double mscsMap::calculate_rmsT() {
	return get_nonMasked("T").rms();
}

//************************************************************************
double mscsMap::get_integralT() {
	return get_nonMasked("T").sum() * fourPI / (double) pixNum();
}

//************************************************************************
// the resolution parameter defines how wide is the binning of theta range, thus defines the number of points that constitute C_th
// angles are given in degrees.
// calculates the correlation function on a map with requested resolution
mscsCorrelationFunction mscsMap::calculate_C_th(double theta_min,
		double theta_max, double resolution) {
	theta_min *= PI180;
	theta_max *= PI180;
	resolution *= PI180;
	
	long int i, j, corr_i;
	double ang;
//	long point_num_C_th = (int) (ceil((theta_max - theta_min) / resolution));
	
//	cpedsList<double> separation_number;
	mscsCorrelationFunction Cth;
	
//	Cth.setPointsNum(point_num_C_th);
//	separation_number.makeLength(point_num_C_th);
	
	msgs->say(
			"calculating C_th from theta=" + msgs->toStr(theta_min)
					+ ", to theta= " + msgs->toStr(theta_max), High);
	if (!coordLoaded()) set_map_coord();
	/* for (i=0;i<point_num_C_th;i++) { C_th[i][0] = C_th[i][1] = 0; separation_number[i] = 0;} // zeroing tables */

	if (maskLoaded()==false) {
		makekill_space_manager("make","m");
		m()=1;
	}
	
	
	
/*
	// corralation function calculation
	long pix_num = pixNum();
	for (i = 0; i < pix_num; i++) {
		if (get_m(i)!=0) {
			for (j = i; j < pix_num; j++) {
				if (get_m(j) != 0) {
					
					ang = get_C(i).angle(get_C(j)); //cpeds_ang_n1n2(map->n[i].b,map->n[i].l,map->n[j].b,map->n[j].l);
					
					if ((ang >= theta_min) && (ang <= theta_max)) {
						corr_i = (long) round((ang - theta_min) / resolution);
						separation_number[corr_i]++; // this stores the number of given separations on a sky ( for normalization purposes)
						Cth[corr_i].rx() += ang;
						Cth[corr_i].ry() += get_T(i) * get_T(j); // the mask check is done in the condition above
					}
				}
			}
		}
		
		printf("calculating: %li of %li\r", i, pix_num);
	}

	// normalization and averaging

	for (i = 0; i < point_num_C_th; i++) {
		Cth[i].rx() /= separation_number[i]; // th  -- this gives the average angle over all angles that fall into this range (bin) limited by the resolution parameter
		Cth[i].ry() /= separation_number[i]; // normalization of C(th)
	}
*/

	
	
	cpedsDirectionSet ds;
	long pix_num = pixNum();
	for (i = 0; i < pix_num; i++) {
		if (get_m(i)!=0) {
			ds.append(get_C(i));
			ds.last().setVal(get_T(i));
		}
	}
	Cth=cpeds_calculate_angular_correlation_fn(ds,theta_min, theta_max, resolution);
	
	return Cth;
}

/* ******************************************************************************************** */
mscsCorrelationFunction mscsMap::calculate_Sth(double theta_min, double theta_max, double resolution) {
	theta_min *= PI180;
	theta_max *= PI180;
	resolution *= PI180;
	
	long int i, j, corr_i;
	double ang;
//	long point_num_C_th = (int) (ceil((theta_max - theta_min) / resolution));
	
//	cpedsList<double> separation_number;
	mscsCorrelationFunction Cth;
	
//	Cth.setPointsNum(point_num_C_th);
//	separation_number.makeLength(point_num_C_th);
	
	msgs->say(
			"calculating C_th from theta=" + msgs->toStr(theta_min)
					+ ", to theta= " + msgs->toStr(theta_max), High);
	if (!coordLoaded()) set_map_coord();
	/* for (i=0;i<point_num_C_th;i++) { C_th[i][0] = C_th[i][1] = 0; separation_number[i] = 0;} // zeroing tables */

	if (maskLoaded()==false) {
		makekill_space_manager("make","m");
		m()=1;
	}
	
	
	cpedsDirectionSet ds;
	long pix_num = pixNum();
	for (i = 0; i < pix_num; i++) {
		if (get_m(i)!=0) {
			ds.append(get_C(i));
			ds.last().setVal(get_T(i));
		}
	}
	Cth=cpeds_calculate_angular_correlation_fn(ds,theta_min, theta_max, resolution);
	
	return Cth;
	
}

//************************************************************************


/* //\************************************************************************ */
/* // this routine derives the full sky angular power spectrum estimator from: */
/* // GIVEN: */
/* // lmax - must be set to indicate the maximum l for C_l reconstruction */
/* // pseudo C_l of the observations loaded on C_l structure */
/* // Blsq = Bl^2*pl^2 - filtering applied to the data due to beams, pixel tf. smoothing etc. -- to be multiplied with M_ll' */
/* // pseudoN - pseudo power spectrum of the noise derived from MC simulations for the experiment (the map doesn't have to be loaded) */
/* // pseudoW - pseudo power spectrum of the sky cut derived from the sky mask (doesn't have to be loaded) */
/* // (optionally) M_ll^-1 - address to the inverse matrix of the mode to mode coupling; if size is 1x1 then it will be calculated from pseudoW */
/* // the power spectra above are not to be flattened (no l(l+1)/2PI factor) */
/* // C - matrix with the simulated (cross)power spectra for calculation of the covariance matrix and error analysis */
/* //  */
/* // the routine  */
/* // RETURNS: */
/* // the address to the estimated power spectrum of the full sky  */
/* // the address to the full covariance matrix of the estimated power spectrum */
/* // the address to the estimated full sky noise power spectrum */

/* power_spectrum * mscsMap::extract_C_l(long lmax_loc, power_spectrum * Blsq, power_spectrum * pseudoN, power_spectrum * pseudoW, matrix <double>* M, matrix <double>* C, power_spectrum * Nfs, string Mllinfo) { */
/*   long l1,l2,l3,  k1,k2,k3,  Nsim; */
/*   double W3J, tmpd,d1,d2,d3,*D,*DD; */
/*   matrix <double> *N, *Mi, *Cl, *Bl; */
/*   power_spectrum *Cfs; */
/*   long stat; */
/*   int i; */
/*   long k,l,m; */
/*   gsl_sf_result res; */
/*   double * W3jt; */
/*   bool noMll;  if (M->RowNo()==1) noMll = true; else  noMll = false;  */
/*   filenamestr tmpch; */

/*   printf("|%s> * extracting C_l from pseudo C^_l and mask window function\n",object_name.c_str()); */
/* /\*   if (map_loaded == 0) { printf("|%s>  -- ERROR: The map is not loaded. Load the map file first \n",object_name.c_str()); return; } *\/ */
/* /\*   if (mask_loaded == 0) { printf("|%s>  -- WARNING: The mask is not loaded. Will assume no mask at all \n",object_name.c_str()); return; } *\/ */
/* /\*   printf("|%s>  -- computing power spectrum of the window function\n",object_name.c_str()); *\/ */
/* /\*   set_alms_lmax(lmax_loc); *\/ */
/* /\*   pseudoCl = new power_spectrum(lmax);   for (l=0;l<=lmax;l++) { (*pseudoCl).set_Cl(l,C_l[l][1]); } // this is to be changed when map class will be updated for more objects *\/ */
/* /\*   calculate_transformF(default_fourier_method,1,"nosmooth",-1,0,"",""); *\/ */
/* /\*   calculate_C_l(0,lmax,1); // this is the power spectrum of the window function *\/ */

/*   // here, check the sanity of the input data */

/*   // ****************************************** */
/*   // declare some usefull matrices */
/*   // ****************************************** */

/*   //if (noMll) { Mi = new matrix <double>;   (*Mi).SetSize((lmax_loc+1),(lmax_loc+1)); } */
/*   Mi = new matrix <double>;   (*Mi).SetSize((lmax_loc+1),(lmax_loc+1)); */
/*   N  = new matrix <double>;   (*N).SetSize((lmax_loc+1),1);     for (l1=0;l1<=lmax_loc;l1++) { (*N)(l1,0) =pseudoN->get_Cl(l1); } */
/*   Cl = new matrix <double>;   (*Cl).SetSize((lmax_loc+1),1);    for (l1=0;l1<=lmax_loc;l1++) { (*Cl)(l1,0)=C_l[l1][1]; } */
/*   Bl = new matrix <double>;   (*Bl).SetSize((lmax_loc+1),1);    for (l1=0;l1<=lmax_loc;l1++) { (*Bl)(l1,0)=Blsq->get_Cl(l1); } */

/*   // ****************************************** */
/*   // deriving  the mode-to-mode coupling matrix */
/*   // ****************************************** */

/*   if (noMll) {  */
/*     printf("|%s>  -- computing mode-to-mode coupling matrix\n",object_name.c_str()); */

/*     (*M).SetSize((lmax_loc+1),(lmax_loc+1)); */
/*     W3jt =  new double[2*(lmax_loc+1)]; */

/*     W3J=0;     */
/*     for (l1=0;l1<=lmax_loc;l1++) { */
/* /\*       k1=2*l1; *\/ */
/*       for (l2=0;l2<=lmax_loc;l2++) { */
/* /\* 	k2=2*l2; *\/ */
/* 	(*M)(l1,l2) = 0;  */
/* 	tmpd = 0; */

/* 	d1=0; d2=(double)lmax_loc; k2=2*(lmax_loc+1); */
/* 	ThreeJSymbolJ((double)l1,(double)l2,(double)0,(double)0,d1,d2,W3jt,k2,i); */
/* 	//if (i != 0) cout << " dupa --- cos nie tak\n"; */

/* /\* 	for (l3=0;l3<=lmax_loc;l3++) { printf("%li %lE,  ",l3,W3jt[l3]); } *\/ */
/* 	k1=0; */
/* 	for (l3=0;l3<=lmax_loc;l3++) { */
/* 	  //gsl implementation */
/* /\* 	  k3=2*l3; *\/ */
/* /\* 	  W3J=gsl_sf_coupling_3j(l1,l2,l3,0,0,0); W3J=W3J*W3J; *\/ */
/* /\* 	  stat = gsl_sf_coupling_3j_e(k1,k2,k3,0,0,0,&res); W3J=0; *\/ */
/* /\* 	  if (stat == 0) { W3J=res.val; } W3J=W3J*W3J; *\/ */

/* 	  // cpeds + lidia implementation -- this is way too slow + there are some problem with header files between lidia and matpack */
/* /\* 	  if (l1 < 400) W3J=cpeds_W3jm0(l1,l2,l3); else W3J=cpeds_W3jm0_lidia(l1,l2,l3);  *\/ */

/* 	  // cpeds implementation */
/* /\* 	  W3J=cpeds_W3jm0(l1,l2,l3); *\/ */
/* /\* 	  W3J=W3J*W3J; *\/ */

/* 	  // matpack implementation */
/* /\* 	  if (cpeds_W3jm0_non_zero(l1,l2,l3)) {W3J=W3jt[k1]*W3jt[k1]; k1++; } *\/ */
/* 	  W3J=0; */
/* 	  if (cpeds_W3jm0_lcond(l1,l2,l3,0)) {W3J=W3jt[k1]*W3jt[k1]; k1++; } */

/* /\* 	  W3J=W3jt[l3]*W3jt[l3]; *\/ */
/* /\* 	  if (cpeds_W3jm0_non_zero) W3J++; *\/ */

/* // test */
/* /\* 	  if (cpeds_W3jm0_non_zero(l1,l2,l3)) { *\/ */
/* /\* 	    d1=cpeds_W3jm0(l1,l2,l3);  *\/ */
/* /\* 	    d2=cpeds_W3jm0_lidia(l1,l2,l3); *\/ */
/* /\* 	    d3=W3jt[k1-1]; *\/ */
/* /\* 	    if (abs(d1-d2)>1e-5 || abs(d1-d3)>1e-5 || abs(d2-d3)>1e-5) printf("----- roznica: cpeds: %lE lidia: %lE: matpack: %lE,  l1 l2 l3 %li %li %li\n",d1,d2,d3,l1,l2,l3); *\/ */
/* /\* 	    //if (abs(d1-d3)>1e-5 ) printf("----- roznica: cpeds: %lE lidia: %lE: matpack: %lE,  l1 l2 l3 %li %li %li\n",d1,d2,d3,l1,l2,l3); *\/ */
/* /\* 	  } *\/ */

/* // test */

/* 	  tmpd+=(double)(2*l3+1)*pseudoW->get_Cl(l3)*W3J; */

/* /\* 	  printf("--------> l1: %li l2: %li, l3: %li    tmpch: %lE W3J^2: %lE\n",l1,l2,l3,tmpd, W3J); *\/ */
/* /\* 	if (l1 == 186 && l2 == 186) printf("--------> l1: %li l2: %li, l3: %li    tmpch: %lE W3J^2: %lE\n",l1,l2,l3,tmpd, W3J); *\/ */
/* /\* 	if (l1 == 400 && l2 >= 286) printf("--------> l1: %li l2: %li, l3: %li    tmpch: %lE W3J^2: %lE\n",l1,l2,l3,tmpd, W3J); *\/ */
/* /\* 	if (l1 >= 250 && l2 >= 540) printf("--------> l1: %li l2: %li, l3: %li    tmpch: %lE W3J^2: %lE\n",l1,l2,l3,tmpd, W3J); *\/ */

/* 	} */
/* 	tmpd*=(2*(double)l2+1)/(4*PI);	(*M)(l1,l2)=tmpd; */
/* /\* 	printf("l1=%li l2=%li\r",l1,l2); *\/ */
/*       } */
/*       printf("l1=%li\r",l1); */
/*     } */

/*     delete [] W3jt; */
/*     printf("|%s>  -- DONE\n",object_name.c_str());  */
/*     sprintf(tmpch,"%s--lmax_%li-Mll.mat",Mllinfo.c_str(),lmax_loc);   printf("|%s>  -- saving window kernel matrix to temporary file: %s\n",object_name.c_str(),tmpch); */
/*     Mscs_matrix_save(M,tmpch); */

/*     printf("|%s>  -- inverting mode to mode coupling matrix...\n",object_name.c_str()); */
/*     *Mi = !(*M);  */

/*     sprintf(tmpch,"%s--lmax_%li-invMll.mat",Mllinfo.c_str(),lmax_loc); printf("|%s>  -- saving inverse window kernel matrix to temporary file: %s (you probably don't need this)\n",object_name.c_str(),tmpch); */
/*     Mscs_matrix_save(Mi,tmpch); */
/*     printf("|%s>  -- DONE\n",object_name.c_str());     */
/*   } */
/*   //else { (*Mi)=(*M); } */

/*   // ****************************************** */
/*   // deriving the full inverse coupling matrix */
/*   // ****************************************** */
/*   printf("|%s>  -- deriving the full coupling matrix\n",object_name.c_str()); */
/*   for (l1=0;l1<=lmax_loc;l1++) {  */
/*     for (l2=0;l2<=lmax_loc;l2++) { (*Mi)(l2,l1) = (*M)(l2,l1)*(*Bl)(l2,0); }} */

/* /\*   sprintf(tmpch,"%s--lmax_%li-Kll.mat",Mllinfo.c_str(),lmax_loc); printf("|%s>  -- saving kernel matrix to temporary file: %s\n",object_name.c_str(),tmpch);  *\/ */
/* /\*   Mscs_matrix_save(Mi,tmpch); *\/ */
/*   printf("|%s>  -- inverting the full coupling matrix\n",object_name.c_str()); */
/*   (*Mi)=!(*Mi); */
/* /\*   sprintf(tmpch,"%s--lmax_%li-invKll.mat",Mllinfo.c_str(),lmax_loc);   printf("|%s>  -- saving kernel matrix to temporary file: %s\n",object_name.c_str(),tmpch);  *\/ */
/* /\*   Mscs_matrix_save(Mi,tmpch); *\/ */

/*   // ****************************************** */
/*   // deriving the full sky C_l estimator */
/*   // ****************************************** */

/*   printf("|%s>  -- deriving the full sky power spectrum estimator\n",object_name.c_str()); */
/*   Cfs = new power_spectrum(lmax_loc+1);  */
/*   *Cl-=*N;  // subtract estimated pseudo noise */
/*   (*Cl)=(*Mi)*(*Cl); */

/*   // copy onto Clfs */
/*   for (l1=0;l1<=lmax_loc;l1++) { Cfs->set_Cl(l1, (*Cl)(l1,0) );   Cfs->set_l(l1,(double)l1); } */

/*   // ****************************************** */
/*   // setting out the inverse mode-mode coupling matrix for output */
/*   // ****************************************** */
/*   (*M)=(*Mi); */

/*   // ****************************************** */
/*   // deriving the full sky noise estimator */
/*   // ****************************************** */
/*   printf("|%s>  -- deriving the full sky noise power spectrum\n",object_name.c_str()); */
/*   (*N)=(*Mi)*(*N);   */
/*   // copy onto Nfs */
/*   for (l1=0;l1<=lmax_loc;l1++) { Nfs->set_Cl(l1, (*N)(l1,0) );   Nfs->set_l(l1,(double)l1); } */

/*   // ****************************************** */
/*   // deriving the estimated C_l covariance matrix  */
/*   // ****************************************** */
/*   printf("|%s>  -- deriving the covariance matrix of the full sky power spectrum\n",object_name.c_str()); */
/*   Nsim = (long)(C->ColNo());  */
/*   printf("|%s>  -- data matrix has %li cols and %li rows\n",object_name.c_str(),(long)C->ColNo(),(long)C->RowNo()); */
/*   D = new double[Nsim*(lmax_loc+1)]; */
/*   m=0; */
/*   for (k=0;k<Nsim;k++) { */
/*     for (l=0;l<=lmax_loc;l++) { D[m]=(*C)(l,k); m++;}} */

/*   DD = cpeds_calculate_covariance_matrix(D,lmax_loc+1,Nsim); */

/*   // copy the cov matrix onto the C matrix structure */
/*   (*C).SetSize(lmax_loc+1,lmax_loc+1); */
/*   m=0; */
/*   for (k=0;k<=lmax_loc;k++) { */
/*     for (l=0;l<=lmax_loc;l++) { (*C)(l,k)=DD[m]; m++; }} */
/*   printf("|%s>  -- data matrix has now %li cols and %li rows\n",object_name.c_str(),(long)C->ColNo(),(long)C->RowNo()); */

/*   // freeing space */

/*   if (noMll) delete Mi; */
/*   delete N; delete Cl; delete Bl; delete [] D; delete [] DD;  */

/*   return Cfs; */
/* } */

//************************************************************************
void mscsMap::shift_mean_to(double val, bool calc_stats) {
	// long i;
	if (calc_stats) calculate_map_stats(0);
	T() += -mapInfo.meanT + val;
}

/* **************************************************************************************************** */
void mscsMap::norm_by_maxT() {
	long i;
	long pix_num = pixNum();
	double max;
	
	msgs->say("normalizing the temperature map by maximal value", Medium);
	max = calculate_maxT();
	
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) T(i) /= max;
		}
	} else {
		T() /= max;
	}
	
	calculate_map_stats(1);
}

/* **************************************************************************************************** */
void mscsMap::norm_by_mean() { // divides the map by it's mean
	long i;
	long pix_num = pixNum();
	double mean;
	
	msgs->say("normalizing the temperature map by mean", Medium);
	mean = calculate_meanT();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) T(i) /= mean;
		}
	} else {
		T() /= mean;
	}
	
	calculate_map_stats(1);
}

/* **************************************************************************************************** */
void mscsMap::norm_by_variance() { // divides the map by it's variance
	long i;
	long pix_num = pixNum();
	double var;
	//  double var=varianceT; // commented out on 2009-01-15

	msgs->say("normalizing the temperature map by variance", Medium);
	var = calculate_varianceT();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) T(i) /= var;
		}
	} else {
		T() /= var;
	}
	
	calculate_map_stats(1);
}

/* **************************************************************************************************** */
void mscsMap::norm_by_stddev() { // divides the map by it's standard deviation
	long i;
	long pix_num = pixNum();
	double stdev;
	
	msgs->say("normalizing the temperature map by standard deviation", Medium);
	stdev = sqrt(calculate_varianceT());
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) T(i) /= stdev;
		}
	} else {
		T() /= stdev;
	}
	calculate_map_stats(1);
}

/* **************************************************************************************************** */
void mscsMap::norm_by_rms() { // divides the map by it's rms
	long i;
	long pix_num = pixNum();
	double rms = calculate_rmsT();
	
	msgs->say("normalizing the temperature map by rms", Medium);
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) T(i) /= rms;
		}
	} else {
		T() /= rms;
	}
	calculate_map_stats(1);
}

/* **************************************************************************************************** */

bool mscsMap::isRing() {
	if (mapInfo.ordering == ring) return true;
	return false;
}
bool mscsMap::isNested() {
	if (mapInfo.ordering == nested) return true;
	return false;
}
/* bool mscsMap::isHealpix() { if (pix_system==1) return true; else return false; } */

/* **************************************************************************************************** */
void mscsMap::mk_abs() {
	long i;
	long pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) {
				T(i) = fabs(get_T(i));
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			T(i) = fabs(T(i));
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::mk_sqrtAbs() {
	long i;
	long pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) {
				T(i) = sqrt(fabs(get_T(i)));
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			T(i) = sqrt(fabs(get_T(i)));
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::invertT() {
	long i;
	long pix_num = pixNum();
	
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) != 0.0 && !isMasked(i)) T(i) = 1.0 / get_T(i);
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) != 0.0) T(i) = 1.0 / get_T(i);
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::logarithmT(double base) {
	long i;
	long pix_num = pixNum();
	double logBase = log(base);
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) != 0.0 && !isMasked(i)) T(i) = log(fabs(T(i)))
					/ logBase;
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			if (get_T(i) != 0.0) T(i) = log(fabs(get_T(i))) / logBase;
		}
	}
}

/* **************************************************************************************************** */
void mscsMap::powerT(double exponent) {
	long i;
	long pix_num = pixNum();
	if (maskLoaded()) {
		for (i = 0; i < pix_num; i++) {
			if (!isMasked(i)) {
				T(i) = pow(get_T(i), exponent);
			}
		}
	} else {
		for (i = 0; i < pix_num; i++) {
			T(i) = pow(get_T(i), exponent);
		}
	}
}

//************************************************************************
// fills the map coordinates within the actuall coordinates pixelizations scheme and sets the coord_loaded status
// it's done with regard to the current settings of nside and pix_system


void mscsMap::set_map_coord(double db, double dl) {
	long int i;
	double th, phi;
	string coord_file;
	makekill_space_manager("make", "C", 1);
	
	
	/* msgs->say("setting actual coordinates for the map for the pixelization system: %i in ordering: %i\n",object_name.c_str(),pix_system,map_ordering); */
	msgs->say(
			"setting actual coordinates for the map for the pixelization system",
			High);
	/* if (pix_system == 1) { */

	if (isNested()) {
		coord_file = MSCS_DATA_DIR + MSCS_GLOBAL__NESTED_COORDINATES_PREF
				+ msgs->toStr(nside()) + MSCS_GLOBAL__NESTED_COORDINATES_SUFF;
	}
	if (isRing()) {
		coord_file = MSCS_DATA_DIR + MSCS_GLOBAL__RING_COORDINATES_PREF
				+ msgs->toStr(nside()) + MSCS_GLOBAL__RING_COORDINATES_SUFF;
	}
	long ret = loadbinC(coord_file);
	if (ret != 0) { // this works only for the nest ordering for now --- RING IS TO BE IMPLEMENTED !!!
		msgs->say(
				"could not open the coordinates file: " + coord_file
						+ ", will calculate the coords (this is slow)", Medium);
		n().setLength(coordNum());
		for (i = 0; i < coordNum(); i++) {
			cpeds_pix2ang_healpix(nside(), i, &th, &phi,
					(long) mapInfo.ordering);
			n(i).set(phi, (PIsnd - th));
		}
		savebinC(coord_file, "");
	}
	if (db != 0 or dl != 0) msgs->warning(
			"Rotating coordinates is currently not supported. Ignoring supplied arguments.",
			High);
	
	
	/* if (pix_system == 2) {  } */
	loaded.n = true;
}
//************************************************************************
void mscsMap::make_gaussian_map(double m, double s, int method, cpedsRNG *rns) {
	bool rnsWasNull;
	msgs->say(
			"generating random gaussian map with mean: " + msgs->toStr(m)
					+ " and variance: " + msgs->toStr(s), High);
	if (rns == NULL) {
		rns = new cpedsRNG("gaussian", "double");
		rnsWasNull = true;
	} else {
		rnsWasNull = false;
	}
	rns->setMeanVariance(m, s);
	rns->setCentralLimitNumbers(100);
	cpedsList<double> cl;
	T() = cl = rns->getRNs(pixNum());
	if (rnsWasNull) delete rns;
}

//************************************************************************
void mscsMap::generate_uniform_vector(double min, double max, int method,
		cpedsRNG *rns) {
	bool rnsWasNull;
	if (rns == NULL) {
		rns = new cpedsRNG("uniform", "double");
		rnsWasNull = true;
	} else {
		rnsWasNull = false;
	}
	rns->setMinMax(min, max);
	cpedsList<double> cl;
	T() = cl = rns->getRNs(pixNum());
	if (rnsWasNull) delete rns;
}
//************************************************************************
/* double mscsMap::calculate_chisq(double * tab1, double * tab2, long int num) { */
/*   long int i; */
/*   double Xsq = 0; */

/*   for (i=0;i<num;i++) { Xsq = Xsq + pow((tab1[i]-tab2[i]),2); } */
/*   return Xsq; */
/* } */

//************************************************************************
// internal rotation methods
// angles given in radians since these are for internal use and in Mscs CS
cpedsDirection mscsMap::Rx(const cpedsDirection& p, double sinAx, double cosAx) {
	cpedsDirection pp;
	double x, y, z, xx, yy, zz, th, phi;
	
	th = PIsnd - p.b();
	phi = p.l(); // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0, th, p.l());
	y = cpeds_sph2cart(1, th, p.l());
	z = cpeds_sph2cart(2, th, p.l());
	
	xx = x;
	yy = y * cosAx - z * sinAx;
	zz = y * sinAx + z * cosAx;
	
	th = cpeds_cart2sph(0, xx, yy, zz);
	phi = cpeds_cart2sph(1, xx, yy, zz);
	cpeds_check_thphi(&th, &phi);
	pp.setLat(PIsnd - th);
	pp.setLon(phi); // change from (th,phi) CS to (l,b) CS
	return pp;
}
//************************************************************************
cpedsDirection mscsMap::Rx(const cpedsDirection& p, double Ax) {
	double sinAx = sin(Ax), cosAx = cos(Ax);
	return Rx(p, sinAx, cosAx);
}

//************************************************************************
cpedsDirection mscsMap::Ry(const cpedsDirection& p, double sinAy, double cosAy) {
	cpedsDirection pp;
	double x, y, z, xx, yy, zz, th, phi;
	
	th = PIsnd - p.b();
	phi = p.l(); // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0, th, p.l());
	y = cpeds_sph2cart(1, th, p.l());
	z = cpeds_sph2cart(2, th, p.l());
	
	xx = x * cosAy + z * sinAy;
	yy = y;
	zz = -x * sinAy + z * cosAy;
	
	th = cpeds_cart2sph(0, xx, yy, zz);
	phi = cpeds_cart2sph(1, xx, yy, zz);
	cpeds_check_thphi(&th, &phi);
	pp.setLat(PIsnd - th);
	pp.setLon(phi); // change from (th,phi) CS to (l,b) CS
	return pp;
}
//************************************************************************
cpedsDirection mscsMap::Ry(const cpedsDirection& p, double Ay) {
	double sinAy = sin(Ay), cosAy = cos(Ay);
	return Ry(p, sinAy, cosAy);
}

//************************************************************************
cpedsDirection mscsMap::Rz(const cpedsDirection& p, double sinAz, double cosAz) {
	cpedsDirection pp;
	double x, y, z, xx, yy, zz, th, phi;
	
	th = PIsnd - p.b();
	phi = p.l(); // change from (l,b) CS to (th,phi) CS
	x = cpeds_sph2cart(0, th, p.l());
	y = cpeds_sph2cart(1, th, p.l());
	z = cpeds_sph2cart(2, th, p.l());
	
	xx = x * cosAz - y * sinAz;
	yy = x * sinAz + y * cosAz;
	zz = z;
	
	th = cpeds_cart2sph(0, xx, yy, zz);
	phi = cpeds_cart2sph(1, xx, yy, zz);
	cpeds_check_thphi(&th, &phi);
	pp.setLat(PIsnd - th);
	pp.setLon(phi); // change from (th,phi) CS to (l,b) CS
	return pp;
}
//************************************************************************
cpedsDirection mscsMap::Rz(const cpedsDirection& p, double Az) {
	double sinAz = sin(Az), cosAz = cos(Az);
	return Rz(p, sinAz, cosAz);
}
//************************************************************************


/* direction mscsMap::Rz(direction p,double Az) { */
/*   direction pp; */
/*   pp.l = (p.l+Az);  pp.b = p.b; */
/*   cpeds_check_phi(&(pp.l)); */
/*   return pp; */
/* } */

//************************************************************************
cpedsDirection mscsMap::RzRy(const cpedsDirection& p, double Az, double Ay) {
	cpedsDirection pp;
	pp = Ry(p, Ay);
	pp = Rz(pp, Az);
	return pp;
}
//************************************************************************


//************************************************************************
// rotates directions in gal. CS
int mscsMap::rotate_map_directions(cpedsDirectionSet& X, int axis, double ang) {
	long i;
	long pix_num = pixNum();
	double cosang = cos(ang);
	double sinang = sin(ang);
	if (axis == 1) for (i = 0; i < pix_num; i++)
		X[i] = Rx(X[i], sinang, cosang); // rotate around x
	if (axis == 2) for (i = 0; i < pix_num; i++)
		X[i] = Ry(X[i], sinang, cosang); // rotate around y
	if (axis == 3) for (i = 0; i < pix_num; i++)
		X[i] = Rz(X[i], sinang, cosang); // rotate around z

	return 0;
}
//************************************************************************

// optimized math routines for SKregstat project
void mscsMap::get_moments_of_distribution(const cpedsList<double>& tab,
		double* mean, double * variance, double * skewness, double * kurtosis) {
	long size = tab.size();
	double mean_loc = 0, variance_loc = 0, skewness_loc = 0, kurtosis_loc = 0,
			sized = (double) size, sizedleo;
	double *t = new double[size];
	double *t2 = new double[size];
	long i;
	
	if (size > 1) {
		sizedleo = sized - 1.0;
		//mean
		/* for (i=0;i<size;i++) mean_loc+=tab[i]; mean_loc/=sized; *mean = mean_loc; */
		*mean = mean_loc = tab.mean();
		
		
		// shift tab to do everything faster since we do the central moments;
		for (i = 0; i < size; i++)
			t[i] = t2[i] = tab[i] - mean_loc;
		
		
		//variance
		for (i = 0; i < size; i++) {
			t2[i] *= t[i];
			variance_loc += t2[i];
		}
		variance_loc /= sizedleo;
		*variance = variance_loc;
		
		
		//skewness
		for (i = 0; i < size; i++) {
			t2[i] *= t[i];
			skewness_loc += t2[i];
		}
		skewness_loc = skewness_loc / (sized * variance_loc
				* sqrt(variance_loc));
		*skewness = skewness_loc;
		
		
		// kurtosis
		for (i = 0; i < size; i++) {
			t2[i] *= t[i];
			kurtosis_loc += t2[i];
		}
		kurtosis_loc = kurtosis_loc / (sized * variance_loc * variance_loc);
		*kurtosis = kurtosis_loc;
	} else {
		if (size == 1) *mean = tab[0];
		else *mean = 0;
		*skewness = 0;
		*kurtosis = 0;
		*variance = 0;
	}
	delete[] t;
	delete[] t2;
}

/********************************************************************************/
void mscsMap::set_nside(long n) {
	if ((mapLoaded()) && (nside() != n)) {
		msgs->warning(
				"the map with nside=" + msgs->toStr(nside())
						+ " is already loaded", High);
	}
	mapInfo.nside = n;
	mapInfo.pix_num = 12 * n * n;
	mapInfo.coord_num = pixNum();
	msgs->say(
			"Set map nside to: " + msgs->toStr(nside()) + " and pix_num to: "
					+ msgs->toStr(pixNum()) + " and coord_num to: "
					+ msgs->toStr(coordNum()), High);
}

//************************************************************************
void mscsMap::set_pixNum(long n) {
	mapInfo.pix_num = n;
	double ns = sqrt(double(n) / 12.0);
	if (n % 12 != 0) msgs->warning(
			"set_pixNum: the requested number of pixels doesn't convert to a valid nside value ("
					+ msgs->toStr(ns) + ").", High);
	if (sqrt(n / 12) != long(ns)) msgs->warning(
			"set_pixNum: the requested number of pixels doesn't convert to a valid nside value ("
					+ msgs->toStr(ns) + ").", High);
	mapInfo.nside = long(ns);
	msgs->say(
			"Setting map nside to: " + msgs->toStr(nside())
					+ " and pix_num to: " + msgs->toStr(pixNum()), High);
}
//************************************************************************
void mscsMap::set_ordering(mapOrderings ord) {
	if (ord == nested) msgs->say("setting ordering to NESTED", Medium);
	if (ord == ring) msgs->say("setting ordering to RING", Medium);
	mapInfo.ordering = ord;
}

//************************************************************************
void mscsMap::copy_map(cpedsList<double> *from, cpedsList<double> *to,
		cpedsList<long> *conv_tab) {
	long i;
	long N;
	// make sure to have enough space in the destination
	*to = *from;
	// convert using the conversion array
	if (conv_tab != NULL) {
		N = from->size();
		for (i = 0; i < N; i++) {
			(*to)[i] = from->at((*conv_tab)[i]);
		}
	}
}
//************************************************************************
// makes a copy of initiated fields from map_from to map_to of the same size pix_num with the conversion table conv_tab
void mscsMap::copy_map(mscsMap *map_from, mscsMap *map_to,
		cpedsList<long> *conv_tab) {
	long i;
	long N;
	// make sure to hava enough space in the destination map
//	*map_to = *map_from; // this assignment is too much, it 
						 // screws up the map flags and other internal data so don't do that
	// convert using the conversion array
	if (conv_tab != NULL) {
		N = map_from->T().size();
		for (i = 0; i < N; i++) {
			map_to->T(i) = map_from->get_T(conv_tab->value(i));
		}
		N = map_from->Q().size();
		for (i = 0; i < N; i++) {
			map_to->Q(i) = map_from->get_Q(conv_tab->value(i));
		}
		N = map_from->V().size();
		for (i = 0; i < N; i++) {
			map_to->V(i) = map_from->get_V(conv_tab->value(i));
		}
		N = map_from->U().size();
		for (i = 0; i < N; i++) {
			map_to->U(i) = map_from->get_U(conv_tab->value(i));
		}
		N = map_from->N().size();
		for (i = 0; i < N; i++) {
			map_to->N(i) = map_from->get_N(conv_tab->value(i));
		}
		N = map_from->m().size();
		for (i = 0; i < N; i++) {
			map_to->m(i) = map_from->get_m(conv_tab->value(i));
		}
		N = map_from->n().size();
		for (i = 0; i < N; i++) {
			map_to->n(i) = map_from->get_C(conv_tab->value(i));
		}
	}
}

/********************************************************************************/
mscsWindowFunction mscsMap::readPixTf(string fn, cpedsStatusCodes *sc) {
	mscsWindowFunction f;
	if (fn == "") fn = MSCS_DATA_DIR + MSCS_GLOBAL__HEALPIX_PIXTF_PREF
			+ msgs->toStr(nside());
	msgs->say("Loading pixel transfer function from file: " + fn, High);
	
	*sc = f.load(fn);
	if (*sc != cpedsSuccess) {
		msgs->warning("Cannot load the pixel window function: " + fn + ".",
				High);
		return f;
	}
	msgs->say("Success", Low);
	return f;
}

/********************************************************************************/
mscsFunction mscsMap::getRing(long r) {
	long i, j, n;
	mscsFunction f;
	n = cpeds_get_pix_num_above_ring_healpix(nside(), r);
	j = n + cpeds_get_pixnum_in_ring_healpix(nside(), r);
	
	for (i = n; i < j; i++) {
		// if (get_m(i) == 1)
		{ //printf("pix in ring ---> %li\n",i);
			f.newPoint(get_C(i).l(), get_T(i));
		}
	}
	//  printf("\n");
	return f;
}
/********************************************************************************/
void mscsMap::draw_circles(const matrix<double>& lbrv, string dstMap,
		string overplot_region_type, long overplot_region_type_dot_points,
		string operation) {
	long N = lbrv.RowNo();
	for (long k = 0; k < N; k++) {
		make_circle_dot(lbrv(k, 0), lbrv(k, 1), lbrv(k, 2), lbrv(k, 3), "T",
				overplot_region_type, overplot_region_type_dot_points,
				operation);
	}
}
/***************************************************************************************/
void mscsMap::makeDipole(cpedsDirection m) {
	if (loaded.T==false) makekill_space_manager("make","T");
	if (loaded.n==false) set_map_coord();
	double amp=m.val();
	m*=PI180;
	for (long i = 0; i < pixNum(); i++) {
		T(i)=amp*cos(cpeds_ang_n1n2(m.get_direction(),n(i).get_direction()));
	}
	
	
}

void mscsMap::setValue(double x, string what) {
	msgs->criticalError("mscsMap::setValue: This method has not been implemented",Top); exit(-1);
}
/********************************************************************************/
// definition of operators
/********************************************************************************/
/* bool mscsMap::operator==(pixel* side) { */
/*   if (map->T == side->T && map->n == side->n && map->m == side->m && map->N == side->N) return true; else return false; */

//************************************************************************
// FLAT MAP CONFIGURING FUNCTION
//************************************************************************

/* void mscsMap::calculate_flat_map_parameters(int printinfo) { */
/*   long int i; */
/* /\*   gsl_vector *x = gsl_vector_alloc(flat_coord_num); *\/ */
/* /\*   gsl_vector *y = gsl_vector_alloc(flat_coord_num); *\/ */
/*   double *tabx = new double[flat_coord_num]; */
/*   double *taby = new double[flat_coord_num]; */
/*   long int iminx_flat,imaxx_flat,iminy_flat,imaxy_flat; */

/* /\*   for (i=0;i<flat_coord_num;i++) { *\/ */
/* /\*     gsl_vector_set(x, i, flatcoord[i].x); *\/ */
/* /\*     gsl_vector_set(y, i, flatcoord[i].y); *\/ */
/* /\*   } *\/ */

/*   for (i=0;i<flat_coord_num;i++) { */
/*     tabx[i] = flatcoord[i].x;  */
/*     taby[i] = flatcoord[i].y;  */
/*   } */

/* /\*   } *\/ */

/*   cpeds_find_minmax_value(tabx,flat_coord_num,&minx_flat,&maxx_flat,&iminx_flat,&imaxx_flat); */
/*   cpeds_find_minmax_value(taby,flat_coord_num,&miny_flat,&maxy_flat,&iminy_flat,&imaxy_flat); */

/* /\*   minx_flat = gsl_vector_min(x);   maxx_flat = gsl_vector_max(x); *\/ */
/* /\*   miny_flat = gsl_vector_min(y);   maxy_flat = gsl_vector_max(y); *\/ */
/*   x2yratio_flat = (maxx_flat- minx_flat)/(maxy_flat - miny_flat); */
/*   if (printinfo != 0) { */
/*     printf("|%s>  -- minx_flat = %lf, maxx_flat = %lf, miny_flat = %lf, maxy_flat = %lf,  x/y ratio = %lf\n",object_name.c_str(),minx_flat,maxx_flat,miny_flat,maxy_flat,x2yratio_flat); */
/*     printf("|%s>  -- minT = %lf, maxT = %lf \n",object_name.c_str(), minT, maxT); */
/*   } */
/*   //  delete x; delete y; */
/*   delete tabx; delete taby; */
/* } */
//************************************************************************
// this returns the delta C_l/C_l uncertainty related to the cosmic variance
/* double mscsMap::calculate_cosmic_variance(double l) { */
/*   return sqrt(2/(l*(l+1))); */
/* } */
//************************************************************************
// fills the C_l array with the default values for cosmic variance according to l
/* void mscsMap::set_cosmic_variance() { */
/*   long i; */
/*   for (i=0;i<=lmax;i++) { C_l[i][2] = calculate_cosmic_variance(i); } */
/* } */
//************************************************************************
// THESE METHODS WERE BLOCKED DUE TO CONFLICTS IN MULTIPLY INCLUDED HEADERS -- problem to be solved in future maybe.
/* void mscsMap::initiate_dodecahedron() { */
/*   //dodec = new dodecahedron(0,0,this); */
/* } */

// calculates the circle statistics for dodecahedron on a map
// return the pointer for array of double of size 7. 6 circles + 1 statistics for all of them together.

/* double * mscsMap::calculate_dodecahedron_circle_stistic(double l, double b, double g, double a, double s) { */
/*   if (dodec == NULL) { dodec = new dodecahedron(l,b,g,a,s,this); }  */
/*   return dodec.Sstat(l,b,g,a,s); */
/* } */
//************************************************************************
// this is dedicated to check the gaussianity of the map or alms by simple
// chisq test
//commented out during transition into version-1.0
/* void mscsMap::test_gaussianity_chisq(strarg what){ */
/*   long int k,i,j=0,bin_num,mini,maxi; */
/*   double * t; */
/*   double x=1,Px,min,max,bin, bin_num_real; */
/*   printf("testing sample for gausianity with chisq test\n"); */
/*   if (strcmp(what,"alms") == 0) { k = 2*alms_num; } */
/*   if (strcmp(what,"map") == 0) { k = pix_num;  } */
/*   t = new double[k]; */

/*   if (strcmp(what,"alms") == 0) { */
/*     for (i=0;i<k;i++) { t[i] = alm[j].R; i++; t[i] = alm[j].I; j++; } */
/*   } */
/*   if (strcmp(what,"map") == 0) { */
/*     for (i=0;i<k;i++) { t[i] = map->T[i]; } */
/*   } */

/*   cpeds_find_minmax_value(t,k,&min,&max,&mini,&maxi);  */
/*   bin_num = 100; bin = (max-min)/(double)bin_num; */

/*   //x = cpeds_test_chisq(k,t,1,0,1,2,bin,&bin_num_real,NULL); */
/*   //Px = cpeds_chisq_prob(x); */
/*   Px = gsl_ran_chisq_pdf(x,(double)k); */
/*   printf("sample size: %li, reduced chisq: %lE, P(%li,%lE) = %lE\n",k,x,k,x,Px); */
/*   delete t; */
/* } */
// calculates the cross power spectrum between the alms in this object and those given in a2 pbject,
// the cross power is stroreg in the C_l structure of this object.

/* void mscsMap::calculate_cross_C_l(mscsMap & a2, int how) {  // calculates the cross power spectrum on F - i.e. from the alms */
/*   long l,m,lm,lmax_loc; */
/*   double temp = 0,ld; */

/*   lmax_loc = cpeds_get_min(lmax,a2.lmax); */
/*   if (lmax != a2.lmax) { printf("WARNING: the size of alms do not match, will calculate upto: %li",lmax_loc); } */

/*   if ((how == 1) || (how == 2)) { // calculate from alms */
/*     printf("  -- calculating cross C_l from l=%li, to l= %li\n",(long)0,lmax_loc); */
/*     for (l=0;l<=lmax_loc;l++) { */
/*       temp = 0; */
/*       for (m=-l;m<=l;m++) { lm=alm2num(l,m); temp += alm[lm].R*a2.alm[lm].R+alm[lm].I*a2.alm[lm].I; } */
/*       ld = (double)l; */
/*       temp /= (2*ld + 1); // 2l + 1 normalization */
/*       C_l[l][0] = ld; */
/*       C_l[l][1] = temp; */
/*       //printf("l=%i Cl=%lE clcheck1=%lE,   clcheck2=%lE\n",l,temp, calculate_single_C_l(l,0), calculate_single_C_l(l,0,0)); */
/*       if (how == 2) { C_l[l][1] = ld*(ld+1)/(2*PI)*C_l[l][1]; } */
/*     } */
/*     lmax_C_l = lmax_loc; // this assumes that the alms always start from l=0 */
/*     lmin_C_l = 0; */
/*     l_num_C_l = (lmax_C_l-lmin_C_l)+1; */
/*   } */

/* } */

/* //\************************************************************************ */
/* // binns the mscsMap power  spectrum structure and returns the binned power spectrum object. */
/* // binning is only from lmin to lmax, outside this range C_l is copied without change; */
/* // binning is done according to the bintab table is size bintabs that holds the long numbers defining  */
/* // how many ls to bin using weights w=1/bintab[i]  */
/* // or weights from 3rd col of the C_l structure for each l in which case w should be 0: w=0; otherwise it should be w=1 */
/* // the sum of values stored in bintab array should be lmax-lmin+1 */
/* // value in array bintabs - 1 corresponds to no binning; 2 - means that two multipoles are binned and so on. */

/* power_spectrum* mscsMap::bin_C_l(long lmin_loc, long lmax_loc, long *bintabs, long* bintab, double w) { */
/*   long i,j,k,binl_sum=0; */
/*   long Cls, Cbs,cls,cbs; */
/*   power_spectrum* cb; */
/*   long tmp; */
/*   double* effl; */
/*   bool EQweights; */
/*   double wsum; */

/*   printf("|%s> * Binning the power spectrum from %li to %li, together: %li multipoles\n",object_name.c_str(),lmin_loc,lmax_loc,lmax_loc-lmin_loc+1); */
/*   printf("  -- bins are:\n"); */

/*   // check safety conditions */
/*   lmax_loc = cpeds_get_min(lmax_loc,lmax); */
/*   lmin_loc = cpeds_get_max(lmin_loc,0); */
/*   if (*bintabs <= 0 || bintab == NULL) { printf("ERROR: wrong aruments\n"); /\* cb = new power_spectrum(-1); *\/ return NULL; } */
/*   for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; printf("%li ",bintab[i]); } printf("\n  -- The overall sum of the multipoles to bin as implied by binning array is: %li \n",binl_sum);  */
/*   if (binl_sum < lmax_loc-lmin_loc+1) { printf("|%s>  -- WARNING !! the binning array is incompatible with requested binning multipole range. Will bin upto the end of the binning array: i.e. upto l=%li.\n",object_name.c_str(),lmin_loc + binl_sum - 1); } */
/*   if (binl_sum < lmax_loc-lmin_loc+1) { lmax_loc = lmin_loc + binl_sum - 1; } */
/*   if (binl_sum >  lmax_loc-lmin_loc+1) { printf("|%s>  -- WARNING !! the binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning array: i.e. upto l==",object_name.c_str()); } */
/*   if (binl_sum >  lmax_loc-lmin_loc+1) { */
/*     printf("number of bins is: %li.\n",*bintabs); */
/*     while (binl_sum > lmax_loc-lmin_loc+1) { */
/*       (*bintabs)--; */
/*       binl_sum=0; */
/*       for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum); */
/*     } */
/*     printf("number of bins is : %li.\n",*bintabs); */
/*     (*bintabs)++; */
/*     bintab[*bintabs-1]=lmax_loc-lmin_loc+2-binl_sum; */
/*     printf("adding last bin of size: %li and increasing number of bins to: %li.\n",bintab[*bintabs-1], *bintabs); */
/*     binl_sum=0; */
/*     for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum); */
/*     printf("bin sum is now: %li.\n",binl_sum); */
/*   } */

/*   if (*bintabs <= 0 || bintab == NULL) { printf("ERROR: wrong aruments\n"); /\* cb = new power_spectrum(-1); *\/ return NULL; } */

/*   // find out the size of the binned C_b */
/*   cls = lmax+1; // original size of the C_l */
/*   cbs = lmax + 1  - (lmax_loc-lmin_loc+1) + *bintabs; // total size of binned C_l: bs */
/*   Cls = binl_sum; // size of original C_l for binning */
/*   Cbs = *bintabs; // size of binned C_l: bs */
/*   if (w == 0) { EQweights = false; set_cosmic_variance(); } else { EQweights = true; } */

/*   printf("  -- total size of cl vector before binning: %li\n",cls); */
/*   printf("  -- size of Cl vector for binning: %li\n",Cls); */
/*   printf("  -- total size of cb - binned Cl: %li\n",cbs); */
/*   printf("  -- size of Cb binned vector: %li\n",Cbs); */

/*   // define the structures needed */
/*   matrix <double> Cl(Cls,1); // original C_l vector for binning */
/*   matrix <double> Cb(Cbs,1); // binned C_b vector */
/*   matrix <double> M(Cbs,Cls); // binning operator */
/*   cb = new power_spectrum(cbs); */
/*   effl = new double[*bintabs]; */

/*   // copy the power spectra outside the reigon for binning */
/*   for (i=0;i<lmin_loc;i++) { cb->set_l(i,C_l[i][0]); cb->set_Cl(i,C_l[i][1]); } // lower end */
/*   tmp = lmax_loc+1; j=lmin_loc+(*bintabs); */
/*   for (i=tmp;i<=lmax;i++)  { cb->set_l(j,C_l[i][0]); cb->set_Cl(j,C_l[i][1]); j++;} // upper end shifted as to fit the binned Cb in range lmin_loc, lmax_loc */

/*   // prepare the power spectrum vector for binning */
/*   for (j=0;j<Cls;j++) {  Cl(j,0) = C_l[j+lmin_loc][1]; } // copy the power spectra part for binning */

/*   // prepare the binning matrix operator */

/*   for (i=0;i<Cbs;i++)  */
/*     for (j=0;j<Cls;j++)  */
/*       M(i,j) = 0;  // zero to binning matrix */

/*   k=0; */
/*   for (i=0;i<Cbs;i++) { */
/*     //if (i==0) jst = i; else jst = bintab[i]-1; */
/*     if (EQweights) w = 1/(double)bintab[i];  */
/*     effl[i] = 0; wsum=0; */
/*     for (j=k;j<k+bintab[i];j++) { */
/*       if (!EQweights) { w = C_l[j+lmin_loc][2]; wsum+=w; } */
/*       M(i,j) = w;   */
/*       effl[i]+=(double)(j+lmin_loc)*w;       */
/*     } */
/*     if (!EQweights) { for (j=k;j<k+bintab[i];j++) { M(i,j)=M(i,j)/wsum; } effl[i]/=wsum; }     */
/*     printf("effective ls: %lf\n",effl[i]); */
/*     k+=bintab[i]; */
/*   } */

/*   // do the binning */
/*   Cb = M*Cl; */

/*   // rewrite the Cb onto the object */

/*   for (i=0;i<Cbs;i++)  {  */
/*     j=lmin_loc+i; cb->set_l(j,effl[i]); cb->set_Cl(j,Cb(i,0));  */
/*     printf("l %lf Cl %lE            binned Cb: l %lf Cb %lE\n",effl[i],Cb(i,0),cb->get_l(j),cb->get_Cl(j)); */
/*   }  */

/*   delete [] effl; */

/*   return cb; */

/* } */

/* long mscsMap::get_seed_offset() {  return seed_offset; } // commented out on 2009-01-22 -- moved to cpeds and individual applications*/
/* void mscsMap::set_seed_offset(long offset) {  seed_offset = offset; } // commented out on 2009-01-22 -- moved to cpeds and individual applications */

/* //\************************************************************************ */
/* // projection_type: 1 -- Mollweide projection */
/* // projection_type: 2 -- Miller projection */
/* // projection_type: 3 -- Orthographic projection */
/* // dl, db are given in degrees */
/* void mscsMap::calculate_flat_coord(int projection_type, double dl, double db) { // convert the corrdinated data into some XY projection for plotting purposes */
/*   filenamestr projtype,eastwest; */
/*   filenamestr command_str; */
/*   projUV p; */
/*   projPJ pj; */
/*   long i; */

/*   printf("|%s>  -- calculating projected coordinates of the map\n",object_name.c_str()); */
/*   if (coord_loaded == 0) { set_map_coord(0,0); } */

/*   //savetxtC(tmpfile,3); */
/*   if (projection_type == 1) { // doing Mollweide projection */
/*     strcpy(projtype,"moll"); */
/*   } */
/*   if (projection_type == 2) { // doing Miller projection */
/*     strcpy(projtype,"mill"); */
/*   } */
/*   if (projection_type == 3) { // doing orthogonal projection */
/*     strcpy(projtype,"ortho"); */
/*   } */
/*   if (projection_type == 4) { // doing Aitoff  projection */
/*     strcpy(projtype,"aitoff"); */
/*   } */
/*   if (projection_type == 5) { // doing stereographic  projection */
/*     strcpy(projtype,"stere"); */
/*   } */

/*   if (dl > 0) { strcpy(eastwest,"w"); } else { strcpy(eastwest,"e"); } */
/*   sprintf(command_str,"+proj=%s +lat_0=%lf +lon_0=%lf%s +R=1",projtype,db,dl,eastwest); */
/*   flat_coord_num = pix_num; makekill_space_manager("make","C",2); */
/*   pj = pj_init_plus(command_str); */
/*   for (i=0;i<flat_coord_num;i++) { */
/*     p.u = PI-map->n[i].l;    p.v = map->n[i].b; // conversion from phi in (0,2pi) --> east and west longitude from (-pi,pi) -- here we inverse the direction of rising l to be consistent with the standards in publications that the l rises form the center of the figure left-wise, jumps to the right edge and rises from pi->2pi towards the picture center */
/* /\*     p.u = map->n[i].l-PI;    p.v = map->n[i].b; // EXPERIMENTAL !! *\/ */
/*     p = pj_fwd(p,pj); */
/*     flatcoord[i].x = p.u;    flatcoord[i].y = p.v; */
/*   } */
/*   pj_free(pj); */
/*   flat_coord_loaded = 1; flat_coord_ordering = map_ordering; */
/* } */

/* //\************************************************************************ */
/* float*  mscsMap::create_resized_flat_map(int printinfo) { */
/*   double factor; */
/*   double factorx, factory; */
/*   long int i,j,x,y; */
/*   long tmp,k,suppress_output_no=0,suppress_output_lim=3; */

/*   calculate_flat_map_parameters(printinfo); */
/* /\*   xpix_num_flat = map_size_x; *\/ */
/*   xpix_num_flat = (long)round(4*(double)nside); // this is the map resolution parameter - to get rid of the empty fields you can reduce this or use interpolation in the gaps in the map */
/*   //xpix_num_flat = 2000; */
/*   //ypix_num_flat = 2*cpeds_get_ring_num_healpix(nside);  */
/*   ypix_num_flat = (int)(ceil((double)xpix_num_flat/x2yratio_flat)); */

/*   factor = 1.00;//(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed */
/*   factorx = factor*(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed */
/*   factory = factor*(maxy_flat-miny_flat)/(double)ypix_num_flat; /// how many times the map was squeezed */
/* /\*   ypix_num_flat = (int)(ceil((maxy_flat-miny_flat)/factor))+1; *\/ */

/*   // *** EXPERIMENTAL *** // */
/* /\*   xpix_num_flat = cpeds_get_pixnum_in_ring_healpix(nside,2*nside); *\/ */
/* /\*   ypix_num_flat = 4*nside; *\/ */
/*   // *** EXPERIMENTAL *** // */

/*   pix_num_flat = (long int)xpix_num_flat*(long int)ypix_num_flat; */

/*   printf("|%s>  -- xpix_num_flat: %i, ypix_num_flat: %i pix_num_flat: %li\n",object_name.c_str(), xpix_num_flat, ypix_num_flat, pix_num_flat); */

/*   // initiate space for flatmap */
/*   makekill_space_manager("make","T",2);  //float * flattab = new float[pix_num_flat]; */
/* /\*   if (xpix_num_flat == 256) exit(0); *\/ */
/*   // create flatmap */
/*   for (i=0;i<pix_num_flat;i++) { flatmap[i] = -1e10; } */
/*   //for (i=0;i<pix_num_flat;i++) { flatmap[i] = 0; } */

/*   // *** EXPERIMENTAL *** // */
/* /\*   conv_nest2ring(map); k=0; *\/ */
/* /\*   for (i=0;i<ypix_num_flat;i++) { *\/ */
/* /\*     tmp =cpeds_get_pixnum_in_ring_healpix(nside,i); *\/ */
/* /\*     for (j=0;j<tmp;j++) { *\/ */
/* /\*       flatmap[xpix_num_flat*i+j] = (float)(map->T[k]); k++; } *\/ */
/* /\*   } *\/ */
/*   // *** EXPERIMENTAL *** // */

/*   for (i=0;i<flat_coord_num;i++) { */
/*     x = (long)round((flatcoord[i].x-minx_flat)/factorx); */
/*     y = (long)round((flatcoord[i].y-miny_flat)/factory); */
/*     //if (mask_loaded == 0) { */
/*     if (x >= xpix_num_flat) { if (printinfo!=0 && suppress_output_no<=suppress_output_lim) { printf("|%s>  -- WARNING !!  x >= xpix_num_flat, setting x = xpix_num_flat-1\n",object_name.c_str()); suppress_output_no++; } x = xpix_num_flat-1; } */
/*     if (y >= ypix_num_flat) { if (printinfo!=0 && suppress_output_no<=suppress_output_lim) { printf("|%s>  -- WARNING !!  y >= xpix_num_flat, setting y = ypix_num_flat-1\n",object_name.c_str()); suppress_output_no++; } y = ypix_num_flat-1; } */
/*     if (x < 0) { if (printinfo!=0 && suppress_output_no<=suppress_output_lim) { printf("|%s>  -- WARNING !!  x < 0, setting x = 0\n",object_name.c_str()); suppress_output_no++; } x = 0; } */
/*     if (y < 0) { if (printinfo!=0 && suppress_output_no<=suppress_output_lim) { printf("|%s>  -- WARNING !!  y < 0, setting y = 0\n",object_name.c_str()); suppress_output_no++; } y = 0; } */

/*     if (xpix_num_flat*y+x > pix_num_flat) { if (printinfo!=0 && suppress_output_no<=suppress_output_lim ) { printf("|%s>  -- WARNING !!   xpix_num_flat*y+x > pix_num_flat\n",object_name.c_str()); suppress_output_no++; }} */
/*     if (suppress_output_no==suppress_output_lim) { printf("|%s>  Further warning messages will be suppressed\n",object_name.c_str());  } */

/*     if (map->m == NULL) { */
/*       if (flatmap_change_zrange) { */
/* 	if ((map->T[i]) > flatmap_minT && (map->T[i]) < flatmap_maxT) { flatmap[xpix_num_flat*y+x] = (float)(map->T[i]); }  */
/* 	else { */
/* 	  if ((map->T[i]) >= flatmap_maxT) flatmap[xpix_num_flat*y+x] = flatmap_maxT;  */
/* 	  else flatmap[xpix_num_flat*y+x] = flatmap_minT;  */
/* 	} */
/*       }  */
/*       else flatmap[xpix_num_flat*y+x] = (float)(map->T[i]);     */
/*     } */
/*     else { */
/*       if (map->m[i] < 1) { */
/* /\* 	flatmap[xpix_num_flat*y+x] = (float)(map->T[i]*map->m[i]); } //printf("maska fm=%E\n",flatmap[xpix_num_flat*y+x]); *\/ */
/* /\*       else { flatmap[xpix_num_flat*y+x] = (float)(map->T[i]); //printf("poza fm=%E\n",flatmap[xpix_num_flat*y+x]);  *\/ */
/* /\*       } *\/ */

/* 	if (flatmap_change_zrange) { */
/* 	  if ((map->T[i]) > flatmap_minT && (map->T[i]) < flatmap_maxT) { flatmap[xpix_num_flat*y+x] = (float)(map->T[i]*map->m[i]); } else { */
/* 	    if ((map->T[i]) >= flatmap_maxT) flatmap[xpix_num_flat*y+x] = flatmap_maxT; else flatmap[xpix_num_flat*y+x] = flatmap_minT; }} else flatmap[xpix_num_flat*y+x] = (float)(map->T[i]*map->m[i]);    } */
/*       else { 	if (flatmap_change_zrange) { */
/* 	  if ((map->T[i]) > flatmap_minT && (map->T[i]) < flatmap_maxT) { flatmap[xpix_num_flat*y+x] = (float)(map->T[i]); } else { */
/* 	    if ((map->T[i]) >= flatmap_maxT) flatmap[xpix_num_flat*y+x] = flatmap_maxT; else flatmap[xpix_num_flat*y+x] = flatmap_minT; }} */
/* 	else flatmap[xpix_num_flat*y+x] = (float)(map->T[i]); */
/*       } */

/*       if (map->m[i] == 0) flatmap[xpix_num_flat*y+x] = -1e10; */
/*     } // here the mask only if loaded also is plotted; multimask is treated as 1 */
/*   } */

/*   // linear interpolation if the gaps resulting from the projections effects that stretch the map onto a plane and some points are missing // this is for visual effect only */
/*   k=0;  for (i=0;i<(int)(0.15*ypix_num_flat);i++) { for (j=0;j<xpix_num_flat;j++) { k++; if ( flatmap[k] == -1e10 && flatmap[k+1] != -1e10 ) flatmap[k]=(flatmap[k-1]+flatmap[k+1])/2; } } */
/*   k=(int)(0.75*ypix_num_flat)*xpix_num_flat; */
/*   for (i=(int)(0.75*ypix_num_flat);i<ypix_num_flat;i++) { for (j=0;j<xpix_num_flat;j++) { k++; if ( flatmap[k] == -1e10 && flatmap[k+1] != -1e10 ) flatmap[k]=(flatmap[k-1]+flatmap[k+1])/2; } } */

/*   flat_map_loaded = 1; */
/*   return flatmap; */
/* } */

//************************************************************************

