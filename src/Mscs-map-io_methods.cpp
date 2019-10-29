//#define _ISOC99_SOURCE
//#include <features.h>
//#include <sys/stat.h>
//#include <cpgplot.h>
//#include <proj_api.h>
#include <errno.h>
#include "Mscs-map.h"
#include "cpeds-fits.h"
//#include <curses.h>
//#include "cpeds-consts.h"
//#include "cpeds-math.h"
//#include "cpeds-templates.h"


using namespace std;
using namespace math;















//  MAP HANDLERS

//************************************************************************
cpedsList<double> mscsMap::get_nonMasked(string what) {
	long i,pix_num=pixNum();
	cpedsList<double> cl;
	if (what == "T") {
		if (maskLoaded()) {
			for (i=0; i<pix_num;i++) { if (!isMasked(i)) { cl.append(get_T(i)); } }
		}
		else {
			cl=get_T();
		}
	}
	else {
		msgs->criticalError("get_nonMasked(string what): 'what!=T' is not implemented yet. Sorry",High);		
	}
	return cl;
}
/***************************************************************************************/
cpedsDirectionSet mscsMap::get_nonMaskedDirs() {
	long i,pix_num=pixNum();
	cpedsDirectionSet dirs;
	if (maskLoaded()) {
		for (i=0; i<pix_num;i++) { if (!isMasked(i)) { dirs.append(get_C(i)); } }
	}
	else {
		return this->get_n();
	}
	return dirs;
}
/* ******************************************************************************************** */
cpedsDirectionSet mscsMap::get_nonMaskedDirsT() {
	long i,pix_num=pixNum();
	cpedsDirectionSet dirs;
	if (maskLoaded()) {
		for (i=0; i<pix_num;i++) { if (!isMasked(i)) { dirs.append(get_C(i)); dirs.last().setVal(get_T(i)); } }
	}
	else {
		return this->get_n();
	}
	return dirs;	
}
/***************************************************************************************/
cpedsDirectionSet mscsMap::get_maskedDirs() {
	long i,pix_num=pixNum();
	cpedsDirectionSet dirs;
	if (maskLoaded()) {
		for (i=0; i<pix_num;i++) { if (isMasked(i)) { dirs.append(get_C(i)); } }
	}
	return dirs;
}

//************************************************************************
// region numbering passed to this handler is from 1 to reg_num of the multimask
void mscsMap::set_Ttoval(double t, double val) {
	long i,pix_num=pixNum();
	for (i=0; i<pix_num;i++) { if (get_T(i) == t) { T(i) = val; } }
}

//************************************************************************
// region numbering passed to this handler is from 1 to reg_num of the multimask
void mscsMap::set_Treg(long reg, double t) {
	long i,pix_num=pixNum();
	for (i=0; i<pix_num;i++) { if (get_m(i) == (double)reg) { T(i) = t; } }
}
//************************************************************************
//************************************************************************
// region numbering passed to this handler is from 1 to reg_num of the multimask
// replaces the temperature if the absolute value given in parameter (T) is greater than the temperature in the map
void mscsMap::set_Treg_if_gr_abs(long reg, double t) {
	long i,pix_num=pixNum();
	for (i=0; i<pix_num;i++) { if (get_m(i) == (double)reg) if (fabs(t) > fabs(get_T(i))) T(i) = t;  }
}
//************************************************************************
// region numbering passed to this handler is from 1 to reg_num of the multimask
// replaces the temperature if the absolute value given in parameter (T) is smaller than the temperature in the map
void mscsMap::set_Treg_if_le_abs(long reg, double t) {
	long i,pix_num=pixNum();
	for (i=0; i<pix_num;i++) { if (m(i) == (double)reg) if (fabs(t) < fabs(get_T(i))) T(i) = t;  }
}

//************************************************************************
// region numbering passed to this handler is from 1 to reg_num of the multimask
// sets the value at whole region towards l,b direction
void mscsMap::set_Treg(double l, double b, double t) {
	long num;
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-b,l,1);
	set_Treg((long)(get_m(num)),t);
}
//************************************************************************
double mscsMap::get_T(const cpedsDirection& nn) const { // returns the temperature from requested direction
	long num;
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-nn.b(),nn.l(),1);
	return get_T(num);
}
//************************************************************************
double mscsMap::get_T(const cpedsDirection& nn, long *i) const { // returns the temperature from requested direction
	long num;
	
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-nn.b(),nn.l(),1);
	(*i)=num;
	return get_T(num);
}

//************************************************************************
void mscsMap::set_T(const cpedsDirection& nn, double t) { // sets the temeperature at requested direction
	long num;
	
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-nn.b(),nn.l(),1);
	T(num) = t;
}
/***************************************************************************************/
void mscsMap::set_T(double l, double b, double val) {
	long num;
	
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-b,l,1);
	T(num) = val;	
}

//************************************************************************
// returns the temperature from direction given in l,b defined as in Mscs package (in rad. with b=90 at north pole)
double mscsMap::get_T(double l, double b) const { // returns the temperature from requested direction
	long int num;
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-b,l,1);
	return get_T(num);
}


//************************************************************************
// this routine returns the coordinates to the hottest/coldest point in the map
// depending on the max parameter; it is assumed that the statistics of the map are
// up to date; the returned coordinates are in radians
cpedsDirection mscsMap::get_TminmaxCoord(bool max) {
	cpedsDirection n;
	double th,phi;
	bool wasRing=false;
	long i;
	
	if (mapInfo.ordering == ring) { conv_ring2nest(); wasRing=true; }
	else { wasRing=false; }
	
	
	if (max)
		cpeds_pix2ang_healpix(nside(),i,&th,&phi,1);
	else
		cpeds_pix2ang_healpix(nside(),i,&th,&phi,1);
	
	n.set(phi, PIsnd-th);
	
	
	if (wasRing) conv_nest2ring();
	return n;
}


//************************************************************************
double mscsMap::get_m(double l, double b) const {
	long num;
	if (!loaded.m) return 1;
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-b,l,1);
	return get_m(num);
}

//************************************************************************
void mscsMap::get_C(long int num, double * l, double * b) const { // returns the coordinates of i'th pixel
	*l = get_C(num).lon(); *b = get_C(num).lat();
}

//************************************************************************
long mscsMap::get_Ci(const cpedsDirection& nn) const { // returns the pix number at given direction
	long num;
	cpeds_ang2pix_healpix(nside(),&num,PIsnd-nn.b(),nn.l(),1);
	return num;
}


// FILES METHODS


//TEMPERATURE


//************************************************************************
long mscsMap::savebinT(string map_file) {
	long res=0;
	msgs->say("Saving temperature map to binary file: "+map_file,High);
	res=T().save(map_file,true,"double");
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadbinT(string  map_file, mapOrderings ordering) {
	long res=0;
	cpedsList<double> t;
	msgs->say("Loading temperature map from binary file: "+map_file,High);
	res=t.load(map_file,true,"double");
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	T()=t;
	loaded.T=true;
	set_pixNum(T().size());
	msgs->say("Success",Low);
	if (T().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(T().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	set_ordering(ordering);
	if (ordering==ring)  conv_ring2nest();
	return res;
}
//************************************************************************
long mscsMap::savetxtT(string map_file) {
	long res=0;
	msgs->say("Saving temperature map to text file: "+map_file,High);
	res=T().save(map_file,false,"double");
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadtxtT(string  map_file, mapOrderings ordering) {
	long res=0;
	cpedsList<double> t;
	msgs->say("Loading temperature map from text file: "+map_file,High);
	res=t.load(map_file,false,"double");
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	T()=t;
	loaded.T=true;
	set_pixNum(T().size());
	msgs->say("Success",Low);
	if (T().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(T().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	set_ordering(ordering);
	if (ordering==ring)  conv_ring2nest();
	return res;
}

//************************************************************************
long mscsMap::savebinm(string fileName) {
	long res=0;
	msgs->say("Saving mask map to binary file: "+fileName,High);
	res=m().save(fileName,true,"double");
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadbinm(string  fileName) {
	long res=0;
	cpedsList<double> t;
	msgs->say("Loading mask from binary file: "+fileName,High);
	res=t.load(fileName,true,"double");
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	m()=t;
	loaded.m=true;
	set_pixNum(m().size());
	msgs->say("Success",Low);
	return res;
}

//************************************************************************
// this is for loading fits healpix system written maps
// long mscsMap::loadfitsT(string  map_file, int how) {
//   loadsave_manager("load","fits","T",how,map_file,-1);
// }
//************************************************************************

//void mscsMap::savefitsT(string  map_file) { // TODO
//}
//************************************************************************


//************************************************************************

//void mscsMap::savefitsT(string  map_file) { // TODO
//}
//************************************************************************
// void mscsMap::make_map_space(map_pixel * map) {
/*   if (map == map_ptr) { */
/*     makekill_space_manager("make","T",1); */
/*     makekill_space_manager("make","C",1); */
/*     makekill_space_manager("make","m",1); */
/*     makekill_space_manager("make","N",1); */
/*   } */
/*   else { */
/*     map->T = new double[pix_num];  */
/*     map->n = new direction[pix_num]; */
/*     map->m = new double[pix_num]; */
/*     map->N = new double[pix_num]; */
/*   } */
// }
//************************************************************************
// frees allocated memory by map structures, NOT changing the flags status
// eg. to be used by some map transformations
// this is only for special applications eg. ^^

// void mscsMap::kill_map_space(map_pixel * map) {
//   if (map->T().size() > 0) { T().clear(); }
//   if (Q().size() > 0) { Q().clear(); }
//   if (U().size() > 0) { U().clear(); }
//   if (V().size() > 0) { V().clear(); }
//   if (N().size() > 0) { N().clear(); }
//   if (n().size() > 0) { n().clear(); }
//   if (m().size() > 0) { m().clear(); }
// }

//************************************************************************
void mscsMap::killAll() {
	if (T().size() > 0) { T().clear(); }
	if (Q().size() > 0) { Q().clear(); }
	if (U().size() > 0) { U().clear(); }
	if (V().size() > 0) { V().clear(); }
	if (N().size() > 0) { N().clear(); }
	if (n().size() > 0) { n().clear(); }
	if (m().size() > 0) { m().clear(); }
}

//************************************************************************
// mscsMap::map_pixel mscsMap::clone_map_space(mscsMap::map_pixel * map, long pix_num) {
//   map_pixel newmap;
//   if (map->T == NULL) newmap.T = NULL; else newmap.T = new double[pix_num];
//   if (map->n == NULL) newmap.n = NULL; else newmap.n = new cpeds_direction[pix_num];
//   if (map->m == NULL) newmap.m = NULL; else newmap.m = new double[pix_num];
//   if (map->N == NULL) newmap.N = NULL; else newmap.N = new double[pix_num];
//   return newmap;
// }
/* void mscsMap::killmapTC() { // kills the themparature map and coordinates */
/*   printf("|%s> * removing the temperature map from memory\n"); */
/*   if (map_loaded == 1) { delete map; } */
/*   map_loaded = 0; */
/*   coord_loaded = 0; */
/*   mask_loaded = 0; */
/* } */
//************************************************************************
//************************************************************************
// how : 1 - copy data as is.
// how : 2 - change the original data resolution first to fit to the sizes in this object
// how : 3 - kills map and initiates  nside as in the source object (CAREFULL!!)
// (this is only applicable to map stuff (not alms flatmap etc.)
void mscsMap::import_map_data(mscsMap &from_here, string what, int how) {
	msgs->say("Importing map data from "+from_here.getName()+" about "+what,Medium);
	if (how == 3) { killAll(); set_nside(from_here.nside()); }
	
	if (how == 2) { from_here.change_map_resolution(nside()); }
	
	if (what=="T") {
		T()=from_here.get_T();
		mapInfo.meanT = from_here.get_meanT();
		mapInfo.varianceT = from_here.get_varianceT();
		mapInfo.skewnessT = from_here.get_skewnessT();
		mapInfo.kurtosisT = from_here.get_kurtosisT();
		loaded.T = true;
		return;
	}
	if (what=="Q") {
		Q()=from_here.get_Q();
		loaded.Q = true;
		return;
	}
	if (what=="U") {
		U()=from_here.get_U();
		loaded.U = true;
		return;
	}
	if (what=="V") {
		V()=from_here.get_V();
		loaded.n = true;
		return;
	}
	if (what=="m") {
		m()=from_here.get_m();
		loaded.m = true;
		mask= from_here.mask;
		return;
	}
	if (what=="n") {
		n()=from_here.get_n();
		loaded.n = true;
		return;
	}
	if (what=="N") {
		N()=from_here.get_N();
		loaded.N = true;
		return;
	}
	msgs->say("!! WARNING: the import has no effect: I don't know what "+what+" is !!",High); //}  //commented out during transition into version-1.0
	
}

//************************************************************************
// how : 1 - copy data as is.
// how : 2 - change the original data resolution first to fit to the sizes in this object
// (this is only applicable to map stuff (not alms flatmap etc.)
// void mscsMap::import_map_data_to(mscsMap &from_here, string what, string where, int how) {
//   long i;
//   double *to,*from;
//   cpeds_direction *ton,*fromn;
//   //a_lm *toalm,*fromalm; //commented out during transition into version-1.0

//   msgs->say("Importing map data from "+from_here.getName()+" about "+what+" to "+where,High);

//   // define where to copy
//   if (how == 2) { from_here.change_map_resolution(nside); }
//   makekill_space_manager("make",where.c_str(),1);
//   if (where == "T") to=map->T; else {
//   if (where == "m") to=map->m; else {
//   if (where == "N") to=map->N; else {
//   if (where == "n") ton=map->n; else {
//     //if (where == "F" || where == "alms") toalm=alm; else {  //commented out during transition into version-1.0
//     msgs->say("!! WARNING: the import has no effect: I don't know what "+what+" is !!",High); return; }}}} //}  //commented out during transition into version-1.0

//   // define from where copy
//   if (what == "T") from=from_here.map->T; else {
//   if (what == "m") from=from_here.map->m; else {
//   if (what == "N") from=from_here.map->N; else {
//   if (what == "n") fromn=from_here.map->n; else {
//     //if (what == "F" || what == "alms") fromalm=from_here.alm; else {  //commented out during transition into version-1.0
//     msgs->say("!! WARNING: the import has no effect: I don't know what "+what+" is !!",High); return; }}}} //}  //commented out during transition into version-1.0
//   //}}} //}  //commented out during transition into version-1.0


//   // make copy
//   if (what == "T" || what == "m" || what == "N") {    for (i=0;i<pix_num;i++) { to[i] = from[i]; }  }
//   if (what == "n") {    for (i=0;i<pix_num;i++) { ton[i] = fromn[i]; }  }
//   //if (what == "F") {    for (i=0;i<pix_num;i++) { toalm[i] = fromalm[i]; } }  //commented out during transition into version-1.0

//   // set operation specific parameters
//   if (where =="T") { map_loaded = 1; if (what == "T") { meanT = from_here.meanT;    varianceT = from_here.varianceT; skewnessT = from_here.skewnessT; kurtosisT = from_here.kurtosisT;} else { calculate_map_stats(0); }}
//   if (where == "m") { mask_loaded = 1;  if (what == "m") multi_mask_reg_num = from_here.multi_mask_reg_num;  }
//   if (where == "N") { Nobs_loaded = 1;  }
//   if (where == "n") { coord_loaded = 1;  }
//   //if (what == "alms") {    alms_loaded = 1;  }  //commented out during transition into version-1.0

// }



void mscsMap::printtxtT() {
	long int i;
	msgs->say("printing temeratures to screen:",High);
	
	for (i=0;i<pixNum();i++) {
		msgs->say("T[ "+msgs->toStr(i)+" ] = "+msgs->toStr(get_T(i)),Medium);
	}
}






//COORDINATES





//************************************************************************

long mscsMap::savebinC(string  fileName, string how) {
	long res=0;
	msgs->say("Saving coordinates to binary file: "+fileName,High);
	res=n().save(fileName,"binary"+how,false);
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadbinC(string fileName, string how) {
	long res=0;
	cpedsDirectionSet c;
	msgs->say("Loading coordinates from binary file: "+fileName,High);
	n().clear();
	res=n().load(fileName,"binary"+how,false);
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	msgs->say("Success",Low);
	if (mapLoaded() and n().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(n().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	loaded.n=true;
	set_pixNum(n().size());
	return res;
}
//************************************************************************
long mscsMap::savetxtC(string  fileName) {
	long res=0;
	msgs->say("Saving coordinates to text file: "+fileName,High);
	res=n().save(fileName,"",false);
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadtxtC(string  fileName) {
	long res=0;
	cpedsDirectionSet c;
	msgs->say("Loading coordinates from text file: "+fileName,High);
	res=c.load(fileName,"",false);
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	n()=c;
	msgs->say("Success",Low);
	if (n().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(n().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	loaded.n=true;
	set_pixNum(n().size());
	return res;
}
//************************************************************************
// prints the coordinates from the database as requested
// what = g - 'galactic' coordinates eg. from the healpix map (from the map - data)
// what = f - flat coordinates from the flatcoord data - this is eg. after making a projection
// what = g proj - prints the galactic coordinates but in a format acceptable by 'proj' program:
//                i.e. the lattidutes are calculated from equator at 0 to + and - 90 deg at the poles
void mscsMap::printtxtC(string what) {
	long int i;
	double l,b;
	
	if (what=="g") { msgs->say("printing galactic (map native) coordinates to screen ",High); }
	//if (strcmp(what,"f") == 0) { msgs->say("flat (eg. from proj program) "); }
	if (what=="g proj") { msgs->say("printing the proj program undestandable galactic (map native but in deg and shifted in lattitude by -180 deg) coordinates to screen",High); }
	
	if (what == "g") n().print();
	
	for (i=0;i<pixNum();i++) {
		//if (strcmp(what,"f") == 0) { x = flatcoord[i].x; y = flatcoord[i].y; msgs->say("i = %li, x = %lE, y = %lE\n",i,x,y);}
		if (what=="g proj") { l = get_C(i).l(); b = get_C(i).b(); msgs->say(msgs->toStrf(180/PI*l)+" "+msgs->toStrf(180/PI*b),Medium);}
	}
}
//************************************************************************

// NUMBER OF OBSERVATIONS

//************************************************************************
long mscsMap::savebinN(string map_file) {
	long res=0;
	msgs->say("Saving pixel observations number to binary file: "+map_file,High);
	res=N().save(map_file,true,"double");
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadbinN(string  map_file) {
	long res=0;
	cpedsList<double> t;
	msgs->say("Loading pixel observations number from binary file: "+map_file,High);
	res=t.load(map_file,true,"double");
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	N()=t;
	msgs->say("Success",Low);
	loaded.N=true;
	set_pixNum(N().size());
	if (N().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(N().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	return res;
}
//************************************************************************
long mscsMap::savetxtN(string map_file) {
	long res=0;
	msgs->say("Saving pixel observations number to text file: "+map_file,High);
	res=N().save(map_file,false,"double");
	if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
	msgs->say("Success",Low);
	return res;
}
//************************************************************************
long mscsMap::loadtxtN(string  map_file) {
	long res=0;
	cpedsList<double> t;
	msgs->say("Loading pixel observations number from binary file: "+map_file,High);
	res=t.load(map_file,false,"double");
	if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
	N()=t;
	msgs->say("Success",Low);
	loaded.N=true;
	set_pixNum(n().size());
	if (N().size()!=pixNum()) { msgs->warning("The size of the loaded data ("+msgs->toStr(N().size())+") is inconsistent with the current pixel number ("+msgs->toStr(pixNum())+").",High); }
	return res;
}
//************************************************************************
























// FILE SAVING HELPER

//************************************************************************
//int  mscsMap::makekill_space_manager(string whattodo, string fileformat, string for_what, string what_file, int how) {
void  mscsMap::makekill_space_manager(string whattodo, string for_what, long how) {
	
	if (for_what=="T") { // allocating memory for tempearature map
		if ((whattodo=="make") || (whattodo=="load")) {
			msgs->say("allocating memory for temperature map",High);
			T().makeLength(pixNum());
		}
		if (whattodo=="kill") {
			msgs->say("deleting memory for temperature map",High);
			T().makeLength(0);
			loaded.T=false;
		}
	}
	
	if (for_what=="m") { // allocating memory for mask
		if ((whattodo=="make") || (whattodo=="load")) {
			msgs->say("allocating memory for mask",High);
			m().makeLength(pixNum());
		}
		if (whattodo=="kill") {
			msgs->say("deleting memory for mask",High);
			m().makeLength(0);
			resetMaskInfo();
			loaded.m=false; // BLmodification (Sep 24, 2011, 1:39:52 AM): inserted when using draw_maps_new with HP multimasks generated on fly
		}
	}
	
	if (for_what=="C") { // allocating memory for mask
		if ((whattodo=="make") || (whattodo=="load")) {
			msgs->say("allocating memory for coordinates",High);
			n().setLength(pixNum());
			mapInfo.coord_num=pixNum();
		}
		if (whattodo=="kill") {
			msgs->say("deleting memory for coordinates",High);
			n().setLength(0);
			mapInfo.coord_num=0;
			loaded.n=false;
		}
	}
	
	
	if (for_what=="N") { // allocating memory for mask
		if ((whattodo=="make") || (whattodo=="load")) {
			msgs->say("allocating memory for pixel observation number",High);
			N().makeLength(pixNum());
		}
		if (whattodo=="kill") {
			msgs->say("deleting memory for pixel observation number",High);
			N().makeLength(0);
			loaded.N = false;
		}
	}
	
	if (for_what=="rot") { // allocating memory for mask
		if ((whattodo=="make") || (whattodo=="load")) {
			msgs->say("allocating memory for rotation map",High);
			rotation_map.makeLength(pixNum());
		}
		if (whattodo=="kill") {
			msgs->say("deleting memory for rotation map",High);
			rotation_map.makeLength(0);
			loaded.rotation = false;
		}
	}
	
}





















//************************************************************************
// void mscsMap::read_binmap_parameters(string  map_file) {
//   double tmp;
//   //long int n = 0;
//   long pix_num_loc,nside_loc;
//   //printf("%s",*map_file);
//   //filenamestr mapfile;
//   struct stat64 info;
//   //strcpy(mapfile,map_file.c_str());
// //  strcat(mapfile,"-T-bin");


//   stat64(map_file.c_str(),&info);
//   pix_num_loc = (long)(info.st_size)/(long)sizeof(tmp);
//   nside_loc = (long int)sqrt((double)pix_num_loc/12);
//   msgs->say("current nside is: "+msgs->toStr(nside),Medium);
//   msgs->say("the file nside appears to be: "+msgs->toStr(nside_loc),Medium);

//   if (nside != 0 && nside_loc != nside) { // nside==0 is here a check whether something is loaded - may it be a temperature map or coordinates of number of observations etc.
//     msgs->say("WARNING: the intended file to load is of incompatible size: Will change resolution of the structures to fit the new size. The different structure sizes within the same object are not yet allowed (and probably don't make sense).",High);
//     change_map_resolution(nside_loc);
//   }
//   nside=nside_loc;
//   pix_num=pix_num_loc;
//   rows_num = cpeds_get_ring_num_healpix(nside);

//   msgs->say("map file parameters: file size (bytes) "+msgs->toStr((long)info.st_size),Medium);

// /*   f = fopen(mapfile,"rb"); */
// /*   while (!feof(f)) { fread(&tmp,sizeof(tmp),1,f); n++;} */
// /*   fclose(f); */
// /*   pix_num = n-1; */

//   coord_num = pix_num;

//   msgs->say("map file parameters: nside = "+msgs->toStr(nside)+", pix_num = "+msgs->toStr(pix_num),Medium);
//   //printf("|%s>  -- map file size: bytes = %li\n",object_name.c_str(),pix_num*(long int)(sizeof(tmp)));
//   msgs->say("map file name: "+map_file,Medium);
//   msgs->say("appriximated pixel size [deg]: "+msgs->toStr(PI180inv * cpeds_pix_size_healpix(nside)),Medium);

// }
//************************************************************************
// void mscsMap::read_txtmap_parameters(string map_file) {
//   double tmp;
//   long int n = 0;
//   long pix_num_loc,nside_loc;
//   //filenamestr mapfile;
//   //strcpy(mapfile,map_file.c_str());
// //  strcat(mapfile,"-T-txt");

//   f = fopen(map_file.c_str(),"r");
//   while (!feof(f)) {  fscanf(f,"%lE",&tmp); n++; }
//   fclose(f);

//   pix_num_loc = n-1;
//   nside_loc = (long int)sqrt((double)pix_num_loc/12);
//   msgs->say("current nside is: "+msgs->toStr(nside),Medium);
//   msgs->say("the file nside appears to be: "+msgs->toStr(nside_loc),Medium);

//   if (nside != 0 && nside_loc != nside) { // nside==0 is here a check whether something is loaded - may it be a temperature map or coordinates of number of observations etc.
//     msgs->say("WARNING: the intended file to load is of incompatible size: Will change resolution of the structures to fit the new size. The different structure sizes within the same object are not yet allowed (and probably don't make sense).",High);
//     change_map_resolution(nside_loc);
//   }
//   nside=nside_loc;
//   pix_num=pix_num_loc;
//   rows_num = cpeds_get_ring_num_healpix(nside);

//   coord_num = pix_num;

//   msgs->say("map file parameters: nside = "+msgs->toStr(nside)+", pix_num = "+msgs->toStr(pix_num),Medium);
//   msgs->say("map file size: bytes = "+msgs->toStr(pix_num*(long int)sizeof(tmp)),Medium);
//   msgs->say("map file name: "+map_file,Medium);
//   msgs->say("appriximated pixel size [deg]: "+msgs->toStr(PI180inv * cpeds_pix_size_healpix(nside)),Medium);
// }
//************************************************************************
// void mscsMap::read_bincoord_parameters(string  coord_file) {
//   long int n = 0;
//   float tmp;
//   filenamestr coordfile;
//   strcpy(coordfile,coord_file.c_str());
// //  strcat(coordfile,"-C-bin");


//   f = fopen(coordfile,"rb");
//   while (!feof(f)) { fread(&tmp,sizeof(tmp),1,f); fread(&tmp,sizeof(tmp),1,f); n++;}
//   fclose(f);

//   //flat_coord_num = n-1; //commented out during transition into version-1.0
//   coord_num = n-1;
//   pix_num = n-1;
//   msgs->say("coordinates file parameters: coord_num = "+msgs->toStr(coord_num),Medium);
//   msgs->say("coordinates file size in memory: bytes = "+msgs->toStr(2*(n-1)*(long int)sizeof(tmp)),Medium);
//   msgs->say("coordinates file name: "+msgs->toStr(coordfile),Medium);
// }
// //************************************************************************
// void mscsMap::read_txtcoord_parameters(string  coord_file) {
//   long int n = 0;
//   double tmp;
//   filenamestr coordfile;
//   strcpy(coordfile,coord_file.c_str());
// //  strcat(coordfile,"-C-txt");

//   f = fopen(coordfile,"r");
//   while (!feof(f)) { fscanf(f,"%*f %*f"); n++; }
//   fclose(f);

//   //flat_coord_num = n-1; //commented out during transition into version-1.0
//   coord_num = n-1;
//   pix_num = n-1;
//   msgs->say("coordinates file parameters: coord_num = "+msgs->toStr(coord_num),Medium);
//   msgs->say("coordinates file size in memory: bytes = "+msgs->toStr(2*(n-1)*(long int)sizeof(tmp)),Medium);
//   msgs->say("coordinates file name: "+msgs->toStr(coordfile),Medium);
// }
long mscsMap::loadfits(const string& fileName, const QList<string>& colNames, const QList<string>& dstStructures) {
	// long mscsMap::loadfits(string fileName, string data) {
	fitsfile *fptr;
	int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
	int nkeys;
	
	//
	// list fits file herader
	//
	list_fits_map_header(fileName);
	//
	// READING MAP PARAMETERS
	//
	fits_open_file(&fptr, fileName.c_str(), READONLY, &status);
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	long ns;
	ns=get_fits_map_nside(&fptr);
	if (ns>0) set_nside(ns);
	mapOrderings ord;
	ord=get_fits_map_ordering(&fptr);
	if (ord>=0) set_ordering(ord);
	
	
	//
	// READING MAP DATA
	//
	
	for (long i=0;i<dstStructures.size();i++) {
		if (dstStructures[i]=="T") { T()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.T=true; set_pixNum(T().size()); }
		if (dstStructures[i]=="Q") { Q()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.Q=true; set_pixNum(Q().size()); }
		if (dstStructures[i]=="V") { V()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.V=true; set_pixNum(V().size()); }
		if (dstStructures[i]=="U") { U()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.U=true; set_pixNum(U().size()); }
		if (dstStructures[i]=="m") { m()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.m=true; set_pixNum(m().size()); }
		//    if (dstStructures[i]=="n") { n()=get_fits_map_data_column(&fptr, colNames.at(i),&status); }
		if (dstStructures[i]=="N") { N()=get_fits_map_data_column(&fptr, colNames.at(i),&status); loaded.N=true; set_pixNum(N().size()); }
	}
	
	
	// if (data=="mask") {
	//   m()=get_fits_map_data_column(&fptr, "N_OBS",&status);
	//   if (status == 0) msgs->error("N_OBS column has been succesfuly read",High);
	//   else msgs->error("could not read in the mask",High);
	// }
	// if (data=="TEMPERATURE") {
	//   T()=get_fits_map_data_column(&fptr, data,&status);
	//   if (status == 0) msgs->error(data+" column has been succesfuly read",High);
	//   else msgs->error("could not read in the "+data+" column",High);
	// }
	// if (data=="TEMPERATURE") {
	//   T()=get_fits_map_data_column(&fptr, data,&status);
	//   if (status == 0) msgs->error(data+" column has been succesfuly read",High);
	//   else msgs->error("could not read in the "+data+" column",High);
	// }
	
	//       if (how == 1) { makekill_space_manager("make","T",1); }// loading temperature map from I map
	//       if (how == 12 || how == 120) { makekill_space_manager("make","T",1); makekill_space_manager("make","N",1); }// loading temperature map from I map
	//       if (how == 14) { makekill_space_manager("make","T",1); makekill_space_manager("make","N",1); }// loading temperature map from IQU map
	//       if (how == 2) makekill_space_manager("make","m",1); // loading temperature map
	// /*       makekill_space_manager("make","T",1); // loading temperature map */
	//       //if (how == 2) { makekill_space_manager("make","mask",1); } // loading mask map
	//       msgs->say("loading temperature vector map from column 1",Medium);
	//       // reading column 2 - temeprature for the WMAP file
	//       fits_open_file(&fptr, whereto.c_str(), READONLY, &status);
	//       if ( fits_get_hdu_num(fptr, &hdunum) == 1 )
	// 	/* This is the primary array;  try to move to the */
	// 	/* first extension and see if it is a table */
	// 	fits_movabs_hdu(fptr, 2, &hdutype, &status);
	//       else
	// 	fits_get_hdu_type(fptr, &hdutype, &status); /* Get the HDU type */
	
	//       if (hdutype == IMAGE_HDU)
	// 	msgs->say("Error: this program only displays tables, not images",High);
	//       else
	//       {
	// 	fits_get_num_rows(fptr, &nrows, &status);
	//         fits_get_num_cols(fptr, &ncols, &status);
	// 	array = new double[nrows];
	// 			   //val = value;
	// 	//msgs->say("lastcol = %li rows=%li",ncols,nrows);
	// 	if (how == 12 || how==14 ) fits_read_col(fptr,TDOUBLE,1,1,1,nrows, NULL,array, NULL, &status); if (status != 0) { exit(0); }
	// 	if (how == 120 ) fits_read_col(fptr,TFLOAT,1,1,1,nrows, NULL,array, NULL, &status); if (status != 0) { exit(0); }
	// 	if (how == 12 || how==14 || how == 120) for (jj=0;jj<nrows;jj++) { map->T[jj]=(double)(array[jj])/1000.0; /* map->m[jj]=1.0; */}//msgs->say("array[%li]=%lE\n" ,jj,array[jj]); } // the WMAP data was released in mK thermodynamic temperature
	// 	//if (how == 2) for (jj=0;jj<nrows;jj++) { map[jj].m=array[jj]; }
	// 	if ((how == 12)||(how==2)) fits_read_col(fptr,TDOUBLE,2,1,1,nrows, NULL,array, NULL, &status); if (status != 0) { exit(0); }
	// 	if (how == 120) fits_read_col(fptr,TFLOAT,2,1,1,nrows, NULL,array, NULL, &status); if (status != 0) { exit(0); }
	// 	if (how == 14) fits_read_col(fptr,TDOUBLE,4,1,1,nrows, NULL,array, NULL, &status); if (status != 0) { exit(0); }
	// 	if (how == 12 || how==14 || how ==120) for (jj=0;jj<nrows;jj++) { map->N[jj]=(double)(array[jj]); }//msgs->say("array[%li]=%lE\n" ,jj,array[jj]); }
	// 	if (how == 2) for (jj=0;jj<nrows;jj++) { map->m[jj]=array[jj]; }
	// 	fits_close_file(fptr, &status);
	// 	msgs->say("done",Medium);
	// 	delete [] array;
	// 	if (how == 12 || how==14 || how==120) { map_loaded = 1; Nobs_loaded = 1; calculate_map_stats(1); tmps=whereto; setLoadedMapFileName(tmps); }
	// 	if (how == 2) { mask_loaded = 1; tmps=whereto; setLoadedMaskFileName(tmps); }
	// 	check_mask();
	//       }
	//     //--------------------------------------------------------------------------
	//     }
	
	return long(status);
	
}

long mscsMap::loadfits(const string& fileName, const string colName, const string dstStructure) {
	fitsfile *fptr;
	int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
	int nkeys;
	
	msgs->say("Loading fits file: "+fileName,High);
	
	//
	// list fits file herader
	//
	list_fits_map_header(fileName);
	//
	// READING MAP PARAMETERS
	//
	fits_open_file(&fptr, fileName.c_str(), READONLY, &status);
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	long ns;
	ns=get_fits_map_nside(&fptr);
	msgs->say("nside is: "+msgs->toStr(ns),High);
	
	if (ns>0) set_nside(ns);
	mapOrderings ord;
	ord=get_fits_map_ordering(&fptr);
	if (ord>=0) set_ordering(ord); else {   msgs->warning("The file doesn't indicate supported (or any) ordering.",High);  }
	
	
	if (dstStructure=="T") { T()=get_fits_map_data_column(&fptr, colName,&status); loaded.T=true; set_pixNum(T().size()); }
	if (dstStructure=="Q") { Q()=get_fits_map_data_column(&fptr, colName,&status); loaded.Q=true; set_pixNum(Q().size()); }
	if (dstStructure=="V") { V()=get_fits_map_data_column(&fptr, colName,&status); loaded.V=true; set_pixNum(V().size()); }
	if (dstStructure=="U") { U()=get_fits_map_data_column(&fptr, colName,&status); loaded.U=true; set_pixNum(U().size()); }
	if (dstStructure=="m") { m()=get_fits_map_data_column(&fptr, colName,&status); loaded.m=true; set_pixNum(m().size()); }
	//  if (dstStructure=="n") { n()=get_fits_map_data_column(&fptr, colName,&status); }
	if (dstStructure=="N") { N()=get_fits_map_data_column(&fptr, colName,&status); loaded.N=true; set_pixNum(N().size()); }
	
	return long(status);
}
/***************************************************************************************/
long mscsMap::loadPLANCK_temp_fits(const string& fileName, string hduName) {
	cpedsFits ff;
	ff.openFile(fileName,"read");
//	ff.selectHDU("FULL SKY MAP","binTable");
	long n,cols;
	double* data;
	ff.readBinTableColumn(hduName,1,&data,&n,&cols);
	msgs->say("read %li rows",n,Medium);
	set_pixNum(n);
//	set_nside(sqrt(n)/12);
	makekill_space_manager("make","T");
	for (long i = 0; i < n; i++) {
		T(i)=data[i];
	}
	loaded.T=true;

}

void mscsMap::list_fits_map_header(string fileName) {
	
	fitsfile *fptr;
	int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */
	int nkeys;
	char card[FLEN_CARD];
	char keyname[FLEN_KEYWORD];
	string str;
	char colname[FLEN_VALUE], coltype[FLEN_VALUE], colunit[FLEN_VALUE];
	int  hdutype, ncols;//, anynul, dispwidth[1000];
	int i;
	long  nrows;
	long naxes[10];
	
	int hdupos, bitpix, naxis;
	//double temp;
	//filenamestr dupa1,dupa2;
	
	
	msgs->say("reading fits file",High);
	msgs->say("file contents:",Medium);
	
	msgs->println("Main header",Medium);
	fits_open_file(&fptr, fileName.c_str(), READONLY, &status);
	// listing stuff
	fits_get_hdrspace(fptr, &nkeys, NULL, &status);
	
	for (i=1; i<=nkeys; i++)  {
		fits_read_record(fptr, i, card, &status); /* read keyword */
		msgs->println(msgs->toStr(card),Medium);
	}
	msgs->println("END",High);  /* terminate listing with END */
	msgs->println("",High);
	fits_close_file(fptr, &status);
	
	msgs->say("reading other headers",High);
	// listing stuff
	fits_open_file(&fptr, fileName.c_str(), READONLY, &status);
	
	// for (; !status; hdupos++) {  /* Main loop for each HDU */
	hdupos=0;
	while (status==0) {
		fits_get_hdu_type(fptr, &hdutype, &status);  /* Get the HDU type */
		
		msgs->println("HDU #"+msgs->toStr(hdupos),Medium);
		if (hdutype == IMAGE_HDU) {  /* primary array or image HDU */
			fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);
			
			msgs->println("Array:  NAXIS = "+msgs->toStr((long)naxis)+"  BITPIX = "+msgs->toStr((long)bitpix),Medium);
			for (i=0;i<naxis;i++)
				msgs->println("   NAXIS"+msgs->toStr((long)(i+1))+" = "+msgs->toStr(naxes[i]),Medium);
		}
		else {  /* a table HDU */
			fits_get_num_rows(fptr, &nrows, &status);
			fits_get_num_cols(fptr, &ncols, &status);
			
			if (hdutype == ASCII_TBL)
				msgs->println("ASCII Table:  ",Medium);
			else
				msgs->println("Binary Table:  ",Medium);
			
			msgs->println(msgs->toStr((long)ncols)+" columns x "+msgs->toStr((long)nrows)+" rows",Medium);
			msgs->println(" COL NAME             FORMAT",Medium);
			
			char tmpch[100];
			
			for (i=1;i<=ncols;i++)
			{
				str="TTYPE"+msgs->toStr(i); //str.copy(keyname,str.length());
				strcpy(tmpch,str.c_str());
				fits_read_key(fptr, TSTRING, tmpch, colname, NULL, &status);
				// fits_make_keyn("TTYPE", i, keyname, &status); /* make keyword */
				// fits_read_key(fptr, TSTRING, keyname, colname, NULL, &status);
				
				// str="TFORM"+msgs->toStr(i); str.copy(keyname,str.length());
				str="TFORM"+msgs->toStr(i); //str.copy(keyname,str.length());
				strcpy(tmpch,str.c_str());
				fits_read_key(fptr, TSTRING, tmpch, coltype, NULL, &status);
				// fits_make_keyn("TFORM", i, keyname, &status); /* make keyword */
				// fits_read_key(fptr, TSTRING, keyname, coltype, NULL, &status);
				
				str="TUNIT"+msgs->toStr(i); //str.copy(keyname,str.length());
				strcpy(tmpch,str.c_str());
				fits_read_key(fptr, TSTRING, tmpch, colunit, NULL, &status);
				// str="TUNIT"+msgs->toStr(i); str.copy(keyname,str.length());
				// fits_make_keyn("TUNIT", i, keyname, &status); /* make keyword */
				// fits_read_key(fptr, TSTRING, keyname, colunit, NULL, &status);
				
				printf(" %3d %-16s %-16s %-16s\n", i, colname, coltype, colunit);
			}
		}
		
		
		fits_movrel_hdu(fptr, 1, NULL, &status);  /* try move to next ext */
	}
	
	if (status == END_OF_FILE) status = 0; /* Reset normal error */
	fits_close_file(fptr, &status);
	
	
}



cpedsList<double> mscsMap::get_fits_map_data_column(fitsfile **fptr, string colName, int * result) {
	int hdunum, hdutype;
	int status = 0;
	cpedsList<double> cl;
	if ( fits_get_hdu_num(*fptr, &hdunum) == 1 )
		/* This is the primary array;  try to move to the */
		/* first extension and see if it is a table */
		fits_movabs_hdu(*fptr, 2, &hdutype, &status);
	else
		fits_get_hdu_type(*fptr, &hdutype, &status); /* Get the HDU type */
	
	if (hdutype == IMAGE_HDU) {
		printf("Error: this program only displays tables, not images\n");
		*result=-1;
		return cl;
	}
	
	
	int ncols;
	char colname[FLEN_VALUE];
	char tmpch[100];
	// char keyname[FLEN_KEYWORD];
	string keyname;
	msgs->say("reading data column by name:"+colName,Medium);
	fits_get_num_cols(*fptr, &ncols, &status);
	msgs->say("the column number is: "+msgs->toStr(ncols)+ " status is:" + msgs->toStr(status),Medium);
	// skim through columns to find the one with the right name
	for (long i=1;i<=ncols;i++) {
		keyname="TTYPE"+msgs->toStr(i);
		strcpy(tmpch,keyname.c_str());
		fits_read_key(*fptr, TSTRING, tmpch, colname, NULL, &status);
		string s=colname;
		msgs->say("get_fits_map_data_column: colName is: '"+s+"'",High);
		if (s==colName) {
			*result=status;
			return get_fits_map_data_column(fptr,i,result);
		}
	}
	msgs->warning("The file has no column by name: "+msgs->toStr(colname)+". Not reading.",High);
	*result=status;
	return cl;
}

cpedsList<double> mscsMap::get_fits_map_data_column(fitsfile **fptr, long colno, int* result) {
	int hdunum, hdutype;
	int status = 0;
	if ( fits_get_hdu_num(*fptr, &hdunum) == 1 )
		fits_movabs_hdu(*fptr, 2, &hdutype, &status);
	else
		fits_get_hdu_type(*fptr, &hdutype, &status);
	
	if (hdutype == IMAGE_HDU) {
		printf("Error: this program only displays tables, not images\n");
		*result=-1;
		cpedsList<double> cl;
		return cl;
	}
	
	long nrows;
	
	msgs->say("reading data column by index"+msgs->toStr(colno),Medium);
	fits_get_num_rows(*fptr, &nrows, &status);
	double* array = new double[nrows];
	
	fits_read_col(*fptr,TDOUBLE,colno,1,1,nrows, NULL,array, NULL, &status);
	cpedsList<double> cl(array,nrows);
	
	delete [] array;
	
	
	// check units
	char keyname[FLEN_KEYWORD];
	char colunit[FLEN_VALUE];
	string str;
	str="TUNIT"+msgs->toStr(colno); str.copy(keyname,str.length());
	fits_read_key(*fptr, TSTRING, keyname, colunit, NULL, &status);
	if (status!=0) msgs->warning("the column has no unit key",High);
	else {
		string u=colunit;
		if (u=="mK") { cl*=double(1e-3); } else
			if (u=="uK") { cl*=double(1e-6); } else
				if (u!= "K") {   msgs->warning("unknown unit",High); }
		
	}
	
	*result=status;
	return cl;
}

long mscsMap::get_fits_map_nside(fitsfile **fptr) {
	int hdutype;
	int status = 0;
	fits_movabs_hdu(*fptr, 1, &hdutype, &status);
	char nsidec[10];
	char comment[80];
	long ns;
	char keyname[FLEN_KEYWORD]="NSIDE";
	
	fits_read_keyword(*fptr, keyname,nsidec,comment,&status);
	if (status==0) { errno=0; ns=strtol(nsidec,NULL,10);
	if (errno!=0) { msgs->error("wrong nside number",High); return -1; }
	else {
		msgs->say("nside comment is: "+msgs->toStr(comment),High);
	}
	}
	else {
		msgs->error("NSIDE keyword not found, not reading",High);
		return -1;
	}
	return ns;
}

mscsMap::mapOrderings mscsMap::get_fits_map_ordering(fitsfile **fptr) {
	int hdutype;
	int status = 0;
	fits_movabs_hdu(*fptr, 1, &hdutype, &status);
	char ordc[20];
	char comment[80];
	string ord;
	char keyname[FLEN_KEYWORD]="ORDERING";
	
	fits_read_keyword(*fptr, keyname, ordc,comment,&status);
	if (status==0) {
		ord=ordc;
		msgs->say("ORDERING keyword is:"+ord,High);
		if (ord.find("NESTED")!=string::npos) return nested;
		if (ord.find("RING")!=string::npos) return ring;
	}
	else {
		msgs->error("ORDERING keyword not found, assuming NESTED ordering",High);
		return nested;
	}
	return unknownMapOrdering;
}



long mscsMap::savefits(string fileName) {
	
	return -1;
}







