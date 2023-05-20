#include <stdlib.h>
#include <stdio.h>
//#ifdef withPGPLOT
#include <cpgplot.h>
//#endif
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <tclap/CmdLine.h>
//#include <tclap/ValuesConstraint.h>
// #include "cpeds-consts.h"
#include "cpeds-math.h"
#include "cpeds-msgs.h"
#include "cpeds-direction.h"
#include "Mscs-global-defs.h"
// #include "Mscs-map.h"
#include "Mscs-map_plot.h"
#include "Mscs-alms.h"
#include "Mscs-function.h"

#include "Mscsgrid.h"

#ifdef GOMPI
#include "mpi.h"
#endif

// #define _REENTRANT
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

//declaration of the global variables
void parseOptions(int argc, char** argv);
vector<string>	_circT_files,_circm_files;
vector<long> ff;
double _eclipticVal, _equatorVal;
filenamestr _input_file_my;
string _input_file,_output_file,_mask_file,_proj,_ft,_ord,_color_scheme,_multi_mask,_cmf[5],_cmf_U,_ctf[5],_ctf_U,_overplot_region_type,_save_overplot_as,_mytitle,_bgcolor,_fgcolor,_colorbar,_colorbar_label,_labels_format, _plot_taken_vals, _BGColor;
double _Zrange[2],_Zrangelog[2],_zrange,_zrange_min,_zrange_max,_delta_latitude,_delta_longitude;
double _camera[3], _mm_Rl,_mm_Rth,_mm_Rphi, _Rl,_Rth,_Rphi,_ml[100],_mb[100],_ms[100],_tl[100],_tb[100],_ts[100],_tv[100],_cmf_r,_cmf_r0,_ctf_r,_ctf_r0,_ctf_v, _reim_ratio,_PGwindow_size, _find_circles_r;
double _rotGmphi,_rotGmth, _rotGsphi,_rotGsth, _l0,_b0;
long _size,_sizeX,_sizeY,_mlevels,_nside=0,_nsidedef=512, _renside=0, _color_num,_lreg_num=0,_breg_num=0,_mmnside,_resX,_resY,_cmf_size,_ctf_size,_ctf_cols,_dpi, _Yl, _Ym,_colorbar_fontsize,_meridians_fontsize;
bool _col_ord,_cyclic_point,_reverse_l,_nopyth,_pyth,_show_ecliptic, _show_equator,_negate_r, _dont_save_mask;
bool _default_Z=true,_default_Zlog=true, default_cam=true,masked=false,_mink=false, _not_from_file=false,multi_masked=false,_randomize,_random_orientation,_mm_random_orientation,_logscale,_abs,_norm_by_max,_mask_from_here,_plotNS,_dont_plot_zero, _meridians, _plot_taken_valsB=false,_save_rot,_load_rot,_do_special_block,_find_circles, _mk_gal_symm,_mk_gal_asymm,_mk_ew_symm, _mk_ew_asymm, _mk_pt_symm,_mk_pt_asymm,_mk_ring_av,_rot_short,_ctf_add_vals,_invert_mask, _mk_thres_map, _area,_trecl, _treq;
bool _isSetRth,_isSetRphi,_isSetRl, _random_orientation_sph, _rotate_random_gauss_lb;
cpeds_queue<double> _th,_trl,_trb,_trs,_trv;
cpeds_queue<long> *cq;
cpeds_queue<long> *circle_points;
long _overplot_region_type_dot_points,_ravns, _seed_offset;
cpeds_VerbosityLevel _verbosity=CPEDS_defaultVerbosityLevel;
//-----------------------------------


int main(int argc, char **argv) {
  cpedsMsgs msgs("draw_maps");
  msgs.setVerbosity(_verbosity);
  // float * tab;
  double l,b,dl,db,l0,*tmpmap,T,Tth;
  cpedsDirection n;
  string file2,tmpch,rotfile;
  filenamestr commandstr;
  long i,j,k,m,nside_orig,reg_num,ret_code;
  int output=0,meridians=0,plotNS=0;
  string fname,histname_file,operation;
  FILE *f,*f1,*f2;
  bool onrlist;
  long p,*tabint;//* ptab;
  double pth,psum,th,phi;
  cpeds_point* cppt;
  mscsMap::mapOrderings mapordering=mscsMap::nested;
  // typedef struct cmfT{
  //   double l,b; // galactic coordinates
  //   double r; // circular mask radius
  //   cmfT *Pnext;
  // };
  // cmfT *cmfB,*cmfn;

  // typedef struct ctfT{
  //   double l,b; // galactic coordinates
  //   double r; // circular mask radius
  //   double v; // value in the temperature map
  //   ctfT *Pnext;
  // };
  // ctfT *ctfB,*ctfn;


#ifdef GOMPI
  MPI_Status mpi_stat;
  int ierr,node_num,world_id=0;
  long Ndirpn;
#endif


  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  Mscs_initiate_global_variables();
  mscsMapPlot map("map"); 
  mscsMap mask("mask");
  mscsAlms a("alms");
  parseOptions(argc,argv);
  map.setVerbosityLevel(High);
  mask.setVerbosityLevel(High);
  //----------------------------------------------------------------------------------------------------
  srand((long)time(NULL)-cpeds_seed_offset);

/*   strcpy(mappref, ARGV[1]); */
/*   mapordering = (int)(strtol(ARGV[2],NULL,10)); */
/*   projection = (int)(strtol(ARGV[3],NULL,10)); */

  //if ((ARGC != 10) && (ARGC != 5)) { print_program_usage(); exit(1);}
/*  if (ARGC == 5) { */
/*     strcpy(txtbin,ARGV[4]); */
/*     //if (strcmp(txtbin,"") == 0) { strcpy(txtbin,"bin"); } */
/*   } */

/*   if (ARGC == 11) { */
/*     start = (int)(strtol(ARGV[4],NULL,10)); */
/*     stop = (int)(strtol(ARGV[5],NULL,10)); */
/*     step = (int)(strtol(ARGV[6],NULL,10)); */
/*     minka = (int)(strtol(ARGV[7],NULL,10)); */
/*     minkc = (int)(strtol(ARGV[8],NULL,10)); */
/*     minks = (int)(strtol(ARGV[9],NULL,10)); */
/*     strcpy(txtbin,ARGV[10]); */
/*   } */
/*   else { start = 1; stop = 1; step = 1; minka = minkc = minks = 1; } */

/*   sim1.minkowski_level_num_circ = 100; */
/*   for (i=start;i<=stop;i=i+step) { */
  //sprintf(out,"simulated_map_sum_lmax%i",i);
/*     if (ARGC == 11) { sprintf(out,"%s%i",mappref,i); } else  { sprintf(out,"%s",mappref); } */
/*     if (strcmp(txtbin,"bin") == 0) { sim1.loadbinT(out,mapordering); } */
/*     if (strcmp(txtbin,"txt") == 0) { sim1.loadtxtT(out,mapordering); } */
  if (_ord == "nest") { mapordering = mscsMap::nested; }
  if (_ord == "ring") { mapordering = mscsMap::ring; }

/*   if (_not_from_file) { // we do not plot any particular file */
/*     if (multi_masked) { // plotting requested multi mask */
/*       map.set_map_nside(256); map.makekill_space_manager("make","T",1); map.clear_map();    */
/*       if (_multi_mask == "LB") { */
/* 	map.make_multi_maskLB(_lreg_num,_breg_num,"test"); */
/* 	reg_num = 4*_lreg_num * _breg_num;  */
/*       } */
/*       if (_multi_mask == "HP") { */
/* 	map.make_multi_maskHP(_mmnside,"test"); */
/* 	reg_num = 12*_mmnside*_mmnside; */
/*       } */
/*       map.mask2map();  // make temperature map from mask */
/*       //for (i=1;i<=reg_num;i++) map.set_Treg(i,(double)(reg_num)); */
/*       for (i=0;i<ffsize;i++) { map.set_Treg(ff[i],(double)(2*reg_num)); } */
/*       map.makekill_space_manager("kill","m",1);  */
/*     } */
/*     else { printf("  -- no map to draw\n"); exit(0); } // we do not plot anything */
/*   } */
/*   else { // plotting map from file */
/*     if (_ft == "bin") { map.loadbinT(_input_file.c_str(),mapordering); } */
/*     if (_ft == "txt") { map.loadtxtT(_input_file.c_str(),mapordering); } */
/*   } */

  if (_not_from_file==false) { // we do  plot from  file
    if (_ft == "bin") { map.loadbinT(_input_file.c_str(),mapordering); }
    if (_ft == "txt") { map.loadtxtT(_input_file.c_str(),mapordering); }
    if (_ft == "fitsWMAP") { map.loadfits(_input_file.c_str(),"TEMPERATURE","T"); }
    if (_ft == "fitsPL") { map.loadfits(_input_file.c_str(),"I_Stokes","T"); }
  }


/* ******** EXPERIMENTAL - testing the ang2pix bug 2008-01-11 - BEGIN **************** */

/*   map.set_map_coord(0,0); */
/*   n.l=1.5707963268; n.b=-1.2823667021 ;  map.set_T(n,0);   printf("in: l:%lE b:%lE    |||   ",n.l*PI180inv,n.b*PI180inv);   // RESULT : THIS WORKS FINE */
/*   n.l=90.0*PI180; n.b=-72.95887684680076938*PI180;  map.set_T(n,0);   printf("in: l:%lE b:%lE    |||   ",n.l*PI180inv,n.b*PI180inv);   // RESULT : THIS THIS DOESN'T WORK */
/*   n.l=PI/2.0; n.b=-72.95887684680076938*PI180;  map.set_T(n,0);   printf("in: l:%lE b:%lE    |||   ",n.l*PI180inv,n.b*PI180inv);   // RESULT : THIS THIS DOESN'T WORK */

/*   n.l=PIsnd; n.b=-72.95887684680076938*PI180;  map.set_T(n,0);   printf("in: l:%lE b:%lE    |||   ",n.l*PI180inv,n.b*PI180inv);   // RESULT : THIS THIS  WORKS FINE */

/*   map.make_circle_dot(_tl[1],_tb[1],_ts[1],_tv[1],"T",_overplot_region_type,_overplot_region_type_dot_points);  // RESULT : THIS DOESN'T WOR */
/*   n.l=90; n.b=-72.95887684680076938;  map.make_circle_dot(n.l,n.b,_ts[1],_tv[1],"T",_overplot_region_type,_overplot_region_type_dot_points);  // RESULT : THIS DOESN'T WORK */

/*   for (i=0;i<map.pix_num;i++) if (map.get_T(i) == 0) n=map.get_C(i); */
/*   printf("out:   l:%lE b:%lE\n",n.l*PI180inv,n.b*PI180inv); */
/*   exit(0); */


/* *****EXPERIMENTAL - END *********************************** */


/* **************************************** */
  if (_do_special_block) {
    printf("DOING EXPERIMENTAL SPECIAL BLOCK - START\n");

    mscsGrid grid(&map,"grid");
    grid.makeGrid(mscsGrid::grid_ESLonLat, 256L, 256L);
    cpeds_matrix_save(grid.field(),"field.mat","");
    mscsMap map=grid.toMap(64);
    map.savebinT("field-Tn-bin");
    exit(0);
/*     map.set_map_nside(512); */
/*     map.set_map_coord(0,0); */
/*     map.makekill_space_manager("make","T",1); */
/*     for (i=0;i<map.pix_num;i++) { */
/*       n=map.get_C(i); */
/*       map.set_T(i,sin(_nside*(PIsnd-n.b))); */
/*     } */
/*     sprintf(tmpch,"sin-digit-n%li",_nside); */
/*     map.savebinT(tmpch,1); */

/*     for (i=0;i<map.pix_num;i++) { */
/*       n=map.get_C(i); */
/*       map.set_T(i,_nside*cos(_nside*(PIsnd-n.b))+2*_nside); */
/*     } */
/*     sprintf(tmpch,"cos-digit-n%li",_nside); */
/*     map.savebinT(tmpch,1); */


    // map.set_nside(512);
    // map.set_map_coord();
    // map.makekill_space_manager("make","T",1);
    // for (i=0;i<map.pix_num;i++) {
    //   n=map.get_C(i);
    //   map.set_T(i,sin(_nside*(n.l)));
    // }
    // sprintf(tmpch,"sinl-digit-n%li",_nside);
    // map.savebinT(tmpch,1);

    // for (i=0;i<map.pix_num;i++) {
    //   n=map.get_C(i);
    //   map.set_T(i,_nside*cos(_nside*(n.l))+2*_nside);
    // }
    // sprintf(tmpch,"cosl-digit-n%li",_nside);
    // map.savebinT(tmpch,1);

    // printf("DOING EXPERIMENTAL SPECIAL BLOCK - END\n");
    // exit(0);
  }
/* **************************************** */





  //--------------------------------------------------------------------------------
  if (_nside == 0 && _not_from_file) _nside = _nsidedef; //
  if (_not_from_file) { map.set_nside(_nside); map.makekill_space_manager("make","T",1); map.clear_map(); }

  // ----------------------------------------------------------------------------------------------------
  msgs.say("zeroing regions", High);
  if (_negate_r) {
    k=(long)map.getMinT();
    while (k<(long)map.getMaxT()) {
      onrlist=false;
      for (i=0;i<ff.size();i++) { if (k == ff[i] && onrlist==false) onrlist=true;  }  // leaving only -r regions and removing all else
      if (!onrlist) map.set_Ttoval((double)k,0.0);
      k++;
    }
  }
  else for (i=0;i<ff.size();i++) { map.set_Ttoval((double)(ff[i]),0.0); }

  // ----------------------------------------------------------------------------------------------------
  if (multi_masked) { // plotting requested multi mask
    if (_mm_random_orientation) {
      _mm_Rth = cpeds_random_uniform_number(0,180); _mm_Rphi = cpeds_random_uniform_number(0,180);  _mm_Rl = cpeds_random_uniform_number(0,180);
      printf("  -- random orientations requested: mm_Rth = %lE, mm_Rphi = %lE, mm_Rz = %lE\n",_mm_Rth,_mm_Rphi,_mm_Rl);    }
    if (_multi_mask == "LB") {
      map.make_multi_maskLB(_lreg_num,_breg_num,"test",_mm_Rth,_mm_Rphi,_mm_Rl,true);
      //reg_num = _lreg_num * _breg_num;
    }
    if (_multi_mask == "HP") {
      map.make_multi_maskHP(_mmnside,"test",_mm_Rth,_mm_Rphi,_mm_Rl,true);
      //reg_num = 12*_mmnside*_mmnside;
    }
/*     for (i=0;i<ffsize;i++) { map.set_Treg(ff[i],0); } */
    if (_not_from_file) {
      map.mask2map();  // make temperature map from mask
      //for (i=1;i<=reg_num;i++) map.set_Treg(i,(double)(reg_num));
      //printf("******debug mm regnum : %li\n",map.multi_mask_reg_num);
      if (_randomize) {
    	  for (i=1;i<=map.multi_mask_reg_num();i++) { map.set_Treg(i,cpeds_random_uniform_number(1,1000)); }
      }
      printf(" draw_maps: *************************** WARNING LOOK IN THE SOURCE BEFORE USING IT *********************************\n");
      for (i=0;i<ff.size();i++) { map.set_Treg(ff[i],(double)(2*reg_num)); } // this is weired
      map.makekill_space_manager("kill","m",1);
    }
  }

  //--------------------------------------------------------------------------------

  if (_input_file == "Ylm") { // plotting SH multiplied by a_lm given from command line;
    msgs.say("generating map of Spherical Harmonics",Top);
    a.lmax(_Yl);
    a.zero_alms_multipole_range(0,_Yl);
    a.set(_Yl,_Ym,1.0,_reim_ratio);
    a.antysymmetrize();
    if (_Ym==0) a.set(_Yl,_Ym,1.0,_reim_ratio);
/*     printf("C_0 = %lE\n",map.calculate_C_l(0)); */
/*     printf("C_1 = %lE\n",map.calculate_C_l(1));    */
    map.SH_synthesis(a,_Yl);
    map.calculate_map_stats(1);
  }


//--------------------------------------------------------------------------------

  // make the map symmetries if requested
  if (_mk_gal_symm || _mk_gal_asymm) {
    msgs.say("making map symmetrical about galactic plane",Top);
    j=map.pixNum();
    if (!map.coordLoaded()) map.set_map_coord();
    for (i=0;i<j;i++) {
      // make the map symmetrical w.r.t galactic plane
      n=map.get_C(i);
      if ((n.b()) >= 0) {
	n.lat()=-n.b(); if (_mk_gal_symm) map.set_T(n,map.get_T(i)); else  map.set_T(n,-map.get_T(i));
      }
    }
  }

//--------------------------------------------------------------------------------
  if (_mk_ew_symm || _mk_ew_asymm) {
    msgs.say("making map symmetrical about zero meridian",Top);
    j=map.pixNum();
    if (!map.coordLoaded()) map.set_map_coord();
    for (i=0;i<j;i++) {
      // make the map symmetrical w.r.t galactic plane
      n=map.get_C(i);
      if (n.l() < PI) {
	n.lon()=twoPI-n.l();
	if (_mk_ew_symm) map.set_T(n,map.get_T(i)); else  map.set_T(n,-map.get_T(i));
      }
    }
  }

//--------------------------------------------------------------------------------
  if (_mk_pt_symm || _mk_pt_asymm) {
    msgs.say("making point symmetry",Top);
    if (_dont_plot_zero) j=map.pixNum(); else j=map.pixNum()/2;
/*     j=map.pix_num; // temporary change =--- testing cpeds */
    if (!map.coordLoaded()) map.set_map_coord();
    // make the map point-symmetrical
    for (i=0;i<j;i++) {
      n=map.get_C(i); n.lat()=-n.b(); n.lon()+=PI; if (n.l() >= twoPI) n.lon()-=twoPI;
      if (_mk_pt_symm) { if (_dont_plot_zero) { if (map.get_T(i)!=0) map.set_T(n,map.get_T(i)); } else map.set_T(n,map.get_T(i)); }
      else { if (_dont_plot_zero) { if (map.get_T(i)!=0) map.set_T(n,-map.get_T(i)); } else map.set_T(n,-map.get_T(i)); }
    }
  }


//--------------------------------------------------------------------------------

  if (masked) {   // load mask if required
    printf("|DRAW_MAPS> LOADING MASK\n");
    if (!_mask_from_here) file2=MSCS_WMAP_LOCAL_DIR+_mask_file; else file2=_mask_file; fname=file2; // fname=_mask_file;
    if (fname.find(".fits",fname.size()-5) != string::npos) { //fname.erase(fname.size()-7,7);
      mask.loadfits(file2,"N_OBS","m");  } // assume that it's fits if there is no -bin extension
    else
      mask.loadbinm(fname.c_str());    //map2.change_map_resolution(map.nside);
      //if (fname.find(".fits",fname.size()-5) != string::npos) map2.loadfitsT(file2,2);  //map2.change_map_resolution(map.nside);
    if (_invert_mask) mask.invert_mask();
    map.import_map_data(mask,"m",2);
    map.calculate_map_stats(1);
    map.check_mask();    map.mask_map_merge();
    map.calculate_map_stats(1);
    if (!multi_masked && _not_from_file && _input_file!="Ylm") { map.mask2map(); }
    printf("|DRAW_MAPS> LOADING MASK: DONE\n");
  }


//-----------------------------------------------------------------------------------------------------

  if (_mk_ring_av) {
    msgs.say("averaging map in rings",Top);
    long mapres=map.nside();
    map.change_map_resolution(_ravns);
    map.average_map_in_rings(NULL);
    map.change_map_resolution(mapres);
  }

  if (_renside > 0 ) { nside_orig = map.nside(); map.change_map_resolution(_renside); map.change_map_resolution(nside_orig); }
  if (_nside != 0 ) map.change_map_resolution(_nside);
/*   if (_zrange !=0) { map.calculate_map_stats(1); map.flatmap_change_zrange = true; map.flatmap_minT=-_zrange*sqrt(map.varianceT); map.flatmap_maxT=_zrange*sqrt(map.varianceT); printf("**********%lE %lE\n",map.flatmap_minT,map.flatmap_maxT);} //EXPERIMENTAL !!! TO BE CORRECTED BECAUSE THIS IS LAME !!!! do this by implementing a procedure of scrammbling temperatures exceeding some threshold, and conditionally call this procedure. */
  if (_zrange != 0) { map.scramble_over_nsigma(_zrange); }
  if (_default_Z == false) {  map.scramble_over_minmax(_Zrange[0],_Zrange[1]);
/*     printf("defaultZ false min %lE max %lE\n",_Zrange[0],_Zrange[1]);     exit(0);  */
  }

//-----------------------------------------------------------------------------------------------------
// make a map from a circular regions defined in the text file: l,b,r,v
// r value can be overridden by option _cmf_r
//

  msgs.say("marking circular regions in the temperature map",Top);
  matrix<double> circ_lbrv_tmp, circ_lbrv;
  for (i=0;i<_circT_files.size();i++) {
    msgs.say("making a map from circular regions from file: "+_circT_files[i],High);
    circ_lbrv_tmp=cpeds_matrix_load(_circT_files[i]);
    circ_lbrv.SetSize(circ_lbrv_tmp.RowNo(),4);

    long N=circ_lbrv_tmp.RowNo();
    if (circ_lbrv_tmp.ColNo() == 4) { circ_lbrv=circ_lbrv_tmp; }
    // if (circ_lbrv_tmp.ColNo() == 3) { for (long i=0;i<N;i++) { circ_lbrv(i,0)=circ_lbrv_tmp(i,0); circ_lbrv(i,1)=circ_lbrv_tmp(i,1); circ_lbrv(i,2)=circ_lbrv_tmp(i,2); circ_lbrv(i,3)=_ctf_v; } }
    if (circ_lbrv_tmp.ColNo() == 3) { for (long i=0;i<N;i++) { circ_lbrv(i,0)=circ_lbrv_tmp(i,0); circ_lbrv(i,1)=circ_lbrv_tmp(i,1); circ_lbrv(i,2)=_ctf_r; circ_lbrv(i,3)=circ_lbrv_tmp(i,2); } }
    if (circ_lbrv_tmp.ColNo() == 2) { for (long i=0;i<N;i++) { circ_lbrv(i,0)=circ_lbrv_tmp(i,0); circ_lbrv(i,1)=circ_lbrv_tmp(i,1); circ_lbrv(i,2)=_ctf_r0+_ctf_r; circ_lbrv(i,3)=_ctf_v; } }
    circ_lbrv_tmp.SetSize(0,0);
    if (_ctf_U == "m") {       for (long i=0;i<N;i++) { circ_lbrv(i,2)/=60.0; } }
    if (_ctf_U == "s") {       for (long i=0;i<N;i++) { circ_lbrv(i,2)/=3600.0; } } // conversion of the input data to degrees

    msgs.say("read "+msgs.toStr(long(circ_lbrv.RowNo()))+" new map values positions",Medium);
    if (_ctf_add_vals) operation="add";
    map.draw_circles(circ_lbrv,"T",_overplot_region_type,_overplot_region_type_dot_points,operation);
  }



//-----------------------------------------------------------------------------------------------------
// make circular mask on the map from command line parameters
  if (_ml[0] != 0 && _mb[0] != 0) {
    msgs.say("masking circular regions in the temperature map from cmd parameters",Top);
/*     map.make_multi_maskHP(long(map.nside)/8,"test",0,0,0,true); */
/*     for (i=1;i<=(long)(_ml[0]);i++) {       map.set_Treg(_ml[i]*PI/180,_mb[i]*PI/180,2); } */
/*     map.clear_multimask(); */
    for (i=1;i<=(long)(_ml[0]);i++) {       map.make_circle_dot(_ml[i],_mb[i],_ms[i],0,"m",_overplot_region_type,_overplot_region_type_dot_points,""); }
    if (_invert_mask) map.invert_mask();
    if (!multi_masked && _not_from_file) { map.mask2map(); } else map.mask_map_merge(); //map.makekill_space_manager("kill","m",1);
  }

//-----------------------------------------------------------------------------------------------------
// draw circular region on the map from command line parameters
  if (_tl[0] != 0 && _tb[0] != 0) {
    msgs.say("marking circular regions in the temperature map from cmd parameters",Top);
/*     map.make_multi_maskHP(long(map.nside)/8,"test",0,0,0,true); */
/*     for (i=1;i<=(long)(_ml[0]);i++) {       map.set_Treg(_ml[i]*PI/180,_mb[i]*PI/180,2); } */
/*     map.clear_multimask(); */

    for (i=1;i<=(long)(_tl[0]);i++) {       
    	map.make_circle_dot(_tl[i],_tb[i],_ts[i],_tv[i],"T",_overplot_region_type,_overplot_region_type_dot_points,""); 
    }

/*     EXPERIMENTAL BEGIN */

/*     // check of the spot cut off with zeros from map with ones */
/*     printf("EXPERIMENTAL: performing check\n"); */
/*     map.set_map_coord(0,0); */
/*     f=fopen("errors","w"); */
/*     for (i=0;i<map.pix_num;i++) { */
/*       if (map.get_T(i) == 0) { */
/* 	map.get_C(i,&l,&b); */
/* 	if (cpeds_ang_n1n2(PIsnd-PI180*_tb[1],PI180*_tl[1],PIsnd-b,l) > 50.5*PI180) {fprintf(f,"%li %.10lE %.10lE %.10lE %.10lE\n",i,l,b,l*PI180inv,b*PI180inv); } */
/*       } */
/*     } */
/*     fclose(f); */

/*     EXPERIMENTAL END */

/* ******** EXPERIMENTAL - testing the ang2pix bug 2008-01-11 - BEGIN **************** */

/*     map.set_map_coord(0,0); */
/*     //    n.l=1.5707963268; n.b=-1.2823667021 ;  map.set_T(n,0);   printf("in: l:%lE b:%lE    |||   ",n.l*PI180inv,n.b*PI180inv); */
/*     for (i=0;i<map.pix_num;i++) if (map.get_T(i) == 0) n=map.get_C(i); */
/*     printf("out:   l:%lE b:%lE\n",n.l*PI180inv,n.b*PI180inv); */
/*     exit(0); */
/*   // RESULT : THIS WORKS FINE */

/* *****EXPERIMENTAL - END *********************************** */


    //if (!multi_masked && _not_from_file) { map.mask2map(); } else map.mask_map_merge(); //map.makekill_space_manager("kill","m",1);
  }
//-----------------------------------------------------------------------------------------------------
// draw great ring at given position given value and size
  msgs.say("marking great ring in the temperature map",Top);
  for (i=0;i<_trl.get_size();i++) {
    map.mk_ring(_trl(i),_trb(i),_trs(i),_trv(i),"T");
    map.calculate_map_stats(1);
  }




//-----------------------------------------------------------------------------------------------------
// make a mask from a circular regions defined in the text file: l,b,r
// r value can be overridden by option _cmf_r
  for (i=0;i<_circm_files.size();i++) {
    msgs.say("masking circular regions in the temperature map",Top);
    msgs.say("making a mask from circular regions from file: "+_circm_files[i],High);
    circ_lbrv_tmp=cpeds_matrix_load(_circm_files[i]);
    circ_lbrv.SetSize(circ_lbrv_tmp.RowNo(),4);

    long N=circ_lbrv_tmp.RowNo();
    if (circ_lbrv_tmp.ColNo() == 3 or circ_lbrv_tmp.ColNo() == 4) { for (long i=0;i<N;i++) { circ_lbrv(i,0)=circ_lbrv_tmp(i,0); circ_lbrv(i,1)=circ_lbrv_tmp(i,1); circ_lbrv(i,2)=circ_lbrv_tmp(i,2); circ_lbrv(i,3)=0.0; } }
    if (circ_lbrv_tmp.ColNo() == 2) { for (long i=0;i<N;i++) { circ_lbrv(i,0)=circ_lbrv_tmp(i,0); circ_lbrv(i,1)=circ_lbrv_tmp(i,1); circ_lbrv(i,2)=_ctf_r0+_ctf_r; circ_lbrv(i,3)=0.0; } }
    circ_lbrv_tmp.SetSize(0,0);
    if (_ctf_U == "m") {       for (long i=0;i<N;i++) { circ_lbrv(i,2)/=60.0; } }
    if (_ctf_U == "s") {       for (long i=0;i<N;i++) { circ_lbrv(i,2)/=3600.0; } } // conversion of the input data to degrees

    msgs.say("read "+msgs.toStr(long(circ_lbrv.RowNo()))+" new map values positions",Medium);
    map.draw_circles(circ_lbrv,"m",_overplot_region_type,_overplot_region_type_dot_points,"");
  }

  map.calculate_map_stats(1);
  //  if (!multi_masked && _not_from_file && _circm_files.size()==0) { map.mask2map(); } else  map.mask_map_merge(); //map.makekill_space_manager("kill","m",1);
  if (!multi_masked &&  _circm_files.size()==0) { map.mask_map_merge(); } else map.mask2map();  //map.makekill_space_manager("kill","m",1);
  map.calculate_map_stats(1);
  map.check_mask();
  map.calculate_map_stats(1);


//-----------------------------------------------------------------------------------------------------
// making overplots
  if (_show_ecliptic) {
    msgs.say("marking ecliptic with value "+msgs.toStr(_eclipticVal),High);
    // map.markEcliptic(_eclipticVal);
  }
  if (_show_equator) {
    msgs.say("marking equator with value "+msgs.toStr(_equatorVal),High);
    // map.markEcliptic(_equatorVal);
  }

//-----------------------------------------------------------------------------------------------------
  msgs.say("MAP ROTATIONS",Top);

  if (_random_orientation) {
    if (!_isSetRth) _Rth = cpeds_random_uniform_number(0,180);
    if (!_isSetRphi) _Rphi = cpeds_random_uniform_number(0,180);
    if (!_isSetRl) _Rl = cpeds_random_uniform_number(0,180);
  }
  if (_random_orientation_sph) {
    i=(long)cpeds_random_uniform_number(0,cpeds_get_healpix_pix_num(512)-1);
    cpeds_pix2ang_healpix(512,i,&th,&phi,1);
    if (!_isSetRth) _Rth = PI180inv*th;
    if (!_isSetRphi) _Rphi = PI180inv*phi;
    if (!_isSetRl) _Rl = cpeds_random_uniform_number(0.0,180.0);
  }

  if (_rotate_random_gauss_lb) {
    if (!_isSetRth) _Rth = cpeds_random_gauss_number(_rotGmth,_rotGsth,100,2);
    if (!_isSetRphi) _Rphi = cpeds_random_gauss_number(_rotGmphi,_rotGsphi,100,2);
    if (!_isSetRl) _Rl = cpeds_random_uniform_number(0.0,180.0);
  }


  if (_save_rot) {
    tmpmap=new double[map.pixNum()]; i=0;    for (i=0;i<map.pixNum();i++) { tmpmap[i]=map.get_T(i); map.set_T(i,(double)i); }
  }

if (_load_rot) {
    // check if there's a correct file for rotation
    //yes
    if (_mask_from_here) rotfile=msgs.toStr(map.nside())+"-Rth_"+msgs.toStrf(_Rth)+"-Rphi_"+msgs.toStr(_Rphi)+"-Rl_"+msgs.toStr(_Rl);
    else rotfile=MSCS_DATA_DIR+msgs.toStr(map.nside())+"-Rth_"+msgs.toStrf(_Rth)+"-Rphi_"+msgs.toStr(_Rphi)+"-Rl_"+msgs.toStr(_Rl);
    cpedsList<double> nums;
    // f = fopen(rotfile.c_str(),"r");
    // if (f != NULL) {
    if (nums.load(rotfile)==0) {
      cpedsList<double> T=map.T();
      msgs.say("will rotate according to the saved file: "+rotfile,Medium);
      // tmpmap=new double[map.pixNum()]; i=0;
      //for (i=0;i<map.pixNum();i++) tmpmap[i]=map.get_T(i);
      for (i=0;i<map.pixNum();i++) map.set_T(i,T[nums[i]]);
      // while (fscanf(f,"%li",&i) != EOF) { map.set_T(j,tmpmap[i]); j++; }
      // fclose(f);
      // delete [] tmpmap;
      T.clear(); nums.clear();
    }
    else {
      msgs.say("there's not correct file for this rotation, will do it in usual way",Medium);
      map.rotate_map(_Rth,_Rphi,_Rl,_rot_short,"T");
    }
  }
  else
    map.rotate_map(_Rth,_Rphi,_Rl,_rot_short,"T");

  if (_save_rot) {
    if (_mask_from_here) rotfile=msgs.toStr(map.nside())+"-Rth_"+msgs.toStrf(_Rth)+"-Rphi_"+msgs.toStr(_Rphi)+"-Rl_"+msgs.toStr(_Rl);
    else rotfile=MSCS_DATA_DIR+msgs.toStr(map.nside())+"-Rth_"+msgs.toStrf(_Rth)+"-Rphi_"+msgs.toStr(_Rphi)+"-Rl_"+msgs.toStr(_Rl);
    msgs.say("will save the rotation information to file: "+rotfile,Medium);
    // save file
    f = fopen(rotfile.c_str(),"w"); for (i=0;i<map.pixNum();i++) fprintf(f,"%.0lf\n",map.get_T(i)); fclose(f);
    // rotate the map
    for (i=0;i<map.pixNum();i++) map.set_T(i,tmpmap[(long)map.get_T(i)]);
  }

//-----------------------------------------------------------------------------------------------------
  if (_abs)   map.mk_abs();
  if (_logscale)   map.logarithmT(10);
  if (_logscale && _default_Zlog == false) {  map.scramble_over_minmax(_Zrangelog[0],_Zrangelog[1]); }

  if (_norm_by_max) { map.norm_by_maxT(); }

//-----------------------------------------------------------------------------------------------------
  if (_mk_thres_map) {
    msgs.say("making thresholded map",Top);

/*     map.conv_nest2ring(map.map); */
/*     for (j=0;j<map.pix_num;j++) { map.map->T[j]=(double)j; } */
/*     map.conv_ring2nest(map.map); */

    map.calculate_map_stats(1);
    map.makekill_space_manager("make","N");
    _th.sort(12);
    msgs.say("sorting map",High);
    mscsFunction ptab=map.get_value_sorted_pixel_numbers(21);

    //dp=4*PI/(double)map.pix_num;
    msgs.say("masked pix num: "+msgs.toStr(map.maskedPixNum()),Medium);
    psum = map.get_meanT()*(double)(map.pixNum()-map.maskedPixNum()); //get_integralT()/dp;
    printf("psum= %lE \n",psum);
    for (i=0;i<map.pixNum();i++) { /* map.map->T[i] /= psum; */ map.set_N(i,0); }
    map.calculate_map_stats(1);


    for (j=0;j<_th.get_size();j++) {
      p=0; pth=0;  th=_th(j)*psum;
      printf("threshold: %lE corresponding integral: %lE\n",_th(j),th);
      if (_th(j) > 1) { printf("|draw_maps> ERROR: the given conficende threshold is larger than 1: %lE\nShould be in 0..1 range\n",_th(j)); exit(0); }

      while (pth <= th) { pth+=map.get_T(long(ptab.f(p)));  /* printf("pth: %lE pthint: %lE; p=%li ptab[p]=%li\n",th,pth,p,ptab[p]); */ if (p < map.pixNum()) p++; else {printf("p>pix_num\n"); exit(0); }}
      printf("th= %lE p=%li\n",_th(j),p);

      for (i=0;i<p;i++) {
	if (map.get_N(long(ptab.f(i))) == 0) map.set_N(long(ptab.f(i)),_th(j));
      }
    }

    // set the rest as outside of the threshold
    for (i=0;i<map.pixNum();i++) {
      if (map.get_N(i) == 0 && map.get_m(i)!= 0) map.set_N(i,1.0);
    }


    //    delete [] ptab;
    ptab.clearFunction();
    // copy the thresholded map on to the temperature map
    for (i=0;i<map.pixNum();i++) { map.set_T(i,map.get_N(i)); }

  }


// ----------------------------------------------------------------------------------------------------
  if (_find_circles) {
    msgs.say("finding circles in the map",High);

    long ST,EN;

    printf("LOOKING FOR THE HOTTEST POINT IN CIRCLES IN THE SKY \n");

#ifdef GOMPI
      ierr = MPI_Init ( &argc, &argv ); // initiate MPI
      ierr = MPI_Comm_size ( MPI_COMM_WORLD, &node_num ); // Get the number of processes
      ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &world_id ); // Get the individual process ID
#endif

    circle_points = new cpeds_queue<long>;
    k=map.pixNum();
#ifndef GOMPI
      ST=0; EN=k;
#endif
#ifdef GOMPI
      Ndirpn=k/node_num;
      ST=world_id*Ndirpn;
      EN=(world_id+1)*Ndirpn;
      if ( world_id == node_num-1 ) { EN=k; } // make sure all directions will be processed
#endif
/*       printf("%li %li\n",ST,EN); */
/*       exit(0); */
      if (!map.maskLoaded()) { map.makekill_space_manager("make","m",1); map.clean_mask(); }

      for (i=ST;i<EN;i++) {
	if (!map.isMasked(i)) {
	  cpeds_pix2ang_healpix(map.nside(),i,&n.lat(),&n.lon(),1);	  n.lat()=PIsnd-n.b(); // get this pixel coordinates and convert to galactic CS
	  n.lon()=n.l()*PI180inv; n.lat()=n.b()*PI180inv;
	  cq = map.get_circle(n.l(),n.b(),_find_circles_r,0.5,_overplot_region_type,_overplot_region_type_dot_points); // get the set of pixels around the current pixel in radius

	  if (cq->get_size() > 1) {
/* 	    cq->printq_long(); */
	    tabint = cq->export_array();
	    cppt = new cpeds_point[cq->get_size()];
	    m=cq->get_size();
	    for (j=0;j<m;j++) { cppt[j].x = map.get_T(cq->getq(j)); cppt[j].y=(double)tabint[j]; }
	    cpeds_sort_data(m,cppt,21);
	    for (j=0;j<m;j++) { tabint[j]=(long)cppt[j].y; }
	    circle_points->add_array(tabint,1); // add to the common queue
	    delete cppt;
	    delete tabint;
	  }
	  delete cq;

	}
#ifdef GOMPI
	if (world_id==0) printf("node: %i> points num: %li / %li points tot %li\n",world_id,i,Ndirpn,circle_points->get_size());
#else
	printf("points num: %li / %li points tot %li\n",i,k,circle_points->get_size());
#endif

      }



#ifdef GOMPI
      // collect data from other nodes
      if ( world_id != 0 ) { // if I'm slave then send the info to master
	m=circle_points->get_size();
	MPI_Send(&m, 1,MPI_LONG,0,1,MPI_COMM_WORLD);
	MPI_Send(circle_points->export_array(), circle_points->get_size(),MPI_LONG,0,2,MPI_COMM_WORLD);
      }
      else {
	if (node_num>1) {
	  // get the points from other nodes
	  for (i=1;i<node_num;i++) {
	    MPI_Recv(&m,1,MPI_LONG,i,1,MPI_COMM_WORLD,&mpi_stat);
	    tabint = new long[m];
	    MPI_Recv(tabint,m,MPI_LONG,i,2,MPI_COMM_WORLD,&mpi_stat);
	    circle_points->add_array(tabint,m);
	    delete tabint;
	  }
	}
      }
#endif

      // mark the points in the map;
      map.clear_map();
      j=circle_points->get_size();
      for (i=0;i<j;i++) {
	k=circle_points->getq(i);
	map.set_T(k,map.get_T(k)+1.0);
      }

      delete circle_points;

  }

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//                           VARIOUS OUTPUTS FROM THE PROGRAM
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//   if (_output_file == "xyz") {  // text file with 3 columns with x,y and z values  of the plot points
//     map.calculate_map_stats(1);
//     if (map.flat_coord_loaded == 0) map.calculate_flat_coord(projection,-180,0);
//     tab = map.create_resized_flat_map(0);
//     fname = _input_file+".xyz";
//     f = fopen(fname.c_str(),"w");
//     k=0;
// /*     for (j=0;j<map.ypix_num_flat;j++) { */
// /*       for (i=0;i<map.xpix_num_flat;i++) { fprintf(f,"%lE ",tab[i*map.xpix_num_flat+j]); k++; } */
// /*       fprintf(f,"\n"); */
// /*     } */

//     for (j=0;j<map.pix_num;j++) {
//       //map.get_C(j,&l,&b);
//       x=map.flatcoord[j].x;       y=map.flatcoord[j].y;
//       fprintf(f,"%lE %lE %lE\n",x,y,map.get_T(j));    }

//     fclose(f);
//     exit(0);
//   }

//-----------------------------------------------------------------------------------------------------
  if (_output_file == "pyth" || _output_file == "lonlat") {  // saves the map on regular longitude lattitude grid into 3 files
	  printf("|DRAW_MAPS> WILL OUTPUT AS A PYTHON PLOT\n");
	  map.calculate_map_stats(1); Tth=map.getMinT()*10; if (Tth >= 0.0) Tth=-1;
	  
	  fname = _input_file+".lonlatT";    f = fopen(fname.c_str(),"w");
	  fname = _input_file+".lon";    f1 = fopen(fname.c_str(),"w");
	  fname = _input_file+".lat";    f2 = fopen(fname.c_str(),"w");
	  
	  //_lreg_num = 2*(_lreg_num+1);    _breg_num = 2*(_breg_num+1);
	  dl=twoPI/(double)_resX;     l0=0; //20*PI/180;
	  db=PI/(double)_resY; 
	  b=PIsnd-db/2;
	  b=-PIsnd+db/2;
	  for (i=0;i<_resY;i++) { 
		  l=dl/2+l0;
		  for (j=0;j<_resX;j++) { 
			  if (_reverse_l) n.lon()=twoPI-l;/* +dl/2; */ else n.lon() = l; 
			  n.lat() = b;
			  T=map.get_T(n,&k);
			  if (_dont_plot_zero && T == 0.0) T=Tth;
			  if (masked) { if (map.get_m(k) == 0) T=Tth; }
			  fprintf(f,"%lE ",T); // save the temperature
			  if (i==0) { fprintf(f1,"%lE\n", 180/PI*l); } // save the equispaced longitude
			  l+=dl;
		  }
		  if (_cyclic_point) {
//			  l=dl/2+l0;
			  n.lon()= dl/2+l0; if (_reverse_l) n.lon()=twoPI-l;/* +dl/2; */ else n.lon() = l;  
			  n.lat() = b;
			  T=map.get_T(n,&k);
			  if (_dont_plot_zero && T == 0) T=Tth;
			  if (masked) { if (map.get_m(k) == 0) T=Tth; }
			  fprintf(f,"%lE\n",T);
			  if (i==0) { fprintf(f1,"%lE\n", PI180inv*l); }
		  }
		  else fprintf(f,"\n");
		  fprintf(f2,"%lE\n", PI180inv*b); // save the equispaced longitude
		  b+=db;
	  }
	  
    fclose(f);    fclose(f1);    fclose(f2);

    if (_output_file == "lonlat") { exit(0); }

    if (_meridians) meridians = 1; else meridians=0;
    if (_plotNS) plotNS = 1; else plotNS=0;

    //sprintf(commandstr,"draw_maps_pyth2.py %s %li",_input_file.c_str(),_color_num);
    //sprintf(commandstr,"draw_maps_pyth2.py %s %li %lE %lE",_input_file.c_str(),_color_num,map.minT,map.maxT);
    //sprintf(commandstr,"draw_maps_pyth3.py %s %li %lE %lE",_input_file.c_str(),_color_num,map.minT,map.maxT);

/*     sprintf(commandstr,"draw_maps_pyth4.py %s %li %lE %lE %li %s %s %s %s %i %i %lf %lf %s %s %li %li %s",_input_file.c_str(),_color_num,map.minT,map.maxT,_dpi,_bgcolor.c_str(),_fgcolor.c_str(),_mytitle.c_str(),_colorbar.c_str(),meridians,plotNS,_delta_latitude,_delta_longitude,_labels_format.c_str(),_colorbar_label.c_str(),_colorbar_fontsize,_meridians_fontsize, _plot_taken_vals.c_str()); // commented out on 2009-04-11 */
    if (_Zrange[0]==_Zrange[1]) { _zrange_min=map.getMinT(); _zrange_max=map.getMaxT(); } else { _zrange_min=_Zrange[0]; _zrange_max=_Zrange[1]; }
    sprintf(commandstr,"draw_maps_pyth8-simple.py %s",_input_file.c_str());
/*     sprintf(commandstr,"draw_maps_pyth4.py %s %li %lE %lE %li %s %s %s %s %i %i %lf %lf %s %s %li %li %s",_input_file.c_str(),_color_num,_zrange_min,_zrange_max,_dpi,_bgcolor.c_str(),_fgcolor.c_str(),_mytitle.c_str(),_colorbar.c_str(),meridians,plotNS,_delta_latitude,_delta_longitude,_labels_format.c_str(),_colorbar_label.c_str(),_colorbar_fontsize,_meridians_fontsize, _plot_taken_vals.c_str()); */
    printf("executing: %s\n", commandstr);
    system(commandstr);
    exit(0);
  }

//-----------------------------------------------------------------------------------------------------
//   if (_output_file == "mat") {
//     map.calculate_map_stats(1);
//     if (map.flat_coord_loaded == 0) map.calculate_flat_coord(projection,-180,0);
//     tab = map.create_resized_flat_map(0);
//     fname = _input_file+".mat";
//     f = fopen(fname.c_str(),"w");
//     fprintf(f,"%i %i\n",map.xpix_num_flat,map.ypix_num_flat);

//     k=0;
//     for (j=0;j<map.xpix_num_flat;j++) {
//       for (i=0;i<map.ypix_num_flat;i++) { fprintf(f,"%lE ",tab[i*map.xpix_num_flat+j]); k++; }
//       fprintf(f,"\n");
//     }
// /*     for (i=0;i<map.flat_coord_num;i++) { */
// /*       fprintf(f,"%lE ",tab[map.xpix_num_flat+y+x]); */
// /*       x = (int)round((map.flatcoord[i].x-minx_flat)/factor);  */
// /*       y = (int)round((map.flatcoord[i].y-miny_flat)/factor);  */
// /*       //if (mask_loaded == 0) { */
// /*       flattab[xpix_num_flat*y+x] = map[i].T*map[i].m; // here the mask only if loaded also is plotted */
// /*   } */

//     fclose(f);
//     exit(0);
//   }
//-----------------------------------------------------------------------------------------------------
  if (_output_file == "lb") {
	  mscsFunction lb;
	  map.set_map_coord(0,0);
	  if (_dont_plot_zero) {
		  for (long i = 0; i < map.pixNum(); i++) {
			  if (map.T(i)!=0) lb.newPoint(map.get_n()[i].toQPointF());
		  }
	  }
	  else {
		  for (long i = 0; i < map.pixNum(); i++) {
			  lb.newPoint(map.get_n()[i].toQPointF());
		  }
	  }
	  lb.multiply(PI180inv);
	  lb.scaleX(PI180inv);
	  lb.save("lb");
	  exit(0);
  }
  if (_output_file == "lbT") {
	  cpedsPointSet3D lbT;
	  map.set_map_coord(0,0);
	  if (_dont_plot_zero) {
		  for (long i = 0; i < map.pixNum(); i++) {
			  if (map.T(i)!=0) lbT.append(cpedsPoint3D(map.get_n()[i].lon()*PI180inv,map.get_n()[i].lat()*PI180inv,map.T(i)));
		  }
	  }
	  else {
		  for (long i = 0; i < map.pixNum(); i++) {
			  lbT.append(cpedsPoint3D(map.get_n()[i].lon()*PI180inv,map.get_n()[i].lat()*PI180inv,map.T(i)));		  
		  }
	  }
	  lbT.save("lbT");
	  exit(0);
  }

  if (_output_file == "x11") { output = 1; }
  // if (_output_file == "eps") { output = 2; sprintf(map.flatmap_file_name,"%s.eps",_input_file.c_str()); }
  // if (_output_file == "gif") { output = 3; sprintf(map.flatmap_file_name,"%s.gif",_input_file.c_str()); }
  // if (_output_file == "jpg") { output = 3; sprintf(map.flatmap_file_name,"%s.jpg",_input_file.c_str()); }


//-----------------------------------------------------------------------------------------------------
  if (_output_file.find("bin",_output_file.size()-3) != string::npos) { map.savebinT(_output_file); return 0;  }
  if (_output_file.find("txt",_output_file.size()-3) != string::npos) { 
	  if (_dont_save_mask) { 
		  cpedsList<double> cl=map.get_nonMasked("T"); cl.save(_output_file);   
	  } 
	  else { map.savetxtT(_output_file);  } return 0; }
  // if (_output_file.find("-TCn-txt",_output_file.size()-8) != string::npos) { _output_file.erase(_output_file.size()-8,8); sprintf(file2,"%s",_output_file.c_str()); strcat(file2,"-TCn-txt");  map.set_map_coord(0,0); map.savetxtT(file2,11); return 0; }
  // if (_output_file.find("-Tr-txt",_output_file.size()-7) != string::npos) { _output_file.erase(_output_file.size()-7,7); sprintf(file2,"%s",_output_file.c_str()); map.conv_nest2ring(map.map); map.savetxtT(file2,0); exit(0); }
  // if (_output_file.find("-Tr-bin",_output_file.size()-7) != string::npos) { _output_file.erase(_output_file.size()-7,7); sprintf(file2,"%s",_output_file.c_str()); map.conv_nest2ring(map.map); map.savebinT(file2,0); exit(0); }
//-----------------------------------------------------------------------------------------------------
  msgs.say("PLOTTING THE MAP",Top);
  ret_code = map.plot(_sizeX,_sizeY,_proj,_l0,_b0,0,0,_color_num,_color_scheme,_input_file.c_str(),_BGColor, _PGwindow_size, output);

  histname_file = _input_file+"-temperature-histogram.txt";
  // if (ret_code == 0) map.plot_temperature_histogram(histname_file);

  // if (_mink) {
  //   map.calculate_minkowski_area(map.sigma0_Q1,_mlevels,map.minT,map.maxT);
  //   map.plot_minkowski_area(1);
  //   //sim1.savetxtM(&out,1);
  //   map.calculate_minkowski_circ(map.sigma0_Q1,_mlevels,map.minT,map.maxT);
  //   map.plot_minkowski_circ(1);
  // }

  // map.kill_map(); map2.kill_map();


  return 0;
}


void parseOptions(int argc, char** argv) {
  long i;
  string::size_type j;
  ValuesConstraint<string>* allowedStrNew;
  
	try {

	  CmdLine cmd("draw_maps: heaplix maps plotter\n", ' ', Mscs_version.c_str() );

	//
	// Define arguments
	//


	UnlabeledValueArg<string> input_file("input_file","file to plot",false,"nofile","string");	cmd.add( input_file );
	ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
	SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
	SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );

	vector<string> allowedStr;
	allowedStr.push_back("moll");
	allowedStr.push_back("aitof");
	allowedStr.push_back("robin");
	allowedStr.push_back("merc");
	allowedStr.push_back("cc");
	allowedStr.push_back("mill");
	allowedStr.push_back("cea");
	allowedStr.push_back("ortho");
	allowedStr.push_back("stere");
	
	allowedStrNew = new ValuesConstraint<string>(allowedStr);	
	ValueArg<string> proj("p","proj","projection to use [moll - Mollweide/ aitoff etc] (default: moll)",false,"moll", allowedStrNew); cmd.add(proj);
	delete allowedStrNew;
	
	ValueArg<string> ft("f","ft","input file type [bin, txt,fitsWMAP, fitsPL] (default: bin)",false,"bin","string"); cmd.add(ft);
	// ValueArg<double> loffset("", "dl", "shift the map in galactic longitude in degr.", false,0,"double");	cmd.add( loffset );
	ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
	SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
	ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );

	ValueArg<string> output("o","out","outfile matrix or eps file; allowed values [x11, mat - matrix, eps - eps file, xyz - text file with coord x,y and temp.Tn-bin .Tr-bin, .Tr-txt, .Tn-txt, -TCn-txt; jpg,gif] default: x11 NOT FULLY IMPLEMENTED",false,"x11","string"); cmd.add(output);
	ValueArg<long> sizeX("","sx","X size of the plot ",false,0,"long"); cmd.add(sizeX);
	ValueArg<long> sizeY("","sy","Y size of the plot ",false,0,"long"); cmd.add(sizeY);
	ValueArg<double> l0("","l0","central meridian (default 0)",false,0,"double"); cmd.add(l0);
	ValueArg<double> b0("","b0","central parallel (default 0)",false,0,"double"); cmd.add(b0);
	ValueArg<long> nside("n", "nside", "nside in which the map is to be plotted (default 512)", false,0,"long");	cmd.add( nside );
	ValueArg<long> renside("", "renside", "change resolution to nside and back before plotting", false,0,"long");	cmd.add( renside );

	allowedStr.clear();
	allowedStr.push_back("spectralBGR");
	allowedStr.push_back("spectralRGB");
	allowedStr.push_back("spectral");
	allowedStr.push_back("color");
	allowedStr.push_back("grayscale");
	allowedStrNew = new ValuesConstraint<string>(allowedStr);	
	ValueArg<string> color_scheme("", "colors", "set map colors scheme [grayscale,color,color-spots,color-spotl,color-spoth,rainbow]", false,"spectralBGR",allowedStrNew);	cmd.add( color_scheme );
	delete allowedStrNew;

	ValueArg<long> color_num("", "color_num", "set number of colors for map (default 30)", false,30,"long");	cmd.add( color_num );
	SwitchArg logscale("","logscale", "draw the log_10 of the map (T<0 are dropped)", false);	cmd.add( logscale );
	SwitchArg abs("","abs", "plot the absolute values of the map ", false);	cmd.add( abs );
	SwitchArg norm_by_max("","norm_by_max", "normalize the map by the maximal value in the map ", false);	cmd.add( norm_by_max );

	ValueArg<string> multi_mask("","mm", "wheather to mask the map for calculations with multi_mask[LB,HP]; number of regions in l and b given by lreg and breg options for thphi mask or mmnside parameter for healpix mask (default: OFF - no multimask)", false,"","string");	cmd.add( multi_mask );
/* 	SwitchArg multi_mask("","mm", "wheather to mask the map for calculations with multi_mask; number of regions in l and b given by lreg and breg options (default: OFF)", false);	cmd.add( multi_mask ); */
	ValueArg<long> lreg_num("","lreg","number of divisions of a hepisphere in l direction (default: 4)",false,4,"long"); cmd.add(lreg_num);
	ValueArg<long> breg_num("","breg","number of divisions of a hepisphere in b direction (default: 4)",false,4,"long"); cmd.add(breg_num);
	ValueArg<long> mmnside("","mmnside","nside for healix multi mask (default: 4)",false,4,"long"); cmd.add(mmnside);
	MultiArg<long>mmfield("r","reg","number of region to mark when plotting multimask",false,"long"); cmd.add(mmfield);
	SwitchArg negate_r("","negate_r", "negates the meaning of the -r list, all regions not on list will be removed from analysis  (default false)", false);	cmd.add( negate_r );

	ValueArg<double> mm_Rl("","mm_Rl","rotation of multimask in l coordinate (first)",false,0,"double"); cmd.add(mm_Rl);
	ValueArg<double> mm_Rth ("","mm_Rth","rotation of multimask to (th,phi) orientation (given in gal.CS) (second)",false,0,"double"); cmd.add(mm_Rth);
	ValueArg<double> mm_Rphi("","mm_Rphi","2rotation of multimask to (th,phi) orientation (given in gal.CS) (second)",false,0,"double"); cmd.add(mm_Rphi);
	SwitchArg randomize("","randomize", "whether or not to randomize the values of regions in the multimask (for visualization purposes)", false);	cmd.add( randomize );
	SwitchArg random_orientation("","random_orientation", "whether or not to make the random orientation of the map", false);	cmd.add( random_orientation );
	SwitchArg mm_random_orientation("","mm_random_orientation", "whether or not to make the random orientation of the multimask", false);	cmd.add( mm_random_orientation );

	SwitchArg cyclic("","cyclic", "whether or not include the cyclic point in the 'lonlat' output", false);	cmd.add( cyclic );
	SwitchArg reverse_l("","reverse_l", "whether or not reverse the direction of in the 'lonlat' output", false);	cmd.add( reverse_l );
	ValueArg<long> resX("","resX","resolution in phi (default: 600)",false,600,"long"); cmd.add(resX);
	ValueArg<long> resY("","resY","resolution in theta (default: 600)",false,600,"long"); cmd.add(resY);
	ValueArg<long> dpi("","dpi","resolution in dpi for lonlat output (default: 300)",false,300,"long"); cmd.add(dpi);
	ValueArg<double> delta_latitude("","dlat","separation between latitude meridians in lonlat plot in [deg] (default: 15)",false,15,"double"); cmd.add(delta_latitude);
	ValueArg<double> delta_longitude("","dlon","separation between longitude meridians in lonlat plot in [deg] (default: 45)",false,45,"double"); cmd.add(delta_longitude);
	SwitchArg plot_taken_vals("","plot_taken_vals", "whether or not to plot info stored when taking values from map with T (default: false)", false);	cmd.add( plot_taken_vals );


	SwitchArg nopyth("","nopyth", "whether or not plot the data with python script", false);	cmd.add( nopyth );

	ValueArg<double> zrange("","zrange","number of standard deviations to plots (above and under, performed before Zmin/max)",false,0,"double"); cmd.add(zrange);
	ValueArg<double> Zmin("", "min", "minimal z value in the map to plot; rest will be scramblesd; done before logscrambling;", false,0,"double");	cmd.add( Zmin );
	ValueArg<double> Zmax("", "max", "maximal z value in the map to plot; rest will be scramblesd; done before logscrambling", false,0,"double");	cmd.add( Zmax );
	ValueArg<double> Zminl("", "minl", "minimal z value in logscale in the map to plot; rest will be scramblesd; done only in logscale mode", false,0,"double");	cmd.add( Zminl );
	ValueArg<double> Zmaxl("", "maxl", "maximal z value in logscale in the map to plot; rest will be scramblesd; done only in logscale mode", false,0,"double");	cmd.add( Zmaxl );

	ValueArg<string> mytitle("", "title", "set plot title (default: none)", false,"notitle","string");	cmd.add( mytitle );
	ValueArg<string> bgcolor("", "bgcolor", "set plot background color (default: none)", false,"white","string");	cmd.add( bgcolor );
	ValueArg<string> fgcolor("", "fgcolor", "set plot foroground color (default: none)", false,"white","string");	cmd.add( fgcolor );
	ValueArg<string> colorbar("", "colorbar", "whether or not put the colorbar in the plot with output lonlat [moll,N,S,all, NS] (default: nobar)", false,"nobar","string");	cmd.add( colorbar );
/* 	SwitchArg colorbar("","colorbar", "whether or not put the colorbar in the plot with output lonlat", false);	cmd.add( colorbar ); */
	SwitchArg meridians("","meridians", "whether or not put the meridians in the plot with output lonlat", false);	cmd.add( meridians );
	SwitchArg plotNS("","plotNS", "whether or not to plot north and south stereographic projection too with output lonlat", false);	cmd.add( plotNS );
	SwitchArg dont_plot_zero("","dont_plot_zero", "whether or not to plot values in the map that have zero value", false);	cmd.add( dont_plot_zero );
	ValueArg<string> labels_format("", "labels_format", "format of the labels in the colorbar [E/f] (default: E)", false,"E","string");	cmd.add( labels_format );
	ValueArg<string> colorbar_label("", "colorbar_label", "xlabel for the color bar [] (default: [K])", false,"[K]","string");	cmd.add( colorbar_label );
	ValueArg<long> colorbar_fontsize("", "colorbar_fontsize", "font size of the color bar labels  (default: [14])", false,14,"long");	cmd.add( colorbar_fontsize );
	ValueArg<long> meridians_fontsize("", "meridians_fontsize", "font size of the meridians labels bar  (default: [14])", false,14,"long");	cmd.add( meridians_fontsize );
	ValueArg<double> PGwindow_size("", "PGws", "size of the PGplot window in inches (eg. 5 - small, 14- big) (default 8)", false,8,"double");	cmd.add( PGwindow_size );


	// do the special block
	SwitchArg do_special_block("","do_special_block", "whether or not go through the special block in the program (special use only)", false);	cmd.add( do_special_block );
	SwitchArg find_circles("","find_circles", "for each point in the map selects a ring of radius 140 deg and stores N hottest points; The resulting map is plotted or saved", false);	cmd.add( find_circles );
	ValueArg<long> find_circles_r("", "find_circles_r", "radius within which to look for pixles -- option to the find_circles command  (default: [140] deg)", false,140,"double");	cmd.add( find_circles_r );

	// inducing symmetries
	SwitchArg mk_gal_symm("","mk_gal_symm", "whether or not force the symmetry T(b,l))=T(-b,l) ", false);	cmd.add( mk_gal_symm);
	SwitchArg mk_gal_asymm("","mk_gal_asymm", "whether or not force the symmetry T(b,l))=-T(-b,l) ", false);	cmd.add( mk_gal_asymm);
	SwitchArg mk_ew_symm("","mk_ew_symm", "whether or not force the symmetry T(b,l))=T(b,-l) ", false);	cmd.add( mk_ew_symm);
	SwitchArg mk_ew_asymm("","mk_ew_asymm", "whether or not force the symmetry T(b,l))=-T(b,-l) ", false);	cmd.add( mk_ew_asymm);
	SwitchArg mk_pt_symm("","mk_pt_symm", "whether or not force point symmetry T(n))=T(-n) ", false);	cmd.add( mk_pt_symm);
	SwitchArg mk_pt_asymm("","mk_pt_asymm", "whether or not force point symmetry T(n))=-T(-n) ", false);	cmd.add( mk_pt_asymm);
	SwitchArg mk_ring_av("","mk_ring_av", "whether average the map in rings before plotting ", false);	cmd.add( mk_ring_av);
	ValueArg<long> ravns("", "ravns", "ns in which perform the ring average  (default: 512)", false,512,"long");	cmd.add( ravns );

	// making thresholded map
	SwitchArg mk_thres_map("","mk_thres_map", "whether or not to threshold the map with given thresholds in th parameter (default: false) ", false);	cmd.add( mk_thres_map);
	SwitchArg area("","area", "whether or not to threat the thresholds given in th parameter as a percentage of maps area (useful for plotting areas of eg 2,3 sigma)  (default: false) ", false);	cmd.add( area);
	MultiArg<double> th("","th","thresholds levels for thresholded map  ",false,"double"); cmd.add(th);

	// rotations of the map
	ValueArg<double> Rl("","Rl","rotation of the map in l coordinate (first)",false,0,"double"); cmd.add(Rl);
	ValueArg<double> Rth ("","Rth","rotation of the map in phi to (th,phi) orientation (given in gal.CS) (second)",false,0,"double"); cmd.add(Rth);
	ValueArg<double> Rphi("","Rphi","rotation of the map in theta to (th,phi) orientation (given in gal.CS) (second)",false,0,"double"); cmd.add(Rphi);
	SwitchArg save_rot("","save_rot", "whether or not to save the rotation information for future quick use with load_rot", false);	cmd.add( save_rot );
	SwitchArg load_rot("","load_rot", "whether or not to load the rotation information to perform the rotation faster", false);	cmd.add( load_rot );
	SwitchArg rot_short("","rot_short", "whether or not to rot_short (default: true)", true);	cmd.add( rot_short );
	SwitchArg random_orientation_sph("","random_orientation_sph", "chooses a random direction from sphere and rotates the map there (default: false)", false);	cmd.add( random_orientation_sph );
	SwitchArg rotate_random_gauss_lb("","rotate_random_gauss_lb", "chooses a random direction in l and b from gaussian distribution centered at (rotGmphi, rotGmth) with standard deviation (rotGsphi,rotGsth) and rotates the map there. IMPORTANT !!! if you use this thing for many runs then make sure there's a time lag between two subsequent runs larger than a 1s to avoid having the same random seeds; or else reprogram the thing for fast runs (default: false)", false);	cmd.add( rotate_random_gauss_lb );
	ValueArg<double> rotGmphi("","rotGmphi","mean l value for random rotation",false,0,"double"); cmd.add(rotGmphi);
	ValueArg<double> rotGmth("","rotGmth","mean b value for random rotation",false,0,"double"); cmd.add(rotGmth);
	ValueArg<double> rotGsphi("","rotGsphi","st.dev. of l value for random rotation",false,0,"double"); cmd.add(rotGsphi);
	ValueArg<double> rotGsth("","rotGsth","st.dev. of b value for random rotation",false,0,"double"); cmd.add(rotGsth);


	// making marks on map
	MultiArg<double> ml("","ml","mask region on the mask at galactic longitude  (upto 100 values accepted) ",false,"double"); cmd.add(ml);
	MultiArg<double> mb("","mb","mask region on the mask at galactic lattitude (upto 100 values accepted)",false,"double"); cmd.add(mb);
	MultiArg<double> ms("","ms","mask region of size (to be connected with ml and mb parameters) (upto 100 values accepted)",false,"double"); cmd.add(ms);
	SwitchArg invert_mask("","invert_mask", "whether or not to invert the mask (default: faslse)", false);	cmd.add( invert_mask );
	MultiArg<string> cmf("","cmf","circular mask text file with galactic coordinates of the position of the circular mask and its sizes (upto 5 files accepted) ",false,"string"); cmd.add(cmf);
	ValueArg<double> cmf_r("","cmf_r","circular mask size to apply to all sources in the file, regardless of the values in the file [deg]",false,0,"double"); cmd.add(cmf_r);
	ValueArg<double> cmf_r0("","cmf_r0","circular mask size 0 offset value to apply to all r data in the file [deg]  (default: 0)",false,0,"double"); cmd.add(cmf_r0);
	ValueArg<string> cmf_U("","cmf_U","Units of the circular mask size in text (default: [deg], can be: d - deg, m - arcmin, s - arcsec",false,"d","string"); cmd.add(cmf_U);

	MultiArg<double> tl("","tl","mark region on the map at galactic longitude  (upto 100 values accepted) ",false,"double"); cmd.add(tl);
	MultiArg<double> tb("","tb","mark region on the map at galactic lattitude (upto 100 values accepted)",false,"double"); cmd.add(tb);
	MultiArg<double> ts("","ts","mark region of size (to be connected with tl and tb parameters) (upto 100 values accepted)",false,"double"); cmd.add(ts);
	MultiArg<double> tv("","tv","mark region with value (to be connected with tl and tb parameters) (upto 100 values accepted)",false,"double"); cmd.add(tv);

	MultiArg<double> trl("","trl","mark great ring with pole at this longitude [deg] ",false,"double"); cmd.add(trl);
	MultiArg<double> trb("","trb","mark great ring with pole at this latitude [deg] ",false,"double"); cmd.add(trb);
	MultiArg<double> trv("","trv","mark great ring with value (default: 0) ",false,"double"); cmd.add(trv);
	MultiArg<double> trs("","trs","mark great ring with ring thickness in [deg] (default: 1deg) ",false,"double"); cmd.add(trs);
	SwitchArg trecl("","trecl", "whether or not to plot the ecliptic with this method (default: false)", false);	cmd.add( trecl );
	SwitchArg treq("","treq", "whether or not to plot the equator with this method (default: false)", false);	cmd.add( treq );
	ValueArg<double> eclVal("","eclVal","value with which to mark ecliptic (default 0)",false,0,"double"); cmd.add(eclVal);
	ValueArg<double> eqVal("","eqVal","value with which to mark equator (default 0)",false,0,"double"); cmd.add(eqVal);

	MultiArg<string> ctf("","ctf","circular map text file with galactic coordinates of the position of the circular map and its sizes and values",false,"string"); cmd.add(ctf);
	ValueArg<double> ctf_r("","ctf_r","circular region size to apply to all sources in the file, regardless of the values in the file [deg]",false,0,"double"); cmd.add(ctf_r);
	ValueArg<double> ctf_v("","ctf_v","circular region value to apply to all sources in the file, regardless of the values in the file [deg]",false,0,"double"); cmd.add(ctf_v);
	ValueArg<double> ctf_r0("","ctf_r0","circular region size 0 offset value to apply to all r data in the file [deg] (default: 0)",false,0,"double"); cmd.add(ctf_r0);
	ValueArg<string> ctf_U("","ctf_U","Units modified for the circular reigon sizes (default: [deg], can be: d - deg, m - arcmin, s - arcsec",false,"d","string"); cmd.add(ctf_U);
	ValueArg<long> ctf_cols("","ctf_cols","number of columns in the ctf file (default: 4) 2 - l,b; 3 - l,b,v; 4 - l,b,r,v",false,4,"long"); cmd.add(ctf_cols);
	SwitchArg ctf_add_vals("","ctf_add_vals", "adds up to the underlying values instead of overplotting (default: false)", false);	cmd.add( ctf_add_vals );


	ValueArg<string> overplot_region_type("","plot_reg_type","region type for plotting on map (dot: default, emptydot) (to be connected with ml and mb parameters) ",false,"dot","string"); cmd.add(overplot_region_type);
	ValueArg<long> overplot_region_type_dot_points("","plot_reg_type_circle_points","numer of points to plot in the width of circle (default: 1) (to be connected with ml and mb parameters) ",false,1,"long"); cmd.add(overplot_region_type_dot_points);

/* 	ValueArg<string> save_overplot_as("","save_over_plot_as","saves the overplots into a txt file with pixel numbers instead of l,b ",false,"dot","string"); cmd.add(save_overplot_as); */

	SwitchArg show_ecliptic("","show_ecliptic", "overplots ecliptic and poles from data/512-poles_ecliptic file", false);	cmd.add( show_ecliptic );
	SwitchArg show_equator("","show_equator", "overplots equator and poles from data/512-poles_equator file", false);	cmd.add( show_equator );

	ValueArg<long> Yl("","Yl","l'th harmonic to plot (default: 4)",false,4,"long"); cmd.add(Yl);
	ValueArg<long> Ym("","Ym","m'th harmonoc to plot (default: 0)",false,0,"long"); cmd.add(Ym);
	ValueArg<double> reim_ratio("","reim_ratio","ratio of re to im part in requested Ylm to plot (default: 0)",false,0,"double"); cmd.add(reim_ratio);

	ValueArg<long> seed_offset("", "seed_offset", "offset for the random seed generator used for many parallel runs.  (default: 0)", false,0,"long");	cmd.add( seed_offset );


	//SwitchArg pyth("","pyth", "whether or not plot the data with python script", false);	cmd.add( pyth );
	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);

	//
	// Set variables
	//
	_input_file = input_file.getValue(); if ( _input_file == "nofile" || _input_file == "Ylm" ) { _not_from_file = true; } else { _not_from_file = false; }
	_mask_file = mask.getValue(); if (_mask_file.length() == 0) { masked = false; } else {masked = true; }
	_mask_from_here = mask_from_here.getValue();
	_invert_mask = invert_mask.getValue();
	_dont_save_mask = dont_save_mask.getValue();
	_output_file = output.getValue();
	_proj = proj.getValue();
	_ft = ft.getValue();
	_ord = ord.getValue();
	// _loffset = loffset.getValue();
	_mink = mink.getValue();
	_mlevels = mlevels.getValue();
	_nside = nside.getValue();
	_renside = renside.getValue();
	_color_scheme = color_scheme.getValue();
	_color_num = color_num.getValue();
	_logscale = logscale.getValue();
	_abs = abs.getValue();
	_norm_by_max = norm_by_max.getValue();

	_multi_mask = multi_mask.getValue(); if (_multi_mask == "LB" || _multi_mask == "HP") multi_masked = true; else multi_masked = false;
	_lreg_num = lreg_num.getValue();
	_breg_num = breg_num.getValue();
	_mmnside = mmnside.getValue();
	_mm_Rl = mm_Rl.getValue();
	_mm_Rth = mm_Rth.getValue();
	_mm_Rphi = mm_Rphi.getValue();
	_save_rot = save_rot.getValue();
	_load_rot = load_rot.getValue();
	_rot_short = rot_short.getValue();
	_randomize = randomize.getValue();
	_random_orientation = random_orientation.getValue();
	_random_orientation_sph = random_orientation_sph.getValue();
	_mm_random_orientation = mm_random_orientation.getValue();
	_rotate_random_gauss_lb = rotate_random_gauss_lb.getValue();
	_rotGmphi = rotGmphi.getValue();
	_rotGmth = rotGmth.getValue();
	_rotGsphi = rotGsphi.getValue();
	_rotGsth = rotGsth.getValue();

	_Rl = Rl.getValue(); _isSetRl=Rl.isSet();
	_Rth = Rth.getValue(); _isSetRth=Rth.isSet();
	_Rphi = Rphi.getValue(); _isSetRphi=Rphi.isSet();

	_zrange = zrange.getValue();
	_Zrange[0] = Zmin.getValue();
	_Zrange[1] = Zmax.getValue();
	_Zrangelog[0] = Zminl.getValue();
	_Zrangelog[1] = Zmaxl.getValue();
	_mytitle = mytitle.getValue();
	_bgcolor = bgcolor.getValue();
	_fgcolor = fgcolor.getValue();
	_dpi = dpi.getValue();
	_meridians = meridians.getValue();
	_meridians_fontsize = meridians_fontsize.getValue();
	_PGwindow_size = PGwindow_size.getValue();
	_plotNS = plotNS.getValue();
	_dont_plot_zero = dont_plot_zero.getValue();
	_delta_latitude = delta_latitude.getValue();
	_delta_longitude = delta_longitude.getValue();
	_plot_taken_valsB = plot_taken_vals.getValue(); if (_plot_taken_valsB) _plot_taken_vals = "true"; else _plot_taken_vals = "false";
	_colorbar = colorbar.getValue();
	_colorbar_label = colorbar_label.getValue();
	_colorbar_fontsize = colorbar_fontsize.getValue();
	_labels_format = labels_format.getValue();

	_cyclic_point = cyclic.getValue();
	_reverse_l = reverse_l.getValue(); //if (_output_file == "pyth") { _reverse_l=true; _cyclic_point=true; }
	_nopyth = nopyth.getValue();
	_resX = resX.getValue();
	_resY = resY.getValue();


	_Yl = Yl.getValue();
	_Ym = Ym.getValue();
	_reim_ratio = reim_ratio.getValue();

	//_pyth = pyth.getValue();

	_l0 = l0.getValue();
	_b0 = b0.getValue();
	_sizeX = sizeX.getValue(); printf("  - size is: %li\n",_sizeX);
	_sizeY = sizeY.getValue(); printf("  - size is: %li\n",_sizeY);
	if (_Zrange[0] == _Zrange[1] &&  _Zrange[1]==0) { _default_Z = true; } else { _default_Z = false;  } //printf("defaultZ false min %lE max %lE\n",_Zrange[0],_Zrange[1]);
	if (_Zrangelog[0] == _Zrangelog[1] && _Zrangelog[1]==0) { _default_Zlog = true; } else { _default_Zlog = false; }

	_negate_r = negate_r.getValue();
	ff = mmfield.getValue();
	// ffsize=v2.size();
	// if (v2.size() > 0) {  ff = new long[ffsize];	for (i=0;i<v2.size();i++) { ff[i] = v2[i]; printf("%li\n",ff[i]); } }
/* 	vector<double> v1 = camera.getValue(); */
/* 	for ( i = 0; (unsigned int)i < v1.size(); i++ ) { _camera[i] = v1[i]; printf("camera pozitio is: \n   - x=%lE\n  - y=%lE\n  -z=%lE\n",v1[0],v1[1],v1[2]); default_cam = false;} */


	vector<double>	v3;
	v3 = ml.getValue(); if (v3.size() > 0) { _ml[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _ml[i] = v3[i-1]; } }
	v3 = mb.getValue(); if (v3.size() > 0) { _mb[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _mb[i] = v3[i-1]; } }
	v3 = ms.getValue(); if (v3.size() > 0) { _ms[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _ms[i] = v3[i-1]; } }

	v3 = tl.getValue(); if (v3.size() > 0) { _tl[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _tl[i] = v3[i-1]; } }
	v3 = tb.getValue(); if (v3.size() > 0) { _tb[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _tb[i] = v3[i-1]; } }
	v3 = ts.getValue(); if (v3.size() > 0) { _ts[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _ts[i] = v3[i-1]; } }
	v3 = tv.getValue(); if (v3.size() > 0) { _tv[0] = (double)v3.size(); for (i=1;i<=(long)v3.size();i++) { _tv[i] = v3[i-1]; } }


	v3 = trl.getValue(); if (v3.size() > 0) { for (i=0;i<(long)v3.size();i++) { _trl.addq(v3[i]); } }
	v3 = trb.getValue(); if (v3.size() > 0) { for (i=0;i<(long)v3.size();i++) { _trb.addq(v3[i]); } }
	v3 = trs.getValue(); if (v3.size() > 0) { for (i=0;i<(long)v3.size();i++) { _trs.addq(v3[i]); } }
	v3 = trv.getValue(); if (v3.size() > 0) { for (i=0;i<(long)v3.size();i++) { _trv.addq(v3[i]); } }
	while (_trs.get_size()< _trl.get_size()) { _trs.addq(1.0); }
	while (_trv.get_size()< _trl.get_size()) { _trv.addq(0.0); }
	_trecl = trecl.getValue(); if (_trecl) { _trl.addq(CPEDS_ECLIPTIC_NPOLE_GAL_L); _trb.addq(CPEDS_ECLIPTIC_NPOLE_GAL_B); }
	_treq = treq.getValue();   if (_treq)  { _trl.addq(CPEDS_EQUATOR_NPOLE_GAL_L);  _trb.addq(CPEDS_EQUATOR_NPOLE_GAL_B); }
	while (_trs.get_size()< _trl.get_size()) { _trs.addq(1.0); }
	while (_trv.get_size()< _trl.get_size()) { _trv.addq(0.0); }

	_circm_files = cmf.getValue(); //if (v4.size() > 0) { _cmf_size = (long)v4.size(); for (i=0;i<(long)v4.size();i++) { _cmf[i] = v4[i]; } }
	_cmf_U = cmf_U.getValue();
	_cmf_r = cmf_r.getValue();
	_cmf_r0 = cmf_r0.getValue();

	_circT_files = ctf.getValue(); // v4 = ctf.getValue(); if (v4.size() > 0) { _ctf_size = (long)v4.size(); for (i=0;i<(long)v4.size();i++) { _ctf[i] = v4[i]; } }
	_ctf_U = ctf_U.getValue();
	_ctf_r = ctf_r.getValue();
	_ctf_v = ctf_v.getValue();
	_ctf_r0 = ctf_r0.getValue();
	// _ctf_cols = ctf_cols.getValue();
	_ctf_add_vals = ctf_add_vals.getValue();

	_overplot_region_type = overplot_region_type.getValue();
	_overplot_region_type_dot_points = overplot_region_type_dot_points.getValue();
	_show_ecliptic = show_ecliptic.getValue();
	_show_equator = show_equator.getValue();

	//j=_input_file.find("-Tn-bin",0); if (j != string::npos) { _input_file.erase(j,(string::size_type)7); }


	_mk_gal_symm = mk_gal_symm.getValue();
	_mk_gal_asymm = mk_gal_asymm.getValue();
	_mk_ew_symm = mk_ew_symm.getValue();
	_mk_ew_asymm = mk_ew_asymm.getValue();
	_mk_pt_symm = mk_pt_symm.getValue();
	_mk_pt_asymm = mk_pt_asymm.getValue();
	_mk_ring_av = mk_ring_av.getValue();
	_ravns = ravns.getValue();

	_mk_thres_map = mk_thres_map.getValue();
	_area = area.getValue();
	v3 = th.getValue(); if (v3.size() > 0) { for (i=1;i<=(long)v3.size();i++) { _th.addq(v3[i-1]); } }

	_seed_offset = seed_offset.getValue(); cpeds_seed_offset=_seed_offset;
	_do_special_block = do_special_block.getValue();
	_find_circles = find_circles.getValue();
	_find_circles_r = find_circles_r.getValue();


/* 	_save_overplot_as = save_overplot_as.getValue(); */
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}




