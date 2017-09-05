
#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"
#ifdef GOMPI
#include "mpi.h"
#endif


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


//typedef struct { 
//    double R; // real part
//    double I; // unreal part
//} a_lm;


typedef struct {
	long alms_num,lmax;
	mscsAlm::mscs_a_lm *a;
	mscsAlm::mscs_a_lm *aa;
	double *r; // |almn|^2/(2l+1)/Cl
	double *rr; // // |almnn|^2/(2l+1)/Cl
	cpeds_direction *n; // keps the direction that maximizes the r
	cpeds_direction *nn; // keps the direction that minimizes the rr
	//  long *i;
} Almn_type;

typedef struct {
	long lmax;
	long *m;
	long *mm;
	double *r;
	double *rr;
	cpeds_direction *n;
	cpeds_direction *nn;
	//  long *i;
} Dln_type;

typedef struct {
	long lmax;
	double *v;
	double *V;
} statl_type;


//declaration of the global variables
void parseOptions(int argc, char** argv);
void copy_alms(mscsAlms &a0j,Almn_type *Aj); 
void calculate_rs(Almn_type *Aj,mscsAlms &a0j); 
void update_Almn(Almn_type *Aj,mscsAlms &map, cpedsDirection n);
void update_Almn(Almn_type *Aj,mscsAlms &map, cpedsDirection n, FILE** frmap);

void init_A(Almn_type **A,long Nsim);
void initiate_Aj(Almn_type *Aj, long alms_num);
void init_D(Dln_type **D,long Nsim);
void initiate_Dj(Dln_type *Dj, long lmax);
void init_stat(statl_type **stat,long Nsim);
void initiate_statj(statl_type *statj, long lmax);

void kill_A(Almn_type **A);
void kill_Aj(Almn_type *A);
void kill_D(Dln_type **D);
void save_Almn(Almn_type *A, filenamestr fname);
void read_Almn(Almn_type **A, filenamestr fname);
void save_Dln(Dln_type *D, filenamestr fname);
void read_Dln(Dln_type **D, filenamestr fname);
void save_stat(statl_type *stat, filenamestr fname, long lmin, long lmax);
void save_chisq_vals(double * chisq,long st, long en,filenamestr fname);
void maximize_Almn(Almn_type *Aj,Dln_type *Dj);
double get_angle(cpedsDirection n1, cpedsDirection n2);
long int alm2num(long l, long m, long lmax);
void num2alm(long int num, long *l, long *m, long lmax);
void num2alm_norm(long int num, long *l, long *m);
long get_almsnumpos(long alms_num);
long get_almsnumtot(long alms_num);
long get_almsnumneg(long alms_num);
long numpos2lmax(long numpos);
bool _something_else();
void calculate_covariance_matrix(Almn_type *A,matrix<double> *Covl,long lmin, long lmax, long NsimC, long ist, string how);
double dot_product(Almn_type *Aj, long l, long m, long lp, long mp);
void copy_nl_to_vec(Almn_type *Aj,matrix<double> *nl, long LMIN, long lmax);
long get_nlmlpmp_vec_size(long lmin,long lmax);
void copy_nlmlpmp_to_vec(Almn_type *Aj, matrix<double> *nlmlpmp, long lmin, long lmax);

long _sim_st,_sim_en,_Nsim,_search_ns,_ns,_l,_LMIN,_LMAX,_Yl,_Ym;
string _data_dir,_output_file, _sim_pref,_sim_suff;
bool _find,_calc_Dstats,_calc_Astats,_calc_Covstats_inter_mode,_calc_Covstats_test,_mapYlm;

// global flag variables
bool _find_done,_calc_Dstats_done,_calc_Astats_done,_calc_Covstats_inter_mode_done,_calc_Covstats_test_done;

//-----------------------------------



int main(int argc, char **argv) {
	Almn_type *A=NULL;  // this keeps alms in orientations n, i that maximize r for each simulation
	Dln_type *D=NULL; // this keeps the number of m mode maximizing r in each multipole, each simulation, and corresponding direcion n,i
	statl_type *stat=NULL; // this keeps the estimator values for all maps and all multipoles
	
	long alms_num,lmax,LMIN;
	cpedsDirection n; 
	long i,j,k,l,m,ln,dir_st,dir_en;
	long Ndir, Nsim,NsimMCPDF,NsimC;
	double phi,th,ang, *chisq,tmpd;
	string tmps;
	filenamestr tmpch;
	
	
#ifdef GOMPI
	MPI_Status mpi_stat;
	int ierr,node_num,world_id=0;
	long Ndirpn;
	double from_slave, *from_slave1,*from_slave3,*from_slave5,*from_slave6,*from_slave7,*from_slave9;
#endif  
	FILE **frmap;//,*f2,*f3,*f4,*f5,*f6;
	
	matrix <double> *nl,*Covl,*ml,*check,*iCovl, *nlmlpmp, *mlmlpmp;
	long dimnl;
	
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	Mscs_initiate_global_variables();
	mscsAlms *aini;
	mscsMap *mini;
	mscsAngularPowerSpectrum *Clini;
	mscsMap map("map"),dirsmap("dirsmap");
	mscsAlms a("alms");
	parseOptions(argc,argv);
	//----------------------------------------------------------------------------------------------------
	
	
	
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// FIND BLOCK -- performs search for almn
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	if (_find) {
		init_A(&A,_Nsim);
		aini = new mscsAlms[_Nsim];
		Clini = new mscsAngularPowerSpectrum[_Nsim];
		mini = new mscsMap[_Nsim];


		
#ifdef GOMPI
		ierr = MPI_Init ( &argc, &argv ); // initiate MPI
		ierr = MPI_Comm_size ( MPI_COMM_WORLD, &node_num ); // Get the number of processes
		ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &world_id ); // Get the individual process ID
#endif
		
		//
		// read in the list of alm files onto the aini structure; initiate variables and data
		//
		for (j=0;j<_Nsim;j++) {
			
			//read in the alms file
			sprintf(tmpch,"alms sim %li",j+_sim_st);    aini[j].setName(tmpch);
#ifdef GOMPI
			sprintf(tmpch,"node:%i, alms sim %li",world_id,j+_sim_st);    aini[j].setName(tmpch);
			sprintf(tmpch,"node:%i, map",world_id);    map.setName(tmpch);
			sprintf(tmpch,"node:%i, dirmap",world_id);    map.setName(tmpch);
#endif
			
			if (_sim_st==-1) { sprintf(tmpch,"%s%s",_data_dir.c_str(),_sim_pref.c_str()); } 
			else sprintf(tmpch,"%s%s%li%s",_data_dir.c_str(),_sim_pref.c_str(),j+_sim_st,_sim_suff.c_str());
			
			aini[j].lmax(_l); aini[j].loadbinAlms(tmpch,"RI"); aini[j].lmax(_LMAX); exit(0); // BLcomment (Jan 25, 2011, 2:22:07 PM): this is nonsens - this was worngly ported from the previous version. the change of lmax does not preserve the ordering of the alms on the list as it erases then from the beginning.
			Clini[j]=aini[j].get_Cl(0,aini[j].lmax(),1); _l=_LMAX;
			exit(0);
			// make map Tj
			mini[j].set_nside(_ns);
			mini[j].SH_synthesis(aini[j],_l);
			
			// initiate the Almnj structure
			alms_num = aini[j].almsNum(); lmax=aini[j].lmax();
			printf("\n\n mpref> INITIATING Aj STRUCTURE FOR %li ALMS CORRESPONDING TO almspos: %li and LMAX: %li\n\n",alms_num,get_almsnumpos(alms_num),lmax);
			initiate_Aj(&A[j],get_almsnumpos(alms_num));
			
			// set the Almnj structure
			copy_alms(aini[j],&A[j]); // copies all alms table from the object to A[j].a structure
			calculate_rs(&A[j],aini[j]); // calculates the |alm^2|/((2l+1) C_l) for all lms in A[j] using info about Cl from object aini[j]
			n.lon()=0; n.lat()=0;
			for (k=0;k<A[0].alms_num;k++) { A[j].n[k]=A[j].nn[k]=n.get_direction(); /* A[j].i[k]=0; */ } // zero the n and i tables
		}
		
		
		//
		// initiate map and coordmap for search; and prepare map for various operations
		//
		map.set_nside(_ns); map.makekill_space_manager("make","T",1);
		a.lmax(lmax); 
		printf("\n\nmpref> the search ns is: %li\n\n",_search_ns);
		dirsmap.set_nside(_search_ns); dirsmap.makekill_space_manager("make","C",1);  dirsmap.set_map_coord(0,0);
		dirsmap.conv_nest2ring();
		
		
#ifndef GOMPI
		// initiate the mapping R_lm files
		if (_mapYlm) { 
			frmap = new FILE*[alms_num];
			for (k=0;k<A[0].alms_num;k++) { 
				aini[0].num2alm(get_almsnumneg(aini[0].almsNum())+k,&l,&m);
				sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li--rmapY%li_%li-Tr-txt",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Almnr",_Nsim,_sim_st,_sim_en,_ns,_search_ns,l,m);
				frmap[k] = fopen(tmpch,"w");
			}
		}
#endif
		
		
		
		//***********************************************
		//***********************************************
		// main loop over the all directions in ns=64 / 2
		//***********************************************
		//***********************************************
		Ndir = 12*_search_ns*_search_ns / 2; // we don't need to go over the full sphere; half is enough
		/*   Ndir = 12*_search_ns*_search_ns; */
		
		
		/*  PARALLEL IMPLEMENTATION -- BEGIN */
		
		// HERE THE PARALLEL IMPLEMENTATION WILL PARARELIZE THE LOOP SO THAT EACH PROCESS
		// WILL TAKE CARE OF DIFFERENT SET OF DIRECTIONS
		
#ifdef GOMPI
		
		if ( world_id == 0 ) {
			printf ( "MPI note>  ************** The number of processes is %i. ****************\n", node_num );
		}
		
		// initiate the runs
		Ndirpn=Ndir/node_num;
		dir_st=world_id*Ndirpn;
		dir_en=(world_id+1)*Ndirpn;
		if ( world_id == node_num-1 ) { dir_en=Ndir; } // make sure all directions will be processed
		printf("MPI note>  node nr: %i will process %li rotations, from %li to %li.\n",world_id,dir_en-dir_st,dir_st, dir_en);
		
		for (i=dir_st;i<dir_en;i++) {
			
			// prepare and store rotation
			n=dirsmap.get_C(i); phi=n.l()*PI180inv; th=90.0-n.b()*PI180inv;
			map.prepare_rotation(th,phi,0);
			
			//
			// main loop over the simulations
			//
			for (j=0;j<_Nsim;j++) {
				map.import_map_data(mini[j],"T",1); //aini[j].makekill_space_manager("kill","T",1);
				
				// rotate map
				map.rotate_quick("T",false);
				
				// calculate rotated alms and power spectrum for better accuracy
//				map.calculate_transformF(1,1,"",-1,0,"",""); map.calculate_C_l(0,map.lmax,1);
				a=map.SH_analysis(lmax); //map.calculate_C_l(0,map.lmax,1);
				
				// update the Almn structure: replace the modes that consumed more power in this orientation than the one stored
				printf("\n\nnode:%i mpref> updating the Almn database for rotation %li/%li\n\n",world_id,i,Ndir);
				update_Almn(&A[j],a,n); 
			}
		}
		
		// need to gather the information from all runs here into A structure of the master process
		// loop over all elements of the A structure
		
		/*   // maximalize the r results over the slaves: approach I  - sending/receiving packets of size one double -- THIS COMMUNICATION IS SLOW */
		/*   if (node_num > 1) { */
		/*     for (j=0;j<_Nsim;j++) { */
		/*       for (i=0;i<A[j].alms_num;i++) { */
		
		/* 	if ( world_id != 0 ) { // if I'm slave then send the info to master */
		/* 	  MPI_Send(&(A[j].a[i].R),1,MPI_DOUBLE,0,1,MPI_COMM_WORLD); // a.R */
		/* 	  MPI_Send(&(A[j].a[i].I),1,MPI_DOUBLE,0,2,MPI_COMM_WORLD); // a.I */
		/* 	  MPI_Send(&(A[j].aa[i].R),1,MPI_DOUBLE,0,3,MPI_COMM_WORLD); // aa.R */
		/* 	  MPI_Send(&(A[j].aa[i].I),1,MPI_DOUBLE,0,4,MPI_COMM_WORLD); // aa.I */
		/* 	  MPI_Send(&(A[j].r[i]),1,MPI_DOUBLE,0,5,MPI_COMM_WORLD); // r */
		/* 	  MPI_Send(&(A[j].rr[i]),1,MPI_DOUBLE,0,6,MPI_COMM_WORLD); // rr */
		/* 	  MPI_Send(&(A[j].n[i].l),1,MPI_DOUBLE,0,7,MPI_COMM_WORLD); // n.l */
		/* 	  MPI_Send(&(A[j].n[i].b),1,MPI_DOUBLE,0,8,MPI_COMM_WORLD); // n.b */
		/* 	  MPI_Send(&(A[j].nn[i].l),1,MPI_DOUBLE,0,9,MPI_COMM_WORLD); // nn.l */
		/* 	  MPI_Send(&(A[j].nn[i].b),1,MPI_DOUBLE,0,10,MPI_COMM_WORLD); // nn.b */
		/* 	} */
		/* 	else { // if I'm master then receive the info from slaves */
		/* 	  for (k=1;k<node_num;k++) { */
		/* 	    MPI_Recv(&from_slave,1,MPI_DOUBLE,k,5,MPI_COMM_WORLD,&mpi_stat); */
		/* 	    if (from_slave > A[j].r[i]) { A[j].r[i]=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,1,MPI_COMM_WORLD,&mpi_stat); A[j].a[i].R=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,2,MPI_COMM_WORLD,&mpi_stat); A[j].a[i].I=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,7,MPI_COMM_WORLD,&mpi_stat); A[j].n[i].l=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,8,MPI_COMM_WORLD,&mpi_stat); A[j].n[i].b=from_slave; */
		/* 	    } */
		/* 	    else { */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,1,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,2,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,7,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,8,MPI_COMM_WORLD,&mpi_stat); */
		/* 	    } */
		/* 	    MPI_Recv(&from_slave,1,MPI_DOUBLE,k,6,MPI_COMM_WORLD,&mpi_stat); */
		/* 	    if (from_slave < A[j].rr[i]) { A[j].rr[i]=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,3,MPI_COMM_WORLD,&mpi_stat);  A[j].aa[i].R=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,4,MPI_COMM_WORLD,&mpi_stat);  A[j].aa[i].I=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,9,MPI_COMM_WORLD,&mpi_stat);  A[j].nn[i].l=from_slave; */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,10,MPI_COMM_WORLD,&mpi_stat); A[j].nn[i].b=from_slave; */
		/* 	    } */
		/* 	    else { */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,3,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,4,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,9,MPI_COMM_WORLD,&mpi_stat); */
		/* 	      MPI_Recv(&from_slave,1,MPI_DOUBLE,k,10,MPI_COMM_WORLD,&mpi_stat); */
		/* 	    } */
		
		/* 	  } */
		/* 	} */
		
		/*       } */
		/*     } */
		/*   } */
		
		
		
		// maximalize the r results over the slaves: approach II -- try to reduce number of calls -- THIS IS FAST
		if (node_num > 1) {
			if (world_id == 0) {
				from_slave1 = new double[A[0].alms_num*2];
				from_slave3 = new double[A[0].alms_num*2];
				from_slave5 = new double[A[0].alms_num];
				from_slave6 = new double[A[0].alms_num];
				from_slave7 = new double[A[0].alms_num*2];
				from_slave9 = new double[A[0].alms_num*2];
			}
			
			for (j=0;j<_Nsim;j++) {
				
				if ( world_id != 0 ) { // if I'm slave then send the info to master
					MPI_Send(&(A[j].a[0].R), A[j].alms_num*2,MPI_DOUBLE,0,1,MPI_COMM_WORLD); // a.R,I
					MPI_Send(&(A[j].aa[0].R),A[j].alms_num*2,MPI_DOUBLE,0,3,MPI_COMM_WORLD); // aa.R,I
					MPI_Send(&(A[j].r[0]), A[j].alms_num,MPI_DOUBLE,0,5,MPI_COMM_WORLD); // r
					MPI_Send(&(A[j].rr[0]),A[j].alms_num,MPI_DOUBLE,0,6,MPI_COMM_WORLD); // rr
					MPI_Send(&(A[j].n[0].l), A[j].alms_num*2,MPI_DOUBLE,0,7,MPI_COMM_WORLD); // n.l
					MPI_Send(&(A[j].nn[0].l),A[j].alms_num*2,MPI_DOUBLE,0,9,MPI_COMM_WORLD); // nn.l
					
				}
				else { // if I'm master then receive the info from slaves
					for (k=1;k<node_num;k++) {
						MPI_Recv(&from_slave1[0],A[j].alms_num*2,MPI_DOUBLE,k,1,MPI_COMM_WORLD,&mpi_stat);
						MPI_Recv(&from_slave5[0],A[j].alms_num,MPI_DOUBLE,k,5,MPI_COMM_WORLD,&mpi_stat);
						MPI_Recv(&from_slave7[0],A[j].alms_num*2,MPI_DOUBLE,k,7,MPI_COMM_WORLD,&mpi_stat);
						
						for (i=0;i<A[j].alms_num;i++) {
							if (from_slave5[i] > A[j].r[i]) { 
								A[j].r[i]=from_slave5[i];
								A[j].a[i].R=from_slave1[2*i];
								A[j].a[i].I=from_slave1[2*i+1];
								A[j].n[i].l=from_slave7[2*i];
								A[j].n[i].b=from_slave7[2*i+1];
							}
						}
						
						MPI_Recv(&from_slave3[0],A[j].alms_num*2,MPI_DOUBLE,k,3,MPI_COMM_WORLD,&mpi_stat);
						MPI_Recv(&from_slave6[0],A[j].alms_num,MPI_DOUBLE,k,6,MPI_COMM_WORLD,&mpi_stat);
						MPI_Recv(&from_slave9[0],A[j].alms_num*2,MPI_DOUBLE,k,9,MPI_COMM_WORLD,&mpi_stat);
						
						for (i=0;i<A[j].alms_num;i++) {
							if (from_slave6[i] < A[j].rr[i]) { 
								A[j].rr[i]=from_slave6[i];
								A[j].aa[i].R=from_slave3[2*i];
								A[j].aa[i].I=from_slave3[2*i+1];
								A[j].nn[i].l=from_slave9[2*i];
								A[j].nn[i].b=from_slave9[2*i+1];
							}
						}
						
					}
				}
				
			}
			if (world_id == 0) {
				delete [] from_slave1;
				delete [] from_slave3;
				delete [] from_slave5;
				delete [] from_slave6;
				delete [] from_slave7;
				delete [] from_slave9;
			}
		}
		
		
		
		if ( world_id == 0) {
			
			//
			// save Almn
			//
			printf("\n\nmpref> SEARCHING PART DONE\n\n");
			if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
			sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Almnr",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
			save_Almn(A,tmpch);
			
			//***********************************************
			//***********************************************
			// maximizing over m in each multipole and each simulation
			//***********************************************
			//***********************************************
			printf("\n\nmpref> STARTING THE MAXIMIZING PART\n\n");
			
			// init the D structure
			init_D(&D,_Nsim);
			
			// set lmax
			lmax = A[0].lmax;
			printf("\n\nmpref> lmax is: %li\n\n",lmax);
			
			for (j=0;j<_Nsim;j++) {
				// initiate the Dlnj sub-structures
				initiate_Dj(&D[j],lmax);
				
				// find maximal r and put onto D structure
				maximize_Almn(&A[j],&D[j]);
				//kill_Aj(&A[j]);
			}
			
			// save the maximized Dln
			if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
			sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Dlnr",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
			save_Dln(D,tmpch);
		}
		
		
		//
		// clean up
		//
		_find_done=true;
		if (!_something_else()) { kill_A(&A); kill_D(&D); }
		
		ierr = MPI_Finalize ( ); // Shut down MPI.
		
#endif
		/*  PARALLEL IMPLEMENTATION -- END */
		
		
		
		
		
		
		// SINGLE PROCESSOR IMPLEMENTATION - START
		// this part is kept only in case a different compiler is used to compile the same source code, which is general
		// now for both c++ and mpic++ compilations
#ifndef GOMPI
		for (i=0;i<Ndir;i++) {
			/*   for (i=0;i<1;i++) { */
			
			// prepare and store rotation
			n=dirsmap.get_C(i); phi=n.l()*PI180inv; th=90.0-n.b()*PI180inv;
			map.prepare_rotation(th,phi,0);
			
			//
			// main loop over the simulations
			//
			for (j=0;j<_Nsim;j++) {
				map.import_map_data(mini[j],"T",1); //aini[j].makekill_space_manager("kill","T",1);
				
				// rotate map
				map.rotate_quick("T",false);
				
				// calculate rotated alms and power spectrum for better accuracy
				a=map.SH_analysis(lmax); //map.calculate_C_l(0,map.lmax,1);
				
				// update the Almn structure: replace the modes that consumed more power in this orientation than the one stored
				printf("\n\nmpref> updating the Almn database for rotation %li/%li\n\n",i,Ndir);
				if (!_mapYlm) update_Almn(&A[j],a,n);
				else update_Almn(&A[j],a,n,frmap); 
			}
		}
		
		//
		// save Almn
		//
		printf("\n\nmpref> SEARCHING PART DONE\n\n");
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Almnr",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		save_Almn(A,tmpch);
		
		if (_mapYlm) {
			for (j=0;j<A[0].alms_num;j++) { 
				for (i=0;i<Ndir;i++) { fprintf(frmap[j],"%lE\n",0.0); } // save the other half of the map
				fclose(frmap[j]);
			}
		}
		
		
		//***********************************************
		//***********************************************
		// maximizing over m in each multipole and each simulation
		//***********************************************
		//***********************************************
		printf("\n\nmpref> STARTING THE MAXIMIZING PART\n\n");
		
		// init the D structure
		init_D(&D,_Nsim);
		
		// set lmax
		lmax = A[0].lmax;
		printf("\n\nmpref> lmax is: %li\n\n",lmax);
		
		for (j=0;j<_Nsim;j++) {
			// initiate the Dlnj sub-structures
			initiate_Dj(&D[j],lmax);
			
			// find maximal r and put onto D structure
			maximize_Almn(&A[j],&D[j]);
			//kill_Aj(&A[j]);
		}
		
		// save the maximized Dln
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Dlnr",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		save_Dln(D,tmpch);
		
		
		//
		// clean up
		//
		_find_done=true;
		if (!_something_else()) { kill_A(&A); kill_D(&D); }
		
#endif    
		// SINGLE PROCESSOR IMPLEMENTATION - END
		
		
	}
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// END OF FIND BLOCK
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates angles between maxmized r's in multipoles
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	if (_calc_Dstats) {
		
		//***********************************************
		//***********************************************
		// calculate angles between the multipoles in the maximized data
		//***********************************************
		//***********************************************
		printf("\n\nmpref> STARTING THE ESTIMATOR CALCULATION PART\n\n");
		
		// check for the Dlmn data if the _find section was not done
		if (!_find) {
			// read the maximized Dln
			sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"Dlnr",_Nsim,_sim_st,_sim_en);
			read_Dln(&D,tmpch);
		}
		
		// init the stat structure;
		init_stat(&stat,_Nsim);
		// set lmax
		lmax=D[0].lmax;
		printf("\n\nmpref> lmax is: %li\n\n",lmax);
		
		
		//
		// loop over maps
		//
		
		// calculate estimators:  alpha_l^1 and Alpha_l^1 = weighted alpha_l^1
		// alpha_l^i = ^n_l*^n_l+i -i.e. dot product of the unit n_lm "axes" not vectors --  ang in [0..90]
		// Alpha_l^i = n_l*n_l+1/r_l/r_l+i -i.e. an r weighted dot product of the unit n_lm "axes" not vectors --  ang in [0..90]
		for (j=0;j<_Nsim;j++) {
			// initiate the statlj sub-structures
			initiate_statj(&stat[j],lmax);
			
			lmax--;
			for (l=0;l<=lmax;l++) {
				ang = PI180inv * get_angle(D[j].n[l],D[j].n[l+1]);
				stat[j].v[l]=ang;
				stat[j].V[l]=ang/D[j].r[l]/D[j].r[l+1];
			}
			lmax++;
		}
		// save the stat
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"a1A1jl",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"a1A1jl",_Nsim,_sim_st,_sim_en);
		save_stat(stat,tmpch,2,lmax-1); // save from 2..20
		
		// calculate estimators:  alpha_l^2 and Alpha_l^2 = weighted alpha_l^2
		lmax-=2;
		for (j=0;j<_Nsim;j++) {
			for (l=0;l<=lmax;l++) {
				ang=get_angle(D[j].n[l],D[j].n[l+2]);
				stat[j].v[l]=PI180inv*ang;
				stat[j].V[l]=PI180inv*ang/D[j].r[l]/D[j].r[l+2];
			}
		}
		lmax+=2; // birng back the initial lmax value
		// save the stat
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"a2A2jl",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"a2A2jl",_Nsim,_sim_st,_sim_en);
		save_stat(stat,tmpch,2,lmax-2); // save from 2..19
		
		// calculate estimators:  omega_l^N, N={3,4,5,6} and  Omega_l^N, N={3,4,5,6} = weighted omega_l^N
		
		
		//
		// clean up
		//
		_calc_Dstats_done=true;
		if (!_something_else()) { kill_A(&A); kill_D(&D); }
		
	}
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates angles between maxmized r's in multipoles END
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates statistics on Almn data; does the full cov statistics to measure probability of getting given n_lm value
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	// * full cov statistics measuring probability of n_lm values in the data
	//   ------------------------------------------------------
	//  C_ll' = 1/(N-1) Sum_sims ( 1/( (lmax+1)(lmax'+1) ) * Sum_m Sum_m' n_lm*n_l'm'
	//  where n_lm = ^n_lm * r_lm ; r_lm = |a_lm|^2/((2l+1) C_l)
	//  chisqn = n_l * C_ll'^-1 n_l'
	//  where n_l = 1/(l+1) Sum_m n_lm
	//  where n_lm is a vector found by the _find block
	
	
	// The output will return a set of chisq values for each sim from each of these tests
	// simID chisq_n chisq_nn chisq_nnn
	//----------------------------------------------------------
	
	if (_calc_Covstats_test) {
		printf("\n\nmpref> STARTING THE FULL COV. STATISTICS PART\n\n");
		
		//
		// chceck for the Almn structure first
		//
		
		if (!_find) {
			sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"Almnr",_Nsim,_sim_st,_sim_en);
			read_Almn(&A,tmpch);
		}
		
		//
		// initiate relevant variables
		//
		printf("\n\nmpref> * initiating variables\n\n");
		LMIN=_LMIN; // minimal value of the multipole that will contribute to the chisq values
		lmax=_LMAX;
		NsimC=50; // number of simulations used for Cov calculation
		NsimMCPDF=_Nsim-NsimC-1; // number of sims used for MC chisq PDF probing; -1 is from the fact that the 0'th slot is reserved for data
		printf("\n\nmpref>  -- LMIN: %li, lmax: %li, NsimC: %li, NsimMCPDF: %li\n\n", LMIN,lmax,NsimC, NsimMCPDF);
		chisq = new double[NsimMCPDF+1]; // the first value will be reserved for the data value; +1 to store also the data value
		
		//
		// initiate matrices
		//
		dimnl=lmax-LMIN+1;
		printf("\n\nmpref>  -- initiating vectors and matrices: size=%li\n\n",dimnl);
		nl = new matrix <double>; (*nl).SetSize(dimnl,3);
		ml = new matrix <double>; (*ml).SetSize(dimnl,3);
		Covl = new matrix <double>; (*Covl).SetSize(dimnl,dimnl);
		iCovl = new matrix <double>; (*iCovl).SetSize(dimnl,dimnl);
		check = new matrix <double>; (*check).SetSize(dimnl,dimnl);
		
		
		//
		// calculate Cov. from NsimC <= _Nsim simulations  : this assumes that the first entry of the A = A[0] is restricted for the data and is not used for Cov calculations; first NsimC following j'th sim. sims will be used for Cov calculation
		//
		printf("\n\nmpref> * calculating cov\n\n");
		j=_Nsim-NsimMCPDF;
		calculate_covariance_matrix(A,Covl,LMIN,lmax,NsimC,j,"Cllp");
		
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"Covll_test",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"Covll_test",_Nsim,_sim_st,_sim_en);
		cpeds_matrix_save(*Covl,tmpch);
		
		//
		// invert Covl
		//
		printf("\n\nmpref> * inverting cov\n\n");
		// add checking for the condition numbers (eg. with octave)
		
		if (!(*Covl).IsSingular()) { (*iCovl)=!(*Covl); } else { printf("******* THE COV MATRIX IS SINGULAR\n"); exit(0); }
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"iCovll_test",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"iCovlll_test",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		cpeds_matrix_save(*iCovl,tmpch);
		(*check)=(*Covl)*(*iCovl);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"CoviCovll_check_test",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"CoviCovll_check_test",_Nsim,_sim_st,_sim_en);
		cpeds_matrix_save(*check,tmpch);
		
		
		
		//
		//    derive chisq for i'th simulation (data)
		//
		printf("\n\nmpref> * deriving chisq values\n\n");
		for (i=0;i<=NsimMCPDF;i++) {
			/*       printf("doing sim:%li\n",i); */
			// set the data vector: copy the data from A[j] to data matrix
			copy_nl_to_vec(&A[i],nl,LMIN,lmax);
			
			// derive ml=C^-1 * nl
			(*ml) = (*iCovl)*(*nl);
			
			// derive chisq=nl*ml
			chisq[i]=0;
			for (l=LMIN;l<=lmax;l++) {  
				ln=l-LMIN;  tmpd = (*nl)(ln,0)*(*ml)(ln,0) + (*nl)(ln,1)*(*ml)(ln,1) + (*nl)(ln,2)*(*ml)(ln,2); chisq[i] +=tmpd; 
				/* 	printf("chisq per l: %li in sim %li:%lE\n",l,i,tmpd);      */
			}
			
			printf("chisq froml: %li sim. is:%lE\n",i,chisq[i]);     
		}  
		
		
		
		//
		// save the chisq values to file from i'th to j'th value in the tab
		//
		i=1; j=NsimMCPDF;
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"chisq_sim_test",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"chisq_sim_test",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		save_chisq_vals(chisq,i,j,tmpch);
		i=0; j=0;
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"chisq_data_test",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"chisq_data_test",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		save_chisq_vals(chisq,i,j,tmpch);
		
		
		// 
		// clean up
		//
		
		delete [] chisq;
		delete Covl; delete iCovl; delete check; delete nl; delete ml;
		_calc_Covstats_test_done=true;
		if (!_something_else()) { kill_A(&A); kill_D(&D); }
	}
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates statistics on Almn data; does the full cov statistics END
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates statistics on Almn data; does the full cov statistics. Measures inter-mode alignments n_lm n_l'm'
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	// * full cov statistics measuring all inter-n_lm-mode correlations (not only inter-l correlations)
	//   ------------------------------------------------------
	
	//  n_lm = ^n_lm * r_lm ; r_lm = |a_lm|^2/((2l+1) C_l)
	//  define a column vector n_lml'm' = n_i = n_lm*n_l'm'  
	//  where n_lm is a vector found by the _find block
	//  where i numerates all combinations of l,m,l',m' from LMIN to lmax in ordering eg.
	//  
	//  2 0 3 0, 2 0 3 1, 2 0 3 2, 2 0 3 3, 2 1 3 0, ... 2 1 3 3,..., 2 2 3 3 ... , 3 0 4 0 ..., 3 3 4 4 , ...,..., 
	//
	//  then define the correlation matrix of such vectors from simulations
	//  C_ij = C_lml'm',LML'M' = 1/(N-1) Sum_sims ( n_lml'm' * n_LML'M' ) = 1/(N-1) Sum_sims ( n_i * n_j )
	//
	//  define chisq value as:
	//  chisq = n_i * C_ij^-1 n_j
	
	
	// The output will return a set of chisq values for each simulation including data
	// simID chisq
	//----------------------------------------------------------
	
	if (_calc_Covstats_inter_mode) {
		printf("\n\nmpref> STARTING THE FULL COV. STATISTICS PART: MEASURING INTER-n_lm-MODES CORRELATIONS\n\n");
		
		//
		// chceck for the Almn structure first
		//
		
		if (!_find) {
			sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"Almnr",_Nsim,_sim_st,_sim_en);
			read_Almn(&A,tmpch);
		}
		
		//
		// initiate relevant variables
		//
		printf("\n\nmpref> * initiating variables\n\n");
		LMIN=_LMIN; // minimal value of the multipole that will contribute to the chisq values
		lmax=_LMAX;
		NsimC=100; // number of simulations used for Cov calculation
		NsimMCPDF=_Nsim-NsimC-1; // number of sims used for MC chisq PDF probing; -1 is from the fact that the 0'th slot is reserved for data
		printf("\n\nmpref>  -- LMIN: %li, lmax: %li, NsimC: %li, NsimMCPDF: %li\n\n", LMIN,lmax,NsimC, NsimMCPDF);
		chisq = new double[NsimMCPDF+1]; // the first value will be reserved for the data value; +1 to store also the data value
		
		//
		// initiate matrices
		//
		dimnl=get_nlmlpmp_vec_size(LMIN,lmax);
		printf("\n\nmpref>  -- initiating vectors and matrices: size=%li\n\n",dimnl);
		nlmlpmp = new matrix <double>; (*nlmlpmp).SetSize(dimnl,1);
		mlmlpmp = new matrix <double>; (*mlmlpmp).SetSize(dimnl,1);
		Covl = new matrix <double>; (*Covl).SetSize(dimnl,dimnl);
		iCovl = new matrix <double>; (*iCovl).SetSize(dimnl,dimnl);
		check = new matrix <double>; (*check).SetSize(dimnl,dimnl);
		
		//
		// calculate Cov. from NsimC <= _Nsim simulations  : this assumes that the first entry of the A = A[0] is restricted for the data and is not used for Cov calculations; first NsimC following j'th sim. sims will be used for Cov calculation
		//
		printf("\n\nmpref> * calculating cov\n\n");
		j=_Nsim-NsimMCPDF;
		calculate_covariance_matrix(A,Covl,LMIN,lmax,NsimC,j,"ClmlpmpLMLPMP");
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"CovlmlpmpLMLPMP",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"ClmlpmpLMLPMP",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		cpeds_matrix_save(*Covl,tmpch);
		
		//
		// invert Covl
		//
		printf("\n\nmpref> * inverting cov\n\n");
		// add checking for the condition numbers (eg. with octave)
		
		if (!(*Covl).IsSingular()) { (*iCovl)=!(*Covl); } else { printf("******* THE COV MATRIX IS SINGULAR\n"); exit(0); }
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"iCovlmlpmpLMLPMP",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"iCovlmlpmpLMLPMP",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		cpeds_matrix_save(*iCovl,tmpch);
		(*check)=(*Covl)*(*iCovl);
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"CovlmlpmpLMLPMP_check",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"CovlmlpmpLMLPMP",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		cpeds_matrix_save(*check,tmpch);
		
		
		
		//
		//    derive chisq for i'th simulation (data)
		//
		printf("\n\nmpref> * deriving chisq values\n\n");
		for (i=0;i<=NsimMCPDF;i++) {
			/*       printf("doing sim:%li\n",i); */
			// set the data vector: copy the data from A[j] to data matrix
			copy_nlmlpmp_to_vec(&A[i],nlmlpmp,LMIN,lmax);
			// derive ml=C^-1 * nl
			
			(*mlmlpmp) = (*iCovl)*(*nlmlpmp);
			
			// derive chisq=nl*ml
			(*mlmlpmp)=~(*nlmlpmp)*(*mlmlpmp);
			
			chisq[i]=(*mlmlpmp)(0,0);
			
			
		}
		
		
		
		//
		// save the chisq values to file from i'th to j'th value in the tab
		//
		i=1; j=NsimMCPDF;
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"chisq_lmlpmp",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"chisq_lmlpmp",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		save_chisq_vals(chisq,i,j,tmpch);
		i=0; j=0;
		//sprintf(tmpch,"%s-%s-Nsim_%li-st_%li-en_%li",_output_file.c_str(),"chisq_lmlpmp_data",_Nsim,_sim_st,_sim_en);
		if (_sim_st == -1) tmps=""; else tmps=_sim_suff;
		sprintf(tmpch,"%s--%sxxx%s--%s-Nsim_%li-st_%li-en_%li--ns%li--sns%li",_output_file.c_str(),_sim_pref.c_str(),tmps.c_str(),"chisq_lmlpmp_data",_Nsim,_sim_st,_sim_en,_ns,_search_ns);
		
		save_chisq_vals(chisq,i,j,tmpch);
		
		
		// 
		// clean up
		//
		
		delete [] chisq;
		delete Covl; delete nl; delete ml;
		_calc_Covstats_inter_mode_done=true;
		if (!_something_else()) { kill_A(&A); kill_D(&D); }
	}
	
	
	// * variance tests of the maximized/minimized maps
	//----------------------------------------------------
	// derives: sigman_l = (2l+1)/4pi * Sum_m |r_lm|^2 and sigmann_l = (2l+1)/4pi * Sum_m |rr_lm|^2
	// if there's no m-prefference then there should be no extra power in any multipole w.r.t. sims
	// the input C_l anomalies do not bias this statistic
	
	// output: simID l sigman_l sigmann_l
	//----------------------------------------------------
	
	
	
	// Similarly for minimized power nn_lm's. -- THIS ISN'T USEFULL SINCE HAVEING ZEROED r for ONE MODE DOESN'T MEAN MUCH ABOUT ALL THE REST m's SO IT'S NOT IMPLEMENTED
	//  CC_ll' = 1/(N-1) Sum_sims ( 1/( (lmax+1)(lmax'+1) ) * Sum_m Sum_m' nn_lm*nn_l'm'
	//  where n_lm = ^nn_lm * rr_lm ; rr_lm = (|aa_lm|^2/((2l+1) C_l))^-1
	//  chisqnn = nn_l * CC_ll'^-1 nn_l'
	//  where nn_l = 1/(l+1) Sum_m nn_lm
	//  where nn_lm is a vector found by the _find block
	
	// Similarly for the combined statistics for min/max data ---- NOT SURE YET THIS MAKES ANY SENSE SO NOT IMPLEMENTED
	//  CCC_ll' = 1/(N-1) Sum_sims ( 1/( 4(lmax+1)(lmax'+1) ) * Sum_m Sum_m' (n_lm+nn_lm)*(n_l'm'+nn_l'm')
	//  where n_lm = ^nn_lm * rr_lm ; rr_lm = (|aa_lm|^2/((2l+1) C_l))^-1
	//  chisq = nn_l * CC_ll'^-1 nn_l'
	//  where nn_l = 1/(l+1) Sum_m n_lm+nn_lm
	//  
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates statistics on Almn data; does the full cov statistics. Measures inter-mode alignments n_lm n_l'm'  END
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	// STATISTICS BLOCK -- calculates statistics on Almn data ; 
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	//***********************************************
	
	
	
	
	// * test of Gaussianity in the power quilized-maximized maps --- NOT SURE THIS MAP SHOULD BE GAUSSIAN
	// * r_lm --> Tmax: Tmax should still be a gaussian map that can be tested with Mink. funct.
	// * rr_lm --> Tmin: Tmin should still be a gaussian map that can be tested with Mink. funct.
	
	// output: Tmax, Tmin maps for each simulation for further tests
	
	
	
	if (_calc_Astats) {
		
		
		
	}
	
	
	
}


// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************



// *************************************************************************************************
// WARNINIG !! the ii defines the direction used for rotation, NOT the axis that maximize the power in mode
// the n too but is converted so that it point the preferred axis
void update_Almn(Almn_type *Aj,mscsAlms &a, cpedsDirection n) {
	long i,l,m,ist,ileist;
	long alms_num=a.almsNum();
	long lmax=(*Aj).lmax;
	double r;
	
	ist=alms_num-get_almsnumpos(alms_num);
	n.lon()+=PI; if (n.l() >= twoPI) n.lon()-=twoPI; // convert to get the axis position
	
	/*   n.b=PIsnd-n.b; // convert to get the axis position */
	if (n.b()<0) { n.lat()=-n.b(); n.lon()+=PI; if (n.l() >= twoPI) n.lon()-=twoPI; } // keep the axes in the NH
	
	for (i=ist;i<alms_num;i++) { 
		a.num2alm(i,&l,&m);
		r = (a.get(i).R() * a.get(i).R()   +   a.get(i).I() * a.get(i).I()) / ((2.0*(double)l+1.0)*a.calculate_single_C_l(l)); // using the recalculated Cl from numerical stabilities reasons even though Cl should be rotationally invariant
		if (m!=0) r*=2.0;
		
		ileist=i-ist;
		if (r>(*Aj).r[ileist]) { (*Aj).a[ileist]=a.get(i).get_a_lm(); (*Aj).r[ileist] =r; (*Aj).n[ileist] =n.get_direction(); /* (*Aj).i[i]=ii; */ } // maximize the r
		if (r<(*Aj).rr[ileist]) { (*Aj).aa[ileist]=a.get(i).get_a_lm(); (*Aj).rr[ileist]=r; (*Aj).nn[ileist]=n.get_direction(); } // minimize the rr
	}
}

void update_Almn(Almn_type *Aj,mscsAlms &a, cpedsDirection n, FILE** frmap) {
	long i,l,m,ist,ileist;
	long alms_num=a.almsNum();
	long lmax=(*Aj).lmax;
	double r;
	
	ist=alms_num-get_almsnumpos(alms_num);
	n.lon()+=PI; if (n.l() >= twoPI) n.lon()-=twoPI; // convert to get the axis position
	
	/*   n.b=PIsnd-n.b; // convert to get the axis position */
	if (n.b()<0) { n.lat()=-n.b(); n.lon()+=PI; if (n.l() >= twoPI) n.lon()-=twoPI; } // keep the axes in the NH
	
	for (i=ist;i<alms_num;i++) { 
		a.num2alm(i,&l,&m);
		r = (a.get(i).R() * a.get(i).R()   +   a.get(i).I() * a.get(i).I()) / ((2.0*(double)l+1.0)*a.calculate_single_C_l(l)); // using the recalculated Cl from numerical statilities reasons even though Cl should be rotationally invariant
		if (m!=0) r*=2.0;
		
		ileist=i-ist;
		if (r>(*Aj).r[ileist]) { (*Aj).a[ileist]=a.get(i).get_a_lm(); (*Aj).r[ileist] =r; (*Aj).n[ileist] =n.get_direction(); /* (*Aj).i[i]=ii; */ } // maximize the r
		if (r<(*Aj).rr[ileist]) { (*Aj).aa[ileist]=a.get(i).get_a_lm(); (*Aj).rr[ileist]=r; (*Aj).nn[ileist]=n.get_direction(); } // minimize the rr
		fprintf(frmap[ileist],"%lE\n",r);
	}
	
}

// *************************************************************************************************
void maximize_Almn(Almn_type *Aj,Dln_type *Dj) {
	long l,lmax,m,mmax,k,kmax,mmin,kmin,s;
	double r,rr;
	
	printf("\n\nmpref> MAXIMIZING Almn's OVER M MODES\n\n");
	lmax = (*Aj).lmax;
	
	s=get_almsnumtot((*Aj).alms_num)-(*Aj).alms_num; // number of all alm with m<0
	printf("\n\nmpref>  -- total number of alms: %li,  m>=0: %li,  m<0: %li, lmax=%li\n\n",get_almsnumtot((*Aj).alms_num),(*Aj).alms_num,s,lmax);
	
	for (l=0;l<=lmax;l++) {
		r=-1; rr=2;
		// find maximal/minimal r in l'th multipole
		for (m=0;m<=l;m++) {      k=alm2num(l,m,lmax)-s;      
		if ((*Aj).r[k]  > r ) { kmax=k; mmax=m; r =(*Aj).r[k]; } 
		if ((*Aj).rr[k] < rr) { kmin=k; mmin=m; rr=(*Aj).rr[k]; } 
		/*       printf("k: %li, l:%li m:%li r:%lE rr: %lE, Ajrr[k]:%lE kmin: %li kmax:%li\n",k,l,m,r,rr, (*Aj).rr[k], kmin,kmax); */
		} 
		
		// store the maximal values on Dj
		(*Dj).m[l] =mmax;    (*Dj).r[l] =(*Aj).r[kmax];     (*Dj).n[l]=(*Aj).n[kmax];    //(*Dj).i[l]=(*Aj).i[kmax];    
		(*Dj).mm[l]=mmin;    (*Dj).rr[l]=(*Aj).rr[kmin];    (*Dj).nn[l]=(*Aj).nn[kmin];  //  (*Dj).i[l]=(*Aj).i[kmax];    
	}
	printf("\n\nmpref> MAXIMIZING Almn's OVER M MODES DONE\n\n");
}

// *************************************************************************************************
double get_angle(cpedsDirection n1, cpedsDirection n2) {
	double ang=cpeds_ang_n1n2(PIsnd-n1.b(),n1.l(),PIsnd-n2.b(),n2.l());
	if (ang > PIsnd) ang=PI-ang;  // the angle between cpedsDirections cannot be larger than 90 deg.
	return ang;
}

// *************************************************************************************************
// calculates the covariance matrix from A vector using NsimC simulations starting from j'th 
// and minimal multipole number LMIN --> the size of the matrix will be lmax-LMIN+1
void calculate_covariance_matrix(Almn_type *A,matrix<double> *Covl,long lmin, long lmax, long NsimC, long ist, string how) {
	long i,iC,ien,l,m,lp,mp,lC,lpC;
	long j,jp,lpo;
	//long L,LP,M,MP,LC,MC,LPC,MPC;
	long dimC;
	
	matrix <double> * Cd;
	
	printf("\n\nmpref> CALCULATING THE COVARIANCE MATRIX\n\n");
	ien=ist+NsimC;
	
	if ( how == "Cllp" ) {
		
		for (l=lmin;l<=lmax;l++) {
			for (lp=lmin;lp<=l;lp++) { // calculating the lower triangular part of the matrix
				lC=l-lmin; lpC=lp-lmin;      (*Covl)(lC,lpC)=0;
				
				for (i=ist;i<ien;i++) {	
					for (m=0;m<=l;m++) {
						for (mp=0;mp<=lp;mp++) { (*Covl)(lC,lpC)+=dot_product(&A[i],l,m,lp,mp); } } 
				}
				(*Covl)(lC,lpC) /= ( (double)((l+1)*(lp+1)*(NsimC-1)) );
				
				(*Covl)(lpC,lC)=(*Covl)(lC,lpC); // set the upper triangular part of the matrix
				
			}
		}
		
	}
	
	if ( how == "ClmlpmpLMLPMP" ) {
		dimC=get_nlmlpmp_vec_size(lmin,lmax);
		Cd = new matrix<double>; (*Cd).SetSize(dimC,NsimC);
		
		// seed the Cd matrix with vectors from A structure
		printf("mpref>  -- making a sim vectors structure\n");
		for (i=ist;i<ien;i++) {	
			iC=i-ist;
			j=0;
			
			for (l=lmin;l<lmax;l++) {
				for (m=0;m<=l;m++) {
					lpo=l+1;
					
					for (lp=lpo;lp<=lmax;lp++) {
						for (mp=0;mp<=lp;mp++) {
							(*Cd)(j,iC)=dot_product(&A[i],l,m,lp,mp);
							j++;
						}
					}
					
				}
			}
			
		}
		cpeds_matrix_save(*Cd,"Cd");
		
		// calculate Cov from Cd matrix
		printf("mpref>  -- calculating Cov.\n");
		for (j=0;j<dimC;j++) {	
			for (jp=0;jp<=j;jp++) {  // calculating the lower triangular part of the matrix
				(*Covl)(j,jp)=0;
				for (i=ist;i<ien;i++) {	
					iC=i-ist;
					(*Covl)(j,jp)+=(*Cd)(j,iC)*(*Cd)(jp,iC);
				}
				(*Covl)(j,jp) /= ( (double)(NsimC-1) );
				(*Covl)(jp,j)=(*Covl)(j,jp); // set the upper triangular part of the matrix
			}
		}
		
		delete Cd;
	}
	
	printf("\n\nmpref> CALCULATING THE COVARIANCE MATRIX DONE\n\n");
}
// *************************************************************************************************
// calculates the dot product between the vectors n_lm and n_lpmp in the j'th simulation of the A structure
double dot_product(Almn_type *Aj, long l, long m, long lp, long mp) {
	double x1,x2,y1,y2,z1,z2,r1,r2,th,ang;
	long i,almsnumneg=(*Aj).alms_num-(*Aj).lmax-1;
	cpedsDirection n1,n2;
	
	/*   i=alm2num(l,m,(*Aj).lmax)-almsnumneg; */
	/*   th=PIsnd-(*Aj).n[i].b; */
	/*   x1=cpeds_sph2cart(0,th,(*Aj).n[i].l); */
	/*   y1=cpeds_sph2cart(1,th,(*Aj).n[i].l); */
	/*   z1=cpeds_sph2cart(2,th,(*Aj).n[i].l); */
	/*   r1=(*Aj).r[i]; */
	
	/*   i=alm2num(lp,mp,(*Aj).lmax)-almsnumneg; */
	/*   th=PIsnd-(*Aj).n[i].b; */
	/*   x2=cpeds_sph2cart(0,th,(*Aj).n[i].l); */
	/*   y2=cpeds_sph2cart(1,th,(*Aj).n[i].l); */
	/*   z2=cpeds_sph2cart(2,th,(*Aj).n[i].l); */
	/*   r2=(*Aj).r[i]; */
	
	/*   return (x1*x2+y1*y2+z1*z2)*r1*r2;   */
	
	i=alm2num(l,m,(*Aj).lmax)-almsnumneg;
	r1=(*Aj).r[i];
	n1=(*Aj).n[i];
	
	i=alm2num(lp,mp,(*Aj).lmax)-almsnumneg;
	r2=(*Aj).r[i];
	n2=(*Aj).n[i];
	
	ang=get_angle(n1,n2);
	return r1*r1*cos(ang);
}

// *************************************************************************************************
// copies the data from the j'th simulation in the A structure to nl vector: nl=sum_m n_lm
// converting from spherical coordinates to cartesian and summing over m's
void copy_nl_to_vec(Almn_type *Aj,matrix<double> *nl, long LMIN, long lmax) {
	long l,m,ln,i,almsnumneg=(*Aj).alms_num-(*Aj).lmax-1;
	double th;
	
	for (l=LMIN;l<=lmax;l++) {
		ln=l-LMIN; (*nl)(ln,0) = (*nl)(ln,1) = (*nl)(ln,2) = 0; 
		for (m=0;m<=l;m++) { 
			i=alm2num(l,m,(*Aj).lmax)-almsnumneg;
			th = PIsnd-(*Aj).n[i].b;
			(*nl)(ln,0) += cpeds_sph2cart(0,th,(*Aj).n[i].l) * (*Aj).r[i]; // x
			(*nl)(ln,1) += cpeds_sph2cart(1,th,(*Aj).n[i].l) * (*Aj).r[i]; // y
			(*nl)(ln,2) += cpeds_sph2cart(1,th,(*Aj).n[i].l) * (*Aj).r[i]; // z
		}
		(*nl)(ln,0) /= (double)(l+1);
		(*nl)(ln,1) /= (double)(l+1);
		(*nl)(ln,2) /= (double)(l+1);
	}
	
}
// *************************************************************************************************
// calculates the size of the 
long get_nlmlpmp_vec_size(long lmin,long lmax) {
	long l,m,lp,mp,i=0,lpo;
	
	for (l=lmin;l<lmax;l++) {
		for (m=0;m<=l;m++) {
			lpo=l+1;
			for (lp=lpo;lp<=lmax;lp++) {
				for (mp=0;mp<=lp;mp++) i++;
			}
		}
	}
	
	return i;
}
// *************************************************************************************************
void copy_nlmlpmp_to_vec(Almn_type *Aj, matrix<double> *nlmlpmp, long lmin, long lmax) {
	long l,m,lp,mp,lpo;
	long j;
	//long L,LP,M,MP,LC,MC,LPC,MPC;
	
	j=0;
	
	for (l=lmin;l<lmax;l++) {
		for (m=0;m<=l;m++) {
			lpo=l+1;
			
			for (lp=lpo;lp<=lmax;lp++) {
				for (mp=0;mp<=lp;mp++) {
					(*nlmlpmp)(j,0)=dot_product(Aj,l,m,lp,mp);
					j++;
				}
			}
			
		}
	}
	
	
}
// *************************************************************************************************








// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************

// IO functions

// *************************************************************************************************
// *************************************************************************************************
// *************************************************************************************************



void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd("finds the axes that maximize the power accumulated in an individual m-modes \n The parallel implementation of this program makes sens only with the --find option, with other options a parallel run will not result in any branching. Only the --find option is CPU intensive. ", ' ', Mscs_version.c_str() );
		
		// 
		// Define arguments
		//
		
		ValueArg<string> data_dir("","data_dir","directory where the alms files in the initial orientation are located",false,"","string"); cmd.add(data_dir);
		
		ValueArg<string> sim_pref("","sim_pref","prefix of the input files  (default: )",false,"","string"); cmd.add(sim_pref);
		ValueArg<string> sim_suff("","sim_suff","suffix of the input files  (default: "")",false,"","string"); cmd.add(sim_suff);
		ValueArg<long> sim_st("","sim_st","start id. if -1 is given then will use sim_pref as the name of the single file to perform search on it   (default: -1)",false,-1,"long"); cmd.add(sim_st);
		ValueArg<long> sim_en("","sim_en","end id   (default: 1000)",false,1000,"long"); cmd.add(sim_en);
		ValueArg<long> search_ns("","sns","ns number defining the density of directions to be probed   (default: 64)",false,64,"long"); cmd.add(search_ns);
		ValueArg<long> ns("","ns","ns of the maps to be generated   (default: 128)",false,128,"long"); cmd.add(ns);
		ValueArg<long> l("l","","lmax multipole to define the input files (ONLY)  (default: 21)",false,21,"long"); cmd.add(l);
		ValueArg<long> LMAX("","LMAX","lmax multipole to define range of statistics performed  (default: l)",false,-1,"long"); cmd.add(LMAX);
		ValueArg<long> LMIN("","LMIN","lmin multipole to define range of statistics performed  (default: 2)",false,2,"long"); cmd.add(LMIN);
		/* 	ValueArg<long> Yl("","Yl","l value for the mapping (see option mapYlm)  (default: 2)",false,2,"long"); cmd.add(Yl); */
		/* 	ValueArg<long> Ym("","Ym","m value for the mapping (see option mapYlm)  (default: 0)",false,0,"long"); cmd.add(Ym); */
		
		SwitchArg find("","find", "perform a search for preferred m's (default:false", false);	cmd.add( find );
		SwitchArg mapYlm("","mapYlm", "saves the values of amount of power accumulated at each probed orientation into a file for a requested l,m values set up by parameters Yl and Ym; results will be stored in a separate file;  option can be used only in ""find"" mode and only works in serial implementation (default:false", false);	cmd.add( mapYlm );
		SwitchArg calc_Dstats("","calc_Dstats", "calculate angles between dirs. of maximized r multipoles, and other stats on them (default:false", false);	cmd.add( calc_Dstats );
		SwitchArg calc_Astats("","calc_Astats", "calculate statistics on Almn data (default:false", false);	cmd.add( calc_Astats );
		SwitchArg calc_Covstats_inter_mode("","calc_Covstats_inter_mode", "calculate statistics on Almn data with full Cov: measures the inter-mode correlations in n_lm*n_l'm' (default:false", false);	cmd.add( calc_Covstats_inter_mode );
		SwitchArg calc_Covstats_test("","calc_Covstats_test", "calculate statistics on Almn data with full Cov: measures consistency with sims of given r_lm value in given direction (default:false", false);	cmd.add( calc_Covstats_test );
		
		ValueArg<string> output("o","out","outfile prefix. Various suffixes will be added by the program automatically",false,"AlmnMAX","string"); cmd.add(output);
		
		
		
		//SwitchArg pyth("","pyth", "whether or not plot the data with python script", false);	cmd.add( pyth );
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		
		_data_dir = data_dir.getValue(); 
		
		_sim_pref = sim_pref.getValue(); 
		_sim_suff = sim_suff.getValue(); 
		_sim_st = sim_st.getValue(); 
		_sim_en = sim_en.getValue(); 
		_Nsim = _sim_en-_sim_st+1; if (_sim_st == -1) { _Nsim = 1; _sim_en=-1; }
		_search_ns = search_ns.getValue();
		_ns = ns.getValue();
		_l = l.getValue();
		_LMAX = LMAX.getValue(); if (_LMAX == -1) _LMAX=_l;
		_LMIN = LMIN.getValue();
		/* 	_Yl = Yl.getValue(); */
		/* 	_Ym = Ym.getValue(); */
		_mapYlm = mapYlm.getValue();
		
		_find = find.getValue(); _find_done=false;
		_calc_Dstats = calc_Dstats.getValue(); _calc_Dstats_done=false;
		_calc_Astats = calc_Astats.getValue(); _calc_Astats_done=false;
		_calc_Covstats_inter_mode = calc_Covstats_inter_mode.getValue(); _calc_Covstats_inter_mode_done=false;
		_calc_Covstats_test = calc_Covstats_test.getValue(); _calc_Covstats_test_done=false;
		
		_output_file = output.getValue(); 
		
		
		
		
		/* 	_save_overplot_as = save_overplot_as.getValue(); */
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}


// *************************************************************************************************
// copies all alms table from the object to A[j].a structure
void copy_alms(mscsAlms &a0j,Almn_type *Aj) {
	long i,ist;
	long alms_num=a0j.almsNum();
	
	ist=alms_num-get_almsnumpos(alms_num);
	for (i=ist;i<alms_num;i++) { (*Aj).a[i-ist]=a0j.get(i).get_a_lm(); (*Aj).aa[i-ist]=a0j.get(i).get_a_lm(); }
}
// *************************************************************************************************

/* double calculate_r(mscsMap a0j, long i) { */
/*   long l,m; */
/*   a0j.num2alm(i,&l,&m); */
/*   return (a0j.alm[i].R*a0j.alm[i].R+a0j.alm[i].I*a0j.alm[i].I)/((2*l+1)*a0j.C_l[l][1]);  */
/* } */

// *************************************************************************************************
void calculate_rs(Almn_type *Aj,mscsAlms &a0j) {
	long i,l,m,ist;
	long alms_num=a0j.almsNum(); // number of alms in the object
	double r,rr;
	
	ist=alms_num-get_almsnumpos(alms_num); // number of negative m alms, and index of the a00
	for (i=ist;i<alms_num;i++) {     
		a0j.num2alm(i,&l,&m);
		r = ( a0j.get(i).R() * a0j.get(i).R()   +   a0j.get(i).I() * a0j.get(i).I() ) / ((2.0*(double)l+1.0)*a0j.calculate_single_C_l(l)); 
		if (m!=0) r*=2.0;
		(*Aj).r[i-ist] = (*Aj).rr[i-ist] = r;
	}
}

// *************************************************************************************************
void init_A(Almn_type **A,long Nsim) {
	long j;
	(*A) = new Almn_type[Nsim]; _Nsim=Nsim;
	for (j=0;j<_Nsim;j++) {
		(*A)[j].a=NULL; (*A)[j].r=NULL; (*A)[j].n=NULL; 
		(*A)[j].aa=NULL; (*A)[j].rr=NULL; (*A)[j].nn=NULL; 
	}
	
}

// *************************************************************************************************
void init_D(Dln_type **D,long Nsim) {
	long j;
	(*D) = new Dln_type[Nsim]; _Nsim=Nsim;
	for (j=0;j<_Nsim;j++) {
		(*D)[j].m=NULL; (*D)[j].r=NULL; (*D)[j].n=NULL; 
		(*D)[j].mm=NULL; (*D)[j].rr=NULL; (*D)[j].nn=NULL; 
	}
}
// *************************************************************************************************
void init_stat(statl_type **stat,long Nsim) {
	(*stat) = new statl_type[Nsim]; _Nsim=Nsim;
}

// *************************************************************************************************
// the input parameter alms_num referes to the number of alms in the Almn structure
// not the non-negative m mode number of alms !!
void initiate_Aj(Almn_type *Aj, long alms_num) { 
	long lmax,mmax,alms_num_tot; 
	// set lmax
	lmax=numpos2lmax(alms_num);
	
	/*   printf("alms_num_pos = %li\n",alms_num); */
	
	// initiate the structures
	if ((*Aj).a  == NULL) {  (*Aj).a  = new mscsAlm::mscs_a_lm[alms_num]; } 	
	if ((*Aj).aa == NULL) {  (*Aj).aa = new mscsAlm::mscs_a_lm[alms_num]; } 
	if ((*Aj).r  == NULL) {  (*Aj).r  = new double[alms_num]; }
	if ((*Aj).rr == NULL) {  (*Aj).rr = new double[alms_num]; }
	if ((*Aj).n  == NULL) {  (*Aj).n  = new cpeds_direction[alms_num]; }
	if ((*Aj).nn == NULL) {  (*Aj).nn = new cpeds_direction[alms_num]; }
	/*   if ((*Aj).i == NULL) {  (*Aj).i = new long[alms_num]; } */
	
	(*Aj).alms_num = alms_num;
	(*Aj).lmax = lmax;
}

// *************************************************************************************************
void initiate_Dj(Dln_type *Dj, long lmax) {
	(*Dj).lmax=lmax;
	lmax++;
	if ((*Dj).m   == NULL) { (*Dj).m  = new long[lmax]; }
	if ((*Dj).mm  == NULL) { (*Dj).mm = new long[lmax]; }
	if ((*Dj).r   == NULL) { (*Dj).r  = new double[lmax]; }
	if ((*Dj).rr  == NULL) { (*Dj).rr = new double[lmax]; }
	if ((*Dj).n   == NULL) { (*Dj).n  = new cpeds_direction[lmax]; }
	if ((*Dj).nn  == NULL) { (*Dj).nn = new cpeds_direction[lmax]; }
	/*   (*Dj).i = new long[lmax]; */
}
// *************************************************************************************************
void initiate_statj(statl_type *statj, long lmax) {
	(*statj).lmax=lmax;
	(*statj).v = new double[lmax];
	(*statj).V = new double[lmax];
}

// *************************************************************************************************
void kill_Aj(Almn_type *Aj) {
	if ((*Aj).a  !=NULL) { delete [] (*Aj).a; (*Aj).a=NULL; }
	if ((*Aj).aa !=NULL) { delete [] (*Aj).aa; (*Aj).a=NULL; }
	if ((*Aj).r  !=NULL) { delete [] (*Aj).r; (*Aj).r=NULL; }
	if ((*Aj).rr !=NULL) { delete [] (*Aj).rr; (*Aj).r=NULL; }
	if ((*Aj).n  !=NULL) { delete [] (*Aj).n; (*Aj).n=NULL; }
	if ((*Aj).nn !=NULL) { delete [] (*Aj).nn; (*Aj).n=NULL; }
	/*   if ((*Aj).i !=NULL) { delete [] (*Aj).i; (*Aj).i=NULL; } */
}
// *************************************************************************************************
void kill_A(Almn_type **A) {
	long j;
	printf("\n\nmpref> KILLING THE A STRUCTURE SPACE: \n\n");
	if ((*A)!=NULL) {
		for (j=0;j<_Nsim;j++) {
			kill_Aj(&(*A)[j]);
		}
		delete [] (*A); (*A)=NULL;
	}
}
// *************************************************************************************************
void kill_D(Dln_type **D) {
	long j;
	printf("\n\nmpref> KILLING THE D STRUCTURE SPACE: \n\n");
	if ((*D)!=NULL) {
		for (j=0;j<_Nsim;j++) {
			if ((*D)[j].m  != NULL) { delete [] (*D)[j].m;  (*D)[j].m = NULL; }
			if ((*D)[j].mm != NULL) { delete [] (*D)[j].mm; (*D)[j].mm= NULL; }
			if ((*D)[j].r  != NULL) { delete [] (*D)[j].r;  (*D)[j].r = NULL; }
			if ((*D)[j].rr != NULL) { delete [] (*D)[j].rr; (*D)[j].rr= NULL; }
			if ((*D)[j].n  != NULL) { delete [] (*D)[j].n;  (*D)[j].n = NULL; }
			if ((*D)[j].nn != NULL) { delete [] (*D)[j].nn; (*D)[j].nn= NULL; }
			
			/*     delete [] (*D)[j].i; */
		}
		delete [] (*D); (*D)=NULL;
	}
}

// *************************************************************************************************
void save_Almn(Almn_type *A, filenamestr fname) {
	FILE *f;
	long l,m,lmax=A[0].lmax;
	long j,k,s=(*A).alms_num-lmax-1;
	
	printf("\n\nmpref> SAVING Almnr TO FILE: %s\n\n",fname);
	
	/*   mscsMap tmp("tmp");  */
	/* num2alm(A[0].alms_num-1,&lmax,&m); tmp.set_alms_lmax(lmax); */
	
	f = fopen(fname,"w");
	
	fprintf(f,"%li %li %li\n",_Nsim,A[0].alms_num, A[0].lmax);
	for (j=0;j<_Nsim;j++) {
		l=m=0;
		for (k=0;k<A[0].alms_num;k++) {
			num2alm(k+s,&l,&m,lmax);
			fprintf(f,"%li %li %li %li %lE %lE %lE %.2lf %.2lf %lE %lE %lE %.2lf %.2lf\n",j,k,l,m,A[j].a[k].R, A[j].a[k].I,   A[j].r[k],   PI180inv*A[j].n[k].l, PI180inv*A[j].n[k].b,  A[j].aa[k].R, A[j].aa[k].I,   A[j].rr[k],   PI180inv*A[j].nn[k].l, PI180inv*A[j].nn[k].b);
		}
	}
	
	
	
	/*   f = fopen(fname,"wb"); */
	
	/*   s=sizeof(_Nsim);  fwrite(&_Nsim,s,1,f); */
	/*   for (j=0;j<_Nsim;j++) { */
	/*     s=sizeof(A[j].alms_num);  fwrite(&(A[j].alms_num),s,1,f); */
	/*     s=sizeof(A[j].a[0]);  fwrite(A[j].a,s,A[j].alms_num,f); */
	/*     s=sizeof(A[j].r[0]);  fwrite(A[j].r,s,A[j].alms_num,f); */
	/*     s=sizeof(A[j].n[0]);  fwrite(A[j].n,s,A[j].alms_num,f); */
	/*     s=sizeof(A[j].i[0]);  fwrite(A[j].i,s,A[j].alms_num,f); */
	/*   } */
	fclose(f);
	printf("\n\nmpref> SAVING Almnr TO FILE: %s DONE\n\n",fname);
}

// *************************************************************************************************
void read_Almn(Almn_type **A, filenamestr fname) {
	FILE *f;
	long j,k,s,alms_num,Nsim,lmax;
	mscsAlm a,aa;
	double r,rr;
	cpedsDirection n,nn;
	
	double are,aim, aare,aaim;
	
	printf("\n\nmpref> READING Almnr FROM FILE: %s\n\n",fname);
	f = fopen(fname,"r");
	
	fscanf(f,"%li %li %li",&Nsim,&alms_num, &lmax);
	printf("\n\nmpref> Sim_num: %li alms_num: %li lmax: %li\n\n",Nsim,alms_num,lmax);
	if ((*A)==NULL) {  init_A(A,Nsim); } else { kill_A(A); init_A(A,Nsim); }
	_Nsim=Nsim;
	
	for (j=0;j<_Nsim;j++) {
		initiate_Aj(&((*A)[j]),alms_num);
		for (k=0;k<alms_num;k++) {
			fscanf(f,"%*li %*li %*li %*li %lE %lE %lE %lf %lf %lE %lE %lE %lf %lf",&are, &aim, &r, &(n.lon()), &(n.lat()),   &aare, &aaim, &rr, &(nn.lon()), &(nn.lat()));
			a.R(are); a.I(aim); aa.R(aare); aa.I(aaim);
			n.lon()*=PI180; nn.lon()*=PI180; n.lat()*=PI180; nn.lat()*=PI180; // convert deg. to rad.
			(*A)[j].a[k] =a.get_a_lm();
			(*A)[j].aa[k]=aa.get_a_lm();
			(*A)[j].n[k] =n.get_direction();
			(*A)[j].nn[k]=nn.get_direction();
			(*A)[j].r[k] =r;
			(*A)[j].rr[k]=rr;           
			/*       printf("sim: %li, k:%li, aI:%lE aR:%lE nl:%lE nb:%lE r:%lE,   aaI:%lE aaR:%lE nnl:%lE nnb:%lE rr:%lE\n",j,k,a.I,a.R,n.l,n.b,r,   aa.I,aa.R,nn.l,nn.b,rr); */
		}
		(*A)[j].lmax=lmax;
		(*A)[j].alms_num=alms_num;
		printf("reading sim: %li of %li\r",j+1,Nsim);
	}
	
	/*   f = fopen(fname,"br"); */
	
	/*   s=sizeof(Nsim);  fread(&Nsim,s,1,f);   */
	/*   if ((*A)==NULL) {  init_A(A,Nsim); } else { kill_A(A); init_A(A,Nsim); } */
	
	/*   for (j=0;j<Nsim;j++) { */
	/*     s=sizeof(alms_num);  fread(&alms_num,s,1,f); */
	/*     initiate_Aj(&(*A)[j],alms_num); */
	/*     s=sizeof((*A)[j].a[0]);  fread((*A)[j].a,s,alms_num,f); */
	/*     s=sizeof((*A)[j].r[0]);  fread((*A)[j].r,s,alms_num,f); */
	/*     s=sizeof((*A)[j].n[0]);  fread((*A)[j].n,s,alms_num,f); */
	/*     s=sizeof((*A)[j].i[0]);  fread((*A)[j].i,s,alms_num,f); */
	/*   } */
	
	fclose(f);
	printf("\n\nmpref> READING Almnr FROM FILE: %s DONE\n\n",fname);
	
}


// *************************************************************************************************
void save_Dln(Dln_type *D, filenamestr fname) {
	FILE *f;
	long j,l;
	long LMIN=0;
	double phi,phi2;
	
	printf("\n\nmpref> SAVING Dlmnr FROM FILE: %s\n\n",fname);
	
	f = fopen(fname,"w");
	fprintf(f,"%li %li\n",_Nsim, D[0].lmax);
	
	for (j=0;j<_Nsim;j++) {
		for (l=LMIN;l<=D[j].lmax;l++) {
			phi= PI180inv*D[j].n[l].l;
			phi2=PI180inv*D[j].nn[l].l;
			if ((long)round(phi*100) ==36000) phi -=360.0; // get rid of the numerics that make 360 deg appear
			if ((long)round(phi2*100)==36000) phi2-=360.0; // get rid of the numerics that make 360 deg appear
			/*       fprintf(f,"%li %li %li %.2lf %.2lf %li %lE\n",j+_sim_st,l,D[j].m[l], phi, PI180inv*D[j].n[l].b, D[j].i[l], D[j].r[l]); */
			fprintf(f,"%li %li %li %.2lf %.2lf %lE %li %.2lf %.2lf %lE\n",j+_sim_st,l,   D[j].m[l], phi, PI180inv*D[j].n[l].b, D[j].r[l],     D[j].mm[l], phi, PI180inv*D[j].nn[l].b, D[j].rr[l]  );
		}
	}
	fclose(f);
	printf("\n\nmpref> SAVING Dlmnr FROM FILE: %s DONE\n\n",fname);
}
// *************************************************************************************************
void read_Dln(Dln_type **D, filenamestr fname) {
	FILE *f;
	long j,l;
	long LMIN=0,lmax,Nsim;
	double phi,phi2,b,b2;
	
	printf("\n\nmpref> LOADING Dlmnr FROM FILE: %s \n\n",fname);
	f = fopen(fname,"r");
	fscanf(f,"%li %li",&Nsim,&lmax);
	printf("\n\nmpref> -- number of sims: %li lmax: %li \n\n",Nsim,lmax);
	
	if ((*D)==NULL) {  init_D(D,Nsim); } else { kill_D(D); init_D(D,Nsim); }
	_Nsim=Nsim;
	
	for (j=0;j<_Nsim;j++) {
		initiate_Dj(&((*D)[j]),lmax);
		for (l=LMIN;l<=lmax;l++) {
			fscanf(f,"%*li %*li %li %lf %lf %lE %li %lf %lf %lE\n",  &((*D)[j].m[l]), &phi, &b, &((*D)[j].r[l]),     &((*D)[j].mm[l]), &phi2, &b2, &((*D)[j].rr[l])  );
			(*D)[j].n[l].l =phi *PI180;  (*D)[j].n[l].b =b *PI180; 
			(*D)[j].nn[l].l=phi2*PI180;  (*D)[j].nn[l].b=b2*PI180;
		}
		printf("reading sim: %li\n",j);
	}
	
	fclose(f);
	printf("\n\nmpref> LOADING Dlmnr FROM FILE: %s DONE\n\n",fname);
}

// *************************************************************************************************
void save_stat(statl_type *stat, filenamestr fname, long lmin, long lmax) {
	FILE *f;
	long j,l;
	
	printf("\n\nmpref> SAVING Dstats TO FILE: %s \n\n",fname);
	f = fopen(fname,"w");
	/*   printf("               dupa: %li\n",lmax); */
	for (j=0;j<_Nsim;j++) {
		for (l=lmin;l<=lmax;l++) {
			/*       printf("               dupa: %li\n",l); */
			fprintf(f,"%li %li %lE %lE\n",j+_sim_st,l,stat[j].v[l],stat[j].V[l]);
		}
	}
	fclose(f);
	printf("\n\nmpref> SAVING Dstats TO FILE: %s DONE\n\n",fname);
}

// *************************************************************************************************
void save_chisq_vals(double * chisq,long ist, long ien,filenamestr fname) {
	long i;
	FILE *f;
	
	printf("\n\nmpref> SAVING chisq values TO FILE: %s \n\n",fname);
	f = fopen(fname,"w");
	for (i=ist;i<=ien;i++) {    fprintf(f,"%li %lE\n",i,chisq[i]); }
	fclose(f);
	printf("\n\nmpref> SAVING chisq values TO FILE: %s DONE\n\n",fname);
}
//**********************************************************************
// conversion taken from the Mscs-map object 
long int alm2num(long l, long m, long lmax) {
	long int x = 0,s=1;
	long lm;
	lm=(lmax*(lmax+1+2*m)+m*(2-abs(m+1)))/2+l;
	return lm;
}
//************************************************************************
// conversion taken from the Mscs-map object 
void num2alm(long int num, long *l, long *m, long lmax) {
	long int x = 0;
	long int ll,lm;
	
	
	for (lm = -lmax;lm<=lmax;lm++) {
		for (ll=abs(lm); ll <= lmax; ll++) { 
			if (x == num) { *l=ll; *m=lm; lm = lmax; ll = lmax; }
			x++;
		}
	}
}
//**********************************************************************
// the alms_num input parameter is the number of alms in the objets files
// returns the number of alms with m>=0
long get_almsnumpos(long alms_num) {
	long lmax,mmax; num2alm_norm(alms_num-1,&lmax,&mmax);
	return (alms_num-lmax-1)/2+lmax+1;
}
//**********************************************************************
// the alms_num input parameter is the number of alms in the objets files
// returns the number of alms with m<0
long get_almsnumneg(long alms_num) {
	return alms_num-get_almsnumpos(alms_num);
}

//**********************************************************************
// the alms_num input parameter is the number of alms with m>=0
long get_almsnumtot(long alms_num) {
	long lmax=numpos2lmax(alms_num);
	return (alms_num-lmax-1)*2+lmax+1;
}

//**********************************************************************
// the numpos input parameter is the number of alms with m>=0
// and the number can only refer to full multipoles
long numpos2lmax(long numpos) {
	return (long)(0.5*(-3+sqrt(1.0+8.0*(double)numpos)));
}

//**********************************************************************
void num2alm_norm(long int num, long *l, long *m) {
	long int ll,lm,llmax=_l;
	long int x=0;
	
	for (ll=0;ll<=llmax;ll++) {
		for (lm=-ll;lm<=ll;lm++) { x++;
		if (x == num+1) {*l = ll; *m = lm; ll=llmax; lm=ll; return;}
		}
	}
	
}
//**********************************************************************
bool _something_else() {
	bool r=false;
	
	if (_find && !_find_done) r=true;
	if (_calc_Dstats && !_calc_Dstats_done) r=true;
	if (_calc_Astats && !_calc_Astats_done) r=true;
	if (_calc_Covstats_inter_mode && !_calc_Covstats_inter_mode_done) r=true;
	if (_calc_Covstats_test && !_calc_Covstats_test_done) r=true;
	
	return r;
}
//**********************************************************************
