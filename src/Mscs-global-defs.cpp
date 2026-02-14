#include <string.h>
#include "Mscs-common.h"
#include "Mscs-global-defs.h"
#include "qtbase.h"
//#include <QtCore/QString>
//#include <QtCore/QStringList>
//#include "cpeds-math.h"

// package version information
string Mscs_version;

/*******************************************************************************/
/*  declaration of shared memory variables that hold running tasks information */
/*******************************************************************************/

Mscs_run_control_t* MSCS_GLOBAL__RUN_SIMULATION_INFO;
char MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY[255];



/******************************************************************/
/* //! declaration of variables that hold naming Mscs conventions */
/******************************************************************/
// heaplix pix system transfer function file prefix
string MSCS_GLOBAL__HEALPIX_PIXTF_PREF;
// heaplix pix system transfer function file suffix
string MSCS_GLOBAL__HEALPIX_PIXTF_SUFF;

// nested galactic coordinates file prefix
string MSCS_GLOBAL__NESTED_COORDINATES_PREF;
// nested galactic coordinates file suffix
string MSCS_GLOBAL__NESTED_COORDINATES_SUFF;
// ring galactic coordinates file prefix
string MSCS_GLOBAL__RING_COORDINATES_PREF;
// ring galactic coordinates file suffix
string MSCS_GLOBAL__RING_COORDINATES_SUFF;

// nested 2 ring conversion array file prefix
string MSCS_GLOBAL__N2R_CONV_TAB_PREF;
// nested 2 ring conversion array file suffix
string MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN;
string MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT;
// ring 2 nested conversion array file prefix
string MSCS_GLOBAL__R2N_CONV_TAB_PREF;
// ring 2 nested conversion array file prefix
string MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN;
string MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT;


/*************************************************************/
/* // definitions of directories where usefull stuff is kept */
/*************************************************************/

string HOME_DIR;// = "/home/blew/cosmo-data/wmap/";
string MSCS_PROGRAM_DIR;// = "/home/blew/cosmo-data/wmap/";
string MSCS_WMAP_DATA_DIR;// = "/home/blew/cosmo-data/wmap/";
string MSCS_WMAP_LOCAL_DIR;// = "/home/blew/programy/Mscs/WMAPstuff/";
string MSCS_DATA_DIR;// =  "/home/blew/programy/Mscs/data/";


/* string DAnames_nums[10];// = {"q1","q2","v1","v2","w1","w2","w3","w4"}; */

/* string DAnames[10];// = {"k","k","q","q","v","v","w","w","w","w"}; */
/* long DAnums[10];// = {1,2,1,2,1,2,1,2,3,4}; */

// definitions of extern files to be used in programs

string GLOBAL_path;
long GLOBAL_nside;
long GLOBAL_pix_num;
long GLOBAL_lmax;
string GLOBAL_file_prefix;
string GLOBAL_DA1;
int GLOBAL_DA1i;
string GLOBAL_DA2;
int GLOBAL_DA2i;
string GLOBAL_mask_file;
string GLOBAL_smoothing_method;
int GLOBAL_fourier_method_num,GLOBAL_fourier_forward_method_num, GLOBAL_fourier_backward_method_num;

string GLOBAL_Plmsfile_forward;
string GLOBAL_Plmsfile_inverse;


void Mscs_initiate_global_variables() {
	
#pragma omp critical ( MSCS_GLOBAL_DIRS )
	{
		Mscs_version =  (string)_Mscs_version_;
		Mscs_initiate_directories();
		/* Mscs_initiate_WMAP_DAs(); */
		
		GLOBAL_fourier_method_num=8;
		GLOBAL_fourier_forward_method_num=8;
		GLOBAL_fourier_backward_method_num=8;
		
		strcpy(MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY,"/MSCS_RUN_SIMULATION_SHM_KEY");
	}
}

void Mscs_initiate_global_variables(long _nside, long _lmax) {
	throw "this part of the code is depreciated and should not be used";
	Mscs_version =  string(_Mscs_version_);
	Mscs_initiate_directories();
	/* Mscs_initiate_WMAP_DAs(); */
	
	GLOBAL_fourier_method_num=7;  GLOBAL_fourier_forward_method_num=7;  GLOBAL_fourier_backward_method_num=7;
	/* Mscs_initiate_fourier_variables(_nside,_lmax); */
	
	//initiation of some other variables (eg. files for fourier transforms )
	/* string __tmp__; */
	/* make_detailed_file_name(&__tmp__,"",_nside,_lmax,"","",0,"",0,"","",GLOBAL_fourier_method_num); */
	
}


/* void Mscs_initiate_fourier_variables(long _nside, long _lmax) { */

/*   sprintf(GLOBAL_Plmsfile_forward,"%sPlms-%li-%li",data_dir,_nside,_lmax); */
/*   if ( GLOBAL_fourier_method_num == 7 ) { sprintf(GLOBAL_Plmsfile_inverse,"%sPlms-%li-%li",data_dir,_nside,_lmax); } */
/*   if (( GLOBAL_fourier_method_num == 8 ) || ( GLOBAL_fourier_method_num == 9 )) { sprintf(GLOBAL_Plmsfile_inverse,"%sPmthls-%li-%li",data_dir,_nside,_lmax); } */

/* } */


void Mscs_print_global_variables() {
	printf("---------------------------------------------------\n");
	printf("MSCS PACKAGE DIRECTORIES\n");
	printf("---------------------------------------------------\n");
	printf("Mscs_version:   %s\n",Mscs_version.c_str());
	printf("HOME_DIR:   %s\n",HOME_DIR.c_str());
	printf("MSCS_PROGRAM_DIR:   %s\n",MSCS_PROGRAM_DIR.c_str());
	printf("MSCS_WMAP_DATA_DIR:   %s\n",MSCS_WMAP_DATA_DIR.c_str());
	printf("MSCS_WMAP_LOCAL_DIR:   %s\n",MSCS_WMAP_LOCAL_DIR.c_str());
	printf("MSCS_DATA_DIR:   %s\n",MSCS_DATA_DIR.c_str());
	printf("---------------------------------------------------\n");
	printf("MSCS PACKAGE CONVERSION ARRAYS NAMING CONVENTIONS\n");
	printf("---------------------------------------------------\n");
	printf("MSCS_GLOBAL__HEALPIX_PIXTF_PREF:   %s\n",MSCS_GLOBAL__HEALPIX_PIXTF_PREF.c_str());
	printf("MSCS_GLOBAL__N2R_CONV_TAB_PREF:   %s\n",MSCS_GLOBAL__N2R_CONV_TAB_PREF.c_str());
	printf("MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT:   %s\n",MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT.c_str());
	printf("MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN:   %s\n",MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN.c_str());
	printf("MSCS_GLOBAL__R2N_CONV_TAB_PREF:   %s\n",MSCS_GLOBAL__R2N_CONV_TAB_PREF.c_str());
	printf("MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT:   %s\n",MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT.c_str());
	printf("MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN:   %s\n",MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN.c_str());
	printf("---------------------------------------------------\n");
	printf("MSCS PACKAGE CONVERSION ARRAYS NAMING CONVENTIONS\n");
	printf("---------------------------------------------------\n");
	printf("MSCS_GLOBAL__NESTED_COORDINATES_PREF:   %s\n",MSCS_GLOBAL__NESTED_COORDINATES_PREF.c_str());
	printf("MSCS_GLOBAL__NESTED_COORDINATES_SUFF:   %s\n",MSCS_GLOBAL__NESTED_COORDINATES_SUFF.c_str());
	printf("MSCS_GLOBAL__RING_COORDINATES_PREF:   %s\n",MSCS_GLOBAL__RING_COORDINATES_PREF.c_str());
	printf("MSCS_GLOBAL__RING_COORDINATES_SUFF:   %s\n",MSCS_GLOBAL__RING_COORDINATES_SUFF.c_str());
	printf("\n\n");
	/* printf(":   %s\n",.c_str()); */
	
}


void  Mscs_initiate_directories() {
	
	//	MSCS_GLOBAL__ENV_PROC_NUM="MSCS_GLOBAL__ENV_PROC_NUM";
	//	MSCS_GLOBAL__ENV_SIMULATION_RUN_FILE_PREF="MSCS_GLOBAL__ENV_SIMULATION_RUN_FILE_PREF";
	//	MSCS_GLOBAL__ENV_SIMULATION_RUN_FILE_SUFF="MSCS_GLOBAL__ENV_SIMULATION_RUN_FILE_SUFF";
	
	
	MSCS_GLOBAL__HEALPIX_PIXTF_PREF="wl-";
	MSCS_GLOBAL__N2R_CONV_TAB_PREF="n2r-convtab-";
	MSCS_GLOBAL__R2N_CONV_TAB_PREF="r2n-convtab-";
	MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT="";
	MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT="";
	MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN=".bin";
	MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN=".bin";
	
	MSCS_GLOBAL__NESTED_COORDINATES_PREF="";
	MSCS_GLOBAL__NESTED_COORDINATES_SUFF="-Cng-bin";
	MSCS_GLOBAL__RING_COORDINATES_PREF="";
	MSCS_GLOBAL__RING_COORDINATES_SUFF="-Crg-bin";
	
	char* cstr;
	cstr=getenv("HOME");
	if (cstr!=NULL)  HOME_DIR=cstr; 
	else HOME_DIR=".";   
	HOME_DIR+="/";
//	cstr=getenv("MSCS_PROGRAM_DIR");
	MSCS_PROGRAM_DIR=HOME_DIR+".cpems/"; 
	MSCS_WMAP_DATA_DIR=HOME_DIR+"cosmo-data/wmap/";
	MSCS_WMAP_LOCAL_DIR=MSCS_PROGRAM_DIR+"WMAPstuff/";
	
	MSCS_DATA_DIR=MSCS_PROGRAM_DIR+"data/";
	
}

/* void Mscs_initiate_WMAP_DAs() { */

/*   strcpy(DAnames_nums[0],"k1"); */
/*   strcpy(DAnames_nums[1],"k2"); */
/*   strcpy(DAnames_nums[2],"q1"); */
/*   strcpy(DAnames_nums[3],"q2"); */
/*   strcpy(DAnames_nums[4],"v1"); */
/*   strcpy(DAnames_nums[5],"v2"); */
/*   strcpy(DAnames_nums[6],"w1"); */
/*   strcpy(DAnames_nums[7],"w2"); */
/*   strcpy(DAnames_nums[8],"w3"); */
/*   strcpy(DAnames_nums[9],"w4"); */

/*   strcpy(DAnames[0],"k"); */
/*   strcpy(DAnames[1],"k"); */
/*   strcpy(DAnames[2],"q"); */
/*   strcpy(DAnames[3],"q"); */
/*   strcpy(DAnames[4],"v"); */
/*   strcpy(DAnames[5],"v"); */
/*   strcpy(DAnames[6],"w"); */
/*   strcpy(DAnames[7],"w"); */
/*   strcpy(DAnames[8],"w"); */
/*   strcpy(DAnames[9],"w"); */

/*   DAnums[0]=1;  */
/*   DAnums[1]=2; */
/*   DAnums[2]=1;  */
/*   DAnums[3]=2; */
/*   DAnums[4]=1; */
/*   DAnums[5]=2; */
/*   DAnums[6]=1; */
/*   DAnums[7]=2; */
/*   DAnums[8]=3; */
/*   DAnums[9]=4; */


/* } */









// definitions of common routines for programs of Mscs package


/* string*  make_detailed_file_name(string *name, */
/* 				    strarg pathP, */
/* 				    long nsideP, */
/* 				    long lmaxP, */
/* 				    strarg file_prefixP, */
/* 				    strarg DA1P, */
/*     				    int DA1iP, */
/* 				    strarg DA2P, */
/* 				    int DA2iP, */
/* 				    strarg mask_fileP, */
/* 				    strarg smoothing_methodP, */
/* 				    int fourier_method_numP) { */

/*   strcpy(GLOBAL_path,pathP); */
/*   GLOBAL_nside = nsideP; */
/*   GLOBAL_pix_num = 12*GLOBAL_nside*GLOBAL_nside; */
/*   GLOBAL_lmax = lmaxP; */
/*   strcpy(GLOBAL_file_prefix,file_prefixP); */
/*   strcpy(GLOBAL_DA1,DA1P); GLOBAL_DA1i = DA1iP; */
/*   strcpy(GLOBAL_DA2,DA2P); GLOBAL_DA2i = DA2iP; */
/*   strcpy(GLOBAL_mask_file,mask_fileP); */
/*   strcpy(GLOBAL_smoothing_method,smoothing_methodP); */
/*   GLOBAL_fourier_method_num = fourier_method_numP; */

/* /\*   sprintf(GLOBAL_Plmsfile_forward,"%sPlms-%li-%li",data_dir,GLOBAL_nside,2*GLOBAL_nside); *\/ */
/*   sprintf(GLOBAL_Plmsfile_forward,"%sPlms-%li-%li",data_dir,GLOBAL_nside,GLOBAL_lmax); */
/*   if ( GLOBAL_fourier_method_num == 7 ) { sprintf(GLOBAL_Plmsfile_inverse,"%sPlms-%li-%li",data_dir,GLOBAL_nside,GLOBAL_lmax); } */
/*   if (( GLOBAL_fourier_method_num == 8 ) || ( GLOBAL_fourier_method_num == 9 )) { sprintf(GLOBAL_Plmsfile_inverse,"%sPmthls-%li-%li",data_dir,GLOBAL_nside,GLOBAL_lmax); } */
/* /\*   strcpy(GLOBAL_Plmsfile_inverse,""); *\/ */
/* /\*   strcpy(GLOBAL_Plmsfile_forward,""); *\/ */

/* /\*   fwhmIN = strtod(ARGV[3],NULL); if (fwhmIN < 0) { strcpy(smoothIN,"nosmooth");} else { strcpy(smoothIN,"smooth"); } //fwhmIN *= PI/180.0; *\/ */
/* /\*   fwhmOUT = strtod(ARGV[4],NULL); if (fwhmIN < 0) { strcpy(smoothOUT,"nosmooth");} else { strcpy(smoothOUT,"smooth"); } //fwhmOUT *= PI/180.0; *\/ */


/*   sprintf(*name,"%s%li-%li-%s-%s%i-%s%i-%s-sm_%s-Fmet_%i",pathP,nsideP,lmaxP,file_prefixP,DA1P,DA1iP,DA2P,DA2iP,mask_fileP,smoothing_methodP,fourier_method_numP); */
/*   return name; */

/* } */

/* void  make_detailed_file_name_comment(string *name, strarg comment) { */
/*   string tmp; */
/*   strcpy(tmp,*name); */
/*   sprintf(*name,"%s-%s",tmp,comment); */
/* } */

/* void make_covariance_matrix_name(string* name, strarg dir, long start_num, long end_num, long sim_num, long th_num, strarg comment) { */
/*   sprintf(*name,"%scovariance_matrix-%s-S%li-E%li-sim%li-th%li",dir,comment,start_num,end_num,sim_num,th_num); */

/* } */

// *************************************************************************************************
/* double get_angle(direction n1, direction n2) { */
/*   return cpeds_ang_n1n2(PIsnd-n1.b,n1.l,PIsnd-n2.b,n2.l); */
/* } */


//--------------------

void Mscs_matrix_print(matrix<double> *M) {
	unsigned long i,j;
	
	for (i=0;i<(*M).RowNo();i++) {
		for (j=0;j<(*M).ColNo();j++)  printf("%lE ",(*M)(i,j));
		printf("\n");
	}
}

/* void Mscs_matrix_save(matrix <double> *M, strarg where) { */
/*   unsigned long i,j; */
/*   FILE * f=NULL; */

/*   f = fopen(where,"w"); */
/*   if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",where);  exit(1); } */

/*   fprintf(f,"%li %li\n",(long)(*M).RowNo(),(long)(*M).ColNo()); */
/*   printf("  -- saving matrix of size: %li rows %li cols to file %s\n",(long)(*M).RowNo(),(long)(*M).ColNo(),where); */
/*   for (i=0;i<(*M).RowNo();i++) { */
/*     for (j=0;j<(*M).ColNo();j++)  fprintf(f,"%.10lE ",(*M)(i,j)); */
/*     fprintf(f,"\n"); */
/*   } */
/*   fclose(f); */
/* } */

/* void Mscs_matrix_save(matrix <double> *M, strarg where, string how) { */
/*   unsigned long i,j; */
/*   FILE * f=NULL; */
/*   double tmpd; */

/*   f = fopen(where,"w"); */
/*   if (f == NULL) { printf("ERROR: couldn't save matrix to: %s\n",where);  exit(1); } */
/*   if (how.find("header",0) != string::npos) { */
/*     fprintf(f,"%li %li\n",(long)(*M).RowNo(),(long)(*M).ColNo()); } */
/*   printf("  -- saving matrix of size: %li rows %li cols to file %s\n",(long)(*M).RowNo(),(long)(*M).ColNo(),where); */
/*   fclose(f); */

/*   if (how.find("binary",0) == string::npos) { */
/*     f = fopen(where,"a"); */
/*     for (i=0;i<(*M).RowNo();i++) { */
/*       for (j=0;j<(*M).ColNo();j++)  { tmpd = (*M)(i,j); fprintf(f,"%.20lE ",tmpd); }   fprintf(f,"\n"); */
/*     } */
/*   } */
/*   else { */
/*     f = fopen(where,"ab"); */
/*     for (i=0;i<(*M).RowNo();i++) { */
/*       for (j=0;j<(*M).ColNo();j++)  { tmpd = (double)((*M)(i,j));  fwrite(&tmpd,sizeof(tmpd),1,f); } */
/*     } */
/*   } */
/*   fclose(f); */
/* } */

/* matrix <double> * Mscs_matrix_load(strarg where) { */
/*   unsigned long i,j; */
/*   FILE * f=NULL; */
/*   long row,col; */
/*   double tmp; */
/*   matrix <double> *M; */

/*   f = fopen(where,"r");  if (f == NULL) { printf("   -- ERROR: no file: %s\n",where); return NULL; } */
/*   fscanf(f,"%li %li",&row,&col);  */
/*   printf("  -- reading matrix of size: %li rows, %li cols from file %s\n",row,col,where); //exit(0); */
/*   M = new matrix <double>; (*M).SetSize(row,col); */
/*   for (i=0;i<(*M).RowNo();i++) { */
/*     for (j=0;j<(*M).ColNo();j++)  { fscanf(f,"%lE ",&tmp); (*M)(i,j) = tmp; }//printf("i=%li, j =%li\n",i,j); } */
/*   } */
/*   fclose(f); */
/*   return M; */
/* } */

/* matrix <double> * Mscs_matrix_load(strarg where, string how) { */
/*   unsigned long i,j; */
/*   FILE * f=NULL; */
/*   long row,col; */
/*   double tmp; */
/*   char c; */
/*   matrix <double> *M; */
/*   cpeds_queue<long> *finfo; */

/*   if (how.find("binary",0) != string::npos) {  f = fopen(where,"r"); } else {f = fopen(where,"rb");}  */
/*   if (f == NULL) { printf("   -- ERROR: no file: %s\n",where); return NULL; } */
/*   if (how.find("header",0) != string::npos) {    fscanf(f,"%li %li",&row,&col); }  */
/*   else { */
/*     finfo=cpeds_get_txt_file_cols_rows(where); */
/*     col=(*finfo)(0); */
/*     row=(*finfo).get_size();     */
/*     printf("|Mscs_matrix_load> * assuming the matrix of size: rows %li cols %li\n",row,col); */
/*   } */
/*   printf("  -- reading matrix of size: %li rows, %li cols from file %s\n",row,col,where); //exit(0); */
/*   M = new matrix <double>; (*M).SetSize(row,col); */

/*   if (how.find("binary",0) == string::npos) { */
/*     for (i=0;i<(*M).RowNo();i++) { */
/*       for (j=0;j<(*M).ColNo();j++)  { fscanf(f,"%lE ",&tmp); (*M)(i,j) = tmp; }//printf("i=%li, j =%li\n",i,j); } */
/*     } */
/*   }  */
/*   else { */
/*     fread(&c,sizeof(c),1,f); // read in the \n character-- 1 byte to get the right values !! this might be system dependent or library dependent !! */
/*     for (i=0;i<(*M).RowNo();i++) { */
/*       for (j=0;j<(*M).ColNo();j++)  { fread(&tmp,sizeof(tmp),1,f); (*M)(i,j) = tmp; }  */
/*     } */
/*   } */

/*   fclose(f); */
/*   return M; */
/* } */

/* void initiate_f_stat(f_stat *fstat, long th_num,long mink_num,long nkbin) { */
/*   long n; */
/*   for (n=0;n<mink_num;n++) { // loop for different covarians matrixes */
/*     fstat[n].k_num = th_num;  fstat[n].k = new double[fstat[n].k_num]; */
/*     fstat[n].fkm = new double[fstat[n].k_num]; fstat[n].fkv = new double[fstat[n].k_num]; fstat[n].fks = new double[fstat[n].k_num]; fstat[n].fkk = new double[fstat[n].k_num]; */
/*     fstat[n].bin_num = nkbin; fstat[n].bin = new double[fstat[n].bin_num]; fstat[n].Nbin = new double[fstat[n].bin_num*th_num]; */
/*   } */
/* } */

/* void delete_f_stat(f_stat *fstat,long mink_num) { */
/*   long n; */
/*   for (n=0;n<mink_num;n++) { // loop for different covarians matrixes */
/*     delete fstat[n].k; */
/*     delete fstat[n].fkm; */
/*     delete fstat[n].fkv; */
/*     delete fstat[n].fks; */
/*     delete fstat[n].fkk; */
/*     delete fstat[n].bin; */
/*     delete fstat[n].Nbin; */
/*   } */
/* } */


// this routine bins the input matrix  by top-hat dx per dy sized bins
// binned matrix pointer is returned.
// work_mode defines the amount of printouts during execution: 0 - no proints 1 - print
matrix<double>* Mscs_matrix_bin(matrix <double> *M, long dx, long dy, long work_mode) {
	matrix <double>* Mb = new matrix<double>;
	long Mx,My,Mbx, Mby;
	long i,j,k,l,p,q,w,lmax,kmax;
	
	
	
	// find out the size of the new matrix
	Mx = M->ColNo(); My = M->RowNo();
	if (Mx % dx == 0) Mbx = Mx/dx; else Mbx = Mx/dx + 1;
	if (My % dy == 0) Mby = My/dy; else Mby = My/dy + 1;
	Mb->SetSize(Mby,Mbx);
	if (work_mode==1) printf(" new matrix X,Y size is: %li %li\n",Mbx,Mby);
	
	q=0; j=0;
	while (j<My) {
		//for (j=0;j<My;j++) {
		p=0; i=0;
		while (i<Mx) {
			//for (i=0;i<Mx;i++) {
			/*       if (j!=0 &&Mbx!=24) { printf("p: %li q: %li\n",p,q); exit(0); } */
			(*Mb)(q,p)=0; w=0;
			
			if (j+dy > My) lmax = My-j; else lmax = dy;
			for (l=0;l<lmax;l++) {
				if (i+dx > Mx) kmax = Mx-i; else kmax = dx;
				for (k=0;k<kmax;k++) {
					w++;
					(*Mb)(q,p)+=(*M)(j+l,i+k);
					//if (j>=800) printf("p %li q %li, i %li j %li, k %li l %li\n",p,q,i,j,k,l);
					
				}
			}
			i+=kmax;
			//printf("shit1\n");
			(*Mb)(q,p)/=(double)w;
			p++;
			
		}
		/*     printf("---- j:%li\n",j); */
		j+=lmax;
		q++;
	}
	/*  if (Mbx!=24) exit(0); */
	//      printf("shit2\n");
	return Mb;
}

/* bool Mscs_get_mask_filename(string in_file, string in_dir,) { */
/* } */


bool Mscs_contains_str(string s, const char* c, long nc) {
	if (s.find(c,0,nc) == string::npos) return false;
	else return true;
}

string Mscs_basename_noext(string file,string ext) {
	string base;
	size_t pos;
	
	// get the basename
	pos = file.rfind("/",string::npos);
	if (pos != string::npos) { file.erase(0,pos+1); }
	
	// remove the requested extension
	if (ext != "") {
		if (file.find(ext,file.size()-ext.size()) != string::npos) { file.erase(file.size()-ext.size(),ext.size()); }
	}
	// remove standard extension
	ext="-Tn-bin";
	if (file.find(ext,file.size()-ext.size()) != string::npos) { file.erase(file.size()-ext.size(),ext.size()); }
	return file;
}

string Mscs_basename(string name) {
	//	 QFileInfo fi(name.c_str());
	QStringList qsl;
	QString s=name.c_str();
	qsl=s.split(".");
	QString base=""; // = fi.baseName();  // base = "archive"	
	if (qsl.size()==1) return name;
	for (unsigned long i = 0; i < qsl.size()-1; i++) {
		base+=qsl[i];
		if (i<qsl.size()-2)
			base+=".";
	}
	return base.toStdString();
	
}

string Mscs_getExtension(string name) {
	QStringList qsl;
	QString s=name.c_str();
	qsl=s.split(".");
	QString ext=""; 
	if (qsl.size()>0) return qsl.last().toStdString();
	return ext.toStdString();
}

vector<string> Mscs_stringToList(string str, string sep) {
	QStringList qsl;
	QString s=str.c_str();
	qsl=s.split(sep.c_str());
	vector<string> v;	
	for (unsigned long i = 0; i < qsl.size(); i++) {
		v.push_back(qsl[i].toStdString());
	}
	return v;
	
}

vector<double> Mscs_stringToDoubleList(string str, string sep) {
	vector<string> sl=Mscs_stringToList(str,sep);
	vector<double> dl;
	QString sv;
	
	for (unsigned long i = 0; i < sl.size(); i++) {
		sv=sl[i].c_str();
		dl.push_back(sv.toDouble());
	}
	return dl;
}
vector<long> Mscs_stringToLongList(string str, string sep) {
	vector<string> sl=Mscs_stringToList(str,sep);
	vector<long> dl;
	QString sv;
	
	for (unsigned long i = 0; i < sl.size(); i++) {
		sv=sl[i].c_str();
		dl.push_back(sv.toLong());
	}
	return dl;	
}

//QList<QString> Mscs_stringToListQ(string str, string sep) {
//	QStringList qsl;
//	QString s=str.c_str();
//	qsl=s.split(sep.c_str());
//	QList<QString> v;	
//	for (unsigned long i = 0; i < qsl.size(); i++) {
//		v.append(qsl[i]);
//	}
//	return v;
//
//}
