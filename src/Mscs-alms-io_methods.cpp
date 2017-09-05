//************************************************************************
void mscsAlms::read_binalm_parameters(strarg  alm_file) {
  double tmp;
  long int n = 0;
  long mmax;
  filenamestr almfile;
  strcpy(almfile,alm_file);
//  strcat(almfile,"-F-bin");

  f = fopen(almfile,"rb");
  while (!feof(f)) { fread(&tmp,sizeof(tmp),1,f); n++;}
  fclose(f);
  //printf("n = %li ",n);
  alms_num = (n-1)/2;
  num2alm_norm(alms_num-1,&lmax,&mmax); // set the lmax and mmax in the file; // -2 because -1 gives the real number of alms and -1 to get the numbering from 0
  //num2alm(n-1,&lmax,&mmax);

  printf("|%s>  -- alm file parameters: lmax = %li, mmax = %li, alms  = %li\n",object_name.c_str(),lmax,mmax,alms_num);
  printf("|%s>  -- alm file size in memory (twice of physical): bytes = %li\n",object_name.c_str(),4*alms_num*(long int)sizeof(tmp)); // this I don't understand (seems to me like this should be n, not n-1)
  printf("|%s>  -- alm file name: %s\n",object_name.c_str(),almfile);

}
//************************************************************************
void mscsAlms::read_txtalm_parameters(strarg  alm_file) {
  double tmp;
  long int n = 0;
  long int mmax;
  filenamestr almfile;
  strcpy(almfile,alm_file);
//  strcat(almfile,"-F-txt");

  f = fopen(almfile,"rb");
  while (!feof(f)) { fscanf(f,"%lE",&tmp); n++;}
  fclose(f);

  alms_num = (n-1)/2;
  num2alm_norm(alms_num-1,&lmax,&mmax); // set the lmax and mmax in the file; // -2 because -1 gives the real number of alms and -1 to get the numbering from 0
  //num2alm(n-1,&lmax,&mmax);

  printf("|%s>  -- alm file parameters: lmax = %li, mmax = %li, alms  = %li\n",object_name.c_str(),lmax,mmax,alms_num);
  printf("|%s>  -- alm file size in memory (twice of physical): bytes = %li\n",object_name.c_str(),4*alms_num*(long int)sizeof(tmp)); // this I don't understand (seems to me like this should be n, not n-1)
  printf("|%s>  -- alm file name: %s\n",object_name.c_str(),almfile);

}






//************************************************************************
// how : 1 - copy data as is.
// how : 3 - kills map and initiates  nside as in the source object (CAREFULL!!)
// (this is only applicable to map stuff (not alms flatmap etc.)
void mscsAlms::import_alms_data(mscsAlms &from_here, string what, int how) {
  long i;

  printf("|%s> * Importing map data from |%s> about %s\n",object_name.c_str(),from_here.object_name.c_str(),what.c_str());
  if (how == 3) { killAlms(); set_alms_lmax(from_here.get_lmax()); }


  if (what == "alms") {
    makekill_space_manager("make","F",1);
    for (i=0;i<alms_num;i++) { set(i,from_here.get(i)); }
    alms_loaded = 1;
  }
  else { printf("|%s> !! WARNING: the import has no effect: I don't know what %s is !!\n",object_name.c_str(),what.c_str()); }}}}}

}

//************************************************************************
void mscsAlms::savebinAlms(strarg alms_file, int how) { // saves alms to a file
  loadsave_manager("save","bin","F",how,alms_file,-1);
}
void mscsAlms::savebinAlms(strarg alms_file, int how, long l) { // saves alms to a file
  loadsave_manager("save","bin","F",how,alms_file,l);
}
//************************************************************************
void mscsAlms::loadbinAlms(strarg  alms_file, int how) { // saves alms to a file
  loadsave_manager("load","bin","F",how,alms_file,-1);
}
void mscsAlms::loadbinAlms(strarg  alms_file, int how, long l) { // saves alms to a file
  loadsave_manager("load","bin","F",how,alms_file,l);
}
//************************************************************************
void mscsAlms::savetxtAlms(strarg  alms_file, int how) { // saves alms to a file
  loadsave_manager("save","txt","F",how,alms_file,-1);
}
void mscsAlms::savetxtAlms(strarg  alms_file, int how, long l) { // saves alms to a file
  loadsave_manager("save","txt","F",how,alms_file,l);
}
//************************************************************************
void mscsAlms::loadtxtAlms(strarg  alms_file, int how) { // saves alms to a file
  loadsave_manager("load","txt","F",how,alms_file,-1);
}
void mscsAlms::loadtxtAlms(strarg  alms_file, int how, long l) { // saves alms to a file
  loadsave_manager("load","txt","F",how,alms_file,l);
}
//************************************************************************
void mscsAlms::printtxtAlms(strarg what) {// prints fourier alms to the screen
  int l,m;
  a_lm tmpalm;

  printf("|%s> * writting ",object_name.c_str());
  //if (strcmp(*what,"g") == 0) { printf("galactic (map native) "); }
  //if (strcmp(*what,"f") == 0) { printf("flat (eg. from proj program) "); }
  printf("alms to screen:\n");

  for (l=0;l<=lmax;l++) {
    for (m=-l;m<=l;m++) {
      tmpalm = alm[alm2num(l,m)];
      printf("l = %i, m = %i, R = %lE, I = %lE, D = %lE, PHI = %lE\n",l,m,tmpalm.R,tmpalm.I,tmpalm.D,tmpalm.P);
    }
  }

}

//************************************************************************
void mscsAlms::killAlms() { // kills the alms from the memory
/*   if (alms_loaded == 1) delete [] alm; */
/*   alms_loaded = 0; */
  makekill_space_manager("kill","F",1);
}


//************************************************************************
void  mscsMap::makekill_space_manager(strarg whattodo, strarg for_what, int how) {
  if (strcmp(for_what,"F") == 0) { // allocating memory for SHT coefficients
    //alms_num = alm2num(lmax,lmax)+1;
    if ((strcmp(whattodo,"make") == 0) || (strcmp(whattodo,"load") == 0)) { if (alm == NULL) { alm = new a_lm[alms_num]; printf("|%s> * allocating memory for alms: alm_num = %li\n",object_name.c_str(),alms_num);  alms_loaded = 0; } }
    if (strcmp(whattodo,"kill") == 0) {   if (alms_loaded == 1) {printf("|%s> * killing memory for alms\n",object_name.c_str()); delete [] alm; alm = NULL; alms_loaded = 0; } }
  }

}

//************************************************************************
void mscsMap::loadsave_manager(strarg whattodo, strarg fileformat, strarg what, int how, strarg where, long whatmultipole) {
//void mscsMap::loadsave_manager(filenamestr whattodo, filenamestr fileformat, filenamestr what, int how, filenamestr where, int whatmultipole) {
  FILE *f = NULL;
  long l,m,tmpl;
  long int i;
  char tmpch[500];
  float lcoord,bcoord;
  double  ld,tmpd,tmpde,theta,dtheta,RE,IM,minkfunc,th1,phi1,th2,phi2,Nobs,*X;
  filenamestr whereto;
  strcpy(whereto,where);
  string tmps;
  //where = &whereto;
 
  //------------- fits variables
  fitsfile *fptr;         
  char card[FLEN_CARD]; 
  int   nkeys;//, ii;
  size_t sofd=sizeof(double);
  bool coords_changed=false;
  //int status = 0
  //int anynul
  //float dupa[100]; 
  //float *dupanull;





  //*************************************************************
  //*************************************************************
  //************************************************************
  // F file format is 2 double numbers (8byte each) for each alm number
  // the alms ordering is the following

/* ALMFILE format - file that keeps the fourier coefficients of map deconvolution (eg. from ccSHTx program) */

/* binary file of double numbers pairs. */
/* each comples alm is written as RE(alm) IM(alm) for each l and -l<m<l */
/* the alms ordering is from left to right and from top to bottom of a  */
/* triangular set of alm numbers eg. */

/*                                      (4) */
/*                                      a00 */
/*                               (2)    (5)     (7) */
/*                              a1-1    a10     a11 */
/*                       (1)     (3)    (6)     (8)    (9) */
/*                      a2-2    a2-1    a20     a21    a22 */
/* and so on */

/*   the each of the following how cases saves or loads only the real (R) and unreal (I) parts of the  */
/*     alm. the modul (D) and the phase (P) are calculated on fly while reading. */

/* how - 1 and 11 reads/saves the whole alm file as for ccSHTx program form l=0 to the very last l; how =11 is for saving to txt files in long format with l,m vals. included */
/* how - 2 reads/saves the full size alms like in 1 but to the sepatare  files for each multipole from 0 to l and with  zeroed all multipoles but one */
/* how - 3 reads/saves the like in 2 but zeroed are the multipoles greater than some given one */
/* how - 4 reads/saves only one given multipole in the file --- !!! this is not compatible with ccSHTx program */ 
  //- if the option whatmultipole = -1 is passed then all multipoles are stored to separate files otherwise the indicated multipole only is stored in the file. this also applies to options how=2 and how=3
  // also in how = 4 the file format is different because what is stored is just one multipole - i.e. the vector. but still this function can be used to load and save the individual multipoles to memory and to combine
  // them with multipoles from some other maps or whatever.
  // the function 4 can save the and read all the fourier space just like the how =1 option but from files that have the information on multipoles stored separatelly.
  // mayebe this is sometimes more convenient.

/* how - 5 read/save the alms in the format according to rising l and within the same l rising m from m=-l to m=l, the file extension will be different in this case */
/*            this method allows for partial file reading upto some requestd multipole which is faster */

  //*************************************************************
  //*************************************************************
  //************************************************************

  if (strcmp(what,"F") == 0) { // doing Fourier stuff

    if (how == 1 || how == 11) { // doing -Falxxx  --- this is not sensitive to value of whatmultipole
      if (how == 1) sprintf(tmpch,"%s-Fal%li",whereto,lmax); 
      if (how == 11) { sprintf(tmpch,"%s-Fal%li_long",whereto,lmax); calculate_C_l(0,lmax,1);}
      if (strcmp(whattodo,"save") == 0) {   printf("writting "); }
      if (strcmp(whattodo,"load") == 0) {   printf("reading "); }
      printf("all a_lms to/from a single ");
      if (strcmp(fileformat,"bin") == 0 )  {   sprintf(whereto,"%s-bin",tmpch); printf("binary "); }
      if (strcmp(fileformat,"txt") == 0 )  {   sprintf(whereto,"%s-txt",tmpch); printf("text "); }
      printf("file: %s\n",whereto);

      if (strcmp(whattodo,"save") == 0) {
        if (strcmp(fileformat,"bin") == 0) {f=fopen(whereto,"wb"); }
        if (strcmp(fileformat,"txt") == 0) {f=fopen(whereto,"w"); }
        //for (l=0;l<lmax;l++) {
	//for (m=-l;m<=l;m++) {

	alms_num = alm2num(lmax,lmax)+1;
	for (i=0;i<alms_num;i++) {
	  RE = alm[i].R; IM = alm[i].I;
	  if (strcmp(fileformat,"bin") == 0) {	    fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	  if (strcmp(fileformat,"txt") == 0) {	    
	    if (how == 1)  fprintf(f,"%lE %lE\n",RE,IM); 
	    if (how == 11) { tmpd=RE*RE+IM*IM; num2alm(i,&l,&m); fprintf(f,"%li %li %lE %lE %lE %lE\n",l,m,RE,IM,sqrt(tmpd),tmpd/((double)(2*l+1)*C_l[l][1])); 	  }}
	}
	
        fclose(f);
      }

      if (strcmp(whattodo,"load") == 0) {
	if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	makekill_space_manager(whattodo,what,how);

        if (strcmp(fileformat,"bin") == 0) {f=fopen(whereto,"rb"); }
        if (strcmp(fileformat,"txt") == 0) {f=fopen(whereto,"r"); }

	//if (strcmp(fileformat,"bin") == 0) {	    fread(alm,16,alms_num,f); }
	alms_num = alm2num(lmax,lmax)+1;
	for (i=0;i<alms_num;i++) {
	  if (strcmp(fileformat,"bin") == 0) {	    fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	  if (strcmp(fileformat,"txt") == 0) 	  fscanf(f,"%lE %lE",&RE,&IM);   
	  alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0.0);
	  //printf("RE=%lf RE=%lf RE=%lf RE=%lf\n",RE,IM,alm[i].D, alm[i].P); // debug thing
         }
        fclose(f);
      }
    }
	//--------------------------------------------------------------------
    if (how == 2) { // doing -F0lxxx
      //long int tmpnum;
      strcpy(tmpch,whereto); strcat(tmpch,"-F0l"); 
      if (strcmp(whattodo,"save") == 0) {   printf("writting "); }
      if (strcmp(whattodo,"load") == 0) {   printf("reading "); }
      printf("a_lms with zeroed all but one to/from ");
      if (strcmp(fileformat,"bin") == 0 )  {   printf("binary "); }
      if (strcmp(fileformat,"txt") == 0 )  {   printf("text "); }
      printf("files: %slxxx\n",tmpch);

      if (strcmp(whattodo,"save") == 0) {
	if (whatmultipole == -1) { 
	  for(l=0;l<=lmax;l++) {
	    if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l); f=fopen(whereto,"wb"); }
	    if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l); f=fopen(whereto,"w"); }

	    for (i=0;i<alms_num;i++) { 
	      num2alm(i,&tmpl,&m);
	      if (tmpl != l) { RE = IM = 0; } else { RE = alm[i].R; IM = alm[i].I; }
	      
	      if (strcmp(fileformat,"bin") == 0) {	fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	      if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	    }
	    fclose(f);
	  }
	}
	else { l=whatmultipole;
	  if (strcmp(fileformat,"bin") == 0) {sprintf(whereto,"%s%li-bin",tmpch,l);  f=fopen(whereto,"wb"); }
	  if (strcmp(fileformat,"txt") == 0) {sprintf(whereto,"%s%li-txt",tmpch,l);  f=fopen(whereto,"w"); }

	  for (i=0;i<alms_num;i++) { 
	    num2alm(i,&tmpl,&m);
	    if (tmpl != l) { RE = IM = 0; } else { RE = alm[i].R; IM = alm[i].I; }

	    if (strcmp(fileformat,"bin") == 0) {	fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	    if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	  }
	  fclose(f);
	}
      }

      if (strcmp(whattodo,"load") == 0) { // this part (loading)  don't make any sens, I keep it here from completness reasons
	                                   //  for now it reads from all zerod files everything overwritting the informations read from the previous file
	if (whatmultipole == -1) {         // this action is performed from l=0 to l = lmax - so the table for alms is overwritten lmax times; the effect would be the 
	  for(l=0;l<=lmax;l++) {           // same just to read the last file and that's it. this option is maybe for developnent purposes.
	    if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);  f=fopen(whereto,"rb"); }
	    if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);  f=fopen(whereto,"r"); }
/* 	    printf("////jeszcze zyje,: %i laduje %s\n",l,whereto);	     */
/* 	    printf("/1//almloaded: %i\n",alms_loaded); */
	    if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    } // allocation of memory
	    if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	    makekill_space_manager(whattodo,what,how);
/* 	    printf("/2//almloaded: %i\n",alms_loaded); */

	    for (i=0;i<alms_num;i++) {
	      if (strcmp(fileformat,"bin") == 0) {	    fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	      if (strcmp(fileformat,"txt") == 0) {     fscanf(f,"%lE %lE",&RE,&IM); }
	      alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0);
	    }
/* 	    printf("/3//almloaded: %i, %i\n",alms_loaded,i); */
	    fclose(f);

	  }
	} else { l=whatmultipole;
	if (strcmp(fileformat,"bin") == 0) {  sprintf(whereto,"%s%li-bin",tmpch,l); f=fopen(whereto,"rb"); }
	if (strcmp(fileformat,"txt") == 0) {  sprintf(whereto,"%s%li-txt",tmpch,l); f=fopen(whereto,"r"); }
	
	if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	makekill_space_manager(whattodo,what,how);

	for (i=0;i<alms_num;i++) {
	  if (strcmp(fileformat,"bin") == 0) {	    fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	  if (strcmp(fileformat,"txt") == 0) {     fscanf(f,"%lE %lE",&RE,&IM);   }
	  alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0);
	}
	fclose(f);
	}
      }
    }
	//--------------------------------------------------------------------
    if (how == 3) { // doing -FSlxxx

      strcpy(tmpch,whereto); strcat(tmpch,"-FSl"); 
      if (strcmp(whattodo,"save") == 0) {   printf("writting "); }
      if (strcmp(whattodo,"load") == 0) {   printf("reading "); }
      printf("a_lms as a sum upto lmax as ");
      if (strcmp(fileformat,"bin") == 0 )  {    printf("binary "); }
      if (strcmp(fileformat,"txt") == 0 )  {    printf("text "); }
      printf("files: %slxxx\n",tmpch);

      if (strcmp(whattodo,"save") == 0) {
	if (whatmultipole == -1) { 
	  for(l=0;l<=lmax;l++) {
	    if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l); f=fopen(whereto,"wb"); }
	    if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l); f=fopen(whereto,"w"); }
	    for (i=0;i<alms_num;i++) {
	      num2alm(i,&tmpl,&m);
	      if (tmpl > l) { RE = IM =0; } else  { RE = alm[i].R; IM = alm[i].I; }
	      if (strcmp(fileformat,"bin") == 0) {	fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	      if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	    }
	    fclose(f);
	  }
	} 
	else { l=whatmultipole;   // this option makes perfect sens - because this essentially reads only once the full set of alms from some file for a given multipole
	  if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l); f=fopen(whereto,"wb"); }
	  if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l); f=fopen(whereto,"w"); }
	  
	  for (i=0;i<alms_num;i++) {
	    num2alm(i,&tmpl,&m);
	    if (tmpl > l) { RE = IM =0; } else  { RE = alm[i].R; IM = alm[i].I; }
	    if (strcmp(fileformat,"bin") == 0) {     fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	    if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	  }
	  fclose(f);
	}
      }

      if (strcmp(whattodo,"load") == 0) { // this doen't make any sens neither; I keep it only from completness and maydy developnemt reasons
	if (whatmultipole == -1) { 
	  for(l=0;l<=lmax;l++) {
	    if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);  f=fopen(whereto,"rb"); }
	    if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);  f=fopen(whereto,"r"); }

	    if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	    if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	    makekill_space_manager(whattodo,what,how);

	    for (i=0;i<alms_num;i++) {
	      if (strcmp(fileformat,"bin") == 0) {	fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	      if (strcmp(fileformat,"txt") == 0) {     fscanf(f,"%lE %lE",&RE,&IM);   }
	      alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0);
	    }
	    fclose(f);
	  }
	}
	else { l=whatmultipole;  // this option makes perfect sens - because this essentially reads only once the full set of alms from some file for a given multipole
	  if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);  f=fopen(whereto,"rb"); }
	  if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);  f=fopen(whereto,"r"); }
	  
	  if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	  if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	  makekill_space_manager(whattodo,what,how);

	  for (i=0;i<alms_num;i++) {
	    if (strcmp(fileformat,"bin") == 0) {     fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	    if (strcmp(fileformat,"txt") == 0) {     fscanf(f,"%lE %lE",&RE,&IM);   }
	    alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0);
	  }
	  fclose(f);
	}
      }
    }

	//--------------------------------------------------------------------
    if (how == 4) { // doing -Fslxxx

      strcpy(tmpch,whereto); strcat(tmpch,"-Fsl"); 
      if (strcmp(whattodo,"save") == 0) {   printf("writting "); }
      if (strcmp(whattodo,"load") == 0) {   printf("reading "); }
      printf("a_lms for a single multipole as ");
      if (strcmp(fileformat,"bin") == 0 )  {    printf("binary "); }
      if (strcmp(fileformat,"txt") == 0 )  {    printf("text "); }
      printf("files: %sxxx\n",tmpch);

      if (strcmp(whattodo,"save") == 0) {
	if (whatmultipole == -1) { 
	  for(l=0;l<=lmax;l++) {
	    if (strcmp(fileformat,"bin") == 0) {sprintf(whereto,"%s%li-bin",tmpch,l);   f=fopen(whereto,"wb"); }
	    if (strcmp(fileformat,"txt") == 0) {sprintf(whereto,"%s%li-txt",tmpch,l);   f=fopen(whereto,"w"); }
	    for (m=-l;m<=l;m++) {
	      RE = alm[alm2num(l,m)].R; IM = alm[alm2num(l,m)].I;
	      if (strcmp(fileformat,"bin") == 0) {	fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	      if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	    }
	    fclose(f);
	  }
	}
	else { l=whatmultipole; 
	  if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);  f=fopen(whereto,"wb"); }
	  if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);  f=fopen(whereto,"w"); }
	  for (m=-l;m<=l;m++) {
	    RE = alm[alm2num(l,m)].R; IM = alm[alm2num(l,m)].I;
	    if (strcmp(fileformat,"bin") == 0) {	fwrite(&RE,sizeof(RE),1,f); fwrite(&IM,sizeof(IM),1,f); }
	    if (strcmp(fileformat,"txt") == 0) {     fprintf(f,"%lE %lE\n",RE,IM);	  }
	  }
	  fclose(f);
	}
      }

      if (strcmp(whattodo,"load") == 0) { // this makes sens but it's the same thing as loading the whole map -- unless of course something has been changed in some particular multipole etc.
	if (whatmultipole == -1) { 
	  for(l=0;l<=lmax;l++) {
	    if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);   f=fopen(whereto,"rb"); }
	    if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);   f=fopen(whereto,"r"); }

	    int tmplmax = lmax; // keep the current lmax. and alms_num it cannot be changed by read_alms_parametera procedures since only one multipole is loaded
	    long int tmpalms_num = alms_num; // remember to set the lmax sufficiently big outside of this procedure so what you load here fits into memory.
	    // space allocation below is done accrding to the lmax value. when loading more then one single multipoles, for higher ls you will 
	    // reach out of memory dedicated to alms sicne the space is allocated only once and it's size cannot be changed when loading some more
	    if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	    if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	    lmax = tmplmax; // restore the old lmax value;
	    alms_num = tmpalms_num; //restore the old alms_num value;
	    makekill_space_manager(whattodo,what,how);

	    for (m=-l;m<=l;m++) {
	      if (strcmp(fileformat,"bin") == 0) {	fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	      if (strcmp(fileformat,"txt") == 0) {     fscanf(f,"%lE %lE",&RE,&IM);   }
	      alm[alm2num(l,m)].R = RE; alm[alm2num(l,m)].I = IM; alm[alm2num(l,m)].D = sqrt(RE*RE+IM*IM); alm[alm2num(l,m)].P = cpeds_cart2sph(1,RE,IM,0);
	    }
	    fclose(f);
	  }
	}
	else { l=whatmultipole; 
	  if (strcmp(fileformat,"bin") == 0) { sprintf(whereto,"%s%li-bin",tmpch,l);   f=fopen(whereto,"rb"); }
	  if (strcmp(fileformat,"txt") == 0) { sprintf(whereto,"%s%li-txt",tmpch,l);   f=fopen(whereto,"r"); }

	  int tmplmax = lmax; // keep the current lmax. it cannot be changed by read_alms_parametera procedures since only one multipole is loaded
	  long int tmpalms_num = alms_num; // remember to set the lmax sufficiently big outside of this procedure so what you load here fits into memory.
	  // space allocation below is done accrding to the lmax value. when loading more then one single multipoles, for higher ls you will 
	  // reach out of memory dedicated to alms sicne the space is allocated only once and it's size cannot be changed when loading some more
	  if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	  if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	  lmax = tmplmax; // restore the old lmax value;
	  alms_num = tmpalms_num; //restore the old alms_num value;
	  makekill_space_manager(whattodo,what,how);

	  for (m=-l;m<=l;m++) {
	    if (strcmp(fileformat,"bin") == 0) {    fread(&RE,sizeof(RE),1,f); fread(&IM,sizeof(IM),1,f);}
	    if (strcmp(fileformat,"txt") == 0) {    fscanf(f,"%lE %lE",&RE,&IM);   }
	    alm[alm2num(l,m)].R = RE; alm[alm2num(l,m)].I = IM; alm[alm2num(l,m)].D = sqrt(RE*RE+IM*IM); alm[alm2num(l,m)].P = cpeds_cart2sph(1,RE,IM,0);
	  }
	  fclose(f);
	}
      }
    }

    //--------------------------------------------------------------------
/*     the file extension for this files will be eg. Fol1024 */
    if (how == 5) { // doing -Folxxx
      if (whatmultipole == -1) whatmultipole=lmax;

      sprintf(tmpch,"%s-Fol%li",whereto,lmax); 
/*       sprintf(tmpch,"%s-Fol%li",whereto,whatmultipole);  */
      if (strcmp(whattodo,"save") == 0) {   printf("writting "); }
      if (strcmp(whattodo,"load") == 0) {   printf("reading "); }
      printf("all a_lms to/from a single ");
      if (strcmp(fileformat,"bin") == 0 )  {   sprintf(whereto,"%s-bin",tmpch); printf("binary "); }
      if (strcmp(fileformat,"txt") == 0 )  {   sprintf(whereto,"%s-txt",tmpch); printf("text "); }
      printf("file: %s\n",whereto);

      if (strcmp(whattodo,"save") == 0) {
        if (strcmp(fileformat,"bin") == 0) { f=fopen(whereto,"wb"); }
        if (strcmp(fileformat,"txt") == 0) { f=fopen(whereto,"w"); }

	for (l=0;l<=whatmultipole;l++) {
	  for (m=-l;m<=l;m++) {
	    i=alm2num(l,m);
	    RE = alm[i].R; IM = alm[i].I;
	    if (strcmp(fileformat,"bin") == 0) { fwrite(&RE,sofd,1,f); fwrite(&IM,sofd,1,f); }
	    if (strcmp(fileformat,"txt") == 0)   fprintf(f,"%lE %lE\n",RE,IM); 
	  }
	}
        fclose(f);
      }

      if (strcmp(whattodo,"load") == 0) {
	if (strcmp(fileformat,"bin") == 0) {      read_binalm_parameters(whereto);    }
	if (strcmp(fileformat,"txt") == 0) {      read_txtalm_parameters(whereto);    }
	
	if (whatmultipole != lmax && whatmultipole > -1) set_alms_lmax(whatmultipole);
	makekill_space_manager(whattodo,what,how);

        if (strcmp(fileformat,"bin") == 0) { f=fopen(whereto,"rb"); }
        if (strcmp(fileformat,"txt") == 0) { f=fopen(whereto,"r"); }

	//if (strcmp(fileformat,"bin") == 0) {	    fread(alm,16,alms_num,f); }
	if (strcmp(fileformat,"bin") == 0) {
	  for (l=0;l<=whatmultipole;l++) {
	    for (m=-l;m<=l;m++) {
	      fread(&RE,sofd,1,f); fread(&IM,sofd,1,f); 
	      i=alm2num(l,m);
	      alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0.0);
	    }
	  }
	}
	if (strcmp(fileformat,"txt") == 0) {
	  for (l=0;l<=whatmultipole;l++) {
	    for (m=-l;m<=l;m++) {
	      fscanf(f,"%lE %lE",&RE,&IM);   
	      i=alm2num(l,m);
	      alm[i].R = RE; alm[i].I = IM; alm[i].D = sqrt(RE*RE+IM*IM); alm[i].P = cpeds_cart2sph(1,RE,IM,0.0);
	    }
	  }
	}

        fclose(f);
      }


    }
    alms_loaded = 1;
  }


}
  //*************************************************************

