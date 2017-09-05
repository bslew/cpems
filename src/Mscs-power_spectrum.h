// this class is depreciated - don't use it; use Mscs-angular_power_spectrum.h instead.
#ifndef MSCS_ANGULAR_POWER_SPECTRUM
#define MSCS_ANGULAR_POWER_SPECTRUM

#include "Mscs-function.h"

/*!
  \class mscsAngularPowerSpectrum
  \brief Encapsulates the power spectrum structure and provides a basic functionality
  \details 

  \note the name of this class has changed from mscsPowerSpectrum to mscsAngularPowerSpectrum to better reflect what it is.
  \date 2009/06/04 23:50:43 
  \author Bartosz Lew
*/
class mscsAngularPowerSpectrum : public mscsFunction {

 public:
  
  typedef struct {
    double l; // effective multipole number
    double C; // power 
    double err; // error
    double cv; // cosmic vatiance
    bool set; // if set flag
  } Cltype;

  Cltype *C_l;

  // power spectrum specifications variables
  long lmax;
  long lmin; // currently this is always assumed to be 0

  double Clmax, Clmin; // minimal and maximal values of Cl
  long l_Clmax, l_Clmin; // l values corresponding to minimal and maximal values of Cl
  string object_name;

  long quietness; // describes how quiet the object should be (0 - no msgs, 1 - msgs on)
  // and so on


  //  mscsPowerSpectrum(long lmax_loc) {
  // the convention is the the Clsize must include the l=0,1 multipoles even if they are not used.


  mscsPowerSpectrum(long Clsize, long quiet) {
    //long i;
    C_l=NULL;
    lmin=0;
    object_name="Cl";
    quietness=quiet;
    initiate_mscsPowerSpectrum(Clsize);
  }

  mscsPowerSpectrum(long Clsize) {
    //long i;
    C_l=NULL;
    lmin=0;
    object_name="Cl";
    quietness=1;
    initiate_mscsPowerSpectrum(Clsize);
  }

  mscsPowerSpectrum() {
    //long i;
    //if (Clsize >= 0) { lmax = Clsize-1; } else { lmax = 0; Clsize=1; }
    //make_Cl();
    //reset_Cl();
    C_l=NULL;
    lmin=0;
    l_Clmin = -1; l_Clmax=-1;
    quietness=1;
    object_name="Cl";
  }

 void initiate_mscsPowerSpectrum(long Clsize) {
    if (Clsize >= 0) { lmax = Clsize-1; } else { lmax = 0; Clsize=1; }
    make_Cl(Clsize);
    reset_Cl();
    if (quietness!=0) printf("|%s> * initiated angular power spectrum of size: %li\n",object_name.c_str(),Clsize);
    //    l_Clmin = 0; l_Clmax=lmax;
  }
/*   mscsPowerSpectrum(long lmin_loc, long lmax_loc) { */
/*     lmax = lmax_loc; */
/*     C_l = new Cltype[lmax+1]; */
/*   } */

  ~mscsPowerSpectrum() {
    kill_Cl();
  }

  // checks for largest nonzero l and rescales the power spectrum accordingly
  // how: 1 - lookin on C; 2 - looking on error; 3 -- looking on cv
  void shrink_in_lmax(long how) {
    long i,lmaxn0=0;
    Cltype * C_lloc=NULL;
    
    // find largest nonzero l
    if (how == 1) for (i=lmin;i<=lmax;i++) { if (get_Cl(i) != 0) lmaxn0=i; }
    if (how == 2) for (i=lmin;i<=lmax;i++) { if (get_err(i) != 0) lmaxn0=i; }
    if (how == 3) for (i=lmin;i<=lmax;i++) { if (get_cv(i) != 0) lmaxn0=i; }

    if (quietness!=0) printf("|%s>  -- shrinking power spectrum to lmax %li\n",object_name.c_str(),lmaxn0);

    if (lmaxn0 != 0) {
      // prepare a new structure with the new size
      C_lloc = new Cltype[lmaxn0+1];
      //copy the data into the new structure
      for (i=lmin;i<=lmaxn0;i++) { C_lloc[i] = C_l[i]; }
      //delete old power spectrum
      kill_Cl();
      // set new power spectrum stuff
      C_l = C_lloc;    
    }

    lmax=lmaxn0;
  }

  void zero_below_threshold(double th) {
    long l;
    for (l=lmin;l<=lmax;l++) { if (C_l[l].C<=th) C_l[l].C=0.0;  }    
  }

  void set_name(string name) {
    printf("|%s> -- changing name of object to: %s\n", object_name.c_str(), name.c_str()); 
    object_name=name;
  }
  string get_name() { return object_name; }
  void set_quietness(long val) {
    quietness=val;
    if (val !=0)
    printf("|%s> -- changing object quietness to: %li\n", object_name.c_str(), quietness); 
  }

  long get_lmin() { return 0; }
  long get_lmax() { return lmax; }
  long get_size() { return lmax+1; }
  double get_Cl(long l) { if (l >= 0 && l <= lmax) return C_l[l].C; else return 0; }
  void set_Cl(long l,double C) { if (l >= 0 && l <= lmax) { C_l[l].C = C; if (C_l[l].set == false) C_l[l].l = l; C_l[l].set = true; }} // here development could be to resize the object if a call outsize the size is given
  void set_C_to(double val) { long i; for (i=lmin;i<=lmax;i++)  C_l[i].C = val;  }
  void reset_Cl() { long i;   for (i=lmin;i<=lmax;i++) { C_l[i].set = false; C_l[i].l = (double)i; C_l[i].C = 0; } }
  void set_l(long l,double dl) { C_l[l].l = dl; C_l[l].set = true;}
  double get_l(long l) {  if (l >= 0 && l <= lmax) return C_l[l].l; else return 0; }

  void set_err(long l,double err) { if (l >= 0 && l <= lmax) { C_l[l].err = err; } }
  void set_cv(long l,double cv) { if (l >= 0 && l <= lmax) { C_l[l].cv = cv; } }
  double get_err(long l) { if (l >= 0 && l <= lmax) { return C_l[l].err; } else return 0; }
  double get_cv(long l) { if (l >= 0 && l <= lmax) { return C_l[l].cv; } else return 0; }

  void compute_cosmic_variance() {
    long i;
    for (i=lmin;i<=lmax;i++) { if (C_l[i].set) { C_l[i].cv = cosmic_variance(i); }}
  }
  void compute_minmax_stats() { // finds the Cl_min/max l_Cl_min/max
    // to be implmented
  } 

  void save(string filename) {
    long l;
    FILE* f;
    f = fopen(filename.c_str(),"w");
    printf("  -- saving the power spectrum to file %s\n",filename.c_str());
    for (l=lmin;l<=lmax;l++) {      fprintf(f,"%lE %lE\n",C_l[l].l,C_l[l].C);     }
    fclose(f);
  }

  void saveint(string filename) {
    long l;
    FILE* f;
    f = fopen(filename.c_str(),"w");
    printf("  -- saving the power spectrum to file %s\n",filename.c_str());
    for (l=lmin;l<=lmax;l++) {      fprintf(f,"%li %lE\n",(long)C_l[l].l,C_l[l].C);     }
    fclose(f);
  }

  void save_all(string filename) {
    long l;
    FILE* f;
    f = fopen(filename.c_str(),"w");
    printf("  -- saving the power spectrum with errors to file %s\n",filename.c_str());
    for (l=lmin;l<=lmax;l++) {      fprintf(f,"%lE %lE %lE %lE\n",C_l[l].l,C_l[l].C,C_l[l].err,C_l[l].cv);     }
    fclose(f);
  }

  void kill_Cl() { delete [] C_l; C_l=NULL; }
  void make_Cl(long size) { if (C_l != NULL) kill_Cl(); C_l = new Cltype[size];    lmin = 0; lmax=size-1;      }
    

  void load(string filename) {
    long l,n=0;
    FILE* f;
    f = fopen(filename.c_str(),"r"); if (f != NULL) {  n=0; while (fscanf(f,"%*E %*E") != EOF) { n++;}         fclose(f); } else { printf("ERROR: no such file %s\n",filename.c_str()); }
    printf("  -- the file to load has %li lines\n",n);
    if (n != lmax+1) { printf("  -- resizing the existing power spectra to size: %li\n",n+1); kill_Cl(); make_Cl(n); reset_Cl(); }
    printf("  -- reading the power spectrum from file %s\n",filename.c_str());
    f = fopen(filename.c_str(),"r"); for (l=lmin;l<=lmax;l++) {      fscanf(f,"%lE %lE",&C_l[l].l,&C_l[l].C);     }
    fclose(f);
  }

  void flush() {
    long l;
    printf("\n");
    for (l=lmin;l<=lmax;l++) {      printf("%lE %lE\n",C_l[l].l,C_l[l].C);     }
    printf("\n");
  }
  
  /*   void set_Cl_multipole_range(long l1, long l2, double val); */
  /*   void Cl_multiply(double factor); // multiplies the C_l power spectrum by a factor */






  /*!
    \brief binns the given power spectrum and returns the binned power spectrum object.
    \details 
    @param lmin_loc - bin from this multipole
    @param lmax_loc - bin to this multipole
    @param bintabs - size of the bintab array
    @param bintab - array containing the bin sizes
    @return returns the binned power spectrum
    
    binning is only from lmin to lmax, outside this range C_l is copied without change;
    binning is done according to the bintab table of size bintabs that holds a set of long type numbers defining 
    how many ls to bin using weights w=1/bintab[i] 
    or weights from the err field (previously the 3rd col) of the C_l structure for each l in which case w should be 0: w=0; otherwise it should be w=1
    the sum of values stored in bintab array should be lmax-lmin+1
    value in array bintabs - 1 corresponds to no binning; 2 - means that two multipoles are binned and so on.
    
    \date 2009/06/04 23:31:24 
    \author Bartosz Lew
  */
  mscsPowerSpectrum* bin_C_l(mscsPowerSpectrum* c_l, long lmin_loc, long lmax_loc, long *bintabs, long* bintab, double w) {
    long i,j,k,binl_sum=0;
    long Cls, Cbs,cls,cbs;
    mscsPowerSpectrum* cb;
    long tmp;
    double* effl;
    bool EQweights;
    double wsum;
  
    printf("|%s> * Binning the power spectrum from %li to %li, together: %li multipoles\n",object_name.c_str(),lmin_loc,lmax_loc,lmax_loc-lmin_loc+1);
    printf("  -- bins are:\n");
    
    // check safety conditions
    lmax_loc = cpeds_get_min(lmax_loc,lmax);
    lmin_loc = cpeds_get_max(lmin_loc,0);
    if (*bintabs <= 0 || bintab == NULL) { printf("ERROR: wrong aruments\n"); /* cb = new power_spectrum(-1); */ return NULL; }
    for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; printf("%li ",bintab[i]); } printf("\n  -- The overall sum of the multipoles to bin as implied by binning array is: %li \n",binl_sum); 
    if (binl_sum < lmax_loc-lmin_loc+1) { printf("|%s>  -- WARNING !! the binning array is incompatible with requested binning multipole range. Will bin upto the end of the binning array: i.e. upto l=%li.\n",object_name.c_str(),lmin_loc + binl_sum - 1); }
    if (binl_sum < lmax_loc-lmin_loc+1) { lmax_loc = lmin_loc + binl_sum - 1; }
    if (binl_sum >  lmax_loc-lmin_loc+1) { printf("|%s>  -- WARNING !! the binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning array: i.e. upto l==",object_name.c_str()); }
    if (binl_sum >  lmax_loc-lmin_loc+1) {
      printf("number of bins is: %li.\n",*bintabs);
      while (binl_sum > lmax_loc-lmin_loc+1) {
	(*bintabs)--;
	binl_sum=0;
	for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum);
      }
      printf("number of bins is : %li.\n",*bintabs);
      (*bintabs)++;
      bintab[*bintabs-1]=lmax_loc-lmin_loc+2-binl_sum;
      printf("adding last bin of size: %li and increasing number of bins to: %li.\n",bintab[*bintabs-1], *bintabs);
      binl_sum=0;
      for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum);
      printf("bin sum is now: %li.\n",binl_sum);
    }
    
    if (*bintabs <= 0 || bintab == NULL) { printf("ERROR: wrong aruments\n"); /* cb = new power_spectrum(-1); */ return NULL; }
    
    // find out the size of the binned C_b
    cls = lmax+1; // original size of the C_l
    cbs = lmax + 1  - (lmax_loc-lmin_loc+1) + *bintabs; // total size of binned C_l: bs
    Cls = binl_sum; // size of original C_l for binning
    Cbs = *bintabs; // size of binned C_l: bs
    if (w == 0) { EQweights = false; set_cosmic_variance(); } else { EQweights = true; }
    
    printf("  -- total size of cl vector before binning: %li\n",cls);
    printf("  -- size of Cl vector for binning: %li\n",Cls);
    printf("  -- total size of cb - binned Cl: %li\n",cbs);
    printf("  -- size of Cb binned vector: %li\n",Cbs);
    
    
    // define the structures needed
    matrix <double> Cl(Cls,1); // original C_l vector for binning
    matrix <double> Cb(Cbs,1); // binned C_b vector
    matrix <double> M(Cbs,Cls); // binning operator
    cb = new mscsPowerSpectrum(cbs);
    effl = new double[*bintabs];
    
    // copy the power spectra outside the reigon for binning
    for (i=0;i<lmin_loc;i++) { cb->set_l(i,c_l->get_l(i)); cb->set_Cl(i,c_l->get_Cl(i)); } // lower end
    tmp = lmax_loc+1; j=lmin_loc+(*bintabs);
    for (i=tmp;i<=lmax;i++)  { cb->set_l(j,c_l->get_l(i)); cb->set_Cl(j,c_l->get_Cl(i)); j++;} // upper end shifted as to fit the binned Cb in range lmin_loc, lmax_loc
    
    // prepare the power spectrum vector for binning
    for (j=0;j<Cls;j++) {  Cl(j,0) = c_l->get_Cl(j+lmin_loc); } // copy the power spectra part for binning
    
    // prepare the binning matrix operator
    
    for (i=0;i<Cbs;i++) 
      for (j=0;j<Cls;j++) 
	M(i,j) = 0;  // zero to binning matrix
    
    k=0;
    for (i=0;i<Cbs;i++) {
      //if (i==0) jst = i; else jst = bintab[i]-1;
      if (EQweights) w = 1/(double)bintab[i]; 
      effl[i] = 0; wsum=0;
      for (j=k;j<k+bintab[i];j++) {
	if (!EQweights) { w = c_l->get_err(j+lmin_loc); wsum+=w; }
	M(i,j) = w;  
	effl[i]+=(double)(j+lmin_loc)*w;      
      }
      if (!EQweights) { for (j=k;j<k+bintab[i];j++) { M(i,j)=M(i,j)/wsum; } effl[i]/=wsum; }    
      printf("effective ls: %lf\n",effl[i]);
      k+=bintab[i];
    }
    
    // do the binning
    Cb = M*Cl;
    
    // rewrite the Cb onto the object
    
    for (i=0;i<Cbs;i++)  { 
      j=lmin_loc+i; cb->set_l(j,effl[i]); cb->set_Cl(j,Cb(i,0)); 
      printf("l %lf Cl %lE            binned Cb: l %lf Cb %lE\n",effl[i],Cb(i,0),cb->get_l(j),cb->get_Cl(j));
    } 
    
    delete [] effl;
    
    return cb;
    
  }




  void divide_llpotwoPI() {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C /= (get_l(l)*(get_l(l)+1)/twoPI); }
  }
  void multiply_llpotwoPI() {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C *= (get_l(l)*(get_l(l)+1)/twoPI); }
  }



  //! sets the cosmic variance field to its default values for cosmic variance according to l
  void set_cosmic_variance() {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].cv = cosmic_variance(l); }
  }



  //------------------------------------------------------------
  // definition of some usefull operators
  //------------------------------------------------------------
  mscsPowerSpectrum & operator += (mscsPowerSpectrum &cl) {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].C+=cl.get_Cl(l); }
    return *this;
  }
  mscsPowerSpectrum & operator -= (mscsPowerSpectrum &cl) {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].C-=cl.get_Cl(l); }
    return *this;
  }
  mscsPowerSpectrum & operator *= (mscsPowerSpectrum &cl) {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].C*=(cl.get_Cl(l)); }
    return *this;
  }
  mscsPowerSpectrum & operator /= (mscsPowerSpectrum &cl) {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C/=cl.get_Cl(l); }
    return *this;
  }


  mscsPowerSpectrum & operator += (double d) {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].C+=d; }
    return *this;
  }
  mscsPowerSpectrum & operator -= (double d) {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C-=d; }
    return *this;
  }
  mscsPowerSpectrum & operator *= (double d) {
    long l;
    for (l=lmin;l<=lmax;l++) { C_l[l].C*=d; }
    return *this;
  }
  mscsPowerSpectrum & operator /= (double d) {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C/=d; }
    return *this;
  }
 

  //------------------------------------------------------------
/*   mscsPowerSpectrum & operator + (double d) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp+=d;        return temp;  } */

/* inline  mscsPowerSpectrum*  operator + (mscsPowerSpectrum cl) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=this;  */
/*     long l; */
/*     for (l=lmin;l<=lmax;l++) { temp.set_l(l,cl.get_l(l)); temp.set_Cl(l,cl.get_Cl(l)); } */
/*     return &temp; */
/*   } */
/*   //------------------------------------------------------------ */
/*   mscsPowerSpectrum & operator - (double d) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp-=d;        return temp;  } */
/*   mscsPowerSpectrum & operator - (mscsPowerSpectrum &cl) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp-=cl;    return temp;  } */
/*   //------------------------------------------------------------ */
/*   mscsPowerSpectrum & operator * (double d) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp*=d;        return temp;  } */

/* /\*   inline mscsPowerSpectrum & operator * (const mscsPowerSpectrum &cl, const mscsPowerSpectrum &dl) { *\/ */
/*   mscsPowerSpectrum & operator * (mscsPowerSpectrum &cl) {  */
/* /\*     mscsPowerSpectrum temp(lmax+1);    temp=*(*this); temp*=cl;    return temp;  } *\/ */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp*=cl;    return temp; } */



  //------------------------------------------------------------
/*   mscsPowerSpectrum & operator / (double d) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp/=d;        return temp;  } */
/*   mscsPowerSpectrum & operator / (mscsPowerSpectrum &cl) { */
/*     mscsPowerSpectrum temp(lmax+1);    temp=*this; temp/=cl;    return temp;  } */


  mscsPowerSpectrum & operator = (mscsPowerSpectrum &cl) {
    long l;
    for (l=lmin;l<=lmax;l++) {  C_l[l].C=cl.get_Cl(l); C_l[l].l=cl.get_l(l);}
    return *this;
  }

 private:
  double cosmic_variance(long l) {
    double ld=(double)l;
    return sqrt(2/(ld*(ld+1)));
  }



};

#endif
