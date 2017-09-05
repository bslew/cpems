
#ifndef MSCS_ANGULAR_POWER_SPECTRUM
#define MSCS_ANGULAR_POWER_SPECTRUM

#include "Mscs-function.h"

/*!
  \class mscsAngularPowerSpectrum
  \brief Encapsulates the angular power spectrum structure and provides a basic functionality
  \details

  \note the name of this class has changed from mscsAngularPowerSpectrum to mscsAngularPowerSpectrum to better reflect what it is.
  \date 2009/06/04 23:50:43
  \author Bartosz Lew
*/
class mscsAngularPowerSpectrum : public mscsFunction {

 public:



/* **************************************************************************************************** */
 mscsAngularPowerSpectrum() : mscsFunction("Cl") {  }
 mscsAngularPowerSpectrum(string name) : mscsFunction(name) {  }

/* **************************************************************************************************** */
  //! this is depreciated - not used at the moment
  void initiate_mscsAngularPowerSpectrum(long Clsize);
/* **************************************************************************************************** */
 ~mscsAngularPowerSpectrum() { }

/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
// some old stuff
/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
  /* mscsAngularPowerSpectrum(long Clsize, long quiet) { */
  /*   //long i; */
  /*   C_l=NULL; */
  /*   lmin=0; */
  /*   object_name="Cl"; */
  /*   quietness=quiet; */
  /*   initiate_mscsAngularPowerSpectrum(Clsize); */
  /* } */

  /* mscsAngularPowerSpectrum(long Clsize) { */
  /*   //long i; */
  /*   C_l=NULL; */
  /*   lmin=0; */
  /*   object_name="Cl"; */
  /*   quietness=1; */
  /*   initiate_mscsAngularPowerSpectrum(Clsize); */
  /* } */


/* **************************************************************************************************** */
  /* // checks for largest nonzero l and rescales the power spectrum accordingly */
  /* // how: 1 - lookin on C; 2 - looking on error; 3 -- looking on cv */
  /* void shrink_in_lmax(long how) { */
  /*   long i,lmaxn0=0; */
  /*   Cltype * C_lloc=NULL; */

  /*   // find largest nonzero l */
  /*   if (how == 1) for (i=lmin;i<=lmax;i++) { if (get_Cl(i) != 0) lmaxn0=i; } */
  /*   if (how == 2) for (i=lmin;i<=lmax;i++) { if (get_err(i) != 0) lmaxn0=i; } */
  /*   if (how == 3) for (i=lmin;i<=lmax;i++) { if (get_cv(i) != 0) lmaxn0=i; } */

  /*   if (quietness!=0) printf("|%s>  -- shrinking power spectrum to lmax %li\n",object_name.c_str(),lmaxn0); */

  /*   if (lmaxn0 != 0) { */
  /*     // prepare a new structure with the new size */
  /*     C_lloc = new Cltype[lmaxn0+1]; */
  /*     //copy the data into the new structure */
  /*     for (i=lmin;i<=lmaxn0;i++) { C_lloc[i] = C_l[i]; } */
  /*     //delete old power spectrum */
  /*     kill_Cl(); */
  /*     // set new power spectrum stuff */
  /*     C_l = C_lloc;     */
  /*   } */

  /*   lmax=lmaxn0; */
  /* } */

/* **************************************************************************************************** */
 //this is depreciated -- check out the setBelow method of the mscsFunction class now
 /* void zero_below_threshold(double th) { */
 /*   long l; */
 /*   for (l=lmin;l<=lmax;l++) { if (C_l[l].C<=th) C_l[l].C=0.0;  }     */
 /* } */

/* **************************************************************************************************** */
  /* void set_name(string name) { */
  /*   printf("|%s> -- changing name of object to: %s\n", object_name.c_str(), name.c_str());  */
  /*   object_name=name; */
  /* } */
  /* string get_name() { return object_name; } */
  /* void set_quietness(long val) { */
  /*   quietness=val; */
  /*   if (val !=0) */
  /*   printf("|%s> -- changing object quietness to: %li\n", object_name.c_str(), quietness);  */
  /* } */

/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */

  /*!
    \brief  computes the returns the cosmic variance function for the allocated multipole range
    \details
    @return

    \note the name of the method changed from calculate_cosmic_variance to get_cosmic_variance
    \date 2010/02/19 14:57:45
    \author Bartosz Lew
  */
  mscsAngularPowerSpectrum get_cosmic_variance(double lmin=0, double lmax=-1);


  const mscsAngularPowerSpectrum& divide_llpotwoPI();
  const mscsAngularPowerSpectrum& multiply_llpotwoPI();
  double cosmic_variance(long l);

  cpedsStatusCodes save(string filename) const;
  cpedsStatusCodes load(string filename, bool commentedFile=false);


  /* void saveint(string filename) { */
  /*   long l; */
  /*   FILE* f; */
  /*   f = fopen(filename.c_str(),"w"); */
  /*   printf("  -- saving the power spectrum to file %s\n",filename.c_str()); */
  /*   for (l=lmin;l<=lmax;l++) {      fprintf(f,"%li %lE\n",(long)C_l[l].l,C_l[l].C);     } */
  /*   fclose(f); */
  /* } */

  /* void save_all(string filename) { */
  /*   long l; */
  /*   FILE* f; */
  /*   f = fopen(filename.c_str(),"w"); */
  /*   printf("  -- saving the power spectrum with errors to file %s\n",filename.c_str()); */
  /*   for (l=lmin;l<=lmax;l++) {      fprintf(f,"%lE %lE %lE %lE\n",C_l[l].l,C_l[l].C,C_l[l].err,C_l[l].cv);     } */
  /*   fclose(f); */
  /* } */


 /********************************************************************************************/
 /* // these are wrappers that work on the mscsFunction methods - for backward compatibility */
 /********************************************************************************************/
  long get_lmin() { return 0; }
  long get_lmax() { return getMaxArg(); }
  long get_size() { return pointsCount(); }
  double get_Cl(long l) const { return Y(l); }
  void set_Cl(long l,double C) { setf(l,C); } // here development could be to resize the object if a call outsize the size is given
  void set_C_to(double val) { setf(val); } //long i; for (i=lmin;i<=lmax;i++)  C_l[i].C = val;  }
  void reset_Cl() { setf(0.0); } //long i;   for (i=lmin;i<=lmax;i++) { C_l[i].set = false; C_l[i].l = (double)i; C_l[i].C = 0; } }
  void set_l(long l,double dl) { setarg(l,dl); } //C_l[l].l = dl; C_l[l].set = true;}
  double get_l(long l) {  return getX(l); } //if (l >= 0 && l <= lmax) return C_l[l].l; else return 0; }

  // these methods are temporary removed - if realy needed they will be reimplemented
  /* void set_err(long l,double err) { if (l >= 0 && l <= lmax) { C_l[l].err = err; } } */
  /* void set_cv(long l,double cv) { if (l >= 0 && l <= lmax) { C_l[l].cv = cv; } } */
  /* double get_err(long l) { if (l >= 0 && l <= lmax) { return C_l[l].err; } else return 0; } */
  /* double get_cv(long l) { if (l >= 0 && l <= lmax) { return C_l[l].cv; } else return 0; } */

  void kill_Cl() { clearFunction(); } //delete [] C_l; C_l=NULL; }
  void set_lmax(long lmax) { setPointsNum(lmax+1); }
  void make_Cl(long size) { setPointsNum(size); }//if (C_l != NULL) kill_Cl(); C_l = new Cltype[size];    lmin = 0; lmax=size-1;      }


  //! print the power spectrum; this is depreciated
  void flush() { print(); }
  void compute_minmax_stats() { checkRanges(); } // finds the Cl_min/max l_Cl_min/max












  /*************/
  /* OPERATORS */
  /*************/
  /* mscsAngularPowerSpectrum & operator += (mscsAngularPowerSpectrum &cl) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) { C_l[l].C+=cl.get_Cl(l); } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator -= (mscsAngularPowerSpectrum &cl) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) { C_l[l].C-=cl.get_Cl(l); } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator *= (mscsAngularPowerSpectrum &cl) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) { C_l[l].C*=(cl.get_Cl(l)); } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator /= (mscsAngularPowerSpectrum &cl) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) {  C_l[l].C/=cl.get_Cl(l); } */
  /*   return *this; */
  /* } */


  /* mscsAngularPowerSpectrum & operator += (double d) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) { C_l[l].C+=d; } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator -= (double d) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) {  C_l[l].C-=d; } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator *= (double d) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) { C_l[l].C*=d; } */
  /*   return *this; */
  /* } */
  /* mscsAngularPowerSpectrum & operator /= (double d) { */
  /*   long l; */
  /*   for (l=lmin;l<=lmax;l++) {  C_l[l].C/=d; } */
  /*   return *this; */
  /* } */



  //------------------------------------------------------------
/*   mscsAngularPowerSpectrum & operator + (double d) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp+=d;        return temp;  } */

/* inline  mscsAngularPowerSpectrum*  operator + (mscsAngularPowerSpectrum cl) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=this;  */
/*     long l; */
/*     for (l=lmin;l<=lmax;l++) { temp.set_l(l,cl.get_l(l)); temp.set_Cl(l,cl.get_Cl(l)); } */
/*     return &temp; */
/*   } */
/*   //------------------------------------------------------------ */
/*   mscsAngularPowerSpectrum & operator - (double d) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp-=d;        return temp;  } */
/*   mscsAngularPowerSpectrum & operator - (mscsAngularPowerSpectrum &cl) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp-=cl;    return temp;  } */
/*   //------------------------------------------------------------ */
/*   mscsAngularPowerSpectrum & operator * (double d) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp*=d;        return temp;  } */

/* /\*   inline mscsAngularPowerSpectrum & operator * (const mscsAngularPowerSpectrum &cl, const mscsAngularPowerSpectrum &dl) { *\/ */
/*   mscsAngularPowerSpectrum & operator * (mscsAngularPowerSpectrum &cl) {  */
/* /\*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*(*this); temp*=cl;    return temp;  } *\/ */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp*=cl;    return temp; } */



  //------------------------------------------------------------
/*   mscsAngularPowerSpectrum & operator / (double d) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp/=d;        return temp;  } */
/*   mscsAngularPowerSpectrum & operator / (mscsAngularPowerSpectrum &cl) { */
/*     mscsAngularPowerSpectrum temp(lmax+1);    temp=*this; temp/=cl;    return temp;  } */


  const mscsAngularPowerSpectrum & operator = (const mscsAngularPowerSpectrum &rhs) {
    if (this!=&rhs) { this->mscsFunction::operator=(rhs); }
    return *this;
  }

 private:



};

#endif
