#include "Mscs-alms.h"

// ****************************************************************************************************
mscsAlms::mscsAlms() : mscsObject("alms", Zero) {}
// ****************************************************************************************************
mscsAlms::mscsAlms(string name, long lmaxl) : mscsObject(name, Zero) { if (lmaxl>=0) lmax(lmaxl); }
// ****************************************************************************************************
mscsAlms::mscsAlms(const cpedsList<mscsAlm>& al) : mscsObject("alms") {
  clear();
  _alms=al;
  set_alms_num(_alms.size());
  long l,m;
  num2alm(almsNum()-1,&l,&m);
  set_lmax(l);
}
// ****************************************************************************************************
mscsAlms::mscsAlms(const mscsAlms& parent) : mscsObject(parent) {
  *this=parent;
}

// ****************************************************************************************************
mscsAlms::mscsAlms(long lmaxl) : mscsObject("alms", Zero) {
  lmax(lmaxl);
}
//************************************************************************
mscsAlm& mscsAlms::a(long l,long m) {  return _alms[alm2num(l,m)]; }
// ****************************************************************************************************
const mscsAlm& mscsAlms::get(long l,long m) const {  return _alms[alm2num(l,m)]; }

// ****************************************************************************************************
const mscsAlm& mscsAlms::get(long i) const { return _alms[i]; }
// ****************************************************************************************************
void mscsAlms::clear() { _alms.clear(); set_lmax(0); set_alms_num(0); }

// ****************************************************************************************************
void mscsAlms::num2alm(long num, long *l, long *m) const {
  long int x = 0;
  long int ll,lm;
  long elmax=lmax();
  for (lm = -elmax;lm<=elmax;lm++) {
    for (ll=abs(lm); ll <= elmax; ll++) {
      if (x == num) { *l=ll; *m=lm; lm = elmax; ll = elmax; }
      x++;
    }
  }
}
// ****************************************************************************************************
long mscsAlms::alm2num(long l, long m) const {
  long elmax=lmax();
  long lm=(elmax*(elmax+1+2*m)+m*(2-abs(m+1)))/2+l;
  return lm;
}
// ****************************************************************************************************
// this function gives the right results only for lmax and mmax. so this is only to check what is the lmax in the unknown file
// this function can also be used for numaration of alms in usulal way. not the one usen in F files and by ccSHTx program.
// this numeration is increasing from lower to higer ls and from lower to higher ms.
// the num argument is the ordering number of alm starting from 0. so a00 has num=0
void mscsAlms::num2alm_norm(long int num, long *l, long *m) const {
  long int ll,lm,llmax=lmax();
  long int x=0;
  msgs->warning("num2alm_norm: this code should be used with care--> please check the sources",High);
  for (ll=0;ll<=llmax;ll++) {
    for (lm=-ll;lm<=ll;lm++) { x++;
    if (x == num+1) {*l = ll; *m = lm; ll=llmax; lm=ll; return;}
    }
  }

}

// ****************************************************************************************************
long mscsAlms::alm2num_norm(long l, long m) const {
  long int ll;//,lm;
  long int x=0,s=1;

  for (ll=0;ll<l;ll++) {
    x = x+s;
    s=s+2;
  }
  x = x + (l+m+1);

  return x-1;
}
// ****************************************************************************************************
void mscsAlms::zero_alms_multipole_range(long l1, long l2) {
  long l,m,num;
  msgs->say("zeroing alms multipole range: from: "+msgs->toStr(l1)+" to: "+msgs->toStr(l2),High);
  for (l=l1;l<=l2;l++) {
    for (m=-l;m<=l;m++) { num = alm2num(l,m); a(num)=0.0;    }
  }
}

// ****************************************************************************************************
// void mscsAlms::set_alms_multipole_range(long l1, long l2,string what, double val) {
//   long l,m,num;
//   for (l=l1;l<=l2;l++) {
//     for (m=-l;m<=l;m++) { num = alm2num(l,m);
//       if (what == "RI") {
// 	alm[num].R = alm[num].I = val;
// 	alm[num].D = sqrt(pow(alm[num].R,2)+pow(alm[num].I,2)); alm[num].P = cpeds_cart2sph(1,alm[num].R,alm[num].I,0); // fill the rest
//       }

//       // rest TO BE IMPLEMENTED ...
//     }
//   }
// }
// ****************************************************************************************************
mscsAlms mscsAlms::tolmOrdering() const {
  long l,m,num;
  long lmaxl=lmax();
  cpedsList<mscsAlm> q;
  for (l=0;l<=lmaxl;l++) {
    for (m=-l;m<=l;m++) { num = alm2num(l,m); q.append(get(num));    }
  }
  return mscsAlms(q);
}
// ****************************************************************************************************
// calculates a single C_l value from the alms table for requested l
mscsAngularPowerSpectrum mscsAlms::get_Cl(long lminl,long lmaxl,long what) const {
  mscsAngularPowerSpectrum cl;
  if (lmaxl==-1) lmaxl=lmax();

  for (long l=lminl;l<=lmaxl;l++) {
    cl.newPoint(l,calculate_single_C_l(l,what));
  }
  return cl;
}

// ****************************************************************************************************
mscsAngularPowerSpectrum mscsAlms::get_cross_Cl(const mscsAlms& a2, long lminl,long lmaxl,long what) {
  double tmp=0;
  mscsAlm tmp2,blm;
  long i;
  mscsAngularPowerSpectrum cl;
  if (lmaxl==-1) lmaxl=lmax();

  for (long l=lminl;l<=lmaxl;l++) {
    tmp=0;
    // for (long m=-l;m<=l;m++) { i=alm2num(l,m); tmp2=a(i)*a2.a(i).conj(); tmp+=tmp2.abs(); }
    for (long m=-l;m<=l;m++) { i=alm2num(l,m); blm=a2.get(i); tmp2=a(i)*blm.conj(); tmp+=tmp2.abs(); }//.conj()
    tmp/= (2.0*double(l)+1.0);
    cl.newPoint(l,tmp);
  }
  return cl;
}
// ****************************************************************************************************
// calculates a single C_l value from the alms table for requested l
double mscsAlms::calculate_single_C_l(long l) const {
  long m;
  double tmp=0,tmp2;

  for (m=-l;m<=l;m++) {  tmp2=abs(get(l,m)); tmp +=tmp2*tmp2; }
  tmp = tmp/(2*(double)l+1);
  return tmp;
}
// ****************************************************************************************************
// the same but you can specify here the power comes from
// what: 0 - all 1 - real part 2 - complex part, 3 - m>=0 part
double mscsAlms::calculate_single_C_l(long l, long what) const {
  long int m;
  double tmp=0,tmp2;
  double num=0;

  if (what == 0) { // calculating whole multipole
    tmp = calculate_single_C_l(l); return tmp; }
  if (what == 1) { // calculating the real part power only
    for (m=-l;m<=l;m++) {  tmp+= pow(get(alm2num(l,m)).real(),2); } num=(2.0*double(l)+1.0); }
  if (what == 2) { // calculating the imaginary part power only
    for (m=-l;m<=l;m++) {  tmp+= pow(get(alm2num(l,m)).imag(),2); } num=(2.0*double(l)+1.0); }
  if (what == 3) { // calculating the m>=0 part power only
    for (m=0;m<=l;m++) {  tmp2=abs(get(l,m)); tmp +=tmp2*tmp2; } num=(double(l)+1.0); }
  if (what == 4) { // calculating the m<0 part power only
    for (m=-l;m<0;m++) {  tmp2=abs(get(l,m)); tmp +=tmp2*tmp2; } num=double(l); }

  tmp/=num;

  return tmp;
}
// ****************************************************************************************************
// the same but you can also specify the convention for the alms ordering -- default is as in ccSHT
// ordering: 0 - default  1 - ordering with increasing ls
double mscsAlms::calculate_single_C_l(long l, long what, long ordering) const {
  int m;
  double tmp=0,tmp2;

  if (ordering == 1) {
    if (what == 0) { // calculating whole multipole
      for (m=-l;m<=l;m++) {  tmp2=abs(get(alm2num_norm(l,m))); tmp+= tmp2*tmp2; }}
    if (what == 1) { // calculating the real part power only
      for (m=-l;m<=l;m++) {  tmp += pow(get(alm2num_norm(l,m)).real(),2); }}
    if (what == 2) { // calculating the complex  part power only
      for (m=-l;m<=l;m++) {  tmp +=pow(get(alm2num_norm(l,m)).imag(),2); }}

    tmp = tmp/(2*(double)l+1);
  } else tmp = calculate_single_C_l(l,what);

  return tmp;
}

// ****************************************************************************************************
//************************************************************************
long mscsAlms::savebinAlms(string alms_file, string how) const { // saves alms to a file
  long res=0;
  msgs->say("Saving alms to binary file: "+alms_file,High);
  if (how == "RI") { res=alms().save(alms_file,true,"complex"); }
  if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
  msgs->say("Success",Low);
  return res;
}
// void mscsAlms::savebinAlms(string alms_file, int how, long l) { // saves alms to a file
//   loadsave_manager("save","bin","F",how,alms_file,l);
// }
//************************************************************************
long mscsAlms::loadbinAlms(string  alms_file, string how) {
  long res=0;
  cpedsList<mscsAlm> al;
  msgs->say("Loading alms from binary file: "+alms_file,High);
  if (how == "RI") { res=al.load(alms_file,true,"complex"); }
  if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
  a()=al;
  msgs->say("Success",Low);

  set_alms_num(alms().size());
  return res;
}
// void mscsAlms::loadbinAlms(string  alms_file, int how, long l) { // saves alms to a file
//   loadsave_manager("load","bin","F",how,alms_file,l);
// }
//************************************************************************
long mscsAlms::savetxtAlms(string  alms_file, string how) const { // saves alms to a file
  long res=0;
  msgs->say("Saving alms to text file: "+alms_file,High);
  if (how == "RI") { res=alms().save(alms_file,false,"complex"); }
  if (how == "lmRI") {
	  FILE* f=fopen(alms_file.c_str(),"w");
	  if (f!=NULL) {
		  for (int l = 0; l < lmax(); ++l) {
			for (int m = -l; m <= l; ++m) {
				mscsAlm alm=get(l,m);
				fprintf(f,"%li %li %lf %lf\n",l,m,alm.R(),alm.I());
			}
		  }
		  fclose(f);
	  }
	  else { msgs->warning("Could not save to file",Medium); return -1; }
  }
  if (res!=0) { msgs->warning("Could not save to file",Medium); return res; }
  msgs->say("Success",Low);
  return res;
}
// void mscsAlms::savetxtAlms(string  alms_file, int how, long l) { // saves alms to a file
//   loadsave_manager("save","txt","F",how,alms_file,l);
// }
//************************************************************************
long mscsAlms::loadtxtAlms(string  alms_file, string how) {
  long res=0;
  cpedsList<mscsAlm> al;
  msgs->say("Loading alms from text file: "+alms_file,High);
  if (how == "RI") { res=al.load(alms_file,false,"mscsAlm");  }
  if (res!=0) { msgs->warning("Could not load from file",Medium); return res; }
  msgs->say("Success",Low);
  a()=al;
  set_alms_num(alms().size());
  return res;
}
// void mscsAlms::loadtxtAlms(string  alms_file, string how, long l) { // saves alms to a file
//   loadsave_manager("load","txt","F",how,alms_file,l);
// }
//************************************************************************
void mscsAlms::printtxtAlms(string what) const {// prints fourier alms to the screen
  int l,m;
  mscsAlm tmpalm;
  long llmax=lmax();
  msgs->say("printing alms:",High);
  for (l=0;l<=llmax;l++) {
    for (m=-l;m<=l;m++) {
      tmpalm = get(alm2num(l,m));
      printf("l = %i, m = %i, R = %lE, I = %lE, D = %lE, PHI = %lE [deg]\n",l,m,tmpalm.real(),tmpalm.imag(),tmpalm.abs(),PI180inv*arg(tmpalm.Z()));
    }
  }

}

// ****************************************************************************************************
void mscsAlms::antysymmetrize() {
  long l,m,num1,num2;
  double sqrt2=sqrt(2.0);
  long llmax=lmax();

  msgs->say("Antysymmetrizing alms",High);
  num1=0;
  for (l=0;l<=llmax;l++) { num2=num1; // remember the number when entering a new multipole
    for (m=-l;m<=l;m++) {
      num1=alm2num(l,m);
      if (m < 0) {
	num2=alm2num(l,-m);
	a(num1)=mscsAlm( pow(-1.0,(double)(-m)) * get(num2).real() ,  -pow(-1.0,(double)(-m)) * get(num2).imag() );      }
      if (m == 0) { a(num1)=mscsAlm(get(num1).real()*sqrt2, 0.0); } // correct the m=0 case for the I_lm=0 - i.e. sqrt(2) factor and update D and phase P
    }
  }

  // num1=0;
  // for (l=0;l<=llmax;l++) { num2=num1; // remember the number when entering a new multipole
  //   for (m=-l;m<=l;m++) {
  //     if (m < 0) {
  // 	num2=alm2num(l,m);
  // 	a(num1)=mscsAlm( pow(-1.0,(double)(-m)) * get(num2).real() ,  -pow(-1.0,(double)(-m)) * get(num2).imag() );
  // 	//alm[num1].P = cpeds_cart2sph(1,alm[num1].R,alm[num1].I,0.0); // fill the rest
  //     }
  //     if (m == 0) { a(num1)=mscsAlm(get(num1).real()*sqrt2, 0.0); // correct the m=0 case for the I_lm=0 - i.e. sqrt(2) factor and update D and phase P
  // 	//	alm[num1].D = fabs(alm[num1].R);  if (alm[num1].R >= 0) { alm[num1].P = 0; } else {alm[num1].P = PI;} // fill the rest
  //     }
  //     num1++;
  //   }
  // }

}
// ****************************************************************************************************
const mscsAlms& mscsAlms::calibrateOn(long ll, double cl) {
   long l,m;

  // if (((method == 5)  || (method == 6)) && (lmax >= lcalib) && (lcalib >= 0)) {
   msgs->say("calibrating alms to power spectrum at l="+msgs->toStr(ll)+", C_"+msgs->toStr(ll)+"="+msgs->toStr(cl),Medium);
   double Clgen=calculate_single_C_l(ll,0);

   double Nfactor = sqrt(cl/Clgen);
   long llmax=lmax();
   for (l=0;l<=llmax;l++) {
     for (m=-l;m<=l;m++) {
       a(l,m)*=Nfactor;
     }
   }

   return *this;
 }
// ****************************************************************************************************
const mscsAlms& mscsAlms::generate_gaussian_alms(long llmax, double fm, double fs, long method, bool half, cpedsRNG* rns) {
  long M;
  bool rnsWasNull;
  msgs->say("generating random gaussian alms with mean: "+msgs->toStr(fm)+" and variance: "+msgs->toStr(fs),High);
  clear();
  lmax(llmax);
  if (rns==NULL) { rns = new cpedsRNG("gaussian","double"); rns->setCentralLimitNumbers(100); rnsWasNull=true; } else { rnsWasNull=false; }
  rns->setMeanStd(fm,fs);


  long N=almsNum();
  if (half) M=N; else M=2*N;
  long j=0;

  cpedsList<double> cl;
  cl=rns->getRNs(M);

  if (half) {
    for (long l=0;l<=llmax;l++) {
      for (long m=0;m<=llmax;m++) {
	if (m==0) { a(l,m)=mscsAlm(cl[j],0.0); j++; }
	else { a(l,m)=mscsAlm(cl[j],cl[j+1]); j+=2; }
      }
    }
  }
  else {
    for (long i=0;i<N;i++) {
      a(i)=mscsAlm(cl[j],cl[j+1]);
      j+=2;
    }
  }

  if (rnsWasNull) { delete rns; rns=NULL; }
  return *this;
}
// ****************************************************************************************************
const mscsAlms& mscsAlms::generate_gaussian_alms(const mscsAngularPowerSpectrum& cl, long llmax, long method, bool half, bool exact, long calibrateon, cpedsRNG* rns) {
  generate_gaussian_alms(llmax,0.0,1.0,method,half,rns);
  {  cpedsList<double> clist=toDoubleList();
	  msgs->say("alms num: "+msgs->toStr(almsNum())+", mean: "+msgs->toStr(double(clist.mean()))+", st.dev: "+msgs->toStr(sqrt(double(clist.variance())))+", skweness: "+msgs->toStr(clist.skewness())+", kurtosis: "+msgs->toStr(clist.kurtosis()),High);
  }
  double TT,ld,Clgen,Clorg,Nfactor;
  long int l,m,k;
  // long int num,num1=0,num2=0,i,lcalib=10;
  // double *vec = NULL, sumDlm2, tmpd1, tmpd2,ld,Clgen=0,Clorg=0,TT=0,C10org,C10gen,Nfactor;
  // double *tab,*tab2;
  // long int istart=0,istart2=0;

  // set_alms_lmax(multipole);
  // lcalib = calibrateon;
  msgs->say("Generating gaussian alms vector up to: lmax = "+msgs->toStr(llmax),High);

  // //------ space allocation in case it's not allocated ----------
  // makekill_space_manager("make","F",-1);
  // //-------------------------------------------------------------

  // if ((method == 4) || (method == 5)) {
    //if (strcmp(whatisgaussian,"RI") == 0) {       tab = cpeds_random_gauss_numbers(fm,fs,2*alms_num,2);    } // this uses the system RNG
    // if (strcmp(whatisgaussian,"RI") == 0) {       tab = cpeds_random_gauss_numbers(fm,fs,2*alms_num,3,cpeds_seed_offset,fastsims);    } // this uses GSL RNG the 'DIEHARD' survivor
    // if (strcmp(whatisgaussian,"D") == 0) {
    //   tab = cpeds_random_gauss_numbers(fm,fs,alms_num,2);
    //   tab2 = cpeds_random_uniform_numbers(0,2*PI,alms_num);
    // }
  // }

  for (l=0;l<=llmax;l++) { // this loop generates the gaussian alms for each l individually but perhaps it's just as good to make one bie gaussian realization of all alms upto a given lmax and redistribute them acuratelly into the alm file hm ? perhaps it's the same thing.
    // printf("multipole number: %li of %li\r",l,lmax);
    // if ( (strcmp(whatkind,"unnormCl") == 0)  || (strcmp(whatkind,"normCl") == 0)) { // read the Cl value
    TT = cl.get_Cl(l); // this is the value of the original power in a given multipole
    if (TT < 0) TT = 0; // safety condition for some crazy Cls
    // fs = 1; // this does not concern the method 4 and 5
    // if (method == 6) fs = sqrt(TT);///C_l[0][1];
    //if (method == 2) { fs = sqrt(2/(2*(double)l+1)); } // this is the cosmic variance used
    //generate_gassuian_distribuant_function((lint)GAUSS_DISTR_TAB_SIZE,fs,fm);
    //fs = 1;
    // }
    // m=2*l+1; // number of m values
    //printf("fm = %lE, fs = %lE n = %i\n",fm,fs,m);

    // if (method == 1) {   vec = cpeds_random_gauss_numbers(fm,fs,m,1); } // using cpeds library for random number generation (iterative)
    // if ((method == 2) || (method == 6)){   vec = cpeds_random_gauss_numbers(fm,fs,m,2); } // using cpeds library for random number generation (by half division)
    // if (method == 3) {   vec = GSL_random_gauss_numbers(fm,fs,m);   } // using gsl   library for random number generation
    // if ((method == 4)||(method == 5)) {   vec = new double[m]; for (i=istart;i<istart+m;i++) { vec[i-istart] = tab[i];  } istart = istart+m; } // this should increase speed as compared to 2 only by eliminating CDF calculations for each multipole
    // sumDlm2=0; tmpd1=tmpd2=0;
    // if (half) { for (k=0;k<=l;k++) { sumDlm2 += calculate_single_C_l(l,3); tmpd1++;}
    // else { for (k=-l;k<=l;k++) { sumDlm2 += calculate_single_C_l(l); tmpd1++;}

// /* EXPERIMENTAL BEGIN */
//     // for (k=1;k<=l;k++) { sumDlm2 += vec[k+l]*vec[k+l]; tmpd1++;} sumDlm2+=vec[l]*vec[l]/2; tmpd1++;// this is only needed by "D" here. for "RI" we don't need to calculate this here !!
// /* EXPERIMENTAL END */

    ld=(double)l;
    //printf("****m= %i, sumSlm2 = %lf tmpd1 = %lf vec[l=0] = %lf\n",m,sumDlm2,tmpd1,vec[0]);

    // if (strcmp(whatkind,"noCl") == 0) { Clgen = 1; Clorg = 1; }
    // if (strcmp(whatkind,"unnormCl") == 0) {
    // if ((method == 5) || (method == 6)) { Clgen = 2; } else { Clgen = sumDlm2/(2*ld+1); } // !!!!!!!!!!!!!!! this is terribly important what goes here.
    if (exact) {
      if (half) { Clgen=calculate_single_C_l(l,3); } else { Clgen=calculate_single_C_l(l); }
    }
    else {
      if (half) Clgen=1.0; else Clgen=2.0;
    }

    Clorg = TT;   // reduction of the original power spectrum normalization. the common normalization with Clgen.
    // }
    // if (strcmp(whatkind,"normCl") == 0) {
    //   if ((method == 5) || (method == 6)) { Clgen = 2; } else { Clgen = sumDlm2/(2*ld+1); } //clgen is not used here for RI
    //   Clorg = TT/(ld*(ld+1)/(2*PI));   // reduction of the original power spectrum normalization (eg. from the cmbfast). the common normalization with Clgen. (but it should be generally done during loading procedure so it's probably not very usefull here
    //   // this reduction is now done during reading or saving so in the program generally the unnormalized power spectrum is kept
    // }
    Nfactor = sqrt(Clorg/Clgen);
    if (cpeds_isnan(Nfactor)) { printf("%li Clorg %lE Clgen %lE Nfact %lE\n",l,Clorg,Clgen,Nfactor); }
    //printf("fact = %lf\n",sqrt(Clorg/Clgen)); exit(0);
    // if (strcmp(whatisgaussian,"D") == 0) { // the modulus should yield the chisq distribution - not the gaussian one !!!
    //   for (k=-l;k<=l;k++)  { alm[num1].D =  vec[k+l] * Nfactor; num1++;} // setting modulus and its normalization of the spectrum to fit to the required power from alm.
    // }
    // if (strcmp(whatisgaussian,"RI") == 0) {
    for (m=-l;m<=l;m++) { a(l,m)*=Nfactor; } // this is the normalization on the spectrum to fit to the required power from alm.
    // }

/*     printf(" - R/D:: l=%li mean: %.3lE skewness: %.3lE kurtosis: %.3lE |",l, cpeds_mean_value(vec,m), cpeds_skewness(vec,m), cpeds_kurtosis(vec,m)); */
    // delete vec;

    //-----------------

    // if (method == 1) { // using cpeds library for random number generation
    //   if (strcmp(whatisgaussian,"D") == 0) { vec = cpeds_random_uniform_numbers(0,2*PI,m); }
    //   if (strcmp(whatisgaussian,"RI") == 0) { vec = cpeds_random_gauss_numbers(fm,fs,m,1);   }
    // }
    // if ((method == 2) || (method == 6)) { // using cpeds library for random number generation
    //   if (strcmp(whatisgaussian,"D") == 0) { vec = cpeds_random_uniform_numbers(0,2*PI,m); }
    //   if (strcmp(whatisgaussian,"RI") == 0) { vec = cpeds_random_gauss_numbers(fm,fs,m,2);   }
    // }
    // if (method == 3) { // using gsl  library for random number generation
    //   if (strcmp(whatisgaussian,"D") == 0)  { vec = GSL_random_uniform_numbers(0,2*PI,m); }
    //   if (strcmp(whatisgaussian,"RI") == 0) { vec = GSL_random_gauss_numbers(fm,fs,m);    }
    // }
    // if ((method == 4) || (method == 5)) {   vec = new double[m];
    //   if (strcmp(whatisgaussian,"D") == 0)  { for (i=istart2;i<istart2+m;i++) { vec[i-istart2] = tab2[i];  } istart2 += m; }
    //   if (strcmp(whatisgaussian,"RI") == 0) { for (i=istart;i<istart+m;i++) {  vec[i-istart] = tab[i];  } istart += m; }
    // }

    // if (strcmp(whatisgaussian,"RI") == 0) { // this is not needed for methods 5 and 6
    //   tmpd1=tmpd2=0;
/* EXPERIMENTAL BEGIN */
// /*       for (k=-l;k<=l;k++) { tmpd1 += vec[k+l]*vec[k+l]; } */
//       for (k=1;k<=l;k++) { tmpd1 += vec[k+l]*vec[k+l]; }  // here summation from 0 will be corrected during anti-symmetrization // this is the current one
/* EXPERIMENTAL END */
//       sumDlm2 += tmpd1;

//       ld=(double)l;  //printf("****m= %i, fact = %lf tmpd1 = %lf vec[l=0] = %lf\n",m,sumDlm2,tmpd1,vec[0]);
//     }
//     //printf("**2)**m= %i, sumSlm2 = %lf tmpd1 = %lf vec[l=0] = %lf\n",m,sumDlm2,tmpd1,vec[0]);
//     if (strcmp(whatkind,"noCl") == 0) { Clgen = 1; Clorg = 1; }
//     if (strcmp(whatkind,"unnormCl") == 0) {
//       if ((method == 5) || (method == 6)) { Clgen = 2; } else { Clgen = 2 * sumDlm2/(2*ld+1); } // here 2 comes from the fact that now I sum sumDlm2 only for m>=0; for m=0 only real part; but used to be without 2 when summed over all the range -l...l
//       Clorg = TT;   // reduction of the original power spectrum normalization. the common normalization with Clgen.
//     }
//     if (strcmp(whatkind,"normCl") == 0)  {
//       if ((method == 5) || (method == 6)) { Clgen = 2; } else { Clgen = 2 * sumDlm2/(2*ld+1); }
//       Clorg = TT/(ld*(ld+1)/(2*PI));   // reduction of the original power spectrum normalization (eg. from the cmbfast). the common normalization with Clgen.
//     }

//     Nfactor = sqrt(Clorg/Clgen);
//     if (method == 6) Nfactor = 1;

//     if (strcmp(whatisgaussian,"D") == 0) {
//       for (k=-l;k<=l;k++) {
// /* 	num = alm2num(l,k); alm[num].P = vec[k+l];  // set the phases */
// /* 	num = alm2num(l,k); alm[num].R = alm[num].D * cos(alm[num].P); alm[num].I = alm[num].D * sin(alm[num].P);  // fill the rest  */
// 	alm[num2].P = vec[k+l];  // set the phases
// 	if (alm[num2].D < 0) { alm[num2].D = -alm[num2].D; alm[num2].P = alm[num2].P + PI; if (alm[num2].P > twoPI) {alm[num2].P = alm[num2].P - twoPI; } }
// 	alm[num2].R = alm[num2].D * cos(alm[num2].P); alm[num2].I = alm[num2].D * sin(alm[num2].P); num2++;  // fill the rest
//       }
//     }
//     if (strcmp(whatisgaussian,"RI") == 0) {
//       for (k=-l;k<=l;k++) {
// /* 	num = alm2num(l,k); alm[num].I = vec[k+l] * sqrt(Clorg/Clgen); alm[num].R = alm[num].R * sqrt(Clorg/Clgen); // set the unreal part and correct the real part */
// /* 	alm[num].D = sqrt(pow(alm[num].R,2)+pow(alm[num].I,2)); alm[num].P = cpeds_cart2sph(1,alm[num].R,alm[num].I,0); // fill the rest */

//         alm[num2].I = vec[k+l] * Nfactor; alm[num2].R *= Nfactor; // set the unreal part and correct the real part (for the R it's done here so it all worked with all methods)
// 	alm[num2].D = sqrt(pow(alm[num2].R,2)+pow(alm[num2].I,2)); alm[num2].P = cpeds_cart2sph(1,alm[num2].R,alm[num2].I,0.0);
// 	//printf("Rlm=%lE Ilm=%lE, Dlm=%lE, philm=%.3lf |||| sqrt(Rlm^2+Ilm^2)=%lE  ||| C_l[][0] = %lE, C_l[][1] = %lE\n",alm[num2].R ,alm[num2].I, alm[num2].D,alm[num2].P,sqrt(alm[num2].R*alm[num2].R+alm[num2].I*alm[num2].I),C_l[l][0],C_l[l][1]);
// 	num2++; // fill the rest
//       }
//     }
// /*     printf(" - I/PHI:: mean: %.3lE skewness: %.3lE  kurtosis: %.3lE \n", cpeds_mean_value(vec,m), cpeds_skewness(vec,m), cpeds_kurtosis(vec,m)); */
//     delete vec;
  }


  // print some statistics
  cpedsList<double> clist=toDoubleList();
  // printf("%i %lE %lE \n",clist.size(),clist[clist.size()-1], clist.mean());
  // printf("dupa\n");
  // exit(0);
  msgs->say("alms num: "+msgs->toStr(almsNum())+", mean: "+msgs->toStr(double(clist.mean()))+", st.dev: "+msgs->toStr(sqrt(double(clist.variance())))+", skweness: "+msgs->toStr(clist.skewness())+", kurtosis: "+msgs->toStr(clist.kurtosis()),High);
//  msgs->say("alms num: "+msgs->toStr(almsNum())+", mean: "+msgs->toStr(double(clist.mean()))+", skweness: "+msgs->toStr(clist.skewness())+", kurtosis: "+msgs->toStr(clist.kurtosis()),High);
  // if ((method == 4) || (method == 5)) {
  //   if (strcmp(whatisgaussian,"RI") == 0) printf("\n\n all alms:: RI :: alms_num=%li mean: %lE   skewness: %lE  kurtosis: %lE \n",2*alms_num, cpeds_mean_value(tab,2*alms_num), cpeds_skewness(tab,2*alms_num), cpeds_kurtosis(tab,2*alms_num));
    // if (strcmp(whatisgaussian,"D") == 0) printf("\n\n all alms:: DP :: alms_num=%li mean: %lE   skewness: %lE  kurtosis: %lE \n",alms_num, cpeds_mean_value(tab,alms_num), cpeds_skewness(tab,alms_num), cpeds_kurtosis(tab,alms_num));

  //   delete tab;
  //   if (strcmp(whatisgaussian,"D") == 0) { delete tab2; }
  // }

  // make the alms antisymmetrical if the method is 5 or 6
  // if (method == 4 || method == 5  || method == 6) { // the expectancy value for the alms with the hermitian symmetry a_l-m = (-1)^m alm* is the same as for the general case -- C_l so this does not chance the power in the modelled Cl
  //   //if ((method == 5)  && (method == 6)) { // this is only to swith this part off
  //   printf("  -- antisymmetrizing the alms\n");
  //   num1=0;
  //   for (l=0;l<=lmax;l++) { num2=num1; // remember the number when entering a new multipole
  //     for (m=-l;m<=l;m++) {
  // 	//printf("a_%li%li.R=%lE  ,   a_%li%li.I=%lE    num=%li, num1=%li, num2 = %li\n",l,m,alm[num1].R,l,m,alm[num1].I,num,num1,num2);
  // 	if (m < 0) {
  // 	  k = 2*l+1; // number of ms in a given multipole
  // 	  //num = num2+k-1-(num1-num2); // number of the a_l+m
  // 	  num = 2*(l+num2)-num1; // number of the a_l+m - positive m counterpart // this is the simplifed form of the expression above
  // 	  alm[num1].R = pow(-1.0,(double)(-m)) * alm[num].R; alm[num1].I = -pow(-1.0,(double)(-m)) * alm[num].I;
  // 	  alm[num1].P = cpeds_cart2sph(1,alm[num1].R,alm[num1].I,0.0); // fill the rest
  // 	  //printf("-:a_%li%li.R=%lE  ,   a_%li%li.I=%lE    num=%li, num1=%li, num2 = %li\n",l,-m,alm[num].R,l,m,alm[num].I,num,num1,num2);
  // 	}
  // 	if (m == 0) { alm[num1].I = 0; if (method==5 || method==6) alm[num1].R *= sqrt(2.0);
  // 	  alm[num1].D = fabs(alm[num1].R);  if (alm[num1].R >= 0) { alm[num1].P = 0; } else {alm[num1].P = PI;} // fill the rest
  // 	} // correct the m=0 case for the I_lm=0 - i.e. sqrt(2) factor and update D and phase P
  //       //printf("+:a_%li%li.R=%lE  ,   a_%li%li.I=%lE    num=%li, num1=%li, num2 = %li\n\n",l,m,alm[num1].R,l,m,alm[num1].I,num,num1,num2);
  // 	num1++;
  //     }
  //     //printf("\n");
  //   }
  // }

  //printf("\n\n all alms:: alms_num=%li mean: %lE   skewness: %lE  kurtosis: %lE \n",2*alms_num, cpeds_mean_value(alm,2*alms_num), cpeds_skewness(alm,2*alms_num), cpeds_kurtosis(alm,2*alms_num));

  // THIS BIT IS NOT NEEDED ANYMORE SINCE THE CORRECT CALCULATION GIVES EX<(a_lm)>_l = C_l
  // here can be added some normalization of the alms according to the power spectrum but no obvious way of how to do that.
  // single multipole calibration might not be good idea because at least at low l's thre is big variance, so the remaining part of the Cl still might not be well fitted
  // anyway here is the single multipole calibration the magic l is l=10;
  if (calibrateon >0) {
    calibrateOn(calibrateon,cl.get_Cl(calibrateon));
  }

  // if (strcmp(lowls,"nol0-1") == 0) {   printf("  -- removing monopole and dipole\n"); for (i=0;i<4;i++) { alm[i].R = alm[i].I = alm[i].D = alm[i].P = 0; } }

  // converting the alms table to the format consistent with the convention in this program and the ccSHT library
  // a_lm* tmpalm = new a_lm[alms_num];  for (i=0;i<alms_num;i++) { tmpalm[i] = alm[i]; }

  // printf("  -- converting alms table... \n ");
  // num=0;
  // for (l=0;l<=lmax;l++) {
  //   for (m=-l;m<=l;m++) { alm[alm2num(l,m)] = tmpalm[num];  num++; }
  //   printf("  -- multipole number: %li of %li\r",l,lmax);
  // }
  // alms_loaded = 1;
  // printf("\n");
  // delete tmpalm;

  return *this;
}
// ****************************************************************************************************
const cpedsList<double> mscsAlms::toDoubleList() const {
  cpedsList<double> cl;
  long N=almsNum();
  long i;
  mscsAlm a;

  for (i=0;i<N;i++) { a=get(i); cl.append(a.real()); cl.append(a.imag()); }
  return cl;
}

// ****************************************************************************************************
void mscsAlms::set_alms_num(long alms_num) {
  long m;
  almsInfo.num=alms_num;
  if (alms_num!=0) {
//    num2alm(alms_num-1,&almsInfo.lmax,&m);
    almsInfo.lmax=sqrt(alms_num)-1;
    if (almsInfo.num!=alm2num(lmax(),lmax())+1) msgs->warning("there's an inconsistency between the number of alms and lmax - Check alms...",High);
  }
  else {
    almsInfo.lmax=-1;
  }
  printf("alms_num %li, lmax=%li\n",alms_num,almsInfo.lmax);
}
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************


  // METHODS FOR COMPUTING ALMS NUMBERS
//**/**********************************************************************

//************************************************************************




//************************************************************************
// whatkind: 1 - dont use the information about the power spectrum in C_l table (noCl)
// whatkind: 2 - DO use that information and generate the gaussian alms but the assumption is that the C_l are not normalized by factor l(l+1)/2pi (default in this program) (unnormCl
// whatkind: 3 - DO use that information and generate the gaussian alms but the assumption is that the C_l are normalized with factor l(l+1)/2pi (normCl)
// whatisgaussian: RI - for making real and imaginary parts drown from gaussian distribution or
// D  - for making modul from gaussian distribution and phase from uniform distribution


/*                        CAUTION !!!  */
/* for the moment the only valid method for the proper generation of gaussian field is through RI NOT DP. This method should follow different distribution function */
/*   The correct distribution for the modulus of the a_lm modes is the Rayleigh distribution: if U is uniform on (0,1] then R=sigma sqrt(-2 ln U) is Rayleigh distributed */

/* ************************************************************************ */
// void mscsAlms::generate_gaussian_alms(string whatisgaussian, string whatkind, string lowls, double fm, double fs, long multipole, int method, long int calibrateon, bool fastsims) { //generates the gaussian alms either with or without the underlying Cl
// }


// void mscsAlms::set_alms_lmax(long l) { //sets the maximal multipole number in the fourier space
//   if ((alms_loaded == 1) && (lmax != l)) { printf("!! WARNING !! - alms table with lmax=%li is already loaded\n",lmax); }
//   lmax = l;
//   alms_num = alm2num(l,l)+1;
//   printf("|%s> * setting lmax = %li alms_num: %li\n",object_name.c_str(),lmax, alms_num);
// }

// void mscsAlms::set_alms_lmax(long l,string how) { //sets the maximal multipole number in the fourier space
//   a_lm *newalms=NULL;
//   long new_num=0,i,ll,mm,alm_oldI,alm_newI,lmaxold=lmax;

//   if ((alms_loaded == 1) && (lmax != l)) {

//     if (how == "copy" || how == "replace") {
//       new_num = alm2num_norm(l,l)+1;
//       newalms = new a_lm[new_num];

//       for (ll=0;ll<=l;ll++) {
// 	for (mm=-ll;mm<=ll;mm++) {
// 	  lmax = lmaxold;
// 	  alm_oldI = alm2num(ll,mm);
// 	  lmax = l;
// 	  alm_newI = alm2num(ll,mm);
// 	  newalms[alm_newI] = alm[alm_oldI];
// 	}
//       }
//       lmax=l;  alms_num = alm2num(l,l)+1;
//       printf("|%s> * setting lmax = %li alms_num: %li\n",object_name.c_str(),lmax, alms_num);

//       if (how == "copy") {
// 	for (i=0;i<new_num;i++) alm[i] = newalms[i];
// 	printf("!! WARNING !! - overriding old alms, to fit the new structure of alms with lmax = %li \n!! WARNING !! - memory is still allocated for old size\n",l);
//       }

//       if (how == "replace") {
// 	makekill_space_manager("kill","F",1);     makekill_space_manager("make","F",1);
// 	for (i=0;i<new_num;i++) alm[i] = newalms[i];
// 	printf("!! NOTE !! -   copying old alms, to fit the new structure of alms with lmax = %li \n!! NOTE !! - memory is  allocated for the new size\n",l);
//       }
//       delete [] newalms;
//     }
// /*     exit(0); */
//   }
//   else { set_alms_lmax(l); }


// }

