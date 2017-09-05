#include "Mscs-angular_power_spectrum.h"

// ****************************************************************************************************
void mscsAngularPowerSpectrum::initiate_mscsAngularPowerSpectrum(long Clsize) {
  // if (Clsize >= 0) { lmax = Clsize-1; } else { lmax = 0; Clsize=1; }
  // make_Cl(Clsize);
  // reset_Cl();
  // if (quietness!=0) printf("|%s> * initiated angular power spectrum of size: %li\n",object_name.c_str(),Clsize);
}

// ****************************************************************************************************
double mscsAngularPowerSpectrum::cosmic_variance(long l) {
  double ld=(double)l;
  return sqrt(2.0/(ld*(ld+1)));
}


// ****************************************************************************************************
const mscsAngularPowerSpectrum& mscsAngularPowerSpectrum::divide_llpotwoPI() {
  long l;
  long N=pointsCount();
  f(0)=double(0);
  for (l=1;l<N;l++) {  f(l)/=(get_l(l)*(get_l(l)+1)/twoPI); }
  return *this;
}
// ****************************************************************************************************
const mscsAngularPowerSpectrum& mscsAngularPowerSpectrum::multiply_llpotwoPI() {
  long l;
  long N=pointsCount();
  for (l=0;l<N;l++) {  f(l)*=(get_l(l)*(get_l(l)+1)/twoPI); }
  return *this;
}
// ****************************************************************************************************
cpedsStatusCodes mscsAngularPowerSpectrum::save(string filename) const {
  msgs->say("Saving power spectrum to file "+filename,High);
  msgs->messagingOff();
  cpedsStatusCodes s=mscsFunction::save(filename);
  msgs->messagingOn();
  return s;
}

// ****************************************************************************************************
cpedsStatusCodes mscsAngularPowerSpectrum::load(string filename, bool commentedFile) {
  msgs->say("Loading power spectrum from file "+filename,High);
  msgs->messagingOff();
  cpedsStatusCodes s=mscsFunction::load(filename, commentedFile,0,1);

  // check if there is a consistency between the position of an element in the array and the multipole number:
  // i.e. the 0'th multipole should sit on 0'th place in the queue and so on

  bool correction=false;
  for (int i = 0; i < getMaxArg(); ++i) {
	  if (getX(i)!=i) {
		  if (getX(i) > i) insertPoint(i,0,i);
		  else deletePoint(i);
		  correction=true;
	  }
  }

  msgs->messagingOn();
  if (correction) {
//	  msgs->warning("the structure of the loaded power spectrum has been corrected to match the index of the elements in the queue with the multipole number. This trigerred the printout.", High);
//	  print();
  }

  return s;
}
// ****************************************************************************************************
mscsAngularPowerSpectrum mscsAngularPowerSpectrum::get_cosmic_variance(double lmin, double lmax) {
  mscsAngularPowerSpectrum cv;
  long i,N=pointsCount();
  if (lmax==-1) lmax=N;
  for (i=lmin;i<=lmax;i++) { cv.newPoint(getX(i), cosmic_variance(i)); }
  return cv;
}
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************
