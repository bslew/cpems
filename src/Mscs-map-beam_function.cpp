#include <math.h>
#include "Mscs-map-beam_function.h"


/* **************************************************************************************************** */
// CONSTRUCTORS
/* **************************************************************************************************** */

mscsBeamFunction::mscsBeamFunction() : mscsFunction() {

}
/* **************************************************************************************************** */
mscsBeamFunction::mscsBeamFunction(string _wf_name) : mscsFunction(_wf_name) {

}
/* **************************************************************************************************** */

mscsBeamFunction::mscsBeamFunction(string _wf_name, double fwhm, long nHWHM, long _points_per_HWHM) : mscsFunction(_wf_name) {
  long pointNum=2*nHWHM*_points_per_HWHM;
  setFWHM(fwhm);
  double range=double(nHWHM)*FWHM();
  double delta=double(range)/pointNum;
  double th=-range/2+delta;
  for (long i=0;i<pointNum;i++) {
    newPoint(th,double(0.0));
    th+=delta;
  }
}

/* **************************************************************************************************** */
// DESTRUCTOR
mscsBeamFunction::~mscsBeamFunction() {
}

/* **************************************************************************************************** */

void mscsBeamFunction::make_gaussian_beam() {
  make_gaussian_beam(FWHM());
}
/* **************************************************************************************************** */
void mscsBeamFunction::make_gaussian_beam(double fwhm) {
  long i;
  setFWHM(fwhm);
  double s = FWHM()/(2*sqrt(2*log(2.0)));
  double s2=s*s;
  double l;

  msgs->say("making gaussian window function kernel. Points:"+msgs->toStr(pointsCount()),Medium);
  long num=pointsCount();
  for (i=0;i<num;i++) {
    setf(i,exp(-getX(i)*getX(i)/(2*s2)));
  }
  printRanges();

}

