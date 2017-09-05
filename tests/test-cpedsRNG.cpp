#include "cpeds-rng.h"
#include "Mscs-function.h"

int main() {

  cpedsRNG rns;
  mscsFunction f;
  double T=1000;
  long N=1000;
  
  rns.setRNsType(rns.gaussian_invcdf);
  
  cpedsList<double> g=rns.getRNs(3000000);
  g.save("gaussian_invcvf-numbers.txt");
  
  exit(0);
  
  f.mkBoltzmann(0,CPEDS_kB*T,CPEDS_kB*T/N,T,0);
  rns.setRNsType(rns.from_array_invCDF);

  rns.setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());

  long n=1000;
  double *t=rns.getRN(n);
//
  cpeds_save_matrix(t,n,1,"cpedsRNS.txt");
  f.clearFunction();
  f.importFunction(rns.getCDFargs(),rns.getCDF(),n);
  f.derivative(true);
  f.save("PDF.txt");
  
  f.clearFunction();
  double T1=1000;
  f.mkBoltzmann(0,4*CPEDS_kB*T1,CPEDS_kB*T1/N,T1,0);
  double c=1.0/f.getMaxArg();
//  f.scaleX(c);
  f/=f.integrate();
  
  printf("expectation value: %lE\n",f.expectationValue()/(CPEDS_kB*T1));
  printf("integral: %lE\n",f.integrate());
  f.save("boltzmann1000.txt");

  f.clearFunction();
  double T2=100;
  f.mkBoltzmann(0,4*CPEDS_kB*T2,CPEDS_kB*T2/N,T2,0);
  f.scaleX(c);
  f.save("boltzmann100.txt");

  return 0;
}
