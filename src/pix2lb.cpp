#include<stdlib.h>
#include<stdio.h>
#include<cpeds-math.h>
#include<cpeds-consts.h>

int main(int argc, char** argv) {
  if (argc!=4) { printf("pix2lb [ns] [ordering r-0,n-1] [pix]\n"); exit(0); }
  long ns=strtol(argv[1],NULL,10);
  long ord=strtol(argv[2],NULL,10);
  long pix=strtol(argv[3],NULL,10);
  double th,phi;
  cpeds_pix2ang_healpix(ns,pix,&th,&phi,ord);
  printf("%lE %lE\n",phi*PI180inv,(PIsnd-th)*PI180inv);
  return 0;
}
