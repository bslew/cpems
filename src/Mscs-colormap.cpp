#include <complex>
#include "Mscs-colormap.h"
#include "cpeds-math.h"

/* **************************************************************************************************** */
mscsColormap::mscsColormap(string colors, long num, double minV, double maxV) {
  msgs=new cpedsMsgs("colormap");
  colorScheme=colors;
  valsCount=num;
  min=minV;
  max=maxV;

  sigma0=0.3;
  sigmaR0=0.60;
  sigmaG0=0.14;
  sigmaB0=0.35;

  gamma=0;
  sigma=sigma0;
  sigmaR=sigmaR0;
  sigmaG=sigmaG0;
  sigmaB=sigmaB0;

  r = new long[valsCount];
  g = new long[valsCount];
  b = new long[valsCount];
  thres = new double[valsCount];

  clear_variables();
  mk_thres();
  generateColors();



}

/* **************************************************************************************************** */
mscsColormap::~mscsColormap() {
  delete [] r;
  delete [] g;
  delete [] b;
  delete [] thres;
  if (msgs!=NULL) { delete msgs; msgs=NULL; }
}

/* **************************************************************************************************** */
void mscsColormap::print_colors() {
  msgs->say("minV: "+msgs->toStr(min)+" maxV: "+msgs->toStr(max)+" colors_num: "+msgs->toStr(valsCount)+" color scheme: "+colorScheme,High);
  msgs->say("thresholds no. th. r g b",Medium);
  for (long i=0;i<valsCount;i++) { 
    msgs->say(msgs->toStr(i)+" "+msgs->toStr(thres[i])+" "+msgs->toStr(r[i])+" "+msgs->toStr(g[i])+" "+msgs->toStr(b[i]),Medium);
  }
}

/* **************************************************************************************************** */
void mscsColormap::mk_thres() {
  long i;
  double d=(max-min)/(double)valsCount;

  thres[0]=min+d/2.0;
  for (i=1;i<valsCount;i++) {    thres[i]=thres[i-1]+d; }

}

/* **************************************************************************************************** */
void mscsColormap::clear_variables() {
  long i;
  for (i=0;i<valsCount;i++) { thres[i]=0.0; r[i]=g[i]=b[i]=0; }
}

/* **************************************************************************************************** */
void mscsColormap::generateColors() {
  long i;
  double z;

  if (colorScheme == "grayscale") { // linear gray scale
    for (i=0;i<valsCount;i++) {
      z=(double)(i+1)/(double)valsCount;
      r[i] = (long)round(z*255);      
      g[i] = (long)round(z*255); 
      b[i] = (long)round(z*255); 
    }
  }

  if (colorScheme == "color") { // linear
    for (i=0;i<valsCount;i++) {
      z=(double)(i+1)/(double)valsCount;
      r[i] = (long)round(z*255.0);      
      if (z<= 0.5) {g[i] = (long)round(2*z*255); } else { g[i] = (long)round((1.0-z)*255.0); }
      b[i] = (long)round((1.0-z)*255);
    }
  }

/*   if (color_scheme == "color-spots") { // color, linear with black and white extremum spots */
/*     l[num-2] = 1; r[num-2] = 0; g[num-2] = 0; b[num-2] = 0; */
/*     dz = 1/(float)numG; z=dz; */
/*     for (i=1;i<=numG;i++) { */
/*       l[i] = z; */
/*       r[i] = z;       */
/*       if (z<= 0.5) {g[i] = 2*z; } else { g[i] = g[i-1]-2*dz; } */
/*       b[i] = 1-z; */
/*       z=z+dz; */
/*     } */
/*     // setting white & black */
/*     l[num-1] = 0; r[num-1] = 1; g[num-1] = 1; b[num-1] = 1; */
/*     l[0] = 0; r[0] = 0; g[0] = 0; b[0] = 0; */
/*   } */

/*   if (color_scheme == "color-spoth") { // color, linear with black and white extremum spots */
/*     l[num-2] = 1; r[num-2] = 0; g[num-2] = 0; b[num-2] = 0; */
/*     dz = 1/(float)numG; z=dz; */
/*     for (i=1;i<=numG;i++) { */
/*       l[i] = z; */
/*       r[i] = z;       */
/*       if (z<= 0.5) {g[i] = 2*z; } else { g[i] = g[i-1]-2*dz; } */
/*       b[i] = 1-z; */
/*       z=z+dz; */
/*     } */
/*     // setting white */
/*     l[num-1] = 0; r[num-1] = 1; g[num-1] = 1; b[num-1] = 1; */
/*     //l[0] = 0; r[0] = 0; g[0] = 0; b[0] = 1; */
/*   } */

/*   if (color_scheme == "color-spotl") { // color, linear with black and white extremum spots */
/*     l[num-2] = 1; r[num-2] = 0; g[num-2] = 0; b[num-2] = 0; */
/*     dz = 1/(float)numG; z=dz; */
/*     for (i=1;i<num;i++) { */
/*       l[i] = z; */
/*       r[i] = z;       */
/*       if (z<= 0.5) {g[i] = 2*z; } else { g[i] = g[i-1]-2*dz; } */
/*       b[i] = 1-z; */
/*       z=z+dz; */
/*     } */
/*     // setting black */
/*     l[0] = 0; r[0] = 0; g[0] = 0; b[0] = 0; */
/*   } */


  if (colorScheme == "spectral") { // same RGB gaussian filters
    for (i=0;i<valsCount;i++) {
      z=(double)(i+1)/(double)valsCount;
      r[i] = (long)round( exp(-z*z/(2.0*sigma*sigma)) * 255.0 );      
      g[i] = (long)round( exp(-( (z-0.5)*(z-0.5) )/(2.0*sigma*sigma)) * 255.0 );      
      b[i] = (long)round( exp(-( (z-1.0)*(z-1.0) )/(2.0*sigma*sigma)) * 255.0 );      
    }
  }

  if (colorScheme == "spectralRGB") { // gaussian filters
    for (i=0;i<valsCount;i++) {
      z=(double)(i+1)/(double)valsCount;
      r[i] = (long)round( exp(-z*z/(2.0*sigmaR*sigmaR)) * 255.0 );      
      g[i] = (long)round( exp(-( (z-0.5)*(z-0.5) )/(2.0*sigmaG*sigmaG)) * 255.0 );      
      b[i] = (long)round( exp(-( (z-1.0)*(z-1.0) )/(2.0*sigmaB*sigmaB)) * 255.0 );      
    }
  }

  if (colorScheme == "spectralBGR") { // gaussian filters
    for (i=0;i<valsCount;i++) {
      z=(double)(i+1)/(double)valsCount;
      r[valsCount-i-1] = (long)round( exp(-z*z/(2.0*sigmaR*sigmaR)) * 255.0 );      
      g[valsCount-i-1] = (long)round( exp(-( (z-0.5)*(z-0.5) )/(2.0*sigmaG*sigmaG)) * 255.0 );      
      b[valsCount-i-1] = (long)round( exp(-( (z-1.0)*(z-1.0) )/(2.0*sigmaB*sigmaB)) * 255.0 );      
    }
  }

/*   if (color_scheme == "spectral") { //  */
/*     l[num-2] = 1; r[num-2] = 0; g[num-2] = 0; b[num-2] = 0; */
/*     dz = 1/(float)numG; z=dz; */
/*     for (i=1;i<=numG;i++) { */
/*       l[i] = z; */
/*       if (z < 0.375) { r[i] = 0; } if ((z >= 0.375) && (z < 0.625)) { r[i] = (z-0.375)/(0.625-0.375); } if ((z >= 0.625) && (z < 0.875)) { r[i] = 1; } if (z >=0.875) { r[i] = 1-(z-0.875)/(1-0.875); } // red */
/*       if (z < 0.125) { g[i] = 0; } if ((z >= 0.125) && (z < 0.375)) { g[i] = (z-0.125)/(0.375-0.125); } if ((z >= 0.375) && (z < 0.625)) { g[i] = 1; } if ((z >=0.625) && (z < 0.875)) { g[i] = 1-(z-0.875)/(1-0.875); } if (z >= 0.875) { g[i] = 0; }// green */
/*       if (z < 0.125) { b[i] = 0.5 + z/(0.125); } if ((z >= 0.125) && (z < 0.375)) { b[i] = 1; } if ((z >= 0.375) && (z < 0.625)) { b[i] = 1-(z-0.375)/(0.625-0.375); } if (z >=0.625) { b[i] = 0; } // blue */

/*       z=z+dz; */
/*     } */
    // setting white & black
    //l[num-1] = 0; r[num-1] = 1; g[num-1] = 1; b[num-1] = 1;
/*   } */



}

/* **************************************************************************************************** */
long mscsColormap::getSchemeID(string colors) {
  long schemeID;
  if (colors=="grayscale") schemeID=0;
  if (colors=="color") schemeID=1;
  if (colors=="spectral") schemeID=2;
  if (colors=="spectralRGB") schemeID=3;
  if (colors=="spectralBGR") schemeID=4;
  return schemeID;
}
/* **************************************************************************************************** */
string mscsColormap::getColorSchemeName(long id) {
  string colors;
  switch (id) {
  case 0: { colors="grayscale"; break; }
  case 1: { colors="color"; break; }
  case 2: { colors="spectral"; break; }
  case 3: { colors="spectralRGB"; break; }
  case 4: { colors="spectralBGR"; break; }
  default: { colors="unknown"; break; }
  }
  return colors;
}
/* **************************************************************************************************** */
string mscsColormap::cycleColorSchemes(string colors) {
  long schemeID=getSchemeID(colors);
  return getColorSchemeName((schemeID+1) % getColorSchemesNumber());
}
/* **************************************************************************************************** */
MscsColor mscsColormap::getC(double v) {
  MscsColor c;
  //  double * tmpd=export_thres();
  long Cidx;

  Cidx=cpeds_find_value(v,thres,valsCount,0,valsCount);
  //  delete [] tmpd;

  c.r=r[Cidx];
  c.g=g[Cidx];
  c.b=b[Cidx];
  
  return c;
}

/* **************************************************************************************************** */
double* mscsColormap::export_thres() {
  double * tmp = new double[colorNum()];
  long i;
  for (i=0;i<colorNum();i++) tmp[i]=thres[i];
  return tmp;
}
/* **************************************************************************************************** */
float* mscsColormap::export_thresf() {
  float * tmp = new float[colorNum()];
  long i;
  for (i=0;i<colorNum();i++) tmp[i]=(float)thres[i];
  return tmp;
}
/* **************************************************************************************************** */
float* mscsColormap::export_thresf_norm() {
  long N=colorNum();
  float * tmp = new float[N];
  long i;
  for (i=0;i<N;i++) tmp[i]=(thres[i]-thres[0])/(thres[N-1]-thres[0]);
  return tmp;
}

/* **************************************************************************************************** */
float* mscsColormap::export_rf() {
  float* R= new float[colorNum()];
  for (long i=0;i<colorNum();i++) { R[i]=float(r[i])/float(255); }
  return R;
}
/* **************************************************************************************************** */
float* mscsColormap::export_gf() {
  float* G= new float[colorNum()];
  for (long i=0;i<colorNum();i++) { G[i]=float(g[i])/float(255); }
  return G;
}
/* **************************************************************************************************** */
float* mscsColormap::export_bf() {
  float* B= new float[colorNum()];
  for (long i=0;i<colorNum();i++) { B[i]=float(b[i])/float(255); }
  return B;
}

/* **************************************************************************************************** */
void mscsColormap::setGamma(double g) {
  gamma=g;
}

/* **************************************************************************************************** */
void mscsColormap::setColorScheme(string c) {
  colorScheme=c;
  generateColors();
}

/* **************************************************************************************************** */
void mscsColormap::setSigma(double s) {
  sigma=s;
}

/* **************************************************************************************************** */
void mscsColormap::setSigma(double sr, double sg, double sb) {
  sigmaR=sr;
  sigmaG=sg;
  sigmaB=sb;
}

/* **************************************************************************************************** */
void mscsColormap::setSigma(MscsColorD c) {
  sigmaR=c.r;
  sigmaG=c.g;
  sigmaB=c.b;
}

/* **************************************************************************************************** */
MscsColorD mscsColormap::getSigmaRGB0() {
  MscsColorD c;
  c.r=sigmaR0;
  c.g=sigmaG0;
  c.b=sigmaB0;
  return c;
}

/* **************************************************************************************************** */
double mscsColormap::getSigma0() {
  return sigma0;
}



