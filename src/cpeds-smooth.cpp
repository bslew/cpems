#include "cpeds-smooth.h"


/* **************************************************************************************************** */
/* cpedsSmooth::cpedsProject(QList<cpedsDirection*> &dirs, string projection) { */
/*     _=dirs; */

/* } */

/* **************************************************************************************************** */
cpedsSmooth::cpedsSmooth() {}
/* **************************************************************************************************** */
cpedsSmooth::cpedsSmooth(cpedsSmooth& p) {
  *this=p;
}
/* **************************************************************************************************** */
cpedsSmooth::cpedsSmooth(const cpedsPointSet3D& points, const cpedsPoint3D& borderSize, const cpedsPoint3D& resolution, double fwhmX, double fwhmY, double alpha) {
// cpedsSmooth::cpedsSmooth(const cpedsPointSet3D &points) {
  clear();
  append(points); // copy the points to be smoothed
  initiate_smoothing_parameters();
  setBorder(borderSize);
  setRes(resolution);
  _fwhmX=fwhmX/res().x();
  _fwhmY=fwhmY/res().y();
  setAlpha(alpha);
}
/* **************************************************************************************************** */
cpedsSmooth::cpedsSmooth(const cpedsPointSet2D &points,double alpha) {
  clear();
  append(cpedsPointSet3D(points)); // copy the points to be smoothed
  setAlpha(alpha);
  initiate_smoothing_parameters();
}
/* **************************************************************************************************** */
cpedsSmooth::cpedsSmooth(const cpedsDirectionSet& d, const cpedsDirection& n, const cpedsDirection& borderSize, const cpedsDirection& resolution, double fwhmLon, double fwhmLat, double alpha) {
  clear();
  d.print();
  append(cpedsPointSet3D( cpedsProject(d).projectOnPlane(n) ));
  setAlpha(alpha);
  resolution.print_direction("res");
  borderSize.print_direction("borderSize");
  cpedsDirection(fwhmLon,fwhmLat).print_direction("smoothing");
  set_smoothing_parameters(resolution,cpedsDirection(fwhmLon,fwhmLat),borderSize);
  projectionPlane()=n;
}

/* **************************************************************************************************** */
void cpedsSmooth::initiate_smoothing_parameters() {

  // // set automatically the smoothing parameters
  //  // set the X and Y direction resolution of the field: this is based on the rms of the points distribution with resolution parameter alpha
   // setAlpha(.01); // this means that the resolution will be 10 times finer than the rms of the input field points
   res().setXY(alpha()*dispersionX(), alpha()*dispersionY());

  // // set the range parameters for the convolution
  //  //
  //  setRange(1,1); // this means that the size of the 
  //  _dx=beta()*_Dx;
  //  _dy=beta()*_Dy;   


   //setFieldRanges();  
}

/* **************************************************************************************************** */
void cpedsSmooth::set_smoothing_parameters(const cpedsDirection& resolution, const cpedsDirection& smoothingLength, const cpedsDirection& fieldBorder) {
  cpedsDirectionSet ds;
  cpedsPointSet3D ps;

  // set resolution parameters on the projective plane
  ds.append(cpedsDirection(-resolution.lon()/2.0,0.0));
  ds.append(cpedsDirection( resolution.lon()/2.0,0.0));
  ds.append(cpedsDirection(0.0,-resolution.lat()/2.0));
  ds.append(cpedsDirection(0.0, resolution.lat()/2.0));
  ps = cpedsProject(ds).projectOnPlane(cpedsDirection(0,0));
  res().setXY( fabs(ps[1].x()-ps[0].x()),  fabs(ps[3].y()-ps[2].y()) );  
  ds.clear();
  
  res().print_point("res");

  // set the smoothing length parameter on the projective plane
  ds.append(cpedsDirection(-smoothingLength.lon()/2.0,0.0));
  ds.append(cpedsDirection( smoothingLength.lon()/2.0,0.0));
  ds.append(cpedsDirection(0.0,-smoothingLength.lat()/2.0));
  ds.append(cpedsDirection(0.0, smoothingLength.lat()/2.0));
  ps = cpedsProject(ds).projectOnPlane(cpedsDirection(0,0));
  setSmoothingLength( fabs(ps[1].x()-ps[0].x()),  fabs(ps[3].y()-ps[2].y()) );  
  ds.clear();

  // set the field border if requested on the projective plane
  if (fieldBorder!=cpedsDirection(0,0)) {
    ds.append(cpedsDirection(-fieldBorder.lon()/2.0,0.0));
    ds.append(cpedsDirection( fieldBorder.lon()/2.0,0.0));
    ds.append(cpedsDirection(0.0,-fieldBorder.lat()/2.0));
    ds.append(cpedsDirection(0.0, fieldBorder.lat()/2.0));
    ps = cpedsProject(ds).projectOnPlane(cpedsDirection(0,0));
    border().setXY( fabs(ps[1].x()-ps[0].x()),  fabs(ps[3].y()-ps[2].y()) );  
    ds.clear();
  }
  else setBorder( cpedsPoint3D(0,0,0) );

  _smoothingLength.print_point("smoothing");
  border().print_point("border");
  res().print_point("res");
  

}
/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
// void cpedsSmooth::setFieldRanges() {
//   getRanges(_minX,_maxX,_minY,_maxY,_Zmin,_Zmax); // get the field ranges
//   // _minX-=border().x();
//   // _maxX+=border().x();
//   // _minY-=border().y();
//   // _maxY+=border().y();
// }

// ****************************************************************************************************

const matrix<double> cpedsSmooth::make_field(double resx, double resy) {
  double xmin, xmax, ymin, ymax, zmin, zmax; 
  if (resx!=0) res().setX(resx);
  if (resy!=0) res().setY(resy);

  getRanges(xmin, xmax, ymin, ymax, zmin, zmax);
  _Lx=xmax-xmin + 2*border().x();
  _Ly=ymax-ymin + 2*border().y();
  cpedsPoint3D borderPixNum=(border()/res()).roundToInteger();
  //borderPixNum.print_point();

  long fieldSizeX=long(round( (xmax-xmin)/res().x()) + 1);
  long fieldSizeY=long(round( (ymax-ymin)/res().y()) + 1);

  // printf("$$$$$$ %li %li\n",fieldSizeX,fieldSizeY);
  _m=exportAsField(fieldSizeX,fieldSizeY,_mask);
  _m=cpeds_pad_matrix(_m,long(borderPixNum.y()),long(borderPixNum.y()),long(borderPixNum.x()), long(borderPixNum.x()));
  // printf("rows %li cols %li\n ",_m.RowNo(),  _m.ColNo());

  return _m;
}
// ****************************************************************************************************
const matrix<double> cpedsSmooth::make_field(long sizeX, long sizeY) {
  double xmin, xmax, ymin, ymax, zmin, zmax; 
  getRanges(xmin, xmax, ymin, ymax, zmin, zmax);
  res().setX((xmax-xmin)/sizeX);
  res().setY((ymax-ymin)/sizeY);
  return make_field(res().x(),res().y());
}


// ****************************************************************************************************
void cpedsSmooth::field_fftY(fftw_complex* t, long vecNum, long vecSize, bool fwd) {
  fftw_plan p;
  long i,j;

  fftw_complex* t2=(fftw_complex*)fftw_malloc(vecSize*sizeof(fftw_complex));  
  int dir; if (fwd) dir=FFTW_FORWARD; else dir=FFTW_BACKWARD;

  p=fftw_plan_dft_1d(vecSize, t2,t2, dir, FFTW_ESTIMATE);
  for (i=0;i<vecNum;i++) {
    for (j=0;j<vecSize;j++) { t2[j][0]=t[i*vecSize+j][0]; t2[j][1]=t[i*vecSize+j][1];}
    fftw_execute(p);
    for (j=0;j<vecSize;j++) { t[i*vecSize+j][0]=t2[j][0]; t[i*vecSize+j][1]=t2[j][1]; }
  }
  fftw_destroy_plan(p);
  fftw_free(t2);

}
// ****************************************************************************************************
void cpedsSmooth::field_fftX(fftw_complex* t, long vecNum, long vecSize, bool fwd) {
  fftw_plan p;
  long i,j;

  fftw_complex* t2=(fftw_complex*)fftw_malloc(vecNum*sizeof(fftw_complex)); 
  int dir; if (fwd) dir=FFTW_FORWARD; else dir=FFTW_BACKWARD;

  p=fftw_plan_dft_1d(vecNum, t2,t2, dir, FFTW_ESTIMATE); 
  for (j=0;j<vecSize;j++) { 
    for (i=0;i<vecNum;i++) { t2[i][0]=t[i*vecSize+j][0]; t2[i][1]=t[i*vecSize+j][1]; }
    fftw_execute(p);
    for (i=0;i<vecNum;i++) { t[i*vecSize+j][0]=t2[i][0]; t[i*vecSize+j][1]=t2[i][1]; }
  }
  fftw_destroy_plan(p);
  fftw_free(t2);
}
// ****************************************************************************************************
void cpedsSmooth::field_fft(fftw_complex* t, long vecNum, long vecSize, bool fwd) {
// void cpedsSmooth::field_fft(fftw_complex* t, long vecNum, long vecSize) {
  fftw_plan p;

  int dir; if (fwd) dir=FFTW_FORWARD; else dir=FFTW_BACKWARD;

  // fftw_complex* t2= new fftw_complex[vecSize];
  // printf("***first row fft is:\n");
  // for (long i=0;i<vecSize;i++) { t2[i][0]=t[i][0]; t2[i][1]=t[i][1]; }
  // p=fftw_plan_dft_1d(vecSize, t2, t2, FFTW_FORWARD, FFTW_ESTIMATE);
  // // for (long i=0;i<vecSize;i++) { printf("before %lE, %lE  ",t2[i][0],t2[i][1]); }
  // // printf("\n***after\n");
  // fftw_execute(p);
  // for (long i=0;i<vecSize;i++) { printf("%lf, %lf  ",t2[i][0],t2[i][1]); }
  // // printf("\n***\n");
  // printf("\n***second row fft is:\n");
  // for (long i=0;i<vecSize;i++) { t2[i][0]=t[vecSize+i][0]; t2[i][1]=t[vecSize+i][1]; }
  // // for (long i=0;i<vecSize;i++) { printf("before %lE, %lE  ",t2[i][0],t2[i][1]); }
  // fftw_execute(p);
  // // printf("\n***after\n");
  // for (long i=0;i<vecSize;i++) { printf("%lf, %lf  ",t2[i][0],t2[i][1]); }
  // printf("\n***third row fft is:\n");
  // for (long i=0;i<vecSize;i++) { t2[i][0]=t[2*vecSize+i][0]; t2[i][1]=t[2*vecSize+i][1]; }
  // fftw_execute(p);
  // for (long i=0;i<vecSize;i++) { printf("%lf, %lf  ",t2[i][0],t2[i][1]); }
  // printf("\n***\n");
  // fftw_destroy_plan(p);

  // fftw_plan p=fftw_plan_dft_r2c_2d(vecNum, vecSize, t, (fftw_complex*)t, FFTW_ESTIMATE);
  // fftw_plan p=fftw_plan_dft_2d(vecNum, vecSize, (fftw_complex*)t, (fftw_complex*)t, FFTW_FORWARD, FFTW_ESTIMATE);

  p=fftw_plan_dft_2d(vecNum, vecSize, t,t, dir, FFTW_ESTIMATE);
  // p=fftw_plan_dft(vecNum, vecSize, (fftw_complex*)t,(fftw_complex*)t, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);

  // delete [] t2;
}

// ****************************************************************************************************
// void cpedsSmooth::field_inv_fft(fftw_complex* t, long vecNum, long vecSize) {
//   // fftw_plan p=fftw_plan_dft_c2r_2d(vecNum, vecSize, (fftw_complex*)t, t, FFTW_ESTIMATE);
//   fftw_plan p=fftw_plan_dft_2d(vecNum, vecSize, (fftw_complex*)t, (fftw_complex*)t, FFTW_BACKWARD, FFTW_ESTIMATE);
//   fftw_execute(p);
//   fftw_destroy_plan(p);
// }

// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGaussX(double fwhmX, matrix<double>* m) {
  _fwhmX=fwhmX/res().x();
  return smoothGaussX(long(_fwhmX),m);
}
// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGaussY(double fwhmY, matrix<double>* m) {
  _fwhmY=fwhmY/res().y();
  return smoothGaussY(long(_fwhmY),m);
}
// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGauss(double fwhmX, double fwhmY, matrix<double>* m) {
  if (fwhmX!=0) { _fwhmX=fwhmX/getRes().x(); }
  if (fwhmY!=0) { _fwhmY=fwhmY/getRes().y(); }
  return smoothGauss(long(_fwhmX),long(_fwhmY),m);
}
// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGauss(long nfwhmX, long nfwhmY, matrix<double>* m) {

  // OLD DEBUG STUFF - BEGIN
  // long paddingSize;
  // long i;

  // make_field();
  // // printf("rows %li cols %li\n ",_m.RowNo(),  _m.ColNo());
  // // add padding at the end of each row
  // if (_m.ColNo()%2==0) paddingSize=2; // if even number of cols we need two padded columns 
  // else paddingSize=1; // otherwise just one

  // // printf("rows %li cols %li  padding: %li\n",_m.RowNo(),  _m.ColNo(),paddingSize);
  // cpeds_print_matrix(_m);
  // cpeds_matrix_save(_m,"matrix.dat", "");

  // //printf("rows %li cols %li before padding in smoothGauss\n",_m.RowNo(),  _m.ColNo());

  // //  _m=cpeds_pad_matrix(_m,0,0,0,paddingSize);

  // //printf("rows %li cols %li\n after padding in smoothGauss",_m.RowNo(),  _m.ColNo());
  // // transform to linear rowMajor array
  // long size;
  // double *t=cpeds_matrix2array(_m,size,true);
  
  // // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // // cpeds_print_matrix(_m);
  // // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // // do fftw

  // fftw_complex* t2 = (fftw_complex*) fftw_malloc(size * sizeof(fftw_complex));
  // { for (long i=0;i<size;i++) { t2[i][0]=t[i]; t2[i][1]=0; } }
  // // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_uniform_number(0,1); t2[i][1]=cpeds_random_uniform_number(0,1); } }
  // // { for (long i=0;i<size;i++) { t2[i][0]=double(i); t2[i][1]=0; } }
  // // printf("initial random array: real part (imaginary is zero)\n");
  // // { for (long i=0;i<size;i++) {printf("%lf ",t2[i][0]); } printf("\n"); }
  

  // // field_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // // field_fft(t2, _m.RowNo(),  _m.ColNo());
  // // field_fft((double*)t2, _m.RowNo(),  _m.ColNo());

  // printf("matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());
  // field_fftY(t2, _m.RowNo(),  _m.ColNo());
  
  // //tmp block begin
  // // double* t2=new double[size/2];			 //
  // // { for (long i=0;i<size;i++) {printf("%lf, %lf  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // // exit(0);
  // // long aa=_m.ColNo();				 //
  // // cpeds_matrix_save(cpeds_array2matrix(t,size, aa,true),"matrixfft.mat","");	 //
  // // exit(0);
  // // { for (long i=0;i<size;i++) {t2[i]=t[2*i]; } }	 //
  // // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // // cpeds_matrix_save(_m,"matrixfft-real.mat","");	 //
  // // {  for (long i=0;i<size;i++) {t2[i]=t[2*i+1]; } }	 //
  // // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // // cpeds_matrix_save(_m,"matrixfft-imaginary.mat",""); //
  // // tmp block end
  
  // // do smoothing
  // field_convolve_Gaussian_kernel(t2,_m.RowNo(),  _m.ColNo());

  // // make the matrix symmetrical
  // make_field_symmetrical(t2,_m.RowNo(),  _m.ColNo());


  // // do inverse fftw
  // // field_inv_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // // field_inv_fft((double*)t2, _m.RowNo(),  _m.ColNo());
  // field_fftY(t2, _m.RowNo(),  _m.ColNo(),false);
  // // field_inv_fft(t, _m.RowNo(),  _m.ColNo());
  // // printf("po inv fft\n");
  // // { for (long i=0;i<size;i++) {printf("%lE,%lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // // exit(0);
  // // convert to matrix
  // for (i=0;i<size;i++) { t[i]=t2[i][0]; }
  // _m=cpeds_array2matrix(t,size, _m.ColNo(),true);

  // // remove the padding
  // //printf("deleting columns from %li to %li\n",_m.ColNo()-paddingSize,_m.ColNo()-1);

  // //_m.delCol(_m.ColNo()-paddingSize,_m.ColNo()-1);

  // // normalize after fft-ing
  // _m/=_m.RowNo() * _m.ColNo();

  // delete [] t;  
  // delete [] t2;  
  // return _m;

  // OLD DEBUG STUFF - END

  // GOOD DEBUG STUFF - BEGIN
  
  // check the smoothing params
  
  if (_fwhmX==0 && _fwhmY!=0) return smoothGaussY(nfwhmY,m);
  if (_fwhmX!=0 && _fwhmY==0) return smoothGaussX(nfwhmX,m);
  if (_fwhmX==0 && _fwhmY==0) {   
    // associate the smoothed matrix with the input point set
    associateMatrix2Field();
    // update the point set  
    _pssm=cpedsPointSet3D().field2set(_m,_referencePointRow,_referencePointCol,at(0),res());
    return _m;
  }

  //////////////////////////////
  // // prepare for smoothing //
  //////////////////////////////

  long i;
  if (m==NULL) make_field(); else _m=(*m);

  print();
#ifdef DEBUG_CPEDS
//  cpeds_print_matrix(_m);
  cpeds_matrix_save(_m,"cpeds-smooth.mat.in");
#endif

  long size;
  double *t=cpeds_matrix2array(_m,size,false); // this is colMajor to be consistent with the convention that the rows are the first dimention (X - although it's vertical) and columns are the last dimention (Y - although it's horizontal)

  fftw_complex* t2 = (fftw_complex*) fftw_malloc(size * sizeof(fftw_complex));
  for (i=0;i<size;i++) { t2[i][0]=t[i]; t2[i][1]=0; }
  // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_gauss_number(0,1,100,2); t2[i][1]=0; } }

  // printf("matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());
  // do fft in one direction
  field_fftX(t2,  _m.ColNo(), _m.RowNo());

  // do smoothing
  if (nfwhmX!=0) { _fwhmX=double(nfwhmX); }
  field_convolve_Gaussian_kernel(t2,  _m.ColNo(),_m.RowNo(),_fwhmX,0);

  // make the matrix symmetrical
  make_field_symmetricalCol(t2,  _m.ColNo(),_m.RowNo());

  // do inverse fftw
  field_fftX(t2,  _m.ColNo(), _m.RowNo(),false);


  // do fft in the other direction
  field_fftY(t2,  _m.ColNo(), _m.RowNo());

  // do smoothing
  if (nfwhmY!=0) { _fwhmY=double(nfwhmY); }
  field_convolve_Gaussian_kernel(t2, _m.ColNo(), _m.RowNo(),0,_fwhmY);

  // make the matrix symmetrical
  make_field_symmetricalRow(t2, _m.ColNo(), _m.RowNo());

  field_fftY(t2,  _m.ColNo(), _m.RowNo(),false);

  // convert to matrix
  for (i=0;i<size;i++) { t[i]=t2[i][0]; }
  _m=cpeds_array2matrix(t,size, _m.RowNo(),false);


  // normalize after fft-ing
  _m/=_m.RowNo() * _m.ColNo();

  delete [] t;  
  delete [] t2;  


  // associate the smoothed matrix with the input point set
  associateMatrix2Field();
  // update the point set  
  _pssm=cpedsPointSet3D().field2set(_m,_referencePointRow,_referencePointCol,at(_referencePointNumber),res());
  

  return _m;
  // GOOD DEBUG STUFF - END

  // EXPERIMENTAL  STUFF - BEGIN

  // long i;
  // make_field();
  // cpeds_matrix_save(_m,"matrix.dat", "");

  // long size;
  // double *t=cpeds_matrix2array(_m,size,false); // this is colMajor to be consistent with the convention that the rows are the first dimention (X - although it's vertical) and columns are the last dimention (Y - although it's horizontal)

  // fftw_complex* t2 = (fftw_complex*) fftw_malloc(size * sizeof(fftw_complex));
  // for (i=0;i<size;i++) { t2[i][0]=t[i]; t2[i][1]=0; }
  // // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_gauss_number(0,1,100,2); t2[i][1]=0; } }

  // // printf("matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());
  // // // do fft in one direction
  // // field_fftX(t2,  _m.ColNo(), _m.RowNo());
  // // // do fft in the other direction
  // // field_fftY(t2,  _m.ColNo(), _m.RowNo());

  // // do fft in both directions
  // field_fft(t2,  _m.ColNo(), _m.RowNo());

  // // // do smoothing
  // // if (fwhmX!=0) { _fwhmX=fwhmX; }
  // // _fwhmY=0;
  // // field_convolve_Gaussian_kernel(t2,  _m.ColNo(),_m.RowNo());
  // // // do smoothing
  // // if (fwhmY!=0) { _fwhmY=fwhmY; }
  // // _fwhmX=0;
  // // field_convolve_Gaussian_kernel(t2, _m.ColNo(), _m.RowNo());
  // // do smoothing
  // if (nfwhmX!=0) { _fwhmX=double(nfwhmX); }
  // if (nfwhmY!=0) { _fwhmY=double(nfwhmY); }
  // // _fwhmX=0;
  // field_convolve_Gaussian_kernel(t2, _m.ColNo(), _m.RowNo());

  // // make the matrix symmetrical
  // make_field_symmetricalCol(t2,  _m.ColNo(),_m.RowNo());
  // // make the matrix symmetrical
  // make_field_symmetricalRow(t2, _m.ColNo(), _m.RowNo());
  // // // make the matrix symmetrical
  // // make_field_symmetrical(t2, _m.ColNo(), _m.RowNo());

  // // do inverse fftw
  // // field_fftY(t2,  _m.ColNo(), _m.RowNo(),false);
  // // field_fftX(t2,  _m.ColNo(), _m.RowNo(),false);
  // field_fft(t2,  _m.ColNo(), _m.RowNo(),false);

  // // convert to matrix
  // for (i=0;i<size;i++) { t[i]=t2[i][0]; }
  // _m=cpeds_array2matrix(t,size, _m.RowNo(),false);


  // // normalize after fft-ing
  // _m/=_m.RowNo() * _m.ColNo();

  // delete [] t;  
  // delete [] t2;  
  // return _m;
  // EXPERIMENTAL  STUFF - END
}

// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGaussX(long nfwhmX, matrix<double>* m) {
  long i;
  if (m==NULL) make_field(); else _m=(*m);
  // printf("rows %li cols %li\n ",_m.RowNo(),  _m.ColNo());
  // add padding at the end of each row
  // if (_m.ColNo()%2==0) paddingSize=2; // if even number of cols we need two padded columns 
  // else paddingSize=1; // otherwise just one

  // printf("rows %li cols %li  padding: %li\n",_m.RowNo(),  _m.ColNo(),paddingSize);
  // cpeds_print_matrix(_m);
  // cpeds_matrix_save(_m,"matrix.dat", "");

  //printf("rows %li cols %li before padding in smoothGauss\n",_m.RowNo(),  _m.ColNo());

  //  _m=cpeds_pad_matrix(_m,0,0,0,paddingSize);

  //printf("rows %li cols %li\n after padding in smoothGauss",_m.RowNo(),  _m.ColNo());
  // transform to linear rowMajor array
  long size;
  double *t=cpeds_matrix2array(_m,size,false); // this is colMajor to be consistent with the convention that the rows are the first dimention (X - although it's vertical) and columns are the last dimention (Y - although it's horizontal)
  
  // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // cpeds_print_matrix(_m);
  // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // do fftw

  fftw_complex* t2 = (fftw_complex*) fftw_malloc(size * sizeof(fftw_complex));
  for (i=0;i<size;i++) { t2[i][0]=t[i]; t2[i][1]=0; }
  // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_uniform_number(0,1); t2[i][1]=cpeds_random_uniform_number(0,1); } }
  // { for (long i=0;i<size;i++) { t2[i][0]=double(i); t2[i][1]=0; } }
  // printf("initial random array: real part (imaginary is zero)\n");
  // { for (long i=0;i<size;i++) {printf("%lf ",t2[i][0]); } printf("\n"); }
  

  // field_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // field_fft(t2, _m.RowNo(),  _m.ColNo());
  // field_fft((double*)t2, _m.RowNo(),  _m.ColNo());

  // printf("matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());
  field_fftX(t2,  _m.ColNo(), _m.RowNo());
  
  //tmp block begin
  // double* t2=new double[size/2];			 //
  // { for (long i=0;i<size;i++) {printf("%lf, %lf  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // exit(0);
  // long aa=_m.ColNo();				 //
  // cpeds_matrix_save(cpeds_array2matrix(t,size, aa,true),"matrixfft.mat","");	 //
  // exit(0);
  // { for (long i=0;i<size;i++) {t2[i]=t[2*i]; } }	 //
  // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // cpeds_matrix_save(_m,"matrixfft-real.mat","");	 //
  // {  for (long i=0;i<size;i++) {t2[i]=t[2*i+1]; } }	 //
  // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // cpeds_matrix_save(_m,"matrixfft-imaginary.mat",""); //
  // tmp block end
  
  // do smoothing
  if (nfwhmX!=0) { _fwhmX=double(nfwhmX); }
  field_convolve_Gaussian_kernel(t2,  _m.ColNo(),_m.RowNo(),_fwhmX,0);

  // make the matrix symmetrical
  make_field_symmetricalCol(t2,  _m.ColNo(),_m.RowNo());


  // do inverse fftw
  // field_inv_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // field_inv_fft((double*)t2, _m.RowNo(),  _m.ColNo());
  field_fftX(t2,  _m.ColNo(), _m.RowNo(),false);
  // field_inv_fft(t, _m.RowNo(),  _m.ColNo());
  // printf("po inv fft\n");
  // { for (long i=0;i<size;i++) {printf("%lE,%lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // exit(0);
  // convert to matrix
  for (i=0;i<size;i++) { t[i]=t2[i][0]; }
  _m=cpeds_array2matrix(t,size, _m.RowNo(),false);

  // remove the padding
  //printf("deleting columns from %li to %li\n",_m.ColNo()-paddingSize,_m.ColNo()-1);

  //_m.delCol(_m.ColNo()-paddingSize,_m.ColNo()-1);

  // normalize after fft-ing
  _m/= _m.ColNo();

  delete [] t;  
  delete [] t2;  

  // associate the smoothed matrix with the input point set
  associateMatrix2Field();
  // update the point set  
  _pssm=cpedsPointSet3D().field2set(_m,_referencePointRow,_referencePointCol,at(0),res());

  return _m;
  
}
// ****************************************************************************************************
const matrix<double> cpedsSmooth::smoothGaussY(long nfwhmY, matrix<double>* m) {
  long i;
  if (m==NULL) make_field(); else _m=(*m);
  // printf("rows %li cols %li\n ",_m.RowNo(),  _m.ColNo());
  // add padding at the end of each row
  // if (_m.ColNo()%2==0) paddingSize=2; // if even number of cols we need two padded columns 
  // else paddingSize=1; // otherwise just one

  // printf("rows %li cols %li  padding: %li\n",_m.RowNo(),  _m.ColNo(),paddingSize);
  //cpeds_print_matrix(_m);
  // cpeds_matrix_save(_m,"matrix.dat", "");

  //printf("rows %li cols %li before padding in smoothGauss\n",_m.RowNo(),  _m.ColNo());

  //  _m=cpeds_pad_matrix(_m,0,0,0,paddingSize);

  //printf("rows %li cols %li\n after padding in smoothGauss",_m.RowNo(),  _m.ColNo());
  // transform to linear rowMajor array
  long size;
  double *t=cpeds_matrix2array(_m,size,false); // this is colMajor to be consistent with the convention that the rows are the first dimention (X - although it's vertical) and columns are the last dimention (Y - although it's horizontal)
  
  // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // cpeds_print_matrix(_m);
  // printf("rows %li cols %li size %li padding: %li\n",_m.RowNo(),  _m.ColNo(),size, paddingSize);
  // do fftw

  fftw_complex* t2 = (fftw_complex*) fftw_malloc(size * sizeof(fftw_complex));
  for (i=0;i<size;i++) { t2[i][0]=t[i]; t2[i][1]=0; }
  // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_uniform_number(0,1); t2[i][1]=cpeds_random_uniform_number(0,1); } }
  // { for (long i=0;i<size;i++) { t2[i][0]=cpeds_random_uniform_number(0,1); t2[i][1]=0; } }
  // { for (long i=0;i<size;i++) { t2[i][0]=double(i); t2[i][1]=0; } }
  // printf("initial random array: real part (imaginary is zero)\n");
  // { for (long i=0;i<size;i++) {printf("%lf ",t2[i][0]); } printf("\n"); }
  

  // field_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // field_fft(t2, _m.RowNo(),  _m.ColNo());
  // field_fft((double*)t2, _m.RowNo(),  _m.ColNo());

  // printf("matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());
  field_fftY(t2,  _m.ColNo(), _m.RowNo());
  
  //tmp block begin
  // double* t2=new double[size/2];			 //
  // { for (long i=0;i<size;i++) {printf("%lf, %lf  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // exit(0);
  // long aa=_m.ColNo();				 //
  // cpeds_matrix_save(cpeds_array2matrix(t,size, aa,true),"matrixfft.mat","");	 //
  // exit(0);
  // { for (long i=0;i<size;i++) {t2[i]=t[2*i]; } }	 //
  // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // cpeds_matrix_save(_m,"matrixfft-real.mat","");	 //
  // {  for (long i=0;i<size;i++) {t2[i]=t[2*i+1]; } }	 //
  // _m=cpeds_array2matrix(t2,size/2, aa/2,true);	 //
  // cpeds_matrix_save(_m,"matrixfft-imaginary.mat",""); //
  // tmp block end
  
  // do smoothing
  if (nfwhmY!=0) { _fwhmY=double(nfwhmY); }
  field_convolve_Gaussian_kernel(t2, _m.ColNo(), _m.RowNo(),0,_fwhmY);

  // make the matrix symmetrical
  make_field_symmetricalRow(t2, _m.ColNo(), _m.RowNo());


  // do inverse fftw
  // field_inv_fft(t, _m.RowNo(),  _m.ColNo()-paddingSize);
  // field_inv_fft((double*)t2, _m.RowNo(),  _m.ColNo());
  field_fftY(t2,  _m.ColNo(), _m.RowNo(),false);
  // field_inv_fft(t, _m.RowNo(),  _m.ColNo());
  // printf("po inv fft\n");
  // { for (long i=0;i<size;i++) {printf("%lE,%lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // exit(0);
  // convert to matrix
  for (i=0;i<size;i++) { t[i]=t2[i][0]; }
  _m=cpeds_array2matrix(t,size, _m.RowNo(),false);
  // printf("end of smoothY: matrix rows: %li cols: %li\n",_m.RowNo(),  _m.ColNo());

  // remove the padding
  //printf("deleting columns from %li to %li\n",_m.ColNo()-paddingSize,_m.ColNo()-1);

  //_m.delCol(_m.ColNo()-paddingSize,_m.ColNo()-1);

  // normalize after fft-ing
  _m/=_m.RowNo();

  delete [] t;  
  delete [] t2;  


  // associate the smoothed matrix with the input point set
  associateMatrix2Field();
  // update the point set  
  _pssm=cpedsPointSet3D().field2set(_m,_referencePointRow,_referencePointCol,at(0),res());


  return _m;

}
// ****************************************************************************************************
// vecNum is the number of rows in the t array (the array is by default in row-major ordering throughout the class
// vecSize if the number of cols in the t array
// perhaps this naming convention is confusing to some and not very general in fact hence these comments might be useful

void cpedsSmooth::field_convolve_Gaussian_kernel(fftw_complex* t, long vecNum, long vecSize, double fwhmX, double fwhmY) {
  double kRow; // this is k-mode along rows - parallel to j direction
  double kCol; // this is k-mode along columns - parallel to i direction
  double sigma2fwhm=2.35482004503; //(2*sqrt(2*log(2.0)));
  double sigmaRow=fwhmY/sigma2fwhm*twoPI;  // this is smoothing dispersion in along rows; QUESTION: WHERE DOES THE 2PI FACTOR COME FROM ???
  double sigmaCol=fwhmX/sigma2fwhm*twoPI;  // this is smoothing dispersion in along cols; QUESTION: WHERE DOES THE 2PI FACTOR COME FROM ???
  double factor;
  long vecSize_snd, vecNum_snd;
  double twoSigmaRow2=sigmaRow*sigmaRow*2;
  double twoSigmaCol2=sigmaCol*sigmaCol*2;
  long i,j;

  vecSize_snd=vecSize/2+vecSize%2; // vecSize should always be even there due to fftw conventions
  vecNum_snd =vecNum/2+vecNum%2; // vecSize should always be even there due to fftw conventions

  // if (_fwhmX!=0 && _fwhmY!=0) {
     // vecSize_snd=vecSize/2+vecSize%2; // vecSize should always be even there due to fftw conventions
     // vecNum_snd =vecNum/2+vecNum%2; // vecSize should always be even there due to fftw conventions
  //   for (long i=0;i<vecNum;i++) {
  //     for (long j=0;j<vecSize_snd;j++) { 
  //   	kx=double(i)/double(vecNum)/_Lx;
  //   	ky=double(j)/double(vecSize_snd)/_Ly;
  //   	factor=exp(-kx*kx*sigmaX2o2-ky*ky*sigmaY2o2);
  //   	t2[i*vecSize+j*2  ]*=factor; // real part
  //   	t2[i*vecSize+j*2+1]*=factor; // imaginary part
  //     }
  //   }


  // fftw_complex* t2=(fftw_complex*)t;
  // printf("vecSize: %li  vecNum %li vecSize_snd: %li  vecNum_snd %li sigmaRow:%lE sigmaCol: %lE _Lx: %lE _Ly: %lE\n",vecSize, vecNum, vecSize_snd, vecNum_snd, sigmaRow, sigmaCol,_Lx,_Ly);
  // printf("before smoothing\n");
  // long size=vecSize*vecNum;
  // { for (i=0;i<size;i++) {printf("%lE, %lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
  // { for (i=0;i<vecNum;i++) {
  //     for (j=0;j<vecSize;j++) { printf("%lE, %lE  ",t[i*vecSize+j][0],t[i*vecSize+j][1]); } printf("\n"); }
  //   printf("\n");
  // }
    
  // if (_fwhmX!=0 && _fwhmY!=0) {
  // long k=0;
    // for (long i=0;i<vecNum;i++) {
  double tmp;
  if (twoSigmaRow2==0) {
    for (i=0;i<vecNum;i++) {
      for (j=0;j<vecSize;j++) { 
	kCol=double(i)/double(vecNum);///double(vecSize);///_Ly; // this is k-mode along rows - parallel to j direction
	tmp=kCol*sigmaCol;
	factor=exp(-0.5 * tmp*tmp );
	t[i*vecSize+j][0]*=factor; // real part
	t[i*vecSize+j][1]*=factor; // imaginary part
      }
    }
  }
  if (twoSigmaCol2==0) {
    for (j=0;j<vecSize_snd;j++) { 
      kRow=double(j)/double(vecSize);///_Ly; // this is k-mode along rows - parallel to j direction
      tmp=kRow*sigmaRow;
      factor=exp(-0.5 * tmp*tmp); // note that this is the same as 2^(-4 j^2/_fwhm^2) - where _fwhm defines here the number of points that corrrespond to fwhm
      // factor= exp( -0.5 * pow(double(j)/double(vecSize),2.0)  *  pow(_fwhmY/double(vecSize)/sigma2fwhm,2.0) );
      // printf("*** tmp: %lE factor: %lE fwhmY: %lE\n", tmp, factor, _fwhmY );
      for (i=0;i<vecNum;i++) {
	t[i*vecSize+j][0]*=factor; // real part
	t[i*vecSize+j][1]*=factor; // imaginary part
      }
    }
  }
  if (twoSigmaCol2!=0 && twoSigmaRow2!=0) {
    double tmp1,tmp2;
    for (i=0;i<vecNum;i++) {
      // for (j=0;j<vecSize_snd;j++) { 
      for (j=0;j<vecSize;j++) { 
	kCol=double(i)/double(vecNum);///_Lx;
	// ky=double(j)/double(vecSize_snd)/_Ly;
	kRow=double(j)/double(vecSize);///_Ly; // this is k-mode along rows - parallel to j direction
	tmp1=kRow*sigmaRow;
	tmp2=kCol*sigmaCol;
	factor=exp(-0.5 * (tmp1*tmp1+tmp2*tmp2) );
	// printf("*** factor %lE tmp: %lE,  kRow: %lE kCol: %lE sigmaRow2o2: %lE , sigmaCol2o2: %lE\n",factor, tmp , kRow,kCol,twoSigmaRow2,twoSigmaCol2);
	t[i*vecSize+j][0]*=factor; // real part
	t[i*vecSize+j][1]*=factor; // imaginary part
	// t2[i*vecSize+j][0]=k;
	// t2[i*vecSize+j][1]=-k;
	// k++;
      }
      // printf("row %li, smoothed 1st half\n",i);
      // { for (long i=0;i<size;i++) {printf("%lE, %lE  ",t[i][0],t[i][1]); } printf("\n"); }
    }
  }
    
}

// ****************************************************************************************************
// WARNING !!!
// something is wrong with this routine -- there's some missing symmetry here or something else screwed up !!!
// WARNING !!!

void cpedsSmooth::make_field_symmetrical(fftw_complex* t, long vecNum, long vecSize) {
  
  long vecSize_snd, vecNum_snd;

  vecSize_snd=vecSize/2+vecSize%2; // vecSize should always be even there due to fftw conventions
  vecNum_snd =vecNum/2+vecNum%2; // vecSize should always be even there due to fftw conventions

  long i,j,ist,jst;
  i=0;
  // make sure the first element in the array is real - this is due to the Hermitian symmetry of course
  t[0][1]=0; 
  // make sure the first column is real due to Y_0==Y_vecSize^* symmetry where * indicates the complex conjugate
  for (i=0;i<vecNum;i++) {
    t[i*vecSize][1]=0; 
  }
  // make sure the first row is real due to Y_0==Y_vecSize^* symmetry where * indicates the complex conjugate
  for (j=0;j<vecSize;j++) {
    t[j][1]=0; 
  }

  // symmetrize the first row
  jst=vecSize_snd;
  for (j=jst;j<vecSize;j++) { 
    t[j][0]= t[vecSize-j][0]; // real part	 
    t[j][1]=-t[vecSize-j][1]; // imaginary part	 
  }

  // symmetrize first column
  ist=vecNum_snd;
  for (i=ist;i<vecNum;i++) {
    t[i*vecSize][0]= t[((vecNum-i)%vecNum)*vecSize][0]; // real part	 
    t[i*vecSize][1]=-t[((vecNum-i)%vecNum)*vecSize][1]; // imaginary part	 
  }

  // symmetrizing the elements that should be real for the even rows number fields due to X_i==X_{vecSize-i} for  even vecSize=2*vecSize_snd we have X_vecSize_snd==X_vecSize-vecSize_snd^* which is the same element and so must be real
  if (vecNum%2==0) {
    ist=vecNum_snd;
    for (j=0;j<vecSize;j++) {
      t[ist*vecSize+j][1]=0;
    }
  }


  // symmetrizing the elements that should be real for the even columns number fields
  for (i=1;i<vecNum;i++) {
    // printf("imposing symmetry\n");
    // imposing Hermitian symmetry
    jst=vecSize_snd;
    if (vecSize%2==0) {
      t[i*vecSize+jst][1]=0;
      jst++;
    }
    
    // //antysymmetrizing the first column
    // if (i>=vecNum_snd) {
    //   t[i*vecSize][0]= t[((vecNum-i)%vecNum)*vecSize][0]; // real part	 
    //   t[i*vecSize][1]=-t[((vecNum-i)%vecNum)*vecSize][1]; // imaginary part	 
    // }
    
    
    // antysymmetrizing the rest of the matrix
    // the weired symmetry comes from the fact that the FFT(f(x)^*)==FFT(f(-x)) and from the periodicity of f_i==f_{vecSize-i} condition
    for (j=jst;j<vecSize;j++) { 
      t[i*vecSize+j][0]= t[((vecNum-i)%vecNum)*vecSize+vecSize-j][0]; // real part	 
      t[i*vecSize+j][1]=-t[((vecNum-i)%vecNum)*vecSize+vecSize-j][1]; // imaginary part	 
    }
    // printf("row %li: anty-symmetrized\n",i);
    // { for (long i=0;i<size;i++) {printf(" %lE, %lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
    
  }



    // BELOW is an old debug stuff

    // antisymmetrizing the first column - it should be antysymmetrical due to the second fft (along the fist direction) because of the reality of the zeroth frquencies in the first column after fft along last dimention
    // for (long i=vecNum_snd;i<vecNum;i++) {
    //   t2[i*vecSize][0]= t2[(vecNum-i)*vecSize][0];
    //   t2[i*vecSize][1]=-t2[(vecNum-i)*vecSize][1];
    // }
    // { for (long i=0;i<size;i++) {printf(" %lE, %lE  ",t2[i][0],t2[i][1]); } printf("\n"); }
   // exit(0);



    // for (long j=0;j<vecSize;j++) { 
    //   t2[j][1]=0;
    //   ist=vecNum_snd;
    //   if (vecNum%2==0) {
    // 	t2[ist*vecSize+j][1]=0;
    // 	ist++;
    //   }

    //   // imposing Hermitian symmetry in x direction
    //   for (long i=ist;i<vecNum_snd;i++) {
    // 	t2[i*vecSize+j][0]= t2[(vecNum-i)*vecSize+j][0]; // real part	 
    // 	t2[i*vecSize+j][1]=-t2[(vecNum-i)*vecSize+j][1]; // imaginary part	 
    //   }
    // }



  // }



}

// ****************************************************************************************************
// enforces the Hermitian symmetry for every row of the t array
void cpedsSmooth::make_field_symmetricalRow(fftw_complex* t, long vecNum, long vecSize) {
  long vecSize_snd, vecNum_snd;
  long i,j,jst;

  vecSize_snd=vecSize/2+vecSize%2; // vecSize should always be even there due to fftw conventions
  vecNum_snd =vecNum/2+vecNum%2; // vecSize should always be even there due to fftw conventions

  // make sure the first column is real due to Y_0==Y_vecSize^* symmetry where * indicates the complex conjugate
  i=0;
  for (i=0;i<vecNum;i++) {
    t[i*vecSize][1]=0; 
  }
  // symmetrize each row
  for (i=0;i<vecNum;i++) {
    jst=vecSize_snd;
    for (j=jst;j<vecSize;j++) { 
      t[i*vecSize+j][0]= t[i*vecSize+vecSize-j][0]; // real part	 
      t[i*vecSize+j][1]=-t[i*vecSize+vecSize-j][1]; // imaginary part	 
    }
  }

  // symmetrizing the elements that should be real for the even columns number fields due to X_i==X_{vecSize-i} for  even vecSize=2*vecSize_snd we have Y_vecSize_snd==Y_vecSize-vecSize_snd^* which is the same element and so must be real
  if (vecSize%2==0) {
    jst=vecSize_snd;
    for (i=0;i<vecNum;i++) {
      t[i*vecSize+jst][1]=0;
    }
  }
    

}

// ****************************************************************************************************
void cpedsSmooth::make_field_symmetricalCol(fftw_complex* t, long vecNum, long vecSize) {
  long vecSize_snd, vecNum_snd;
  long i,j,ist;

  vecSize_snd=vecSize/2+vecSize%2; // vecSize should always be even there due to fftw conventions
  vecNum_snd =vecNum/2+vecNum%2; // vecSize should always be even there due to fftw conventions

  // make sure the first row is real due to Y_0==Y_vecSize^* symmetry where * indicates the complex conjugate
  i=0;
  for (j=0;j<vecSize;j++) {
    t[j][1]=0; 
  }
  // symmetrize each column
  ist=vecNum_snd;
  for (i=ist;i<vecNum;i++) {
    for (j=0;j<vecSize;j++) { 
      t[i*vecSize+j][0]= t[((vecNum-i)%vecNum)*vecSize+j][0]; // real part	 
      t[i*vecSize+j][1]=-t[((vecNum-i)%vecNum)*vecSize+j][1]; // imaginary part	 
    }
  }

  // symmetrizing the elements that should be real for the even rows number fields due to X_i==X_{vecSize-i} for  even vecSize=2*vecSize_snd we have X_vecSize_snd==X_vecSize-vecSize_snd^* which is the same element and so must be real
  if (vecNum%2==0) {
    ist=vecNum_snd;
    for (j=0;j<vecSize;j++) {
      t[ist*vecSize+j][1]=0;
    }
  }
  

}

// ****************************************************************************************************
cpedsSmooth& cpedsSmooth::operator=(const cpedsSmooth& rhs) {
  clearAll();
  append(rhs.getPoints());
  _m=rhs.getSmoothedField();
  res()=rhs.getRes();
  border()=rhs.getBorder();
  setSmoothingLength(rhs.getSmoothingLength().x(), rhs.getSmoothingLength().y());
  _pssm=rhs.getSmoothedPointSet();
  _projectionPlane=rhs.getProjectionPlane();
  _referencePointRow=rhs._referencePointRow;
  _referencePointCol=rhs._referencePointCol;

  // printf("\nrhs\n\n");
  // rhs.printInfo();
  // printf("\nlhs\n\n");
  // printInfo();

  return *this;
}
// ****************************************************************************************************
// const cpedsSmooth& cpedsSmooth::operator=(const cpedsPointSet3D &rhs) {
//   clear();
//   long N=rhs.count();
//   for (long i=0;i<N;i++) { append(rhs.at(i));  }
//   return *this;
// }
// ****************************************************************************************************
 void cpedsSmooth::associateMatrix2Field() {
   double minx,maxx,miny,maxy,minz,maxz;
   _referencePointNumber=0;
   getRanges(minx,maxx,miny,maxy,minz,maxz);
   cpedsPoint3D borderSize=getBorder();

    if (maxx==minx) { _referencePointCol=0; }
    else {
      _referencePointCol=long((at(_referencePointNumber).x()-minx) / (maxx-minx) * _m.ColNo() - 0.5 - borderSize.x());
    }

    if (maxy==miny) { _referencePointRow=0; }
    else { 
      _referencePointRow=long((at(_referencePointNumber).y()-miny) / (maxy-miny) * _m.RowNo() - 0.5 - borderSize.y());
    }

 }
// ****************************************************************************************************
void cpedsSmooth::printInfo() const {

  printf("cpedsSmooth object information:\n\n");
  printf("ON PLANE QUANTITIES\n");
  printf("border size: %lE,  %lE\n",getBorder().x(),getBorder().y());
  printf("resolution: %lE,  %lE\n",getRes().x(),getRes().y());
  printf("smoothing length: %lE,  %lE\n",getSmoothingLength().x(),getSmoothingLength().y());
  printf("smoothing length pix: %lE,  %lE\n",getSmoothingLengthPix().x(),getSmoothingLengthPix().y());
  printf("size of the set to be smooted: %li\n",long(count()));
  printf("size of the smooted set: %li\n",long(_pssm.count()));
  printf("field size: rows %li, cols %li\n",long(_m.RowNo()), long(_m.ColNo()));

  printf("\n\n");

}
// ****************************************************************************************************
