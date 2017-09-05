/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/* IMPLEMENTATION OF THE cpeds_PDF CLASS */

/* the input arrays of the pdf function should be resonably resonable eg. sorter in rising abscissa order etc :P */
/* but itdoesn't need to be normalized I guess */
/* p1 p2 p3 CLs in rising order to the corresponding returned CL ranges */

/* the confidence limits can be derived in a couple of different ways. Here, I derive the CL ranges as the ranges */
/* integrating from the ML point (the mode) saidways untill a requested CL is reached. */

#include "cpeds-pdf.h"

/*****************************************************************************************************************/
cpeds_PDF::cpeds_PDF(double *x, double *y, long n, cpeds_queue_double *cl) {
  msgs=new cpedsMsgs("PDF");
  initiate_CLs_info(long(0));
  initiate_CLs_info(cl);

  N=n;
  if (N==1) { msgs->say("WARNING !! PDF of size 1 doesn't make much sense",High); }
  pdf = new cpeds_point[N];
  // copy PDF onto object structure
  for (i=0;i<n;i++) { pdf[i].x=x[i]; pdf[i].y=y[i]; }
  DX=(pdf[1].x-pdf[0].x);
  pdf_fields_swapped=false;

  EX=0;
  MODE=0;
  norm=0;
  EXdone=false;
  MODEdone=false;
  normdone=false;
  MODES=NULL;
  MODESidx=NULL;

  object_name="PDF";
}

/*****************************************************************************************************************/
cpeds_PDF::cpeds_PDF(double *x, double *y, long n) {
  msgs=new cpedsMsgs("PDF");
  initiate_CLs_info((long)0);

  N=n;
  if (N==1) { msgs->say("WARNING !! PDF of size 1 doesn't make much sense",High); }
  pdf = new cpeds_point[N];
  // copy PDF onto object structure
  for (i=0;i<n;i++) { pdf[i].x=x[i]; pdf[i].y=y[i]; }
  DX=(pdf[1].x-pdf[0].x);
  pdf_fields_swapped=false;

  EX=0;
  MODE=0;
  norm=0;

  object_name="PDF";
}

/*****************************************************************************************************************/
cpeds_PDF::~cpeds_PDF() {
  if (msgs!=NULL) { delete msgs; msgs=NULL; }
  destroy_CLs_info();
  delete [] pdf;
}

/*****************************************************************************************************************/
void cpeds_PDF::initiate_CLs_info(cpeds_queue_double *cl) {
  if (cl==NULL) {    initiate_CLs_info((long)0);  }
  else {
    destroy_CLs_info();
    NCL=cl->get_size();  
    CL          =new double[NCL]; for (i=0;i<NCL;i++) CL[i]=(*cl)(i);
    CL2side_down=new double[NCL]; 
    CL2side_up  =new double[NCL];  
    CLRpts = new cpeds_queue_double[NCL];
  }
}
/*****************************************************************************************************************/
void cpeds_PDF::initiate_CLs_info(long NCLs) {
  cpeds_queue_double tmp; 
 
  if (NCLs==0) { 
    NCL=0;
    CL=NULL;
    CL2side_down=NULL;
    CL2side_up=NULL;
    CLRpts=NULL;
  }
  else { 
    tmp.mk_list(0,NCLs);
    initiate_CLs_info(&tmp);  
  }
}


/*****************************************************************************************************************/
void cpeds_PDF::destroy_CLs_info() {
  if (CL!=NULL) {delete [] CL; CL=NULL; }
  if (CL2side_down!=NULL) { delete [] CL2side_down; CL2side_down=NULL; }
  if (CL2side_up!=NULL) { delete [] CL2side_up; CL2side_up=NULL; }
  if (CLRpts!=NULL) { delete [] CLRpts; CLRpts=NULL; }
}

/*****************************************************************************************************************/
void cpeds_PDF::set_name(string name) { object_name=name; }
string cpeds_PDF::get_name() { return object_name; }

/*****************************************************************************************************************/
void cpeds_PDF::get_info() {
  derive_EX();
  derive_MODE();
  derive_CLs();
}

/*****************************************************************************************************************/
void cpeds_PDF::find_normalization() {
  norm=0;
  for (i=0;i<N;i++) {   norm+=pdf[i].y; }
  //  for (i=0;i<N;i++) {   pdf[i].y/=norm; }
  normdone=true;
}

/*****************************************************************************************************************/
void cpeds_PDF::derive_EX() {
  //////////////////////
  // derive the EX value
  //////////////////////
  find_normalization();
  for (i=0;i<N;i++) {   EX+=pdf[i].x*pdf[i].y; }
  EX/=norm;
  EXdone=true;
}

/*****************************************************************************************************************/
void cpeds_PDF::derive_MODE() {
  /////////////////
  // find MODE value
  /////////////////
  long tmp;
  cpeds_point p;
  switch_pdf_fileds();
  p=cpeds_find_max_value(pdf,N,0,&tmp);
  MODE=p.y;
  switch_pdf_fileds();
  MODEdone=true;
}

/*****************************************************************************************************************/
void cpeds_PDF::switch_pdf_fileds() {
  double tmp;
//  msgs->say("swapping PDF fields",Low);
  for (i=0;i<N;i++) { tmp=pdf[i].x; pdf[i].x=pdf[i].y; pdf[i].y=tmp; }
  if (pdf_fields_swapped) pdf_fields_swapped=false; else pdf_fields_swapped=true;
}

/*****************************************************************************************************************/
void cpeds_PDF::derive_CLs() {
  double A;
  double p;

  find_normalization();
  switch_pdf_fileds();
  //
  // sort the PDF according to decreasing likelihood keeping the arguments information
  //
  cpeds_sort_data(N,pdf,21);
  //
  // find the set of points corresponding to the CL range of the requested CL
  //
  switch_pdf_fileds();
  
  A=0;
  i=0;
  for (j=0;j<NCL;j++) { // loop over the CLs to be resolved
    p=CL[j]*norm;
    while (A<p && i<N) {
      A+=pdf[i].y;
      for (k=j;k<NCL;k++) CLRpts[k].addq(pdf[i].x);
/*       printf("i=%li N=%li norm*CL[j=%li]=%lE A=%lE pdf[i].x=%lE pdf[i].y=%lE\n",i,N,j,p,A,pdf[i].x,pdf[i].y); */
      i++;
    }
  }

  //
  // set the CL ranges assuming an unimodal PDF
  //
  for (j=0;j<NCL;j++) {
    CLRpts[j].sort(12);
    CL2side_down[j]=CLRpts[j](0);
    CL2side_up[j]=CLRpts[j](CLRpts[j].get_size()-1);
  }
}

/*****************************************************************************************************************/
void cpeds_PDF::derive_CLs(cpeds_queue_double *cl) {
  destroy_CLs_info();
  initiate_CLs_info(cl);
  derive_CLs();
}

/*****************************************************************************************************************/
void cpeds_PDF::integrate_to_CLr(double CLr, bool gaps, long nthmode) {
  double A;
  double xi;
  long imode,i;
  double*args;
  //long CLri;
  bool gapsOK,CLr_reached,modesOK, finish, CLRpts_sorted,imode_updated;

  // find the indes of the argument of the PDF up to which we gonna interate
  args=cpeds_point2array(N,pdf,0);
  cpeds_sort_data(N,args,12);
  CLr=args[cpeds_find_value(CLr,args,N,0,N)];
  delete [] args;

  //  printf("will look for %lE",CLr);

  find_normalization();
  switch_pdf_fileds();
  //
  // sort the PDF according to decreasing likelihood keeping the arguments intormation
  //
  cpeds_sort_data(N,pdf,21);
  //
  // find the set of points corresponding to the CL range of the requested CL
  //
  switch_pdf_fileds();

  
  initiate_CLs_info((long)1);
  A=0;
  i=0;
  imode=0;
  CLr_reached=false;
  if (nthmode==0) modesOK=true; else modesOK=false; 
  CLRpts_sorted=false;
  finish=false;
  imode_updated=false;
  if (gaps) gapsOK=true; else gapsOK=false; // check if there are gaps in the integrated part so far if it matters

  do {
    xi=pdf[i].x; // chose next argument
    A+=pdf[i].y; // integrate it

    if ( CLRpts_sorted ) CLRpts[0].insertq_ord(xi,12); // store the integrated argument in a sorted manner
    else CLRpts[0].addq(xi); // store the integrated argument in an unsorted manner

/*     printf("pts queue is: %li large\n",CLRpts[0].get_size()); */

    if (xi==CLr) CLr_reached=true; // watch if we have integrated at least up to the requested CLr value


    if (CLr_reached) { // if CLr is reached then other contidions can be checkek
      if (! modesOK) {

	if (! imode_updated) { imode=MODES->count_common_vals(&CLRpts[0]); imode_updated=true; } // update the imode counter
      	else {
	  if (MODES->val_in_list(xi)) { // cound how many modes we have already integrated over
	    imode++; 
	  }
	}
	if (imode==nthmode) modesOK=true;
      }
      else {
	if (! gapsOK) {
	  if ( ! CLRpts_sorted) { CLRpts[0].sort(12); CLRpts_sorted=true;  }
	  gapsOK=!(CLr_gaps(0)); 
	  if (gapsOK) finish=true;
	}
	else finish=true;	
      }
    }
        



/*     printf("i=%li xi=%lE A=%lE imode=%li CLr = %lE ",i,xi,A,imode,CLr); */
    msgs->say("gaps exist: "+msgs->toStr(gapsOK),Medium);
    msgs->say("gapsOK exist: "+msgs->toStr(gapsOK),Medium);
    msgs->say("CLr_reached: "+msgs->toStr(CLr_reached),Medium);
    msgs->say("modesOK: "+msgs->toStr(modesOK),Medium);
    i++;
    
/*   } while ( ! (CLr_reached && modesOK && gapsOK) );      */
  } while ( ! finish );

  //
  // set the CL ranges assuming an unimodal PDF
  //
  for (j=0;j<NCL;j++) {
    CL[j]=A/norm;
    if (! CLRpts_sorted) CLRpts[j].sort(12);
    CL2side_down[j]=CLRpts[j](0);
    CL2side_up[j]=CLRpts[j](CLRpts[j].get_size()-1);
  }


}
/*****************************************************************************************************************/
void cpeds_PDF::get_PDF_modes(cpeds_queue_double** PDF_modes, cpeds_queue_long** PDF_modes_idx) {
  long i,Nleo=N-1;
  double p1,p3;

  cpeds_queue_double *mv = new cpeds_queue_double();
  cpeds_queue_long *mi = new cpeds_queue_long();

  if (pdf_fields_swapped) {     switch_pdf_fileds(); }
  cpeds_sort_data(N,pdf,12);

  for (i=1;i<Nleo;i++) {
    p1=pdf[i-1].y;
    p3=pdf[i+1].y;

/*     printf("i: %li x: %lE y:%lE \n",i,pdf[i].x,pdf[i].y); */

    if (pdf[i].y-p1 >= 0 && pdf[i].y-p3 > 0) {  // be ware of some stupid PDFs with flat features in it. or if you care, then elaborate on this
      mv->addq(pdf[i].x); mi->addq(i);
/*       printf("adding i: %li x: %lE y:%lE p1: %lE p3: %lE     pdf[i].y-p1: %lE, pdf[i].y-p3: %lE \n",i,pdf[i].x,pdf[i].y, p1,p3,  pdf[i].y-p1,  pdf[i].y-p3  ); */
    }
  }

  if (MODES != NULL) delete MODES;        
  MODES=mv;
  if (MODESidx != NULL) delete MODESidx;  
  MODESidx=mi;
  if (PDF_modes !=NULL) (*PDF_modes)=mv;
  if (PDF_modes_idx !=NULL) (*PDF_modes_idx)=mi;     
}

/*****************************************************************************************************************/
// handlers
double cpeds_PDF::get_EX() {   if (!EXdone) derive_EX(); return EX; }
double cpeds_PDF::get_MODE() { if (!MODEdone) derive_MODE(); return MODE; }
double cpeds_PDF::get_CL_range(bool lower, long i, bool* cont) {
  double r=0;
  if (i>=0 && i<NCL) {
    if (lower) r=CL2side_down[i]; 
    else r=CL2side_up[i]; 
  }
  if (cont!=NULL) (*cont)=CLr_gaps(i);
  return r;
}

double cpeds_PDF::get_CL(long i) { return CL[i]; } // returns the i'th CL from internal CL array
/*****************************************************************************************************************/
bool cpeds_PDF::CLr_gaps(long j) {
  long imax=CLRpts[j].get_size()-1;
  double delta=1.1*DX;

  if (imax < 1) return false;
  for (i=0;i<imax;i++) {
    //    printf("checking gaps: delta=%lE d=%lE\n",delta,CLRpts[j](i+1)-CLRpts[j](i));
    if (CLRpts[j](i+1)-CLRpts[j](i) > delta) return true;
  }
    
  return false;  
}

/*****************************************************************************************************************/
void cpeds_PDF::save_CLr(string f,long n) {
  CLRpts[0].saveq(f,"double");  
}







































/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/
/*****************************************************************************************************************/

/* the input arrays of the pdf function should be resonably resonable eg. sorter in rising abscissa order etc :P */
/* but itdoesn't need to be normalized I guess */
/* p1 p2 p3 CLs in rising order to the corresponding returned CL ranges */

/* the confidence limits can be derived in a couple of different ways. Here, I derive the CL ranges as the ranges */
/* integrating from the ML point (the mode) saidways untill a requested CL is reached. */

/* cpeds_PDF* cpeds_get_PDF_info(double *X,double *Y, long N, cpeds_queue_double *CL) { */
/*   /////////// */
/*   // initiate */
/*   /////////// */
/*   long i,j,iEX,Nleo=N-1; */
/*   double x,y,norm,*XY,A,tmp; */
/*   cpeds_point * pt = new cpeds_point[N]; */
/*   cpeds_PDF_info* pdf=new cpeds_PDF_info; */
 
/*   long iclsort[pdf->NCL]; */
/*   //////////////// */
/*   // normalize PDF; as for the discrete case */
/*   //////////////// */
/* /\*   norm=cpeds_integrate_1d(N,X,Y); *\/ */
/*   norm=0; */
/*   for (i=0;i<N;i++) { norm+=Y[i]; } */
/* /\*   printf("normalization: %lE\n",norm); *\/ */
/*   for (i=0;i<N;i++) { Y[i]/=norm; } */

/*   ///////////////// */
/*   // find MODE value */
/*   ///////////////// */
/*   y=Y[0]; */
/*   for (i=0;i<N;i++) {    if (Y[i]>y) { pdf->MODE=X[i]; y=Y[i]; }  } */

/*   ////////////////////// */
/*   // derive the EX value */
/*   ////////////////////// */
/*   XY=new double[N]; */
/*   for (i=0;i<N;i++) {    XY[i]=X[i]*Y[i]; } */
/*   pdf->EX=cpeds_integrate_1d(N,X,XY); */
/*   iEX=cpeds_find_value(pdf->EX,X,N,0,N); */
/*   delete [] XY; */

/*   ////////////////////////////////////// */
/*   // derive the confidence limits   */
/*   ////////////////////////////////////// */
/* /\*   printf("deriving the CLrs\n"); *\/ */


/*   // copy PDF onto points structure */
/*   for (i=0;i<N;i++) { pt[i].x=Y[i]; pt[i].y=X[i]; } */
 
/*   // */
/*   // sort the PDF according to decreasing likelihood keeping the arguments intormation */
/*   // */
/*   cpeds_sort_data(N,pt,21); */

/*   // */
/*   // find the set of points corresponding to the CL range of the requested CL */
/*   // */
  

/*   A=0; */
 
/*   for (j=0;j<CL->get_size();j++) { // loop over the CLs to be resolved */
/*     while (A<cl(j) && i<Nleo) { */
/*       A+=pt[i].X; */
/*       printf("i=%li A=%lE pt[i].x=%lE pt[i].y=%lE\n",i,A,pt[i].x,pt[i].y); */
/*       i++; */
/*     } */
/*     iclsort[j]=i-1; */
/*   } */



/*   // */
/*   // clean up */
/*   // */
/*   delete [] pt; */
/*   delete [] iclsort; */



/*   if (iEX<=1) { printf("EX value not sufficiently resolved, cant help you; try giving pdf with larger amount of point\n"); } */

/*   // lower CL ranges */
/*   A=0; */
/*   i=iEX-1; */
 
/*   for (j=0;j<CL->get_size();j++) { // loop over the CLs to be resolved */
/*     while (A<cl(j) && i>=0) { */
/*       tmp=(X[i+1]-X[i]) * 0.5*(Y[i+1]+Y[i]); */
/*       A+=tmp; */
/*       printf("i=%li A=%lE dA=%lE X[i]=%lE\n",i,A,tmp,X[i]); */

/*       //A+=(X[i+1]-X[i]) * 0.5*(Y[i+1]+Y[i]); */
/*       i--; */
/*     } */
/*     pdf->CL2side_down[j]=X[i+1]; */
/*   } */
/* /\*   printf("deriving the up CLrs\n"); *\/ */

/*   // upper CL ranges */
/*   A=0; */
/*   i=iEX; */
 
/*   for (j=0;j<CL->get_size();j++) { // loop over the CLs to be resolved */
/*     while (A<cl(j) && i<Nleo) { */
/*       //tmp=(X[i+1]-X[i]) * 0.5*(Y[i+1]+Y[i]); */
/*       //A+=tmp; */
/* /\*       printf("i=%li A=%lE dA=%lE X[i]=%lE\n",i,A,tmp,X[i]); *\/ */
/*       A+=(X[i+1]-X[i]) * 0.5*(Y[i+1]+Y[i]); */
/*       i++; */
/*     } */
/*     pdf->CL2side_up[j]=X[i-1]; */
/*   } */

/*   return pdf; */
/* } */
