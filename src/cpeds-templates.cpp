#include "cpeds-templates.h"
#include "cpeds-math.h"



template <class T> cpeds_queue<T>::cpeds_queue(string name) {
  Q=NULL;
  Qlast=NULL;
  Qsize=0;
  Qmax=Qmin=0;
  iQmax=iQmin=0;
  set=false;
  qname=name;
/*     printf("initiating new queue and resetting parameters\n"); */
}

template <class T> cpeds_queue<T>::cpeds_queue() {
  Q=NULL;
  Qlast=NULL;
  Qsize=0;
  Qmax=Qmin=0;
  iQmax=iQmin=0;
  set=false;
  qname="queue";
/*     printf("initiating new queue and resetting parameters\n"); */
}

template <class T> cpeds_queue<T>::~cpeds_queue() {
  if (Qsize!=0) delete_all_queue();
}

template <class T> void cpeds_queue<T>::delete_all_queue() {
  queue_type *qp=Q;

  /*     printf("deleting queue\n"); */
  while (Q!=NULL) {  Q=(*Q).n; delete qp; qp=Q; }
  Q=NULL;
  Qlast=NULL;
  Qsize=0;
  Qmax=Qmin=0;
  iQmax=iQmin=0;
  set=false;
}

template <class T> long cpeds_queue<T>::get_size() { return Qsize; }

template <class T> void cpeds_queue<T>::addq(cpeds_queue<T> *q) {
  long i,num=q->get_size();
  for (i=0;i<num;i++) {
    addq(q->getq(i));
  }
}


template <class T> void cpeds_queue<T>::addq(T v) {
  queue_type *tmp = new queue_type;
  /*     printf("adding to queue\n"); */
  if (!set) { Q = tmp; Qlast = tmp; }
  tmp->n=NULL;
  tmp->v=v;
  Qsize++;
  if (set) {
    Qlast->n=tmp; Qlast=tmp; 
    if (v < Qmin) { Qmin=v; iQmin=Qsize-1; }
    if (v > Qmax) { Qmax=v; iQmax=Qsize-1; }
  }
  else {
    set = true; 
    Qmin=v; iQmin=Qsize-1;
    Qmax=v; iQmax=Qsize-1;
  }
}

template <class T> void cpeds_queue<T>::setq(long I, T v) { // sets the I'th cell to value v; the quque is extended if needed and filled with zeros upto the requested index I
  long i,N;
  
  if (I >= Qsize) {
    N = I-Qsize+1; 
    for (i=0;i<N;i++) { addq((T)0); }
  }
  setqOK(I,v);
}


template <class T> T cpeds_queue<T>::getq(long I) {  // returns the I'th value in the queue; numbering starts from 0
  long i;
  queue_type* tmp=Q;
  
  if (I >= Qsize || I<0) { return 0; }
  for (i=0;i<I;i++) { tmp=tmp->n; }
  return tmp->v;
}

/* template <class T> queue_type* cpeds_queue<T>::getQ(long I) {  // returns the I'th value in the queue; numbering starts from 0 */
/*   long i; */
/*   queue_type* tmp=Q; */
  
/*   if (I > Qsize || I<0) { return NULL; } */
/*   for (i=0;i<I;i++) { tmp=tmp->n; } */
/*   return tmp; */
/* } */



template <class T> T cpeds_queue<T>::getq_max(long *i) { if (i!=NULL) *i=iQmax; return Qmax; }
template <class T> T cpeds_queue<T>::getq_min(long *i) { if (i!=NULL) *i=iQmin; return Qmin; }

template <class T> void cpeds_queue<T>::delq(long I) { // deletes the I'th cell from the list and merges the rest
  queue_type *prev,*rem,*next;

  if (I>=Qsize) return;
  if (I<0) return;

  prev = getptr(I-1);
  rem = prev->n; if (rem == NULL) return;
  next = rem->n;
  delete rem;
  prev->n=next;    
  Qsize--;
  if (I == iQmin || I == iQmax) { update_Qminmax(); }
}

template <class T> void cpeds_queue<T>::delq_from(long I) { // deletes all cells from I'th cell to the end of the list
  long i,num=get_size();
  for (i=I;i<num;i++) { delq(I); }
}

template <class T> void cpeds_queue<T>::delq() { // deletes the last cell in the list
  delq(Qsize-1);
}


template <class T> void cpeds_queue<T>::update_Qminmax() {
  queue_type *tmp=Q;
  //T v;
  long i;
  
  Qmin=Qmax=tmp->v;
  
  for (i=0;i<Qsize;i++) { 
    if (tmp->v < Qmin) { Qmin=tmp->v; iQmin=i; }
    if (tmp->v > Qmax) { Qmax=tmp->v; iQmax=i; }      
    tmp=tmp->n; 
  }
}

template <class T> void cpeds_queue<T>::printq_long() {
  long i;
  queue_type *tmp=Q;    
  
  printf("|%s> ",qname.c_str()); 
  for (i=0;i<Qsize;i++) { 
    printf("%li ",tmp->v); tmp=tmp->n;
  }
  printf("\n%s","");
}

template <class T> void cpeds_queue<T>::printq_double() {
  long i;
  queue_type *tmp=Q;    
  
  printf("|%s> ",qname.c_str()); 
  for (i=0;i<Qsize;i++) { 
    printf("%lE ", tmp->v); tmp=tmp->n;
  }
  printf("\n");
}

template <class T>  void cpeds_queue<T>::printq_string() {
  long i;
  queue_type *tmp=Q;    
  
  for (i=0;i<Qsize;i++) { 
    printf("%s ",tmp->v.c_str()); tmp=tmp->n;
  }
  printf("\n");
}

template <class T> bool cpeds_queue<T>::all_same() {
  /*     queue_type* tmp=Q; */
  long i;
  T v=getq(0);
  
  for (i=0;i<Qsize;i++) { if (v!=getq(i)) return false; }
  return true;
}

template <class T> T* cpeds_queue<T>::export_array() {
  T* t = new T[get_size()];
  long i;
  for (i=0;i<get_size();i++) { t[i]=getq(i); }
  return t;
}

template <class T> QList<T> cpeds_queue<T>::exportQList() {
  QList<T> q;
  long i;
  for (i=0;i<get_size();i++) { q.append(getq(i));  }
  /* for (i=0;i<get_size();i++) { q.append(0); printf("appending %li\n",i); } */
  return q;
}

template <class T> void cpeds_queue<T>::import_array(T* t) {
  long i;
  for (i=0;i<get_size();i++) { 
    setq(i,t[i]); }
}

template <class T> void cpeds_queue<T>::import_array(T* t,long N) {
  long i;
  delete_all_queue();
  for (i=0;i<N;i++) { 
    addq(t[i]); 
  }
}

template <class T> void cpeds_queue<T>::add_array(T* t,long N) {
  long i;
  for (i=0;i<N;i++) { 
    addq(t[i]); 
  }
}

template <class T> void cpeds_queue<T>::sort(long direction) { // direction: 12 - ascending, 21 - decending
  T *t;
  // get array of values for sorting
  t=export_array();
  cpeds_sort_data(get_size(),t,direction);
  import_array(t);
}
  
template <class T>  void cpeds_queue<T>::mk_list(T from, T to, T step) {
  T i;
  if (step < 0) step=-step;

  if (from <= to) 
    for (i=from;i<=to;i+=step) {    addq(i); }
  else 
    for (i=from;i<=to;i-=step) {    addq(i); }
}

template <class T>  void cpeds_queue<T>::mk_list(T val, long size) {
  long i;
  for (i=0;i<size;i++) addq(val);
}


//****************************************************************************************************************
template <class T>  bool cpeds_queue<T>::val_in_list(T val) {
  long i;
  for (i=0;i<Qsize;i++) { if (getq(i) == val) return true; }
  return false;
}


//****************************************************************************************************************
template <class T>  void cpeds_queue<T>::saveq(string f,string qtype) {
  FILE* F;
  long i;

  F=fopen(f.c_str(),"w");
  if (qtype == "double") { for (i=0;i<Qsize;i++) { fprintf(F,"%.10lE\n",(double)getq(i)); } fclose(F); }
  if (qtype == "long")   { for (i=0;i<Qsize;i++) { fprintf(F,"%.10li\n",(long)getq(i)); } fclose(F); }

}

//****************************************************************************************************************
template <class T>  void cpeds_queue<T>::saveqidx(string f,string qtype) {
  FILE* F;
  long i;

  F=fopen(f.c_str(),"w");
  if (qtype == "double") { for (i=0;i<Qsize;i++) { fprintf(F,"%li %lE\n",i,(double)getq(i)); } fclose(F); }
  if (qtype == "long")   { for (i=0;i<Qsize;i++) { fprintf(F,"%li %li\n",i,(long)getq(i)); } fclose(F); }

}
//****************************************************************************************************************
template <class T>  long cpeds_queue<T>::loadq(string f,string qtype) {
  FILE* F;
  double tmpd;
  long tmpl;

  F=fopen(f.c_str(),"r");
  if ( F == NULL) { printf("ERROR: no file %s\n",f.c_str()); return -1; }
  else {
    if (qtype == "double") { while ( fscanf(F,"%lE",&tmpd) != EOF ) { addq(tmpd); } fclose(F); }
    if (qtype == "long")   { while ( fscanf(F,"%li",&tmpl) != EOF ) { addq(tmpl); } fclose(F); }
  }
  return 0;
}

//****************************************************************************************************************
template <class T>  void cpeds_queue<T>::insertq_ord(T v,long ord) {
  long i,addi=0;

  if (Qsize>0) {
    i=0;
    if (ord==12)  {  while (i<Qsize) { if (getq(i) < v) {  i++; addi=i; } else { addi=i;  i=Qsize;} } }
    else  { while (i<Qsize) { if (getq(i) > v) { i++; addi=i; } else { addi=i;  i=Qsize;} } }

    insertq(v,addi); 
/*     printf("inserting at position %li\n",addi); */
    
/*     printq_double(); */
  }
  else addq(v);
}

//****************************************************************************************************************
// insets a value at the i'th index of the queue pushing all larger  indexes (starting from the current i'th element ) up
template <class T>  void cpeds_queue<T>::insertq(T v, long I) {
  queue_type *tmp = new queue_type;
  queue_type *curr = getQ(I);
  queue_type *prev = getQ(I-1);

  if (I>=Qsize) addq(v); 
  else {
    tmp->n=curr;
    tmp->v=v;
    if (prev!=NULL) prev->n=tmp;
    Qsize++;
    if (I==0) Q=tmp;
    if (set) {
      if (v < Qmin) { Qmin=v; iQmin=Qsize-1; }
      if (v > Qmax) { Qmax=v; iQmax=Qsize-1; }
    }
  }
}

//****************************************************************************************************************
template <class T> long cpeds_queue<T>::count_common_vals(cpeds_queue<T> * q) {
    long i,comm=0;
    
    for (i=0;i<Qsize;i++) {
      if (q->val_in_list(getq(i))) comm++;
    }
    return comm;
  }

//****************************************************************************************************************
template <class T> bool cpeds_queue<T>::all_positive() {
  long i;
  for (i=0;i<Qsize;i++) { if (getq(i) <= 0) return false; }
  return true;
}
//****************************************************************************************************************
template <class T> bool cpeds_queue<T>::all_nonnegative() {
  long i;
  for (i=0;i<Qsize;i++) { if (getq(i) < 0) return false; }
  return true;
}

//****************************************************************************************************************
template <class T> long cpeds_queue<T>::get_first_nonpositive_idx() {
  long i=-1;
  for (i=0;i<Qsize;i++) { if (getq(i) <= 0) return i; }
  return i;
  
}
//****************************************************************************************************************
template <class T> long cpeds_queue<T>::get_first_negative_idx() {
  long i=-1;
  for (i=0;i<Qsize;i++) { if (getq(i) < 0) return i; }
  return i;
  
}

//****************************************************************************************************************
template <class T> void cpeds_queue<T>::zero_all_below(double val) {
  long i;
  for (i=0;i<Qsize;i++) { if (getq(i) < 0) setq(i,val); }
}


//****************************************************************************************************************
template <class T> long cpeds_queue<T>::find_closest(T val) {
  double *tmpd;
  tmpd=export_array();
  return cpeds_find_value(val,tmpd,Qsize,0,Qsize);
}

