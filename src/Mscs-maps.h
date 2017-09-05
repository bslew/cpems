/* this is an attempt to create a list of Mscs_map class pointers but */
/* much more resonable idea would seem to be just to use the existing */
/* class for operating on lists of pointers to objects of the  */
/* QPtrList from the Qt library. (include "qptrlist.h") */

#include "Mscs-map.h"

class Mscs_maps {

 public:

  typedef struct {
    map_class * MP;
    map_class ** next;
  } MscsQtype;


  MscsQtype mq;
  MscsQtype mqBeg,mqEnd;

  map_class* operator () (long i) { return getq(i); }
  map_class* operator () (map_class*, long i) { return setq(i); }


 private:
};

map_class* Mscs_maps::getq(long i) {  // returns the I'th value in the queue; numbering starts from 0
  long i;
  queue_type* tmp=Q;
  
  if (I >= Qsize || I<0) { return 0; }
  for (i=0;i<I;i++) { tmp=tmp->n; }
  return tmp->v;
}


