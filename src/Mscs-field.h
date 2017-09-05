#define _ISOC99_SOURCE
#MSCS_FIELD_MAX_SIZE 1024 // maximal size of the field that can be stored in the memory beyond which it will have to be read from disk for any operations
#include <features.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

class field_class {

 public:

  // data types definitions

  

  // data structures definitions

  long field_size; // this assumes cubed shaped field 
  bool kspace; // this flag indicates in which space the field is generated: true=k-space, false=real-space
  

  // methods definitions


 private:

}
