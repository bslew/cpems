#include <stdio.h>
#include "Mscs-alm.h"

void mscsAlm::print() {


}

mscsAlm::mscs_a_lm mscsAlm::get_a_lm() const { //!< returns the alm in a C form
	mscs_a_lm a;
	a.R=R();
	a.I=I();
	return a;
}

ostream& operator<<(ostream& out, mscsAlm& a) { 
  out.setf ( ios::scientific  );  
//  return out<<a.R()<<" "<<a.I(); 
  return out<<a.R()<<a.I(); 
}
//istream& mscsAlm::operator<<(istream& in, mscsAlm& a) { 
//  in.setf ( ios::scientific  );  
//  return in>>a.R()<<a.I(); 
//}
