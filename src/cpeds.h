/* standard library for CPEDS routines 
this library contains routines that are required by cpeds main program to work however are not very usefull 
in general separate, external use in any other possible applications
*/

//#include "cpeds-defs.h"



//**************************************************************************************
// funkcja zwracajaca skok dz do obliczenia calki ktory zalezy do akualnego miejsca calkowania
// i dokladnosci calkowania

//extern double ZINF;
//extern double pow(double x, double y);

double skok(double z) {
  if (z > 0)   return (ZINF/pow( 10,floor(log10(ZINF/z)) ))/NP; 
  if (z <= 0)  return (1.0/NF);
  
}

//**************************************************************************************

