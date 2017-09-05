/* File created automatically by makeSizedType */

#ifndef _H_sizedType
#define _H_sizedType


#define long4 long long int
#define long8 long long int
#define long4mpiType MPI_INT
#define long8mpiType MPI_LONG
#define flt4mpiType MPI_FLOAT
#define flt8mpiType MPI_DOUBLE

// definitions for long double
#define SHT_FFTW_COMPLEX fftwl_complex
#define SHT_FFTW_PLAN_DFT_1D fftwl_plan_dft_1d
#define SHT_FFTW_PLAN fftwl_plan
#define SHT_FFTW_EXECUTE fftwl_execute
#define SHT_FFTW_DESTROY_PLAN fftwl_destroy_plan
#define flt4 long double
#define flt8 long double
#define SQRT sqrtl
#define LDEXP(a,b) ( ldexpl(a,(int)b) )
//#define FABS fabsl
//#define ABS llabs
#define FREXP(a,b) ( frexpl(a,(int*)b) )
#define SIN sinl
#define COS cosl
#define ATAN2 atan2l
// Largest and samllest exponent in flt16 precision IEEE standard. 
#define MIN_EXP -16382
#define MAX_EXP 16383

// definitions for double
//#define SHT_FFTW_COMPLEX fftw_complex
//#define SHT_FFTW_PLAN_DFT_1D fftw_plan_dft_1d
//#define SHT_FFTW_PLAN fftw_plan
//#define SHT_FFTW_EXECUTE fftw_execute
//#define SHT_FFTW_DESTROY_PLAN fftw_destroy_plan
//#define flt4 double
//#define flt8 double
//#define SQRT sqrt
//#define LDEXP ldexp
#define FABS fabs
#define ABS abs
//#define FREXP frexp
//#define SIN sin
//#define COS cos
//#define ATAN2 atan2
//// Largest and samllest exponent in flt8 precision IEEE standard. 
//#define MIN_EXP -1022
//#define MAX_EXP 1023


#endif


