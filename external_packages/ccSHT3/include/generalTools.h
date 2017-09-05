/*******************************************************************************
 *  generalTools:                                                               *
 *    A set of miscellaneous useful functions:                                  *
 *      errorCheck()                                                            *
 *        Altered version of Julian Borrill's error checking function.  If it   *
 *        is desired that errorCheck call MPI_Finalize() when exiting           *
 *        compile with "-D USE_MPI"                                             *
 *      input_ftod()                                                            *
 *        Julian Borrill's array file reading function.                         *
 *      output_dtof()                                                           *
 *        Julian Borrill's array file writing function.                         *
 *      readTextFlt8()                                                          *
 *        Reads in ASCII text file of flt4ing polong8 numbers.                   *
 *      ccSHT_round()                                                           *
 *        The rounding function (takes a flt8 returns an long8)                 *
 *      findSmallestEl()                                                        *
 *        Returns the smallest element of an array.                             *
 *      myMod()                                                                 *
 *        The way that I think the modulus function should work                 *
 *        (range 0 to n-1).                                                     *
 *                                                                              *
 *  C.M. Cantalupo  10/29/01 last updated 7/11/03                               *
 *******************************************************************************/

/*******************************************************************************
 *   Version 1.03 July 2003                                                     *
 *                                                                              *
 *   Copyright (C) 2003  C.M. Cantalupo                                         *
 *                                                                              *
 *   This program is free software; you can redistribute it and/or modify       *
 *   it under the terms of the GNU General Public License as published by       *
 *   the Free Software Foundation; either version 2 of the License, or          *
 *   (at your option) any later version.                                        *
 *                                                                              *
 *   This program is distributed in the hope that it will be useful,            *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 *   GNU General Public License for more details.                               *
 *                                                                              *
 *   You should have received a copy of the GNU General Public License          *
 *   along8 with this program; if not, write to the Free Software                *
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA  *
 *                                                                              *
 *******************************************************************************/

#ifndef _H_generalTools
#define _H_generalTools

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <math.h>
#include "sizedType.h"

/*******************************************************************************
 *  If it is desired that errorCheck call MPI_Finalize() when exiting           *
 *  compile with "-D USE_MPI"                                                   *
 *******************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif



#define STRLEN 128
#define minimum(a,b) (((a) < (b)) ? (a) : (b))
#define maximum(a,b) (((a) > (b)) ? (a) : (b))
#define PI 3.14159265358979

typedef struct {
		flt8 re,im;
} fftw_complex2;


void errorCheck(int my_pe, const char * op_string, long8 ok, long8 errorFlag);
void input_ftod(FILE *f, flt8 *data, long8 no_data);
void output_dtof(FILE *f, flt8 *data, long8 no_data);
//void readTextFlt8(char *fileName,long8 n, flt8 *output);
long8 ccSHT_round(flt8 a);
flt8 findSmallestEl( flt8 *theVect, long8 nEl);
long8 myMod( long8 a, long8 b);

#endif
