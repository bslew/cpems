/*******************************************************************************
 *  ccSHT:                                                                      *
 *    This is a set of functions which will perform a spherical harmonic        *
 *    transform.  The principal functions  are the following:                   *
 *      forwardSHT():                                                           *
 *        Forward spherical harmonic transform.                                 *
 *      backwardSHT():                                                          *
 *        Backward spherical harmonic transform.                                *
 *      calculateQlm():                                                         *
 *        a function which will calculate normalized associated Legendre        *
 *        polynomials Qlm, where Qlm = sqrt((l-m)!/(l+m)!)*Plm                  *
 *        Note that this function scales the recurrence so that it is accurate  *
 *        for a wide range of input parameters                                  *
 *                                                                              *
 *     See ccSHT(1) man page for documentation.                                 *
 *  C.M. Cantalupo  10/29/01 last updated 7/11/03                               *
 *******************************************************************************/


/*******************************************************************************
 *   Version 1.03 July 2003                                                     *
 *                                                                              *
 *   Copyright (C) 2003  C.M. Cantalupo                                         *
 *                                                                              *
 *   ccSHT is free software; you can redistribute it and/or modify              *
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

#ifndef _H_ccSHT
#define _H_ccSHT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "fftw3.h"
#include "generalTools.h"
#include "fftw_complex_helper.h"

/*******************************************************************************
 *  Qindex is the indexing macro for the scaled associated Legendre             *
 *  polynomials.  This is an l major index which increases with both l and m.   *
 *******************************************************************************/
#define Qindex(l,m,L) ( (L)*(m) - ((m)*((m)-1))/2 + (l) )

/*******************************************************************************
 *  almIndex is the indexing macro for the spherical harmonic coefficients.     *
 *  This is an l major index which increases with both l and m.                 *
 *******************************************************************************/
#define almIndex(l,m,L) ( (((L)*((L) + 2*(m) + 1)) + (m)*(2 - ( ((m) >= -1) ? ((m) + 1) : (-((m) + 1)) ) ))/2 + (l))

/* pown1 is a macro which returns negative one raised to the integer power a. */
#define pown1(a) (((a)%2) ? (-1) : (1))

/* Small number with respect to pixel separation */
#define MY_EPSILON 1e-6


/*******************************************************************************
 *  Name of directory where temporary calculations are stored if the directory  *
 *  exists.                                                                     *
 *******************************************************************************/
//#define FILE_DIR "ccSHT_files"
#define FILE_DIR "/home/blew/programy/ccSHT3/data"

/*******************************************************************************
 *  coordStruct is the data structure for the pixelization information.  See    *
 *  man page for ccSHT  for detailed information.                               *
 *******************************************************************************/
typedef struct 
{
		flt8 pixSize;
		long8 nPix;  
		long8 nThetaVals;
		flt8 *thetaVals;
		long8 *thetaBreaks;
		flt8 *phi0;
		flt8 *deltaPhi;
		long8 nGaps;
		long8 *gaps;
} coordStruct;

/* Library functions */
void forwardSHT( void *map, long8 inputIsComplex, coordStruct coords, long8 lmax, SHT_FFTW_COMPLEX *alm);
void backwardSHT( SHT_FFTW_COMPLEX *alm, coordStruct coords, long8 lmax, void *map, long8 outputIsComplex);
void calculateQlm(long8 lmax, flt8 x, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4, flt8 *Q);
void lmCalculations(long8 lmax, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4);
//coordStruct convertRaDec2Coords(flt8 pixSize, long8 nPix, flt8 *raDecArray, long8 ndra, flt8 *dra);
//void destroyCoords(coordStruct *coords);
void destroyCoordsCPP(coordStruct *coords);
void calculateSmallQlm( long8 lmax, flt8 x, flt8 *calc1, flt8 *calc2, flt8 *calc3, flt8 *calc4, flt8 *Qstart, long8 *lStart);

/* Internal functions */
long* createPlanSizes(coordStruct coords);
//void destroyPlans(coordStruct coords, fftw_plan **plans);


#endif
