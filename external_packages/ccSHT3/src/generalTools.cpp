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
#include <string.h>
#include "generalTools.h"

void errorCheck(int my_pe, const char * op_string, long8 ok, long8 errorFlag)
{ 
	/*******************************************************************************
	 * errorCheck:                                                                  *
	 *                                                                              *
	 *   This function has been adapted from error_check in Julian Borrill's        *
	 *   MADCAP tools library madtools.c.  errorCheck is a error function with      *
	 *   parallel capabilities.  The function takes as input:                       *
	 *                                                                              *
	 *   my_pe:  This is the rank of the processor who calls the function.  If the  *
	 *           function which calls this is serial then -1 should be passed.      *
	 *   op_string:  This string is used in the error message if ok (see below) is  *
	 *               non-zero.                                                      *
	 *   ok:  The value being tested.  If ok is non-zero then the function          *
	 *        displays a message and exits.                                         *
	 *   errorFlag:  A flag which determines the error message to be displayed and  *
	 *               how op_string is to be used as follows:                        *
	 *               errorFlag == 0                                                 *
	 *                 This is a general error.  The text in op_string will be      *
	 *                 displayed.                                                   *
	 *               errorFlag == 1                                                 *
	 *                 This is a memory allocation error.  op_string should be the  *
	 *                 name of the function in which the allocation failed.         *
	 *               errorFlag == 2                                                 *
	 *                 This is a file opening error when the file was to be read.   *
	 *                 op_string should be the name of the file which could not     *
	 *                 be opened.                                                   *
	 *               errorFlag == 3                                                 *
	 *                 This is a file opening error when the file was to be         *
	 *                 written.  op_string should be the name of the file which     *
	 *                 could not be opened.                                         *
	 *                                                                              *
	 *   Note that this function can be either a parallel or a serial error         *
	 *   function.  See generalTools.h for information about USE_MPI.               *
	 *                                                                              *
	 *******************************************************************************/
	
	if (ok==0) 
	{
		if( my_pe != -1 )
			fprintf(stderr, "ERROR: PE %d - error code %d\n", my_pe, errno);
		else
			fprintf(stderr, "ERROR: error code %d\n", errno);
		switch(errorFlag)
		{
			case 1:
				fprintf(stderr, "     Could not allocate memory in function called %s, exiting.\n\n", op_string);
				break;
			case 2:
				fprintf(stderr,"     Could not open file called %s for reading, exiting.\n\n", op_string);
				break;
			case 3:
				fprintf(stderr,"     Could not open file called %s for writing, exiting. \n\n", op_string);
				break;
			default:
				fprintf(stderr, "     %s\n\n", op_string);
				break;
		}
		
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, errno);
#endif
		exit(errno);
	}
}


void input_ftod(FILE *f, flt8 *data, long8 no_data)
{
	/*******************************************************************************
	 * input_ftod:                                                                  *
	 *                                                                              *
	 *   This is a function written by Julian Borrill which will read in a four     *
	 *   byte floating point binary file and promote it to an eight byte vector.    *
	 *   This function takes as input:                                              *
	 *                                                                              *
	 *   f:  A polong8er to a file which contains the vector to be read, and is       *
	 *       positioned at the beginning of the vector.                             *
	 *   data:  The output vector.                                                  *
	 *   no_data:  The number of elements to be read.                               *
	 *                                                                              *
	 *******************************************************************************/
	
	long8 n;
	flt4 *fptr;
	flt8 *dptr;
	
	/* Read in flt4 data */
	
	fread(data, 4, no_data, f);
	
	/* Convert flt4 data to flt8 data */
	
	for (n=no_data, dptr=data+n-1, fptr=(flt4 *)data+n-1; n>0; n--, dptr--, fptr--) *dptr = (flt8)(*fptr);
}


void output_dtof(FILE *f, flt8 *data, long8 no_data)
{
	/*******************************************************************************
	 * output_dtof:                                                                 *
	 *                                                                              *
	 *   This is a function written by Julian Borrill which writes an array of      *
	 *   eight byte flt4ing polong8 numbers to disk as four byte binary numbers.     *
	 *   The input data will be truncated to four byte precision when the function  *
	 *   returns.  This function takes as input:                                    *
	 *                                                                              *
	 *   f:  A polong8er to a file where the vector will be written.                  *
	 *   data:  The input vector of eight byte numbers.                             *
	 *   no_data:  The number of elements to be written to disk.                    *
	 *                                                                              *
	 *******************************************************************************/
	
	long8 n;
	flt4 *fptr;
	flt8 *dptr;
	
	/* Convert flt8 data to flt4 data */
	
	for (n=0, dptr=data, fptr=(flt4 *)data; n<no_data; n++, dptr++, fptr++) *fptr = (flt4)(*dptr);
	
	/* Write out flt4 data */
	
	fwrite(data, 4, no_data, f);
	fflush(f);
	
	/* Convert flt4 data back to flt8 data */
	
	for (n=no_data, dptr=data+n-1, fptr=(flt4 *)data+n-1; n>0; n--, dptr--, fptr--) *dptr = (flt8)(*fptr);
}


//void readTextFlt8(char *fileName, long8 n, flt8 *output)
//{
//	/*******************************************************************************
//	 * readTextFlt8:                                                                *
//	 *                                                                              *
//	 *   Reads in an ascii text file as a vector of floating point numbers.         *
//	 *   This function takes as input:                                              *
//	 *                                                                              *
//	 *   fileName:  A string which contains the name of the file to be read in.     *
//	 *   n:  The number of elements to be read in.                                  *
//	 *   output:  The vector which will be filled with the contents of the file.    *
//	 *                                                                              *
//	 *******************************************************************************/
//	
//	long8 i;
//	flt4 temp;
//	FILE *theFile;
//	
//	theFile = fopen(fileName,"r");
//	errorCheck(-1, fileName, (long8)theFile, 2);
//	
//	for( i = 0; i < n; i++)
//	{
//		if(fscanf(theFile, "%f", &temp) == EOF)
//		{
//			fprintf(stderr,"ERROR:  Fewer than expected entries in file named %s.\n",fileName);
//			fprintf(stderr,"        Zero padding the rest of the array.  \n\n");
//			
//			memset( output + i, 0, sizeof(flt8)*(n-i));
//			i = n;
//		}
//		else
//			output[i] = temp;
//	}
//	fclose(theFile);
//}
//

long8 ccSHT_round(flt8 a)
{
	if( a >= 0 )
		return((long8)(a+.5));
	return((long8)(a-.5));
}


flt8 findSmallestEl( flt8 *theVect, long8 nEl)
{
	long8 i;
	flt8 minEl;
	flt8 absMinEl;
	
	minEl = theVect[0];
	absMinEl = FABS(theVect[0]);
	for( i = 1; i < nEl; i++ )
		if(FABS(theVect[i]) < absMinEl)
		{
			absMinEl = FABS(theVect[i]);
			minEl = theVect[i];
		}
	return(minEl);
}

long8 myMod( long8 a, long8 b)
{
	long8 c;
	b = ABS(b);
	
	c = a%b;
	
	if( c < 0 )
		c += b;
	return(c);
} 

