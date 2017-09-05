/*!
 * \file lssKitFileHandler.h
 *
 *  Created on: Mar 22, 2012
 *      Author: blew
 */

#ifndef LSSKITFILEHANDLER_H_
#define LSSKITFILEHANDLER_H_

#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <string.h>
#include  "cpeds-msgs.h"

#define lssKitPyramidSliceFileHeaderSize 500 // defines the size of the header in bytes

/*!
	\brief defines a header of the binary data file format used to store slices of the lss pyramid
	\details 
	The file format consists of header which is of the fixed size defined by the lssKitPyramidSliceFileHeaderSize

	In that area the below defined structure is stored. The remaining space is not used.
	The second block contains the actual data consisting of Npart quadruples of float 4-byte numbers.
	These are typically x,y,z, and Trho -- density weighted temperature of SPH gas particle

	Currently the actual size of the header is 176 bytes.

	UPDATE
		Feb 28, 2013, 6:28:54 PM - added hubble param to the header file. Now the header size is 176+8 = 184 bytes
	
	\date Mar 2, 2012, 1:18:08 PM
	\author Bartosz Lew
*/
typedef struct {
		long Npart; //!< number of particles in the file
		long scalarColNum; //!< number of scalar columns per particle (default:1 = original 3NN format, 2 = 3NNN format)
		double Sx,Sy,Sz,Sdx,Sdy,Sdz; //!< defines the slice in simulation box coordinates in Mpc
		double SX,SY,SZ; //!< defines the slice in pyramid coordinates in Mpc with Z=0 at pyramid apex 
		double Dx,Dy,Dz,Ddx,Ddy,Ddz; //!< defines the domain size - this is the definition of the actual space ocupied by the data in the file [Mpc]
		double simBoxZ0; //!< deep field z-coordinate of the closer face of the simulation box wrt which the z-coordinates in the file are calculated [Mpc]
		double simID; //!< simulation ID for the simBoxZ0
		long SNx,SNy,SNz; //!< number of pixels in the slice
		long DNx,DNy,DNz; //!< number of pixels in the domain
		double Wb0,mu_e; //!< baryon density today in critical density units (today) and electron number per proton mass (typically mu_e=1.136 for Y=0.24)
		double y0; //!< integral constant: y0 = sigma_T kB rhoC0(h) Wb0_WMAP / ( mu_e_BBN m_p * m_e c^2 ) * CPEDS_MPC  = 1.711*10^-16*h^2 [ K^-1 Mpc^-1 ], so this constant includes the hubble parameter h that is also specified in this header 
		double z; //!< redshift
		double h; //!< unitless hubble param
} lssPyramidSliceFileHeader_t;


typedef struct {
		float x,y,z,T;
} lssPyramidSliceFileQuadruplet_t;
typedef struct {
		float x,y,z,T,hsml;
} lssPyramidSliceFileQuintuple_t;



float* lssKit3NtoN(float* data3N, long sizeN, int coord) {
	float* v=new float[sizeN];
	for (long i = 0; i < sizeN; i++) {
		v[i]=data3N[3*i+coord];
	}
	return v;
}


/*!
	\brief saves eg. positions and a eg temperatures to a binary file
	\details 
	@param data1 - vector of size 3N
	@param data2 - vector of size N
	@param size - N
	@param append - if true then if the fname file exists then the new data will be appended to the end of the file and the header information will be updated
	
	The format of the output file is as follows
	header block consisting of 100 chars
	data block consisting of quadrupoles of float numbers corresponding to eg. x,y,z,T
	
	In the header block the first 8 bytes is a type long number containing the total number of quadruples in the file
	All remaining chars in the header are not used

	\date Mar 2, 2012, 9:35:48 AM
	\author Bartosz Lew
*/
void lssKitSave3NNdata(string fname, float* data1, float* data2, lssPyramidSliceFileHeader_t h, cpedsMsgs* msgs, bool append) {
	long size=h.Npart;
	char head[lssKitPyramidSliceFileHeaderSize];
	bzero(head,lssKitPyramidSliceFileHeaderSize*sizeof(char));
	lssPyramidSliceFileHeader_t header,headerOld;
	string fmt;
	h.scalarColNum=4;
	
	
	if (append) {
		fmt="rb";
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { 
			header=h;
			header.Npart=0;
			msgs->criticalError(string("lssKitSave3NNdata> Appending to 3NN file which doesn't exist (")+fname+"). This should not happen. Please check the code. I will stop now.",Top);
		}
		else {
			fread(&headerOld,sizeof(lssPyramidSliceFileHeader_t),1,f);
			header=headerOld;
			printf("previous particles count: %li\n",headerOld.Npart);
			fclose(f);			
		}
		header.Npart+=size;
		printf("current particles count: %li\n",header.Npart);
		printf("the domain sizes will not be altered - the particles from the neighboring domain should be inside of the the domain centered on the slice.\n");
		printf("the redshift will not be altered for the slice due to adding particles from a different epoch. This is good, because we will calculate"
				"interpolated values for the redshift of the slice and the particles from neighboring domains should only improve the smoothness "
				"of the field and provide neighborhood mass distribution info.\n");
		
//		printf("checking domain ranges\n");
//		float *v;
//		float xmin,xmax,ymin,ymax,zmin,zmax;
//		v=lssKit3NtoN(data1,size,0);
//		cpeds_find_minmax_value(v,size,&xmin,&xmax,NULL,NULL);
//		delete [] v;
//		v=lssKit3NtoN(data1,size,1);
//		cpeds_find_minmax_value(v,size,&ymin,&ymax,NULL,NULL);
//		delete [] v;
//		v=lssKit3NtoN(data1,size,2);
//		cpeds_find_minmax_value(v,size,&zmin,&zmax,NULL,NULL);
//		delete [] v;
//		if (h.Sx < headerOld.Sx) {
//			header.Sx=h.Sx;
//			printf("changing slice Sx range from %lf to %lf\n",headerOld.Sx,header.Sx);
//		}
//		if (h.Sy < headerOld.Sy) {
//			header.Sy=h.Sy;
//			printf("changing slice Sy range from %lf to %lf\n",headerOld.Sy,header.Sy);
//		}
//		if (h.Sz < headerOld.Sz) {
//			header.Sz=h.Sz;
//			printf("changing slice Sz range from %lf to %lf\n",headerOld.Sz,header.Sz);
//		}
		
		
		//
		// saving header
		//
		
		fmt="r+b";
		f=fopen(fname.c_str(),fmt.c_str());
		fseek(f,0,SEEK_SET);
		memcpy(head,&header,sizeof(lssPyramidSliceFileHeader_t));
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		fseek(f,0,SEEK_END);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
		}
		fclose(f);
	}
	else {
		fmt="wb";
		memcpy(head,&h,sizeof(lssPyramidSliceFileHeader_t));
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { msgs->criticalError("could not write to a file",High); }
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
		}
		fclose(f);
	}
}

/*!
	\brief saves eg. positions and a eg temperatures and eg hsml to a binary file
	\details 
	@param data1 - vector of size 3N
	@param data2 - vector of size N
	@param data3 - vector of size N
	@param size - N
	@param append - if true then if the fname file exists then the new data will be appended to the end of the file and the header information will be updated
	
	The format of the output file is as follows
	header block consisting of 100 chars
	data block consisting of quintupoles of float numbers corresponding to eg. x,y,z,T,hsml
	
	In the header block the first 8 bytes is a type long number containing the total number of quadruples in the file
	All remaining chars in the header are not used
	
	This is a stupid temporary implementation until a general save/load functions are implemented for arbitrary number of additional scalar
	columns.

	\date Mar 2, 2012, 9:35:48 AM
	\author Bartosz Lew
*/
void lssKitSave3NNNdata(string fname, float* data1, float* data2, float* data3, lssPyramidSliceFileHeader_t h, cpedsMsgs* msgs, bool append) {
	long size=h.Npart;
	char head[lssKitPyramidSliceFileHeaderSize];
	bzero(head,lssKitPyramidSliceFileHeaderSize*sizeof(char));
	lssPyramidSliceFileHeader_t header,headerOld;
	string fmt;
	h.scalarColNum=5;
	
	
	if (append) {
		fmt="rb";
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { 
			header=h;
			header.Npart=0;
			msgs->criticalError(string("lssKitSave3NNdata> Appending to 3NN file which doesn't exist (")+fname+"). This should not happen. Please check the code. I will stop now.",Top);
		}
		else {
			fread(&headerOld,sizeof(lssPyramidSliceFileHeader_t),1,f);
			header=headerOld;
			printf("previous particles count: %li\n",headerOld.Npart);
			fclose(f);			
		}
		header.Npart+=size;
		printf("current particles count: %li\n",header.Npart);
		printf("the domain sizes will not be altered - the particles from the neighboring domain should be inside of the the domain centered on the slice.\n");
		printf("the redshift will not be altered for the slice due to adding particles from a different epoch. This is good, because we will calculate"
				"interpolated values for the redshift of the slice and the particles from neighboring domains should only improve the smoothness "
				"of the field and provide neighborhood mass distribution info.\n");
		
		//
		// saving header
		//
		
		fmt="r+b";
		f=fopen(fname.c_str(),fmt.c_str());
		fseek(f,0,SEEK_SET);
		memcpy(head,&header,sizeof(lssPyramidSliceFileHeader_t));
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		fseek(f,0,SEEK_END);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
			fwrite(&data3[i],sizeof(float),1,f);
		}
		fclose(f);
	}
	else {
		fmt="wb";
		memcpy(head,&h,sizeof(lssPyramidSliceFileHeader_t));
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { msgs->criticalError("could not write to a file",High); }
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
			fwrite(&data3[i],sizeof(float),1,f);
		}
		fclose(f);
	}
}

void lssKitSave3NNNdata(string fname, float* data1, float* data2, float* data3, float* data4, lssPyramidSliceFileHeader_t h, cpedsMsgs* msgs, bool append) {
	long size=h.Npart;
	char head[lssKitPyramidSliceFileHeaderSize];
	bzero(head,lssKitPyramidSliceFileHeaderSize*sizeof(char));
	lssPyramidSliceFileHeader_t header,headerOld;
	string fmt;
	h.scalarColNum=6;
	
	
	if (append) {
		fmt="rb";
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { 
			header=h;
			header.Npart=0;
			msgs->criticalError(string("lssKitSave3NNdata> Appending to 3NN file which doesn't exist (")+fname+"). This should not happen. Please check the code. I will stop now.",Top);
		}
		else {
			fread(&headerOld,sizeof(lssPyramidSliceFileHeader_t),1,f);
			header=headerOld;
			printf("previous particles count: %li\n",headerOld.Npart);
			fclose(f);			
		}
		header.Npart+=size;
		printf("current particles count: %li\n",header.Npart);
		printf("the domain sizes will not be altered - the particles from the neighboring domain should be inside of the the domain centered on the slice.\n");
		printf("the redshift will not be altered for the slice due to adding particles from a different epoch. This is good, because we will calculate"
				"interpolated values for the redshift of the slice and the particles from neighboring domains should only improve the smoothness "
				"of the field and provide neighborhood mass distribution info.\n");
		
		//
		// saving header
		//
		
		fmt="r+b";
		f=fopen(fname.c_str(),fmt.c_str());
		fseek(f,0,SEEK_SET);
		memcpy(head,&header,sizeof(lssPyramidSliceFileHeader_t));
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		fseek(f,0,SEEK_END);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
			fwrite(&data3[i],sizeof(float),1,f);
			fwrite(&data4[i],sizeof(float),1,f);
		}
		fclose(f);
	}
	else {
		fmt="wb";
		memcpy(head,&h,sizeof(lssPyramidSliceFileHeader_t));
		FILE* f=fopen(fname.c_str(),fmt.c_str());
		if (f==NULL) { msgs->criticalError("could not write to a file",High); }
		fwrite(head,lssKitPyramidSliceFileHeaderSize*sizeof(char),1,f);
		for (long i = 0; i < size; i++) {
			fwrite(&data1[3*i],sizeof(float),3,f);
			fwrite(&data2[i],sizeof(float),1,f);
			fwrite(&data3[i],sizeof(float),1,f);
			fwrite(&data4[i],sizeof(float),1,f);
		}
		fclose(f);
	}
}

void lssKitPrint3NNFileHeader(const lssPyramidSliceFileHeader_t& h) {
	printf("lssPyramidSliceFileHeader:\n");
	printf("Npart: %li\n",h.Npart);
	printf("Number of scalar columns: %li\n",h.scalarColNum);
	printf("Slice information\n");	
	printf("Sx: %lf [Mpc]\n",h.Sx);
	printf("Sy: %lf [Mpc]\n",h.Sy);
	printf("Sz: %lf [Mpc]\n",h.Sz);
	printf("Sdx: %lf [Mpc]\n",h.Sdx);
	printf("Sdy: %lf [Mpc]\n",h.Sdy);
	printf("Sdz: %lf [Mpc]\n",h.Sdz);
	printf("\n");	
	printf("SX: %lf [Mpc]\n",h.SX);
	printf("SY: %lf [Mpc]\n",h.SY);
	printf("SZ: %lf [Mpc]\n",h.SZ);
	printf("\n");	
	printf("Domain information\n");	
	printf("Dx: %lf [Mpc]\n",h.Dx);
	printf("Dy: %lf [Mpc]\n",h.Dy);
	printf("Dz: %lf [Mpc]\n",h.Dz);
	printf("Ddx: %lf [Mpc]\n",h.Ddx);
	printf("Ddy: %lf [Mpc]\n",h.Ddy);
	printf("Ddz: %lf [Mpc]\n",h.Ddz);
	printf("\n");	
	printf("Pixels information\n");	
	printf("SNx: %li \n",h.SNx);
	printf("SNy: %li \n",h.SNy);
	printf("SNz: %li \n",h.SNz);
	printf("DNx: %li \n",h.DNx);
	printf("DNy: %li \n",h.DNy);
	printf("DNz: %li \n",h.DNz);
	printf("\n");	
	printf("Simulation box information\n");	
	printf("simBoxZ0: %lf [Mpc]\n",h.simBoxZ0);
	printf("simBoxID: %.0lf\n",h.simID);
	printf("\n");	
	printf("Cosmology information\n");	
	printf("mu_e: %lE\n",h.mu_e); //  - number of electrons per proton mass
	printf("Wb0: %lE \n",h.Wb0);
	printf("y0 [K^-1 Mpc^-1]: %lE\n",h.y0); // - comptonization parameter integral constant
	printf("z: %lE\n",h.z);
	printf("h: %lE\n",h.h);
	
	
}

lssPyramidSliceFileHeader_t lssKitPrint3NNData(string fname, cpedsMsgs* msgs, bool printAll=false) {
	lssPyramidSliceFileHeader_t h;
	
	FILE* f=fopen(fname.c_str(),"rb");
	if (f==NULL) { msgs->criticalError("could not open file:"+fname,High); }
	fread(&h,sizeof(lssPyramidSliceFileHeader_t),1,f);
	lssKitPrint3NNFileHeader(h);
	fseek(f,lssKitPyramidSliceFileHeaderSize,SEEK_SET);
	if (printAll) {
		size_t Ncols=h.scalarColNum;
		size_t pSize=sizeof(float)*Ncols;
		float data[Ncols];
		for (long i = 0; i < h.Npart; i++) {
			fread(data,pSize,1,f);
			for (long j = 0; j < Ncols; j++) {
				printf("%lE ",data[j]);
			}
			printf("\n");
		}
	}
	
	fclose(f);
	return h;
}

lssPyramidSliceFileQuadruplet_t* lssKitLoad3NNData(string fname,  long* N) {
	lssPyramidSliceFileHeader_t h;
	lssPyramidSliceFileQuadruplet_t *data;
	
	FILE* f=fopen(fname.c_str(),"rb");
	fread(&h,sizeof(lssPyramidSliceFileHeader_t),1,f);
//	lssKitPrint3NNFileHeader(h);
	fseek(f,lssKitPyramidSliceFileHeaderSize,SEEK_SET);
	*N=h.Npart;
	data = new lssPyramidSliceFileQuadruplet_t[*N];
//	for (long i = 0; i < h.Npart; i++) {
	fread(data,sizeof(lssPyramidSliceFileQuadruplet_t),*N,f);
//	}
	
	fclose(f);
	return data;
}


float* lssKitLoad3NNNData(string fname,  lssPyramidSliceFileHeader_t* H) {
	lssPyramidSliceFileHeader_t h;
	float *data=NULL;
	
	FILE* f=fopen(fname.c_str(),"rb");
	fread(&h,sizeof(lssPyramidSliceFileHeader_t),1,f);
	fseek(f,lssKitPyramidSliceFileHeaderSize,SEEK_SET);
	size_t pSize=sizeof(float)*(h.scalarColNum);
	*H=h;
	if (h.Npart==0) return NULL;
	data = new float[h.Npart*(h.scalarColNum)];
//	for (long i = 0; i < h.Npart; i++) {
	fread(data,pSize,h.Npart,f);
//	}
	
	fclose(f);
	return data;
}


lssPyramidSliceFileHeader_t lssKitLoad3NNData(string fname,  lssPyramidSliceFileQuadruplet_t** data=NULL, long* N=NULL) {
	lssPyramidSliceFileHeader_t h;
	
	FILE* f=fopen(fname.c_str(),"rb");
	if (f!=NULL) {
		fread(&h,sizeof(lssPyramidSliceFileHeader_t),1,f);
		
		if (data!=NULL) {
			fseek(f,lssKitPyramidSliceFileHeaderSize,SEEK_SET);
			*N=h.Npart;
			printf("Npart: %li\n",*N);
			*data = new lssPyramidSliceFileQuadruplet_t[*N];
			printf("memory allocated\n");
//			for (long i = 0; i < h.Npart; i++) {
//				printf("reading particle: %li\n",i);
				fread(*data,sizeof(lssPyramidSliceFileQuadruplet_t),*N,f);
//			}
		}

		fclose(f);
	}
	else {
		printf("lssKitLoad3NNData>> no such file %s",fname.c_str());
		exit(0);
	}
	
	return h;
}


#endif /* LSSKITFILEHANDLER_H_ */
