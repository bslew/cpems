/*!
 * \file test-hdf-openMP.cpp
 *
 *  Project: Mscs
 *  Created on: Jan 31, 2016 2:38:34 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef ENABLE_OMP
#include <omp.h>
#endif
#include <string>
#include <hdf5.h>



#ifndef _NO_NAMESPACE
using namespace std;
#define STD std
#else
#define STD
#endif


bool fileExists(string fname);
int setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment);
int saveHDF5(string filename, string datasetName, double *data,int ix, int iy, int iz, int nx, int ny, int nz);
void setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue);


int main(int argc, char **argv) {
	
	long i,k,l=0;
	
	
#ifdef ENABLE_OMP
	omp_lock_t  _lock;
	omp_init_lock(&_lock);
#endif
	
	
#ifdef ENABLE_OMP
#pragma omp parallel for schedule(guided) default(none) shared(k,_lock,l) private(i)  num_threads(100)
#endif
	for (i = 0; i < 10000; ++i) {
		double data[1000];
		char dsetName[100];

		// a time consuming calculation
		for (unsigned long j = 0; j < 100000; j++) {
#pragma omp atomic
			k++;
		}
		
#ifdef ENABLE_OMP
		omp_set_lock(&_lock);
#endif
		sprintf(dsetName,"data_%li",i);
		printf("saving dset %li/1000\n",++l);
#ifdef ENABLE_OMP
		omp_unset_lock(&_lock);
#endif
		// save data set and scalar parameter to the same hdf file from different threads but in a serial manner
		saveHDF5("test.hdf5",dsetName,data,0,0,0,10,10,10);
		setHDF5_scalarStringAttribute("test.hdf5",dsetName,"dupa","kupa z lokiem");
		
	}
#ifdef ENABLE_OMP
	omp_destroy_lock(&_lock);
#endif
	
	return 0;
}


int saveHDF5(string filename, string datasetName, double *data,int ix, int iy, int iz, int nx, int ny, int nz) {
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (nx==0 or ny==0 or nz==0) {			return -1;		}
	
	int rank=3;
	hid_t file, space, dset, dcpl;    /* Handles */
	herr_t status;
	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
	double *wdata;          /* Write buffer */
	
	dims[0]=nx;
	dims[1]=ny;
	dims[2]=nz;
	chunk[0]=nx;
	chunk[1]=ny;
	chunk[2]=1;
	
	/*
	 * Create a new file using the default properties.
	 */
	if (fileExists(filename)) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		
	}
	else {
		file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);    	
	}
	
	/*
	 * Create dataspace.  Setting maximum size to NULL sets the maximum
	 * size to be the current size.
	 */
	space = H5Screate_simple(rank, dims, NULL);
	
	
	/*
	 * Create the dataset creation property list, and set the chunk
	 * size.
	 */
	dcpl = H5Pcreate (H5P_DATASET_CREATE);
	status = H5Pset_chunk(dcpl, rank, chunk);
	if (status!=0) { printf("non-zero HDF status code while setting chunk (%i)\n",status); }
	
	double defalutVal=1;
	H5Pset_fill_value(dcpl,H5T_IEEE_F64LE,&defalutVal);
	
	/*
	 * Create the chunked dataset.
	 */
	//		string linkName="/"+datasetName;
	//		if (hdf5DatasetExists(file,linkName))
	//			H5Ldelete(file,linkName.c_str(),H5P_DEFAULT);
	//		else {
	//			// create groups if needed
	//			QString qs=linkName.c_str();
	//			QStringList qsl=qs.split("/");
	//			if (qsl.size()==3)
	//				hdf5createGroup(file,linkName);
	//		}
	dset = H5Dcreate(file, datasetName.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
	if (dset<0) { printf("non-zero HDF status code while creating dataset (%s)\n",datasetName.c_str()); }
	
	/*
	 * Write the data to the dataset.
	 */
	long Nslice=nx*ny*nz;
	wdata = new double[Nslice];
	
	long ist,ien,jst,jen,kst,ken;
	ist=ix;
	jst=iy;
	kst=iz;
	ien=ist+nx;
	jen=jst+ny;
	ken=kst+nz;
	
	long l=0;
	
	// prepare slice data    	
	for (long i = ist; i < ien; i++) {
		for (long j = jst; j < jen; j++) {
			for (long k = kst; k < ken; k++) {
				wdata[l]=data[k+j*nz+i*nz*ny];
				l++;
			}
		}
	}
	
	/*
	 * Define and select the first part of the hyperslab selection.
	 */
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	block[0]=1;
	block[1]=1;
	block[2]=1;
	stride[0]=1;
	stride[1]=1;
	stride[2]=1;
	count[0]=nx;
	count[1]=ny;
	count[2]=nz;
	
	status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count, block);
	if (status!=0) { printf("non-zero HDF status code while selecting hyperslab (%i)\n",status); }
	status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, wdata);
	if (status!=0) { printf("non-zero HDF status code while writing (%i)\n",status); }
	delete [] wdata;
	
	setHDF5_scalarDoubleAttribute(dset,"x0",0.0, "function domain lower range in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"y0",0.0, "function domain lower range in y direction");
	setHDF5_scalarDoubleAttribute(dset,"z0",0.0, "function domain lower range in z direction");
	setHDF5_scalarDoubleAttribute(dset,"dx",1.0, "function cell size in x direction");
	setHDF5_scalarDoubleAttribute(dset,"dy",1.0, "function cell size in y direction");
	setHDF5_scalarDoubleAttribute(dset,"dz",1.0, "function cell size in z direction");	
	setHDF5_scalarDoubleAttribute(dset,"xMax",10.0, "function domain upper range in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"yMax",10.0, "function domain upper range in y direction");
	setHDF5_scalarDoubleAttribute(dset,"zMax",10.0, "function domain upper range in z direction");
	setHDF5_scalarDoubleAttribute(dset,"dxo2",0.5, "function domain cell center offset in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"dyo2",0.5, "function domain cell center offset in y direction");
	setHDF5_scalarDoubleAttribute(dset,"dzo2",0.5, "function domain cell center offset in z direction");
	
	/*
	 * Close and release resources.
	 */
	status = H5Pclose (dcpl);
	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Fclose (file);
	
	return 0;
}

int setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment) {
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	hid_t aid = H5Screate(H5S_SCALAR); // Create scalar attribute.
	hid_t attr = H5Acreate(dset, attributeName.c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
	hid_t status = H5Awrite(attr, H5T_IEEE_F64LE, &value); // Write scalar attribute.
	status = H5Sclose(aid); // Close attribute dataspace.
	status = H5Aclose(attr); //Close attribute.
	
	//
	// create comment string attribute
	//		
	
	herr_t ret2;
	hid_t atype2 = H5Tcopy(H5T_C_S1);
	string attributeCommentName = attributeName+"_comment";
	hid_t aid2 = H5Screate(H5S_SCALAR); // Create scalar attribute.
	ret2 = H5Tset_size(atype2, attributeComment.size());
	ret2 = H5Tset_strpad(atype2,H5T_STR_NULLTERM);
	hid_t attr2 = H5Acreate(dset, attributeCommentName.c_str(), atype2, aid2, H5P_DEFAULT, H5P_DEFAULT);
	ret2 = H5Awrite(attr2, atype2, attributeComment.c_str()); 
	
	status = H5Tclose(atype2);
	status = H5Sclose(aid2); // Close attribute dataspace.
	status = H5Aclose(attr2); //Close attribute.        	
	return status;
}

bool fileExists(string fname) {
	if (fname=="") return false;
	FILE* f = fopen(fname.c_str(),"r");
	if (f==NULL) return false;
	fclose(f);
	return true;
}

/***************************************************************************************/
void setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue) {
	hid_t dset,file;
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (fileExists(fname)) {
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
		
		//
		// create comment string attribute
		//
		
		herr_t ret2;
		hid_t atype2 = H5Tcopy(H5T_C_S1);
		hid_t aid2 = H5Screate(H5S_SCALAR); // Create scalar attribute.
		ret2 = H5Tset_size(atype2, attributeValue.size());
		ret2 = H5Tset_strpad(atype2,H5T_STR_NULLTERM);
		hid_t attr2 = H5Acreate(dset, attributeName.c_str(), atype2, aid2, H5P_DEFAULT, H5P_DEFAULT);
		ret2 = H5Awrite(attr2, atype2, attributeValue.c_str()); 
		ret2 = H5Sclose(aid2); // Close attribute dataspace.
		ret2 = H5Aclose(attr2); //Close attribute.        	
		ret2 = H5Tclose(atype2);
		
		ret2 = H5Dclose (dset);
		ret2 = H5Fclose (file);
		
	}
}
