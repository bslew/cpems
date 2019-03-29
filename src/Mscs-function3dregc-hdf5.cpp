/*!
 * \file Mscs-function3dregc-hdf5.cpp
 *
 *  This is an add-on to the Mscs-function3dregc class - compiled in optionally to support HDF5 files format.
 *  
 *  Created on: Apr 7, 2012
 *      Author: blew
 */

#ifndef NO_HDF5

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <sys/dir.h>
//#include <QtCore/QStringList>
#include "qtbase.h"
//#include <pthread.h>
//#include <H5Cpp.h>
//#include <fftw3.h>
#include "Mscs-function3dregc.h"
#include "Mscs-map-window_function.h"
#include "cpeds-math.h"
#include "cpeds-point3d.h"


//string _HDF5currentGroup_tmp;
//vector<string> _HDF5_datasets;

/***************************************************************************************/
void mscsFunction3dregc::initiatie_hdf5_params() {
	//#ifdef HAVE_OPENMP
	//#pragma omp critical (mscsFunction3dregc_hdfAttr)	
	//#endif
	//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	//	omp_hdf_lock=NULL;
	/*
	 * Comment: The above code was moved to all public methods without any locking
	 * Usage of HDF5 assumes implicitly that the programmer must serialize access to HDF5 methods
	 * 
	 * 
	 * 
	 * author: blew
	 * date: Dec 19, 2013 2:22:38 PM
	 *
	 */
	
}

/***************************************************************************************/
void mscsFunction3dregc::setHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, double value, string attributeComment) {
	hid_t dset,file;
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (cpeds_fileExists(fname)) {
		// if (omp_hdf_lock!=NULL) omp_set_lock(omp_hdf_lock);
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
		if (dset<0) {printf("setHDF5_scalarDoubleAttribute attribute: could not open dataset %s\n",dsetName.c_str()); return; }

		setHDF5_scalarDoubleAttribute(dset,attributeName,value,attributeComment);

		hid_t status;
		status = H5Dclose (dset);
		status = H5Fclose (file);
	}
}

void mscsFunction3dregc::setHDF5_scalarStringAttribute(hid_t& dset, string attributeName, string attributeValue) {
	//
	// create comment string attribute
	//
	if (attributeValue.size()==0) attributeValue="NONE";

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
	
}

/***************************************************************************************/
/*
 * Comment: This routine does not work correctly in either parallel or serial mode.
 * While the attributes are indeed saved to file, the HDF calls print errors to stderr.
 * 
 * It seems that a better way to save string arguments is at the dataset creation time
 * when the same hid_t dset can be used as the one used for dataset creation.
 * 
 * To shut-up the error messages you can enable H5Eset_auto(H5E_DEFAULT, NULL, NULL); line
 * but generally this routine should be debuged.
 * 
 * author: blew
 * date: Feb 9, 2016 11:25:13 AM
 *
 */

void mscsFunction3dregc::setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue) {

	hid_t dset,file;
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (attributeValue.size()==0) attributeValue="NONE";

	if (cpeds_fileExists(fname)) {
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
		if (dset<0) {printf("setHDF5_scalarStringAttribute attribute: could not open dataset %s\n",dsetName.c_str()); return; }
		
		//
		// create comment string attribute
		//
		setHDF5_scalarStringAttribute(dset,attributeName,attributeValue);
		herr_t ret2;
//		hid_t atype2 = H5Tcopy(H5T_C_S1);
//		hid_t aid2 = H5Screate(H5S_SCALAR); // Create scalar attribute.
//		ret2 = H5Tset_size(atype2, attributeValue.size());
//		ret2 = H5Tset_strpad(atype2,H5T_STR_NULLTERM);
//		hid_t attr2 = H5Acreate(dset, attributeName.c_str(), atype2, aid2, H5P_DEFAULT, H5P_DEFAULT);
//		ret2 = H5Awrite(attr2, atype2, attributeValue.c_str()); 
//		ret2 = H5Sclose(aid2); // Close attribute dataspace.
//		ret2 = H5Aclose(attr2); //Close attribute.        	
//		ret2 = H5Tclose(atype2);
		
		ret2 = H5Dclose (dset);
		ret2 = H5Fclose (file);
			
	}
}
/***************************************************************************************/
string mscsFunction3dregc::getHDF5_stringAttribute(string fname, string dsetName, string attributeName, int* errCode) {
	hid_t dset,file;
	hid_t status;
	string attr;
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (cpeds_fileExists(fname)) {
		// if (omp_hdf_lock!=NULL) omp_set_lock(omp_hdf_lock);
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		attr=getHDF5_stringAttribute(file,dsetName,attributeName,errCode);
		status = H5Fclose (file);
		// if (omp_hdf_lock!=NULL) omp_unset_lock(omp_hdf_lock);
	}
	else {
		if (errCode!=NULL) *errCode=-1;
	}
	return attr;	
}

/***************************************************************************************/
double mscsFunction3dregc::getHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, int* errCode) {
	hid_t dset,file;
	hid_t status;
	double attr=0;
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (cpeds_fileExists(fname)) {
		// if (omp_hdf_lock!=NULL) omp_set_lock(omp_hdf_lock);
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		//        dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
		
		attr=getHDF5_scalarDoubleAttribute(file,dsetName,attributeName,errCode);
		
		//        status = H5Dclose (dset);
		status = H5Fclose (file);
		// if (omp_hdf_lock!=NULL) omp_unset_lock(omp_hdf_lock);
	}
	else {
		if (errCode!=NULL) *errCode=-1;
	}
	return attr;	
}

/***************************************************************************************/
void mscsFunction3dregc::setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment) {
	hid_t aid = H5Screate(H5S_SCALAR); // Create scalar attribute.
	hid_t attr = H5Acreate(dset, attributeName.c_str(), H5T_IEEE_F64LE, aid, H5P_DEFAULT, H5P_DEFAULT);
	hid_t status = H5Awrite(attr, H5T_IEEE_F64LE, &value); // Write scalar attribute.
	status = H5Sclose(aid); // Close attribute dataspace.
	status = H5Aclose(attr); //Close attribute.
	
	//
	// create comment string attribute
	//		

	setHDF5_scalarStringAttribute(dset,attributeName+"_comment",attributeComment);
	
}
/***************************************************************************************/
double mscsFunction3dregc::getHDF5_scalarDoubleAttribute(hid_t& file, string dsetName, string attributeName, int* errCode) {
	double val=0;
	//#pragma omp critical (mscsFunction3dregc_hdfAttr)
	//	{
	/*
	 * Attach to the scalar attribute using attribute name, then read and
	 * display its value.
	 */
	//	printf("attrib name: %s\n",attributeName.c_str());
	htri_t tfe=H5Aexists_by_name(file,dsetName.c_str(), attributeName.c_str(),	H5P_DEFAULT);
	if (tfe>0) {
		hid_t attr = H5Aopen_by_name(file, dsetName.c_str(), attributeName.c_str(),	H5P_DEFAULT, H5P_DEFAULT);
		herr_t errc = H5Aread(attr, H5T_IEEE_F64LE, &val);
		hid_t status = H5Aclose(attr);
		if (errCode!=NULL) *errCode=errc;
	}
	else {
		if (tfe==0)	*errCode=-1;
		else *errCode=tfe;
	}
	return val;
}
/***************************************************************************************/
string mscsFunction3dregc::getHDF5_stringAttribute(hid_t& file, string dsetName, string attributeName, int* errCode) {
	string strAttr="";
	
	//#pragma omp critical (mscsFunction3dregc_hdfAttr)
	//	{
	
	
	/*
	 * Attach to the scalar attribute using attribute name, then read and
	 * display its value.
	 */
	//	printf("attrib name: %s\n",attributeName.c_str());
	hid_t attr = H5Aopen_by_name(file, dsetName.c_str(), attributeName.c_str(),	H5P_DEFAULT, H5P_DEFAULT);
	hid_t filetype = H5Aget_type (attr);
	size_t sdim = H5Tget_size (filetype);
	sdim++;                         /* Make room for null terminator */
	
	
	hsize_t     dims[1];
	hid_t space = H5Aget_space (attr);
	//    int ndims = H5Sget_simple_extent_ndims(space);
	int ndims = H5Sget_simple_extent_dims (space, dims, NULL);
	//    /*
	//     * Allocate array of pointers to rows.
	//     */
	//    char* rdata = (char *) malloc (dims[0] * sizeof (char));
	
	//    /*
	//     * Allocate space for integer data.
	//     */
	//    rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));
	
	//    /*
	//     * Set the rest of the pointers to rows to the correct addresses.
	//     */
	//    for (i=1; i<dims[0]; i++)
	//        rdata[i] = rdata[0] + i * sdim;
	
	//    /*
	//     * Create the memory datatype.
	//     */
	hid_t memtype = H5Tcopy (H5T_C_S1);
	herr_t status = H5Tset_size (memtype, sdim);
	
	/*
	 * Read the data.
	 */
	//    status = H5Aread (attr, memtype, rdata[0]);
//	char attrch[HDF5_stringAttributeMaxLength];

	hsize_t attSize=H5Aget_storage_size(attr);
	char attrch[attSize];
	herr_t errc = H5Aread(attr, memtype, attrch);
	//	herr_t errc = H5Aread(attr, H5T_STR_NULLTERM, attrch);
	
	//    /*
	//     * Output the data to the screen.
	//     */
	//    for (i=0; i<dims[0]; i++)
	//        printf ("%s[%d]: %s\n", attributeName.c_str(), i, rdata[i]);
	
	strAttr=attrch;
	/*
	 * Close and release resources.
	 */
	//    free (rdata[0]);
	//    free (rdata);
	//    H5Aclose (attr);
	H5Tclose (memtype);
	H5Sclose (space);
	//    H5Dclose (attr);
	//    status = H5Tclose (filetype);
	//    status = H5Fclose (file);
	
	
	
	
	
	H5Aclose(attr);
	
	if (errCode!=NULL) *errCode=errc;
	//	}
	return strAttr;
}
/***************************************************************************************/
bool mscsFunction3dregc::hdf5DatasetExists(hid_t& file, string dsetName) {
//	bool tf=false;
//	htri_t tfe=0;
//	QString qs=dsetName.c_str();
//	if (qs[0]!='/') qs="/"+qs; // make sure we operate on absolute path
//	
//	QStringList qsl=qs.split("/");
//	if (qsl.size()==2) {
//		tfe=H5Lexists(file, dsetName.c_str(), H5P_DEFAULT);		
//	} 
//	else {
//		if (qsl.size()==3) {
//			qs=qsl[0]+"/"+qsl[1];
//			if (H5Lexists(file, qs.toStdString().c_str(), H5P_DEFAULT)) {
//				qs=qsl[0]+"/"+qsl[1]+"/"+qsl[2];
//				if (H5Lexists(file, dsetName.c_str(), H5P_DEFAULT)) { return true; }
//				else return false;
//			}
//			else return false;
//		}
//	}
//	
//	if (tfe) return true;
//	return tf;
	bool tf=false;
	htri_t tfe=0;
	QString qs=dsetName.c_str();
	if (qs[0]!='/') qs="/"+qs; // make sure we operate on absolute path
	
	QStringList qsl=qs.split("/");
	if (qsl.size()==2) {
		tfe=H5Lexists(file, dsetName.c_str(), H5P_DEFAULT);		
	} 
	else {
		if (qsl.size()==3) {
			qs=qsl[0]+"/"+qsl[1];
			if (H5Lexists(file, qs.toStdString().c_str(), H5P_DEFAULT)) {
				qs=qsl[0]+"/"+qsl[1]+"/"+qsl[2];
				if (H5Lexists(file, dsetName.c_str(), H5P_DEFAULT)) { return true; }
				else return false;
			}
			else return false;
		}
	}
	
	if (tfe) return true;
	return tf;
	
	
}

/***************************************************************************************/
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::saveHDF5(string filename, string datasetName, int ix, int iy, int iz, int nx, int ny, int nz, int part) {
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (nx==-1) nx=Nx()-ix;
	if (ny==-1) ny=Ny()-iy;
	if (nz==-1) nz=Nz()-iz;
	
	if (size()==0 or nx==0 or ny==0 or nz==0) {
		return cpedsError;
	}
	
	
	
	/***************************************************************************************/
	int rank=3;
	hid_t file, space, dset, dcpl;    /* Handles */
	herr_t status;
	H5D_layout_t layout;
	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
	double *wdata;          /* Write buffer */
	
	dims[0]=nx;
	dims[1]=ny;
	dims[2]=nz;
	chunk[0]=nx;
	chunk[1]=ny;
	chunk[2]=1;
	
	if (cpeds_fileExists(filename)) {
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
	string linkName=datasetName;
	if (hdf5DatasetExists(file,linkName)) {
//		printf("hdf5DatasetExists(file=%i,linkName=%s)\n",file,linkName.c_str());
		H5Ldelete(file,linkName.c_str(),H5P_DEFAULT);
//		exit(-1);
//		printf("deleting dataset\n");
	}
	else {
		// create groups if needed
		QString qs=linkName.c_str();
		QStringList qsl=qs.split("/");
		if (qsl.size()==3)
			hdf5createGroup(file,linkName);
	}
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
	if (part==0) {
		for (long i = ist; i < ien; i++) {
			for (long j = jst; j < jen; j++) {
				for (long k = kst; k < ken; k++) {
					wdata[l]=fRe(i,j,k);
					l++;
				}
			}
		}		
	}
	if (part==1) {
		for (long i = ist; i < ien; i++) {
			for (long j = jst; j < jen; j++) {
				for (long k = kst; k < ken; k++) {
					wdata[l]=fIm(i,j,k);
					l++;
				}
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
	
#ifndef DEBUG_HDF5_ATTRIBUTES
	setHDF5_scalarDoubleAttribute(dset,"x0",_param.x0, "function domain lower range in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"y0",_param.y0, "function domain lower range in y direction");
	setHDF5_scalarDoubleAttribute(dset,"z0",_param.z0, "function domain lower range in z direction");
	setHDF5_scalarDoubleAttribute(dset,"dx",_param.dx, "function cell size in x direction");
	setHDF5_scalarDoubleAttribute(dset,"dy",_param.dy, "function cell size in y direction");
	setHDF5_scalarDoubleAttribute(dset,"dz",_param.dz, "function cell size in z direction");	
	setHDF5_scalarDoubleAttribute(dset,"xMax",_param.xMax, "function domain upper range in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"yMax",_param.yMax, "function domain upper range in y direction");
	setHDF5_scalarDoubleAttribute(dset,"zMax",_param.zMax, "function domain upper range in z direction");
	setHDF5_scalarDoubleAttribute(dset,"dxo2",_param.dxo2, "function domain cell center offset in x direction"); 
	setHDF5_scalarDoubleAttribute(dset,"dyo2",_param.dyo2, "function domain cell center offset in y direction");
	setHDF5_scalarDoubleAttribute(dset,"dzo2",_param.dzo2, "function domain cell center offset in z direction");
#endif
	//	setHDF5_scalarDoubleAttribute(dset,"redshift",);	
	
	/*
	 * Close and release resources.
	 */
	status = H5Pclose (dcpl);
	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Fclose (file);
	
	return cpedsSuccess;
}
/***************************************************************************************/
bool mscsFunction3dregc::hdf5DatasetExists(string filename, string dsetName) {
	bool dsexists=false;
	vector<string> dsets=getHDF5dataSets(filename);
	for (unsigned long i = 0; i < dsets.size(); i++) {
		if (dsets[i]==dsetName) return true;
	}
	return dsexists;
}

/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::loadHDF5(string filename, string datasetName, int part) {
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	int rank=3;
	hid_t file, space, dset, dcpl;    /* Handles */
	herr_t status;
	H5D_layout_t layout;
	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
	double *wdata;          /* Write buffer */
	
	//    chunk[0]=nx;
	//    chunk[1]=ny;
	//    chunk[2]=1;
	
	/*
	 * Open file using the default properties.
	 */
	if (cpeds_fileExists(filename)) {
		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);    	    	
	}
	else {
		msgs->error("cannot open file: "+filename,High);
		return cpedsNoSuchFile;
	}
	/* Open an existing dataset. */
	// if (omp_hdf_lock!=NULL) omp_set_lock(omp_hdf_lock);
	dset = H5Dopen(file, datasetName.c_str(), H5P_DEFAULT);
	// if (omp_hdf_lock!=NULL) omp_unset_lock(omp_hdf_lock);
	
	if (dset<0) {
		msgs->error("dataset '"+datasetName+"' not found in file: "+filename,High);
		return cpedsNoSuchHDF5dataset;
	}
	
	
	/*
	 * Create dataspace.  Setting maximum size to NULL sets the maximum
	 * size to be the current size.
	 */
	//#pragma omp critical (mscsFunction3dregc_hdf)
	//	{
	// if (omp_hdf_lock!=NULL) omp_set_lock(omp_hdf_lock);
	space=H5Dget_space(dset);
	rank=H5Sget_simple_extent_ndims(space);
	if (rank!=3) msgs->criticalError("the hf5 file specified dataset is not of rank 3. Cannot load the data",High);
	H5Sget_simple_extent_dims(space,dims,NULL);
	//	}	
	long nx,ny,nz;
	nx=dims[0];
	ny=dims[1];
	nz=dims[2];
	
	
	
	//    /*
	//     * Create the dataset creation property list, and set the chunk
	//     * size.
	//     */
	//    dcpl = H5Pcreate (H5P_DATASET_CREATE);
	//    status = H5Pset_chunk(dcpl, rank, chunk);
	//    double defalutVal=1;
	//    H5Pset_fill_value(dcpl,H5T_IEEE_F64LE,&defalutVal);
	
	//    /*
	//     * Create the chunked dataset.
	//     */
	//    dset = H5Dcreate(file, datasetName.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
	
	
	
	/*
	 * Read the data from file
	 */
	
	long Nslice=nx*ny*nz;
	wdata = new double[Nslice];
	
	
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
	
	//#pragma omp critical (mscsFunction3dregc_hdf)
	//	{
	
	status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count, block);
	status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, wdata);
	//	}	
	
	double x0=getHDF5_scalarDoubleAttribute(file,datasetName,"x0");
	double y0=getHDF5_scalarDoubleAttribute(file,datasetName,"y0");
	double z0=getHDF5_scalarDoubleAttribute(file,datasetName,"z0");
	double dx=getHDF5_scalarDoubleAttribute(file,datasetName,"dx");
	double dy=getHDF5_scalarDoubleAttribute(file,datasetName,"dy");
	double dz=getHDF5_scalarDoubleAttribute(file,datasetName,"dz");
	
	//	printf("%lf %lf %lf %lf %lf %lf\n",x0,y0,z0,dx,dy,dz);
	/*
	 * Close and release resources.
	 */
	//    status = H5Pclose (dcpl);
	//#pragma omp critical (mscsFunction3dregc_hdf)
	//	{
	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Fclose (file);
	//	}	
	// if (omp_hdf_lock!=NULL) omp_unset_lock(omp_hdf_lock);
	
	/*
	 * copy the data onto the function
	 */    
	
	
	long ist,ien,jst,jen,kst,ken;
	ist=0;
	jst=0;
	kst=0;
	ien=nx;
	jen=ny;
	ken=nz;
	
	long l=0;
	
	setSize(nx,ny,nz,dx,dy,dz,x0,y0,z0);
	allocFunctionSpace();
	//	initiate(true);
	
	// prepare slice data    	
	if (part==0) {
		for (long i = ist; i < ien; i++) {
			for (long j = jst; j < jen; j++) {
				for (long k = kst; k < ken; k++) {
					fRe(i,j,k)=wdata[l];
					fIm(i,j,k)=0;
					l++;
				}
			}
		}		
	}
	if (part==1) {
		for (long i = ist; i < ien; i++) {
			for (long j = jst; j < jen; j++) {
				for (long k = kst; k < ken; k++) {
					fRe(i,j,k)=0;
					fIm(i,j,k)=wdata[l];
					l++;
				}
			}
		}		
	}
	
	delete [] wdata;
	
	return cpedsSuccess;
}



/***************************************************************************************/
long mscsFunction3dregc::getHDF5dsetCount(hid_t& file) {
	H5F_info_t info;
	hsize_t num;
	H5G_info_t group_info;
	//	H5Fget_info(file,&info);
	
//		printf("%li\n",long(H5Fget_obj_count(file,H5F_OBJ_ALL | H5F_OBJ_LOCAL)));
//	exit(0);
//		printf("%li\n",long(H5Gget_num_objs(file,&num)));
	H5Gget_info(file,&group_info);
	//	printf("%li\n",group_info.nlinks);
	//	printf("%li\n",group_info.info.super_ext_size);
	return group_info.nlinks;
	//    return long(H5Fget_obj_count(file,H5F_OBJ_DATASET));
}


/***************************************************************************************/
vector<string> mscsFunction3dregc::getHDF5dataSets(string filename) {
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	vector<string> dsetNames;

	
    hid_t           file;           /* Handle */
    herr_t          status;
    H5O_info_t      infobuf;
    struct opdata   od;

    /*
     * Open file and initialize the operator data structure.
     */
    file = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    status = H5Oget_info (file, &infobuf);
    od.recurs = 0;
    od.prev = NULL;
    od.addr = infobuf.addr;

    /*
     * Print the root group and formatting, begin iteration.
     */
    status = H5Literate (file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func,
                (void *) &od);

    /*
     * Close and release resources.
     */
    status = H5Fclose (file);
	
    
    dsetNames=_HDF5_datasets;
	
    _HDF5_datasets.clear();
    _HDF5currentGroup_tmp="";
	
//	//	int rank=3;
//	hid_t file, dset;    /* Handles */
//	herr_t status;
//	//	H5D_layout_t layout;
//	//	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
//	//    double *wdata;          /* Write buffer */
//	
//	//    chunk[0]=nx;
//	//    chunk[1]=ny;
//	//    chunk[2]=1;
//	
//	/*
//	 * Open file using the default properties.
//	 */
//	if (cpeds_fileExists(filename)) {
//		file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//	}
//	else {
//		msgs->criticalError("cannot open file: "+filename,High);
//	}
//	long dsetCount=getHDF5dsetCount(file);
//	msgs->say("There are %li datasets in the file",dsetCount,Medium);
////	printf("%li\n",dsetCount);
//	H5G_info_t group_info;
//	for (long i = 0; i < dsetCount; i++) {
//		//		H5Gget_info_by_idx(file,"/",H5_INDEX_NAME,H5_ITER_NATIVE,0,&group_info,0);
//		ssize_t objSize=H5Gget_objname_by_idx(file,i,NULL,0);
//		if (objSize>0) {
//			char name[objSize+1];
//			H5Gget_objname_by_idx(file,i,name,objSize+1);
////			printf("name: %s\n",name);
//			dsetNames.push_back(name);
//		}
//	}
//	H5Fclose(file);
	return dsetNames;
}

/***************************************************************************************/
//#ifdef HAVE_OPENMP
//void mscsFunction3dregc::setOMPlock(omp_lock_t **lock) { 
//	omp_hdf_lock=*lock; 
//	omp_init_lock(omp_hdf_lock);
//}
//#endif
/***************************************************************************************/
bool mscsFunction3dregc::hdf5createGroup(hid_t& file,string linkName) {
	bool tf;
	QString qs=linkName.c_str();
	QStringList qsl=qs.split("/");
	qs="/"+qsl[1];
	if (hdf5DatasetExists(file,qs.toStdString())) return true;
	hid_t group1_id = H5Gcreate(file, qs.toStdString().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group1_id);
	if (group1_id!=0) return false;
	return true;
}

#endif


/***************************************************************************************/

