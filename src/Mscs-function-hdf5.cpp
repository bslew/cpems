/*!
 * \file Mscs-function-hdf5.cpp
 *
 *  This is an add-on to the Mscs-function class - compiled in optionally to support HDF5 files format.
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
//#include <pthread.h>
#include <omp.h>
//#include <H5Cpp.h>
//#include <fftw3.h>
#include "Mscs-function.h"
#include "Mscs-map-window_function.h"
#include "cpeds-math.h"
#include <QtCore/QStringList>
#include <QtCore/QString>
//#include "cpeds-point3d.h"

string _HDF5currentGroup_tmp;
vector<string> _HDF5_datasets;


/***************************************************************************************/
void mscsFunction::initiatie_hdf5_params() {
//#ifdef HAVE_OPENMP
//#pragma omp critical
//#endif
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
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
/***************************************************************************************/
bool mscsFunction::hdf5DatasetExists(hid_t& file, string dsetName) {
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
cpedsStatusCodes mscsFunction::saveHDF5(string filename, string datasetName, string xUnit,string yUnit, string xComment,string yComment) {
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	int nx,ny;
	
	nx=2;
	ny=pointsCount();
	if (ny==0) { return cpedsError; }
    /***************************************************************************************/
	int rank=2;
	hid_t file, space, dset, dcpl;    /* Handles */
	herr_t status;
	H5D_layout_t layout;
	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
    double *wdata;          /* Write buffer */

    dims[0]=nx;
    dims[1]=ny;
    chunk[0]=nx;
    chunk[1]=ny;
    
    /*
     * Create a new file using the default properties.
     */
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
    double defalutVal=1;
    H5Pset_fill_value(dcpl,H5T_IEEE_F64LE,&defalutVal);

    /*
     * Create the chunked dataset.
     */
    string linkName=datasetName;
    if (hdf5DatasetExists(file,linkName))
    	H5Ldelete(file,linkName.c_str(),H5P_DEFAULT);
	else {
		// created groups if needed
		QString qs=linkName.c_str();
		QStringList qsl=qs.split("/");
		if (qsl.size()==3)
			hdf5createGroup(file,linkName);
	}
    dset = H5Dcreate(file, datasetName.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);



    /*
     * Write the data to the dataset.
     */
    long Nslice=nx*ny;
    wdata = new double[Nslice];
    
    long ist,ien,jst,jen;
    ist=0;
    jst=0;
    ien=ist+nx;
    jen=jst+ny;
    
	long l=0;
    	
	// prepare slice data    	
	for (long i = ist; i < ien; i++) {
		if (i==0) {
			for (long j = jst; j < jen; j++) {
				wdata[l]=getx(j);
				l++;
			}			
		}
		if (i==1) {
			for (long j = jst; j < jen; j++) {
				wdata[l]=f(j);
				l++;
			}			
		}
	}		
	
	
	/*
	 * Define and select the first part of the hyperslab selection.
	 */
	start[0] = 0;
	start[1] = 0;
	block[0]=1;
	block[1]=1;
	stride[0]=1;
	stride[1]=1;
	count[0]=nx;
	count[1]=ny;
	status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count, block);
	status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, wdata);

	delete [] wdata;

//#ifndef DEBUG_HDF5_ATTRIBUTES
//    setHDF5_scalarStringAttribute(filename,datasetName,"xUnit",xUnit);
//    setHDF5_scalarStringAttribute(filename,datasetName,"yUnit",yUnit);
//    setHDF5_scalarStringAttribute(filename,datasetName,"xComment",xComment);
//    setHDF5_scalarStringAttribute(filename,datasetName,"yComment",yComment);
//#endif
    setHDF5_scalarStringAttribute(dset,"xUnit",xUnit);
    setHDF5_scalarStringAttribute(dset,"yUnit",yUnit);
    setHDF5_scalarStringAttribute(dset,"xComment",xComment);
    setHDF5_scalarStringAttribute(dset,"yComment",yComment);

    /*
     * Close and release resources.
     */
    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
//    printf("closing dset: %li\n",status);
    status = H5Sclose (space);
//    printf("closing space: %li\n",status);
    status = H5Fclose (file);
//    printf("closing file: %li\n",status);

}
/***************************************************************************************/
void mscsFunction::setHDF5_scalarStringAttribute(hid_t& dset, string attributeName, string attributeValue) {
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
void mscsFunction::setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment) {
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
cpedsStatusCodes mscsFunction::loadHDF5(string filename, string datasetName) {
//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	int rank=2;
	hid_t file, space, dset, dcpl;    /* Handles */
	herr_t status;
	H5D_layout_t layout;
	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
    double *wdata;          /* Write buffer */

    
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
    dset = H5Dopen(file, datasetName.c_str(), H5P_DEFAULT);

    if (dset<0) {
		msgs->error("dataset '"+datasetName+"' not found in file: "+filename,High);
    	return cpedsNoSuchHDF5dataset;
    }
    
    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space=H5Dget_space(dset);
    rank=H5Sget_simple_extent_ndims(space);
    if (rank!=2) msgs->criticalError("the hf5 file specified dataset is not of rank 2. Cannot load the data",High);
    H5Sget_simple_extent_dims(space,dims,NULL);
    
    long nx,ny;
    nx=dims[0];
    ny=dims[1];
    


    /*
     * Read the data from file
     */

    long Nslice=nx*ny;
    wdata = new double[Nslice];
    
	
	/*
	 * Define and select the first part of the hyperslab selection.
	 */
	start[0] = 0;
	start[1] = 0;
	block[0]=1;
	block[1]=1;
	stride[0]=1;
	stride[1]=1;
	count[0]=nx;
	count[1]=ny;
	status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count, block);
	status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, space, H5P_DEFAULT, wdata);

    
//	double x0=getHDF5_scalarDoubleAttribute(file,datasetName,"x0");
//	double y0=getHDF5_scalarDoubleAttribute(file,datasetName,"y0");
//	double z0=getHDF5_scalarDoubleAttribute(file,datasetName,"z0");
//	double dx=getHDF5_scalarDoubleAttribute(file,datasetName,"dx");
//	double dy=getHDF5_scalarDoubleAttribute(file,datasetName,"dy");
//	double dz=getHDF5_scalarDoubleAttribute(file,datasetName,"dz");

//	printf("%lf %lf %lf %lf %lf %lf\n",x0,y0,z0,dx,dy,dz);
    /*
     * Close and release resources.
     */
//    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Fclose (file);

    
    /*
     * copy the data onto the function
     */    
    
    
    long ist,ien,jst,jen;
    ist=0;
    jst=0;
    ien=nx;
    jen=ny;
    
	long l=0;
    	
	setPointsNum(ny);

	
	// prepare slice data    	
	for (long i = ist; i < ien; i++) {
		if (i==0) {
			for (long j = jst; j < jen; j++) {
				setarg(j,wdata[l]);
				l++;
			}			
		}
		if (i==1) {
			for (long j = jst; j < jen; j++) {
				setf(j,wdata[l]);
				l++;
			}
		}

	}
	
	delete [] wdata;

	return cpedsSuccess;
}



/***************************************************************************************/
long mscsFunction::getHDF5dsetCount(hid_t& file) {
	H5F_info_t info;
	hsize_t num;
	H5G_info_t group_info;
	H5Gget_info(file,&group_info);
	return group_info.nlinks;
}

//long mscsFunction::getHDF5dsetCount(hid_t& file) {
//	printf("getHDF5dsetCount> THIS IS NOT IMPLEMENTED YET\n");
//	exit(-1);
//	printf("%li\n",long(H5Fget_obj_count(file,H5F_OBJ_ALL)));
//    return long(H5Fget_obj_count(file,H5F_OBJ_DATASET));
//}


/***************************************************************************************/
//vector<string> mscsFunction::getHDF5dataSets(string filename) {
////	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
//
//	vector<string> dsetNames;
////	int rank=3;
//	hid_t file, dset;    /* Handles */
//	herr_t status;
////	H5D_layout_t layout;
////	hsize_t dims[rank],	chunk[rank], start[rank], stride[rank],	count[rank], block[rank];
////    double *wdata;          /* Write buffer */
//
////    chunk[0]=nx;
////    chunk[1]=ny;
////    chunk[2]=1;
//    
//    /*
//     * Open file using the default properties.
//     */
//    if (cpeds_fileExists(filename)) {
//        file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//    }
//    else {
//		msgs->criticalError("cannot open file: "+filename,High);
//    }
//    long dsetCount=getHDF5dsetCount(file);
//    msgs->say("There are %li datasets in the file",dsetCount,Medium);
//    
////    H5Fget_obj_ids(file,H5F_OBJ_DATASET,)
//	exit(-1);
//}
vector<string> mscsFunction::getHDF5dataSets(string filename) {
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
    status = H5Literate (file, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func,(void *) &od);

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
//	
//	H5G_info_t group_info;
//	for (long i = 0; i < dsetCount; i++) {
//		//		H5Gget_info_by_idx(file,"/",H5_INDEX_NAME,H5_ITER_NATIVE,0,&group_info,0);
//		ssize_t objSize=H5Gget_objname_by_idx(file,i,NULL,0);
//		if (objSize>0) {
//			char name[objSize+1];
//			H5Gget_objname_by_idx(file,i,name,objSize+1);
//			dsetNames.push_back(name);
//		}
//	}
//	H5Fclose(file);
	return dsetNames;
}

/***************************************************************************************/
//void mscsFunction::setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue) {
//	hid_t dset,file;
////	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
//
//	if (cpeds_fileExists(fname)) {
//        file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
//        dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
//
//        //
//        // create comment string attribute
//        //
//        
////        string name=attributeName.substr(0,HDF5_stringAttributeMaxLength-1);
////        string val=attributeValue.substr(0,HDF5_stringAttributeMaxLength-1);
//        
//        herr_t ret2;
//        hid_t atype2 = H5Tcopy(H5T_C_S1);
//        hid_t aid2 = H5Screate(H5S_SCALAR); // Create scalar attribute.
////        ret2 = H5Tset_size(atype2, HDF5_stringAttributeMaxLength);
//        ret2 = H5Tset_size(atype2, attributeValue.size());
//        ret2 = H5Tset_strpad(atype2,H5T_STR_NULLTERM);
//        hid_t attr2 = H5Acreate(dset, attributeName.c_str(), atype2, aid2, H5P_DEFAULT, H5P_DEFAULT);
//        ret2 = H5Awrite(attr2, atype2, attributeValue.c_str()); 
//        ret2 = H5Sclose(aid2); // Close attribute dataspace.
//        ret2 = H5Aclose(attr2); //Close attribute.        	
//		ret2 = H5Tclose(atype2);
//   
//        ret2 = H5Dclose (dset);
//        ret2 = H5Fclose (file);
//    }
//}
void mscsFunction::setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue) {

	hid_t dset,file;
	//	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	
	if (attributeValue.size()==0) attributeValue="NONE";
	
	if (cpeds_fileExists(fname)) {
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		if (hdf5DatasetExists(file,dsetName)==false) {printf("setHDF5_scalarStringAttribute attribute: could not open dataset %s\n",dsetName.c_str()); return; }
		dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
		
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
void mscsFunction::setHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, double value, string attributeComment) {
	hid_t dset,file;
	
	if (cpeds_fileExists(fname)) {
		file = H5Fopen(fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);        
		if (hdf5DatasetExists(file,dsetName)==false) {printf("setHDF5_scalarDoubleAttribute attribute: could not open dataset %s\n",dsetName.c_str()); return; }
		dset = H5Dopen(file, dsetName.c_str(),H5P_DEFAULT);
//		if (dset<0) {printf("setHDF5_scalarDoubleAttribute attribute: could not open dataset %s\n",dsetName.c_str()); return; }

		setHDF5_scalarDoubleAttribute(dset,attributeName,value,attributeComment);

		hid_t status;
		status = H5Dclose (dset);
		status = H5Fclose (file);
	}
}



/***************************************************************************************/
bool mscsFunction::hdf5createGroup(hid_t& file,string linkName) {
	bool tf;
	QString qs=linkName.c_str();
	QStringList qsl=qs.split("/");
	qs="/"+qsl[1];
	if (hdf5DatasetExists(file,qs.toStdString())) return true;
//	printf("creating group\n");
	hid_t group1_id = H5Gcreate(file, qs.toStdString().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	printf("res: %li\n",group1_id);
	H5Gclose(group1_id);
	if (group1_id!=0) return false;
	return true;
}

/***************************************************************************************/
/************************************************************

  Operator function.  This function prints the name and type
  of the object passed to it.  If the object is a group, it
  is first checked against other groups in its path using
  the group_check function, then if it is not a duplicate,
  H5Literate is called for that group.  This guarantees that
  the program will not enter infinite recursion due to a
  circular path in the file.

 ************************************************************/
herr_t op_func (hid_t loc_id, const char *name, const H5L_info_t *info,void *operator_data)
{
    herr_t          status, return_val = 0;
    H5O_info_t      infobuf;
    struct opdata   *od = (struct opdata *) operator_data;
                                /* Type conversion */
    unsigned        spaces = 2*(od->recurs+1);
                                /* Number of whitespaces to prepend
                                   to output */

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by
     * the Library.
     */
    status = H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
//    printf ("%*s", spaces, "");     /* Format output */
    switch (infobuf.type) {
        case H5O_TYPE_GROUP:
//            printf ("Group: %s {\n", name);
            _HDF5currentGroup_tmp=name;
            /*
             * Check group address against linked list of operator
             * data structures.  We will always run the check, as the
             * reference count cannot be relied upon if there are
             * symbolic links, and H5Oget_info_by_name always follows
             * symbolic links.  Alternatively we could use H5Lget_info
             * and never recurse on groups discovered by symbolic
             * links, however it could still fail if an object's
             * reference count was manually manipulated with
             * H5Odecr_refcount.
             */
            if ( group_check (od, infobuf.addr) ) {
                printf ("%*s  Warning: Loop detected!\n", spaces, "");
            }
            else {

                /*
                 * Initialize new operator data structure and
                 * begin recursive iteration on the discovered
                 * group.  The new opdata structure is given a
                 * pointer to the current one.
                 */
                struct opdata nextod;
                nextod.recurs = od->recurs + 1;
                nextod.prev = od;
                nextod.addr = infobuf.addr;
                return_val = H5Literate_by_name (loc_id, name, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, (void *) &nextod,H5P_DEFAULT);
            }
//            printf ("%*s}\n", spaces, "");
            break;
        case H5O_TYPE_DATASET:
//            printf ("Dataset: %s\n", name);
            _HDF5_datasets.push_back(_HDF5currentGroup_tmp+"/"+name);
            break;
        case H5O_TYPE_NAMED_DATATYPE:
//            printf ("Datatype: %s\n", name);
            break;
        default:
//            printf ( "Unknown: %s\n", name);
            break;
    }

    return return_val;
}


/************************************************************

  This function recursively searches the linked list of
  opdata structures for one whose address matches
  target_addr.  Returns 1 if a match is found, and 0
  otherwise.

 ************************************************************/
int group_check (struct opdata *od, haddr_t target_addr)
{
    if (od->addr == target_addr)
        return 1;       /* Addresses match */
    else if (!od->recurs)
        return 0;       /* Root group reached with no matches */
    else
        return group_check (od->prev, target_addr);
                        /* Recursively examine the next node */
}





#endif
