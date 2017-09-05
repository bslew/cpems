/*!
 * \file test-mscsfn-hdf5.cpp
 *
 *  Created on: Mar 23, 2012
 *      Author: blew
 */

#include "Mscs-function3dregc.h"

int main(int argc, char** argv) {
	if (argc==1) {
		printf("usage: test-mscsfn-hdf5 datasetName\n");
		exit(0);
	}
	
	mscsFunction3dregc f,g;
	long N=10;
	f.setSize(N,N,N,1,1,1,0,0,0);
//	f.initiate(true);
	f.allocFunctionSpace();
	long l=0;
	
	for (long i = 0; i < N; i++) {
		for (long j = 0; j < N; j++) {
			for (long k = 0; k < N; k++) {
				f.fRe(i,j,k)=l;
				l++;
			}
		}
	}
	
	f.saveHDF5("test.hdf5",argv[1],1,1,1,N-2,N-2,N-2,0);
	f.setHDF5_scalarDoubleAttribute("test.hdf5",argv[1],"redshift",1.234);
	f.setHDF5_scalarStringAttribute("test.hdf5",argv[1],"asd","asdasdads");

	
	g.loadHDF5("test.hdf5",argv[1],0);
//	g.savetxtlin("txtlin",true);
	
	return 0;
}
