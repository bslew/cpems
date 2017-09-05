/*!
 * \file cpedsHDF5.h
 * 
 * UNDER CONSTRUCTION - NOT IMPLEMENTED YET
 *
 *  Created on: Mar 23, 2012
 *      Author: blew
 */

#ifndef CPEDSHDF5_H_
#define CPEDSHDF5_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates standard IO operations on HDF5 files
 \details 
 
 \date Mar 23, 2012, 7:09:56 PM 
 \author Bartosz Lew
 */

class cpedsHDF5 {
	public:
		cpedsHDF5();
		void openFile(string fname);
		void createFile(string fname);
		void addDataset(string dsetName,);
		
		virtual ~cpedsHDF5();
};

#endif /* CPEDSHDF5_H_ */
