/*!
 \file MscsDateTimeData.h - 
 */

#ifndef SRC_MSCSDATETIMEDATA_H_
#define SRC_MSCSDATETIMEDATA_H_

//#include "Mscs-function3dregc.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
//#include <QtCore/QList>
//#include <QtCore/QDateTime>
#include "qtbase.h"
#include "cpeds-math.h"


/*!
 \class MscsDateTimeData
 \brief Encapsulates IO operations on text files containing date,time and other data
 \details 
 
 \date created: Feb 9, 2018, 3:38:47 PM 
 \author Bartosz Lew
 */
class MscsDateTimeData {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		typedef struct {
			QList<double> data;
			QDateTime dt;
		} row_t;
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
				QList<row_t> rows;
				long dtCol;
				QString dtFmt;
				QString fname;
		} DateTimeData_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		MscsDateTimeData();
		~MscsDateTimeData();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		/*!
			\brief load file fname
			\details 
			@param dtCol column containing date and time
			@param dtFmt - date time format
			@return
		
			\date Feb 9, 2018, 3:47:28 PM
		*/
		void load(string fname, long dtCol, string dtFmt);
		
		void printData();
		
		long linesCount() { return _DTdata.rows.size(); }
		
		/*!
			\brief return first function data from column col after given date
			\details 
			@param dt - date string in the same format as used during loading data
			@param col - column from the data block; The column indexes start from 0
			and are consistent with the input file data with the date/time columns removed
			@param startFrom - index referece from which the search should start. When found,
			this reference is set to the index of the row that matches the selection criteria.
			If not found, this is set to the linestCount() value
			@return
		
			\date Feb 9, 2018, 6:51:19 PM
		*/
		double getFirstValueAfter(string dt, long col, long startFrom, long* idx=NULL);
		double getFirstValueAfter(QDateTime qdt, long col, long startFrom, long* idx=NULL);
		double getFirstValueAfter(double JD, long col, long startFrom, long* idx=NULL);
		
		QDateTime getTime(long idx) { return _DTdata.rows[idx].dt; }
		double getTimeJD(long idx);
		double getData(long col, long row) { return _DTdata.rows[row].data[col]; }
		long getDataCols(long row=0) { return _DTdata.rows[row].data.size(); }
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */

		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */

		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */

		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		DateTimeData_t _DTdata;
//		FILE* _f;
};
#endif /* SRC_MSCSDATETIMEDATA_H_ */ 

