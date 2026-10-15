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
#include <ctime>
#include <string>
#include <tuple>
#include "cpeds-math.h"
#include "cpeds-list.h"

struct MscsDateTime {
	std::tm value = {};

	bool operator<(const MscsDateTime& other) const {
		return std::tie(value.tm_year, value.tm_mon, value.tm_mday, value.tm_hour, value.tm_min, value.tm_sec) <
				std::tie(other.value.tm_year, other.value.tm_mon, other.value.tm_mday, other.value.tm_hour, other.value.tm_min, other.value.tm_sec);
	}
};


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
			cpedsList<double> data;
			MscsDateTime dt;
		} row_t;
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
				cpedsList<row_t> rows;
				long dtCol;
				std::string dtFmt;
				std::string fname;
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
		void load(std::string fname, long dtCol, std::string dtFmt);
		
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
		double getFirstValueAfter(std::string dt, long col, long startFrom, long* idx=NULL);
		double getFirstValueAfter(MscsDateTime qdt, long col, long startFrom, long* idx=NULL);
		double getFirstValueAfter(double JD, long col, long startFrom, long* idx=NULL);
		
		MscsDateTime getTime(long idx) { return _DTdata.rows[idx].dt; }
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
