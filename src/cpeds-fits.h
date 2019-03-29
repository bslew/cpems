/*!
  \file Implements support for viewing, saving and loading fits files, field manipulations etc.
 */

#ifndef CPEDSFITS
#define CPEDSFITS

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <stdio.h>
#include <string>
#include <fitsio.h>
//#include <QtCore/QList>
#include "cpeds-common.h"
#include "cpeds-msgs.h"
#include "cpeds-list.h"
#include "Mscs-function.h"
#include "Mscs-function3dregc.h"


/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */
using namespace std;


/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
  \class 
  \brief Encapsulates support for viewing, saving and loading fits files, field manipulations etc.
  \details 
  This is just a wrapper of 
  
  \date 2010/02/24 16:45:19 
  \author Bartosz Lew
 */
class cpedsFits {
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		typedef struct {
			string keyName;
			string valStr;
			double valDouble;
			long valInt;
			string comment;
		} headerRow_t;
		
		/* ------------- */
		/* CLASS FRIENDS */
		/* ------------- */
		
		
		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		cpedsFits(cpeds_VerbosityLevel level=Medium);
		cpedsFits(string fname, string action="new", cpeds_VerbosityLevel level=Medium);
		
		~cpedsFits();
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		cpedsStatusCodes closeFile();
		/*!
			\brief opens the fits file
			\details 
			@param fname - file name
			@param action - "new" for creating/overwritting a new file, "append" for appending to an existing file, "read" - for read only access
			@return
		
			\date Nov 21, 2011, 1:53:14 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes openFile(string fname="", string action="new");
		
		/*!
			\brief add binary table extension 
			\details 
			@param Tname - extension name
			@param rows - number of rows in the table
			@param colName - list of names of the columns
			@param colType - list of columns types; The allowed values are:\n
				long - for 8byte integer values
				double - for 8 byte floats
				float - for 4 byte floats
				longX - for 8byte integer type array of length X
				doubleX - for 8 byte double type array of length X
				floatX - for 4 byte float type array of length X
				stringX - for string type columns where X is a number that defines the length of the string. This is converted internally to a format AX understood by cfitsio
				bool - bool logical type - data insert handlers are not implemented yet for this type.
				PdoubleX - for a variable length array of type double of maximal size X elements  (converted internally to type PD(X))
				PfloatX - for a variable length array of type float of maximal size X elements (converted internally to type PE(X))
				PlongX - for a variable length array of type long of maximal size X elements  (converted internally to type PJ(X))
			
			@return
		
			\date Jul 9, 2012, 1:04:23 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addBinTable(string Tname, long rows, QList<string> &colName, QList<string>  &colType, QList<string>  &colUnit);
		/*!
			\brief generates binary table extension of zero length and without any columns
			\details 
			@param Tname - name of the extension
			@return
		
			\date Nov 20, 2011, 1:46:31 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addBinTable(string Tname);
		/*!
			\brief reads binary table from the selected hdu.
			\details 
			@param Tname - table name (hdu name)
			@param data - pointer the to array pointer where the read data will be placed in column-major order
			@param size - total size of the data array
			@param rows - number of rows in the array
			@return
			
			It is assumed that the data are all convertible to double numbers and that there are no ascii/string type columns in the 
			array.
		
			\date Mar 6, 2013, 1:08:26 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes readBinTable(string Tname, double** data, long* size, long *rows);
		/*!
			\brief reads binary table from the selected hdu.
			\details 
			@param Tname - table name (hdu name)
			@param col - column number to read from the table (first column index is 1)
			@param data - pointer the to array pointer where the read data will be placed in column-major order
			@param size - total size of the data array
			@param cols - not used for one dimensional columns. Set in case of multidimensional columns.
			In this case the data array is in row-major format.
			@return
			
			It is assumed that the data are all convertible to double numbers and that there are no ascii/string type columns in the 
			array.
		
			\date Mar 6, 2013, 1:08:26 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes readBinTableColumn(string Tname, int col, double** data, long* size, long* cols=NULL);
		/*!
			\brief inserts data into the selected column
			\details 
			@param colNo - number of the column to fill with the data. starts with 1
			@param firstRow - index of the row to start filling with (from 1)
			@param data - data to be put into the column
			@return
		
			\date Apr 11, 2012, 11:56:12 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addBinTableData(int colNo, long firstRow, const cpedsList<double>& data);
		cpedsStatusCodes addBinTableData(int colNo, long firstRow, const cpedsList<long>& data);
		cpedsStatusCodes addBinTableData(int colNo, long firstRow, const cpedsList<string>& data);
//		/*!
//			\brief inserts 2d array data into the selected column
//			\details 
//			@param colNo - number of the column to fill with the data. starts with 1
//			@param firstRow - index of the row to start filling with (from 1)
//			@param data - 2D data to be put into the column in a row-major order
//			@param sizeN0 - size of the 0th dimension. This is used to convert the linear array into 2D C-style array
//			@return
//		
//			\date Jul 9, 2012, 1:48:45 PM
//			\author Bartosz Lew
//		*/
//		cpedsStatusCodes addBinTableData2D(int colNo, long firstRow, cpedsList<double>& data, long sizeN0);
//		cpedsStatusCodes addBinTableData(int colNo, long firstRow, cpedsList<bool>& data);

//		cpedsStatusCodes saveFunctionsListAsBinArray(QList<mscsFunction> data, QList<string> colName, QList<string> colType, QList<string> colUnit);
		
		/*!
			\brief save a function to an opened fits file to a new binary table extension
			\details  
			@param Tname - extension name
			@param fn - function to save
			@param argName - 1st column name
			@param valName - 2nd column name
			@param argUnit - 1st column unit name
			@param valUnit - 2nd column unit name
			@return
		
			\date Jul 15, 2013, 11:46:53 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveFunction(const string Tname, const mscsFunction& fn,const string argName, const string valName, const string argUnit, const string valUnit);

		mscsFunction readFunction(string Tname);
		/*!
			\brief save a 2d array to an opened fits file to a new binary table extension
			\details  
			@param Tname - extension name
			@param fn - function to save
			@param colName - list with column names
			@param colUnit - list with column unit names
			@return
			
			The size of the lists must be the same as the number of columns in the input function.
		
			\date Jul 15, 2013, 11:46:53 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveArray(const string Tname, const mscsFunction3dregc& fn,QList<string> colName, QList<string> colUnit);
		
		/*!
			\brief insets column data into colPosition in HDU hduName
			\details 
			@param hduName - name of bin table extension
			@param colPosition - position in which a new column should be placed (eg. 1, to place before all existing columns)
			@param col - column data
			@param colName - column name
			@param colUnit - column unit 
			@return
		
			\date Jul 15, 2013, 11:05:11 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes insertColumn(const string hduName, long colPosition, cpedsList<double> col, string colName, string colUnit);

		
		
		/*!
			\brief adds a keyword to the current HDU
			\details 
			@param
			@return
		
			\date Nov 19, 2011, 10:43:53 PM
			\author Bartosz Lew
		 */
		cpedsStatusCodes addKey(string key, string value, string comment, bool isLong=false);
		cpedsStatusCodes addKeyCH(char* key, char* value, char* comment, long N=100, bool isLong=false);
		cpedsStatusCodes addKey(string key, double value, string comment);
		cpedsStatusCodes addKey(string key, long value, string comment);
		cpedsStatusCodes addKey(string key, bool value, string comment);
		
		/*!
			\brief add a comment key into the current extension
			\details 
			@param comment - your comment
			@param section - if true then the comment will be sandwiched between two comment lines made up of lines --------------- to create
			a better visible section in the fits header
			@return
		
			\date Nov 20, 2011, 2:02:05 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addComment(string comment, bool section=false, string sectionStr="----------------------------------------");

		/*!
			\brief returns keyword values from the requested HDU by w given name
			\details 
			@param key - name of the key
			@param status - pointer to the cfitsio return status code
			@param comment - pointer to existing string if comment was requested
			@param unit - pointer to existing string if unit was requested
			@param HDU - name of the extension to be used, by default the current extension is used.
			@return
		
			\date Nov 19, 2011, 10:44:56 PM
			\author Bartosz Lew
		 */
		string getStringKey(string key, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0);
		string getStringKeyLong(string key, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0);

		double getDoubleKey(string key, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0);
		long getLongKey(string key, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0);
		bool getBoolKey(string key, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0);
		
		/*!
			\brief returns all header rows as a Qlist 
			\details 
			@param
			@return
		
			\date Nov 19, 2011, 10:59:02 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes getKeys(QList<headerRow_t> keys, string HDU="");
		
		long getKeysCount(int HDU=1);
		
		/*!
			\brief moves the file pointer to the requested extension given by name, type and version
			\details 
			@param HDU - extension name
			@param hduType - "image" - for image, "binTable" - for binary table and "asciiTable" - for ascii table
			@param hduVersion - version of the extension. if 0 then ignored
			@return
		
			\date Nov 21, 2011, 9:13:36 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes selectHDU(string HDU, string hduType, int hduVersion=0);
		
		cpedsStatusCodes selectHDU(int HDU);
		
		/*!
			\brief return the number of the current hdu
			\details 
			@return
		
			\date Sep 4, 2012, 6:27:09 PM
			\author Bartosz Lew
		*/
		int getHDUnum();
		cpedsStatusCodes createNewHDU();
		cpedsStatusCodes createPrimaryHDU();
		
		/*!
			\brief return the total count of HDUs in the file
			\details 
			@param
			@return
		
			\date Sep 4, 2012, 6:27:42 PM
			\author Bartosz Lew
		*/
		int getHDUCount();

		
		
		
		
		
		/*!
			\brief save data as image to fits file
			\details 
			@param IMGname - name of the extension
			@param rows - numer of rows
			@param cols - number of cols
			@param data - row major ordered data linear array
			@param quantizeLevel - what's that ?
			@return
		
			\date Aug 7, 2012, 10:45:17 AM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addIMG(string IMGname, long rows, long cols, cpedsList<double>& data, float quantizeLevel=16);
		cpedsStatusCodes addIMG(string IMGname, long rows, long cols, double* data, float quantizeLevel=16);
		/*!
			\brief add picture the fits file from 3D function slice (XY-plane) at coordinate Z
			\details 
			@param IMGname - name of the extension
			@param img - function to store
			@param Zplane - coordinate of the Z plane to save
			@return operation exit status code
			
			currently it is not possible to save other planes than XY-plane
		
			\date Jan 25, 2015, 12:31:05 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes addIMG(string IMGname, const mscsFunction3dregc& img, long Zplane=0);

		cpedsStatusCodes readIMG(string IMGname, mscsFunction3dregc& img);
		
		
		
		
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		void printLastError();
		
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
		void getKey(string key, int keyType, int* status=NULL, string* comment=NULL, string* unit=NULL, string HDU="", string hduType="binTable", int hduVersion=0, bool isLong=false);
		
		
		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		typedef struct {
				string s; //!< string data
				char *chlong; //!< long string data
				char ch[80]; //!< string data
				double d; //!< double data
				long l; //!< long data
				bool b; //!< bool data
				char comment[80]; //!< comment string
				char unit[80]; //!< unit info
				char key[80]; //!< key string
		} key_data_t;
		
		key_data_t _key_data;
		
		fitsfile* _fptr;
		string _fname;
		int _status;
		cpedsMsgs* msgs;
		char* errorText;
		
};
#endif /* CPEDSFITS */ 

