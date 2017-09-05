/*!
 * \file cpeds-fits.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: blew
 */

#include "cpeds-fits.h"
//#include <QtCore/QString>

cpedsFits::cpedsFits(cpeds_VerbosityLevel level) {
	_status=0;
	msgs = new cpedsMsgs("cpedsFits");
	msgs->setVerbosity(level);	
	errorText=new char[30];
}
cpedsFits::cpedsFits(string fname, string action, cpeds_VerbosityLevel level) {
	_status=0;
	msgs = new cpedsMsgs("cpedsFits");
	errorText=new char[30];
	msgs->setVerbosity(level);
	if (fname!="") {
		if (openFile(fname,action)==cpedsError) exit(0);
	}
}

cpedsFits::~cpedsFits() {
	delete msgs;
	delete [] errorText;
}

cpedsStatusCodes cpedsFits::openFile(string fname, string action) {
	if (action=="new" or cpeds_fileExists(fname)==false) {
		fname="!"+fname;
		fits_create_file(&_fptr,fname.c_str(),&_status);
		if (_status!=0) { msgs->error("cpedsFits::openFile>> cannot open new fits file : "+fname,High);	_status=0; return cpedsError; }
	}
	else
	if (action=="append") {
		fits_open_file(&_fptr, fname.c_str(), READWRITE, &_status);		
		if (_status!=0) { msgs->error("cpedsFits::openFile>> cannot open file: "+fname,High);	_status=0; return cpedsError; }
		int hduType,hduNum;
		fits_get_num_hdus(_fptr,&hduNum,&_status);
		if (_status!=0) { msgs->error("cpedsFits::openFile>> cannot obtain number of HDUs ("+msgs->toStr(_status)+") for file "+fname,High);	_status=0; return cpedsError; }
		fits_movabs_hdu(_fptr, hduNum, &hduType, &_status);
		if (_status!=0) { msgs->error("cpedsFits::openFile>> cannot move to the last extension in file: "+fname+" The number of extensions in this file is: "+msgs->toStr(hduNum),High);	_status=0; return cpedsError; }
	}
	else
	if (action=="read") {
		fits_open_file(&_fptr, fname.c_str(), READONLY, &_status);				
		if (_status!=0) { msgs->error("cpedsFits::openFile>> cannot open fits file for reading: "+fname,High);	_status=0; return cpedsError; }
	}
	return cpedsSuccess; 
}

/***************************************************************************************/
cpedsStatusCodes cpedsFits::addBinTable(string Tname, long rows, QList<string> &colName, QList<string>  &colType, QList<string>  &colUnit) {
	long cols=colName.size();
	char *tname[cols];
	char *tform[cols];
	char *tunit[cols];
	char buff[100];
	long w;

	// setup stuff for cfitsio
	
	for (long i = 0; i < cols; i++) {
		tname[i] = new char[8];
		strcpy(tname[i],colName[i].c_str());
		if (colUnit[i]=="") tunit[i]=NULL; else { tunit[i] = new char[20]; strcpy(tunit[i],colUnit[i].c_str()); }
	}

	
	for (long i = 0; i < cols; i++) {
		tform[i] = new char[8];
//		printf("--> %s\n", colType[i].c_str());
		if (colType[i]=="double") strcpy(tform[i],"1D");
		else
			if (colType[i].find("double",0)!=string::npos) { 
				strcpy(buff,colType[i].replace(0,6,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"%liD",w);
			}
		if (colType[i]=="float") strcpy(tform[i],"1E");
		else
			if (colType[i].find("float",0)!=string::npos) { 
				strcpy(buff,colType[i].replace(0,5,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"%liE",w);
			}
		if (colType[i]=="long") strcpy(tform[i],"1J");
		else
			if (colType[i].find("long",0)!=string::npos) { 
				strcpy(buff,colType[i].replace(0,4,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"%liJ",w);
			}
		if (colType[i].find("string",0)!=string::npos) { 
			strcpy(buff,colType[i].replace(0,6,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"%liA",w);
//			printf("adding string field of width: %li\n",w);
		}
		
		if (colType[i]=="bool") strcpy(tform[i],"1L");
		if (colType[i].find("Pfloat",0)!=string::npos) { 
			strcpy(buff,colType[i].replace(0,6,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"PE(%li)",w);
		}
		if (colType[i].find("Pdouble",0)!=string::npos) { 
			strcpy(buff,colType[i].replace(0,7,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"PD(%li)",w);
		}
		if (colType[i].find("Plong",0)!=string::npos) { 
			strcpy(buff,colType[i].replace(0,5,"").c_str()); w=strtol(buff,NULL,10); sprintf(tform[i],"PJ(%li)",w);
		}
//		printf("->> %s\n",tform[i]);
	}
	
	// create binary table
//	int hduType,hduNum;
	
//	fits_movabs_hdu(_fptr, fits_get_num_hdus(_fptr,&hduNum,&_status)+1, &hduType, &_status);
//	if (_status) { msgs->error("addBinTable>> cannot create table: "+Tname,High);	return cpedsError;	}
	fits_create_tbl(_fptr, BINARY_TBL, rows, colName.size(), tname, tform, tunit, Tname.c_str(), &_status);
	
	// free memory
	for (long i = 0; i < cols; i++) {
		delete [] tname[i];
		delete [] tform[i];
		if (tunit[i]!=NULL) delete [] tunit[i];
	}
	if (_status!=0) { msgs->error("addBinTable (error "+msgs->toStr(_status)+") >> cannot create table: "+Tname,High);	_status=0; return cpedsError;	}
	
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addBinTable(string Tname) {
	QList<string> cols;
	return addBinTable(Tname,0,cols,cols,cols);
}

/***************************************************************************************/
cpedsStatusCodes cpedsFits::addBinTableData(int colNo, long firstRow, const cpedsList<double>& data) {
	double* d=data.toCarray();
	fits_write_col(_fptr, TDOUBLE, colNo, firstRow, 1, data.size(),d,&_status);
	delete [] d;
	if (_status!=0) { msgs->error("addBinTableData (error "+msgs->toStr(_status)+") >> cannot add data into column no,: "+msgs->toStr(colNo),High);	printLastError(); _status=0; return cpedsError;	}
	return cpedsSuccess;
}
/***************************************************************************************/
//cpedsStatusCodes cpedsFits::addBinTableData2D(int colNo, long firstRow, cpedsList<double>& data, long sizeN0) {
//	
//}
cpedsStatusCodes cpedsFits::addBinTableData(int colNo, long firstRow, const cpedsList<long>& data) {
	long* d=data.toCarray();
	fits_write_col(_fptr, TLONG, colNo, firstRow, 1, data.size(),d,&_status);
	delete [] d;
	if (_status!=0) { msgs->error("addBinTableData (error "+msgs->toStr(_status)+") >> cannot add data into column no,: "+msgs->toStr(colNo),High);	printLastError(); _status=0; return cpedsError;	}
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addBinTableData(int colNo, long firstRow, const cpedsList<string>& data) {
	char** d = new char*[data.size()];
	for (long i=0;i<data.size();i++) { d[i] = new char[data[i].size()]; strcpy(d[i],data[i].c_str()); } //printf("test: %s\n",d[i]); }
	fits_write_col(_fptr, TSTRING, colNo, firstRow, 1, data.size(),d,&_status);

	for (long i=0;i<data.size();i++) { delete [] d[i]; }	delete [] d;
	if (_status!=0) { msgs->error("addBinTableData (error "+msgs->toStr(_status)+") >> cannot add data into column no,: "+msgs->toStr(colNo),High);	printLastError(); _status=0; return cpedsError;	}
	return cpedsSuccess;
}
//cpedsStatusCodes cpedsFits::addBinTableData(int colNo, long firstRow, cpedsList<bool>& data) {
//	double* d=data.toCarray();
//	fits_write_col(_fptr, TLOGICAL, colNo, firstRow, 1, data.size(),d,&_status);
//	delete [] d;
//	if (_status!=0) { msgs->error("addBinTableData (error "+msgs->toStr(_status)+") >> cannot add data into column no,: "+msgs->toStr(colNo),High);	return cpedsError;	}
//	return cpedsSuccess;
//}


cpedsStatusCodes cpedsFits::addKey(string key, string value, string comment, bool isLong) {
#ifdef DEBUG_FITS
    printf("addKey:: saving string fits key. The key name is: |%s|, value: |%s|, comment: |%s|\n",key.c_str(),value.c_str(),comment.c_str());
#endif
	if (key.length()==0) { msgs->error("addKey >> key length cannot be zero",High); return cpedsError; }
	if (value.length()==0) { msgs->warning("addKey >> value length is zero, will use the default value: ''",High); }
	char keych[key.length()];
	strcpy(keych,key.c_str());

	const long N=80;
//	N=value.length();
//	if (N==0) N=1;
//	char valch[N];
//	bzero(valch,N);
	char valch[N];
	bzero(valch,N);

//	N=comment.length();
//	if (N==0) N=1;
	char comm[N];
	bzero(comm,N);	
	
	if (value.length()==0) { 
//		printf("length is 0\n"); 
		strcpy(valch,""); 
	} 
	else {
		strcpy(valch,value.c_str());  
	} 
//	printf("valch |%s|, value: |%s|, N: %li\n",valch,value.c_str(),value.length()); 
	
	if (comment.length()==0) { strcpy(comm,""); } 
	else {strcpy(comm,comment.c_str()); } 
//	printf("comm |%s|, comment: |%s|, N: %li\n",comm,comment.c_str(),comment.length()); 
	return addKeyCH(keych,valch,comm,N,isLong);
//	fits_write_key(_fptr,TSTRING, key.c_str(), (char*)(buff), comment.c_str(), &_status);
//	fits_write_key(_fptr,TSTRING, key.c_str(), buff, comment.c_str(), &_status);
//	fits_write_key(_fptr,TSTRING, key.c_str(), &c, comment.c_str(), &_status);
}
cpedsStatusCodes cpedsFits::addKeyCH(char* key, char* value, char* comment, long N, bool isLong) {
//	printf("adding string key: |%s|, value |%s|, comment: |%s|, isLong: %i\n",key,value,comment,int(isLong));
	if (isLong) {
		fits_update_key_longstr(_fptr, key, value, comment, &_status);
		if (_status!=0) { msgs->error("addKey>> cannot add long key: "+msgs->toStr(key)+" with value: "+msgs->toStr(value)+" (Error: "+msgs->toStr(_status)+")",High);	_status=0; return cpedsError;	}
	}
	else {
		fits_update_key(_fptr,TSTRING, key, value, comment, &_status);
		if (_status!=0) { msgs->error("addKey>> cannot add key: "+msgs->toStr(key)+" with value: "+msgs->toStr(value)+" (Error: "+msgs->toStr(_status)+")",High);	_status=0; return cpedsError;	}		
	}
	return cpedsSuccess;		
}

cpedsStatusCodes cpedsFits::addKey(string key, double value, string comment) {
	fits_write_key(_fptr,TDOUBLE, key.c_str(), &value, comment.c_str(), &_status);
	if (_status!=0) { msgs->error("addKey>> cannot add key: "+key+" with value: "+msgs->toStr(value)+" (Error: "+msgs->toStr(_status)+")",High);	_status=0; return cpedsError;	}
	return cpedsSuccess;	
}
cpedsStatusCodes cpedsFits::addKey(string key, long value, string comment) {
	fits_write_key(_fptr,TLONG, key.c_str(), &value, comment.c_str(), &_status);
	if (_status!=0) { msgs->error("addKey>> cannot add key: "+key+" with value: "+msgs->toStr(value)+" (Error: "+msgs->toStr(_status)+")",High); _status=0;	return cpedsError;	}
	return cpedsSuccess;	
}
cpedsStatusCodes cpedsFits::addKey(string key, bool value, string comment) {
	fits_write_key(_fptr,TLOGICAL, key.c_str(), &value, comment.c_str(), &_status);
	if (_status!=0) { msgs->error("addKey>> cannot add key: "+key+" with value: "+msgs->toStr(value)+" (Error: "+msgs->toStr(_status)+")",High);	_status=0; return cpedsError;	}
	return cpedsSuccess;	
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addComment(string comment, bool section, string sectionStr) {

	if (section) { 
		fits_write_comment(_fptr,sectionStr.c_str(),&_status);
	}
	fits_write_comment(_fptr,comment.c_str(),&_status);
	if (_status!=0) { msgs->error("addComment>> cannot add comment: "+comment+" (Error: "+msgs->toStr(_status)+")",High);	_status=0; return cpedsError;	}
	if (section) { 
		fits_write_comment(_fptr,sectionStr.c_str(),&_status);
	}
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::closeFile() {
	fits_close_file(_fptr,&_status);
	return cpedsSuccess;	
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::selectHDU(string HDU, string hduType, int hduVersion) {
	int htype;
	if (hduType=="image") htype=IMAGE_HDU;
	if (hduType=="binTable") htype=BINARY_TBL;
	if (hduType=="asciiTable") htype=ASCII_TBL;
	char h[80]; strcpy(h,HDU.c_str());
//	printf("moving to: %s\n",h);
	fits_movnam_hdu(_fptr, htype, h, hduVersion, &_status);	
	if (_status!=0) { msgs->error("selectHDU>> cannot find requested extension name: "+HDU,High); printLastError();	_status=0; return cpedsError;	}
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::selectHDU(int HDU) {
//	printf("moving to: %i\n",HDU);
	int htype;
	fits_movabs_hdu(_fptr, HDU, &htype, &_status);	
	if (_status!=0) { msgs->error("selectHDU>> cannot find requested extension ID: "+msgs->toStr(HDU),High); printLastError();	_status=0; return cpedsError;	}
	return cpedsSuccess;
	
}
/***************************************************************************************/
int cpedsFits::getHDUnum() {
	int i;
	fits_get_hdu_num(_fptr,&i);
	return i;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::createNewHDU() {
	fits_create_hdu(_fptr,&_status);
	if (_status!=0) { msgs->error("selectHDU>> cannot create a new HDU",High); printLastError();	_status=0; return cpedsError;	}
	return cpedsSuccess;
}
cpedsStatusCodes cpedsFits::createPrimaryHDU() {
	fits_create_hdu(_fptr,&_status);
	if (_status!=0) { msgs->error("selectHDU>> cannot create a new HDU",High);	printLastError();	_status=0; return cpedsError;	}
	addKey("SIMPLE",true,"");
	addKey("BITPIX", long(8),"");
	addKey("NAXIS", long(0),"");
	addKey("EXTEND", true,"");
	addKey("BLOCKED", true,"");
//	if (_status!=0) { msgs->error("selectHDU>> cannot create a new HDU",High); printLastError();	_status=0; return cpedsError;	}
	return cpedsSuccess;
}

/***************************************************************************************/
void cpedsFits::getKey(string key, int keyType, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion, bool isLong) {
	if (HDU!="") { selectHDU(HDU,hduType,hduVersion); }
	strcpy(_key_data.key,key.data());
	if (isLong) {
		fits_read_key_longstr(_fptr,  _key_data.key,&_key_data.chlong, _key_data.comment,&_status); 
		
		if (_status==0) { 
			_key_data.s=_key_data.chlong; 
			delete [] _key_data.chlong; _key_data.chlong=NULL; 
		}
	}
	else {
		if (keyType==TSTRING) { fits_read_key(_fptr, keyType, _key_data.key,_key_data.ch, _key_data.comment,&_status); _key_data.s=_key_data.ch; }
		if (keyType==TDOUBLE) fits_read_key(_fptr, keyType, _key_data.key,&_key_data.d, _key_data.comment,&_status);
		if (keyType==TLONG) fits_read_key(_fptr, keyType, _key_data.key,&_key_data.l, _key_data.comment,&_status);
		if (keyType==TLOGICAL) fits_read_key(_fptr, keyType, _key_data.key,&_key_data.b, _key_data.comment,&_status);
	}
	if (status!=NULL) *status=_status;
	if (_status!=0) {
		msgs->error("getKey>> cannot find requested key name ("+key+") in this extension ("+HDU+").",High);	
		_status=0;
	}

	if (comment!=NULL) comment->assign(_key_data.comment);
	
	if (unit!=NULL) {
		fits_read_key_unit(_fptr, _key_data.key, _key_data.unit,&_status);
		unit->assign(_key_data.unit);
	}
}
/***************************************************************************************/
string cpedsFits::getStringKey(string key, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion) {
	getKey(key,TSTRING,status,comment,unit,HDU, hduType, hduVersion);
	return _key_data.s;
}
/***************************************************************************************/
string cpedsFits::getStringKeyLong(string key, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion) {
	getKey(key,TSTRING,status,comment,unit,HDU, hduType, hduVersion,true);
	return _key_data.s;
}
/***************************************************************************************/
double cpedsFits::getDoubleKey(string key, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion) {
	getKey(key,TDOUBLE,status,comment,unit,HDU, hduType, hduVersion);
	return _key_data.d;
}
/***************************************************************************************/
long cpedsFits::getLongKey(string key, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion) {
	getKey(key,TLONG,status,comment,unit,HDU, hduType, hduVersion);
	return _key_data.l;
}
/***************************************************************************************/
bool cpedsFits::getBoolKey(string key, int* status, string* comment, string* unit, string HDU, string hduType, int hduVersion) {
	getKey(key,TLOGICAL,status,comment,unit,HDU, hduType, hduVersion);
	return _key_data.b;
}
/***************************************************************************************/
void cpedsFits::printLastError() {
	fits_get_errstatus(_status,errorText);
	cout << "\n	cfitsio error description: " << errorText << "\n";
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addIMG(string IMGname, long rows, long cols, cpedsList<double>& data, float quantizeLevel) {
	double* d=data.toCarray();
	cpedsStatusCodes s=addIMG(IMGname,rows,cols,d);
	delete [] d;	
	return s;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addIMG(string IMGname, long rows, long cols, double* data, float quantizeLevel) {
	fits_set_quantize_level(_fptr,quantizeLevel,&_status);
	long naxis=2;
	long naxes[2] = { cols, rows }; 
	fits_create_img(_fptr,  DOUBLE_IMG, naxis, naxes, &_status);
	if (_status!=0) { msgs->error("addIMG>> create img",High);	printLastError();	_status=0; return cpedsError;	}
	if (data!=NULL)	
		fits_write_img(_fptr,TDOUBLE, 1,rows*cols,data,&_status);
	if (_status!=0) { msgs->error("addIMG>> write img",High);	printLastError();	_status=0; return cpedsError;	}
	addKey("EXTNAME",IMGname,"name of this binary image extension",false);
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::readIMG(string IMGname, mscsFunction3dregc& img) {
	selectHDU(IMGname,"image");
	
//	fits_set_quantize_level(_fptr,quantizeLevel,&_status);
	long cols=getLongKey("NAXIS1",&_status);
	long rows=getLongKey("NAXIS2",&_status);
	printf("%li %li\n",cols,rows);
//	long naxes[2] = { cols, rows }; 
//	fits_create_img(_fptr,  DOUBLE_IMG, naxis, naxes, &_status);
//	if (_status!=0) { msgs->error("addIMG>> create img",High);	printLastError();	_status=0; return cpedsError;	}
//	if (data!=NULL)	
	long nelem=rows*cols;
	double* data = new double[nelem];
	fits_read_img(_fptr,TDOUBLE, 1,nelem,0,data,0,&_status);
	if (_status!=0) { msgs->error("addIMG>> write img",High);	printLastError();	_status=0; delete [] data; return cpedsError;	}

//	img.setSize(cols,rows);
//	img.setSizeRange(cols,rows,1,0,0,0,cols,rows,0);
	img.setSize(cols,rows,1,1,1,1,0,0,0);
	img.allocFunctionSpace();
	long idx=0;
	for (long j = 0; j < rows; j++) {
		for (long i = 0; i < cols; i++) {
			img.setf(i,j,0,data[idx],0);
			idx++;
		}
	}
	delete [] data;
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::addIMG(string IMGname, const mscsFunction3dregc& img, long Zplane) {
	float quantizeLevel=16;
	fits_set_quantize_level(_fptr,quantizeLevel,&_status);
	long naxis=2;
	long naxes[2] = { img.Nx(), img.Ny() }; 
	fits_create_img(_fptr,  DOUBLE_IMG, naxis, naxes, &_status);
	if (_status!=0) { msgs->error("addIMG>> create img",High);	printLastError();	_status=0; return cpedsError;	}
	double *data=img.exportRe("ZYXmajor");
	if (data!=NULL)	{
		fits_write_img(_fptr,TDOUBLE, 1,img.Nx()*img.Ny(),data,&_status);
		delete [] data;
	}
	if (_status!=0) { msgs->error("addIMG>> write img",High);	printLastError();	_status=0;  return cpedsError;	}
	addKey("EXTNAME",IMGname,"name of this binary image extension",false);
	return cpedsSuccess;
}


/***************************************************************************************/
long cpedsFits::getKeysCount(int HDU) {
	selectHDU(HDU);
	int keysCount;
	int moreKeys;
	fits_get_hdrspace(_fptr,  &keysCount, &moreKeys, &_status); 
	if (_status!=0) { msgs->error("getKeysCount>> ",High);	printLastError();	_status=0; exit(0);	}
	return keysCount;
}

/***************************************************************************************/
int cpedsFits::getHDUCount() {
	int hduCount;
	fits_get_num_hdus(_fptr,&hduCount,&_status);
	if (_status!=0) { msgs->error("getHDUCount>> write img",High);	printLastError();	_status=0; exit(0);	}
	return hduCount;	
}

/***************************************************************************************/
//cpedsStatusCodes cpedsFits::saveFunctionsListAsBinArray(QList<cpedsList> data, QList<string> colName, QList<string> colType, QList<string> colUnit) {
//	
//}
cpedsStatusCodes cpedsFits::readBinTable(string Tname, double** data, long* size, long *rows) {
	selectHDU(Tname,"binTable");
	long Nr, Nc;
	string comment,unit;
	Nc=getLongKey("TFIELDS",&_status,&comment,&unit,Tname,"binTable");
	Nr=getLongKey("NAXIS2",&_status,&comment,&unit,Tname,"binTable");
	long N=Nc*Nr;
//	printf("Nc: %li\n",Nc);
//	printf("Nr: %li\n",Nr);
	*data = new double[N];
	double* nularray=new double[Nr];
	int* anynul=new int[Nr];
	for (long j = 1; j <= Nc; j++) {		
		fits_read_col(_fptr,TDOUBLE,j,1,1,Nr,nularray,&((*data)[(j-1)*Nr]),anynul,&_status);
		if (_status!=0) { msgs->error("getHDUCount>> write img",High);	printLastError();	_status=0; exit(0);	}
	}
	*size=N;
	*rows=Nr;
//	for (long i = 0; i < N; i++) {
//		printf("%li %lf\n",i+1, (*data)[i]);
//	}
	
	delete [] nularray;
	delete [] anynul;
	
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::readBinTableColumn(string Tname, int col, double** data, long* size, long* cols) {
	selectHDU(Tname,"binTable");
	long Nr;
	string comment,unit;
	Nr=getLongKey("NAXIS2",&_status,&comment,&unit,Tname,"binTable");
	string colfmt="TFORM"+msgs->toStr(col);
	string colType=getStringKey(colfmt,&_status,&comment,&unit,Tname,"binTable");
	printf("colType: %s\n",colType.c_str());
	// check what format do we have
	if (colType.size()<=2) { // eg 1D, then we have simple table column
		long N=Nr;
		*data = new double[N];
		double* nularray=new double[Nr];
		int* anynul=new int[Nr];
		fits_read_col(_fptr,TDOUBLE,col,1,1,Nr,nularray,*data,anynul,&_status);
		if (_status!=0) { msgs->error("readBinTableColumn>> fits_read_col",High);	printLastError();	_status=0; exit(0);	}
		*size=N;
		delete [] nularray;
		delete [] anynul;
	}
	else { // this is multidimensional array
		// figure out how many columns there are in this multidimensional array
		// I'm assuming the format is xY, where x is vector size and Y is element type -- assumed to be TDOUBLE, so
		string s=colType.substr(0,colType.size()-1);
		int vSize=strtol(s.c_str(),NULL,10);
		
		long N=Nr*vSize;
		*data = new double[N];
		double* nularray=new double[N];
		int* anynul=new int[N];
		fits_read_col(_fptr,TDOUBLE,col,1,1,N,nularray,*data,anynul,&_status);
		if (_status!=0) { msgs->error("readBinTableColumn>> fits_read_col",High);	printLastError();	_status=0; exit(0);	}
		*size=N;
		*cols=vSize;
		delete [] nularray;
		delete [] anynul;
	}
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::saveFunction(const string Tname, const mscsFunction& fn,const string argName, const string valName, const string argUnit, const string valUnit) {
	QList<string> colName, colType, colUnit;
	colName << argName << valName;
	colType << "double" << "double";
	colUnit << argUnit << valUnit;
	addBinTable(Tname,fn.pointsCount(),colName,colType,colUnit);
	addBinTableData(1,1,fn.toXList());
	addBinTableData(2,1,fn.toYList());
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::saveArray(const string Tname, const mscsFunction3dregc& fn,QList<string> colName, QList<string> colUnit) {
	if (fn.Nz()!=1) return cpedsError;
	QList<string> colType;
	for (unsigned long i = 0; i < colName.size(); i++) {
		colType << "double";
	}
	fn.printInfo();
	addBinTable(Tname,fn.Ny(),colName,colType,colUnit);
	for (unsigned long i = 0; i < fn.Nx(); i++) {
		addBinTableData(i+1,1,fn.getSlice1Dfn(1,0,i,0).toYList());
//		fn.getSlice1Dfn(1,0,i,0).print();
	}
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes cpedsFits::insertColumn(const string hduName, long colPosition, cpedsList<double> col, string colName, string colUnit) {
	selectHDU(hduName,"binTable");
	
	char tname[8];
	char tform[8];
	char tunit[20];
	char buff[100];
	long w;

	bzero(tname,8);
	bzero(tform,8);
	bzero(tunit,20);
	// setup stuff for cfitsio
	
	strcpy(tname,colName.c_str());
	if (colUnit=="") tunit[0]=0; else { strcpy(tunit,colUnit.c_str()); }
	
	strcpy(tform,"1D");

	fficol(_fptr,colPosition,tname,tform,&_status);
	if (_status!=0) { msgs->error("addBinTableData (error "+msgs->toStr(_status)+") >> cannot insert column no,: "+msgs->toStr(colPosition),High);	printLastError(); _status=0; return cpedsError;	}

	addKey("TUNIT"+msgs->toStr(colPosition),colUnit,"physical unit of field",false);
	return addBinTableData(colPosition,1,col);
}

/***************************************************************************************/
mscsFunction cpedsFits::readFunction(string Tname) {
//	selectHDU(Tname,"binTable");
	long nX,nY;
	double *dataX,*dataY;
	readBinTableColumn(Tname,1,&dataX,&nX);
	readBinTableColumn(Tname,2,&dataY,&nY);
	if (nX!=nY) msgs->criticalError("cpedsFits::readFunction, column sizes are different",Top);
	mscsFunction f("",dataX,dataY,nX,msgs->getVerbosity());
	delete [] dataX;
	delete [] dataY;
	return f;
}
