/*!
 * \file test-fits.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: blew
 */

#include "cpeds-fits.h"
#include "cpeds-list.h"
#include "Mscs-function3dregc.h"
#include "string.h"


int main() {
//	string fname="plikTestowy.fits";
//	cpedsFits f(fname,"new");
//	QList<string> colName;
//	colName.append("string");
//	colName.append("double");
//	colName.append("long");
//	colName.append("bool");
//	QList<string> colType;
//	colType.append("string20");
//	colType.append("double");
//	colType.append("long");
//	colType.append("bool");
//	QList<string> colUnit;
//	colUnit.append("");
//	colUnit.append("m/s");
//	colUnit.append("erg/m^3");
//	colUnit.append("");
//	
//	f.addBinTable("scan_information",0,colName,colType,colUnit);
//	f.addComment("telescope information",true);
//	f.addKey("key1",string("text"),"to jest pole string");
//	f.addKey("key2",double(12.4),"[erg*s^-1] to jest pole double");
//	f.addKey("key3",long(12),"to jest pole long");
//	f.addKey("key4",true,"to jest pole bool:true");
//	f.addKey("key5",false,"to jest pole bool:false");
//	f.addKey("info","aja;ijefalj iaweio jawfjoawiejf aoiwfj aowfj aowijf a;iowjf ;awjf a;ij fa;wij fa;wjf ;awifj a;wjf a;wfj a;wfj a;wfj a;w","to jest bardzo dlugie pole",true);
//	
//	cpedsList<double> data;
//	data.append(1);
//	data.append(2);
//	data.append(3);
//	data.append(4);
//	f.addBinTableData(2,1,data);
//	cpedsList<long> data2;
//	data2.append(1);
//	data2.append(2);
//	data2.append(3);
//	data2.append(4);
//	
//	f.addBinTableData(3,1,data2);
//	
//	
//	f.closeFile();
//	f.openFile(fname,"append");
////	f.addBinTable("Naglowek2",0,colName,colType,colUnit);
//	
//	f.addKey("key11",string("text"),"to jest pole string");
//	f.addKey("key22",double(12.4),"to jest pole double");
//	f.addKey("key33",long(12),"to jest pole long");
//	f.addKey("key44",true,"to jest pole bool:true");
//	f.addKey("key55",false,"to jest pole bool:false");
//
//	f.closeFile();
//	f.openFile(fname,"read");
//	f.selectHDU("scan_information","binTable");
////	string mykeyStr=f.getStringKey("key1");
////	printf("read key is: %s\n",mykey.c_str());
//
//	int stat;
//	string kName="key2";
//	string kUnit;
//	string kComment;
//	double mykeyD=f.getDoubleKey("key2",&stat,&kComment,&kUnit);
//	printf("read key name is: %s\n",kName.c_str());
//	printf("read key value is: %lf\n",mykeyD);
//	printf("read key comment is: %s\n",kComment.c_str());
//	printf("read key unit is: %s\n",kUnit.c_str());
//	string mykeyS=f.getStringKeyLong("info",&stat,&kComment,&kUnit);
//	printf("read key value is: %s\n",mykeyS.c_str());
//	printf("read key comment is: %s\n",kComment.c_str());
//	printf("read key unit is: %s\n",kUnit.c_str());
//	
//	f.closeFile();
	
	
	cpedsFits ff;
//	ff.openFile(string("/home/blew/projects/APRICOT/scanning-optimization/grids/gridShapeSmall/gridShapeSmall/location-lon_18.564-lat_53.095/shape-0.500_x_8.000/st_0-ptsOrd_0-skyTracking-H_1-V_1/scan-trajectory.txt-gridSpacing_0.010-tint_0.240-dir_341.fits"),"read");
	mscsFunction3dregc m;
//	ff.readIMG(string("tobs"),m);
//	ff.closeFile();
////	m.printFunction();
//	m.saveSlice(2,0,"readImgtest",0);
	
	printf("write test\n");
	m.setSizeRange(100,50,1,0,0,0,100,50,1);
	m.allocFunctionSpace();
	m.mkBall3D(5);
	m.shift(-40,-20,0);
	m.saveSlice(2,0,"IMGtest.txt",0);
	ff.openFile("fitsIMGtest.fits","new");
	ff.createPrimaryHDU();
	ff.addIMG("test",m);
	double *data=m.exportRe("ZYXmajor");
	ff.addIMG("test2",m.Ny(),m.Nx(),data);
	ff.closeFile();
	printf("read test\n");
	ff.openFile("fitsIMGtest.fits","read");
//	exit(0);
	ff.readIMG(string("test2"),m);	
	
//	exit(0);

	ff.closeFile();
	m.saveSlice(2,0,"IMGtest-read.txt",0);
	
	
	
	return 0;
}
