/*!
 * \file DateTimeData.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: blew
 */
#include <MscsDateTimeData.h>
//#include <QtCore/QStringList>
#include "qtbase.h"

MscsDateTimeData::MscsDateTimeData() {
//	_f=NULL;
}

MscsDateTimeData::~MscsDateTimeData() {
//	if (_f!=NULL) { fclose(_f); _f=NULL; }
}

void MscsDateTimeData::load(string fname, long dtCol, string dtFmt) {
//	cout << "loading\n";
	ifstream is;
	is.open(fname.c_str(), ios::in);
	
	if (!is.is_open()) return;
	
	_DTdata.dtCol=dtCol;
	_DTdata.dtFmt=dtFmt.c_str();
	_DTdata.fname=fname.c_str();
	
	long dtSt,dtEn;
	
	QString qdtFmt=dtFmt.c_str();
	int dtFmt_spaces=qdtFmt.count(' ');
	dtSt=dtCol;
	dtEn=dtSt+dtFmt_spaces;
//	cout << "qdtFmt: " << qdtFmt.toStdString() << "\n";
//	printf("dtSt: %li dtEn: %li\n",dtSt,dtEn);
	
	double entry;
	string line;
	while ( getline (is,line) ) {
		QString qs=line.c_str();
//		cout << qs.toStdString() << "\n";
		if (qs[0]!='#') {
			QStringList qsl=qs.simplified().split(" ");
			
			row_t row;
			QString dt;
			for (long i = 0; i < qsl.size(); i++) {
				if (i<dtSt) row.data.append(qsl[i].toDouble());
				else {
					if (i>=dtSt and i<=dtEn) {
						dt+=qsl[i]+" ";
					}
					else {
						row.data.append(qsl[i].toDouble());
					}
				}
				
			}
	//		cout << dt.toStdString() << "\n";
	//		cout << "dt: " << dt.toStdString() << " qdtFmt: "<< qdtFmt.toStdString() << "\n";
			row.dt=QDateTime::fromString(dt.simplified(),qdtFmt);
	//		cout << "test: " << row.dt.date().year() << " " << row.dt.date().month() << " " << row.dt.date().day() << " "<< row.dt.time().hour() << ":" << row.dt.time().minute() << ":" << row.dt.time().second() <<  "\n";
			_DTdata.rows.append(row);
		}
	}
	is.close();
	
//	fclose(_f);
	
//	cout << "loading done\n";
}

void MscsDateTimeData::printData() {
	for (long i = 0; i < _DTdata.rows.size(); ++i) {
		cout <<	_DTdata.rows[i].dt.toString("yyyy-MM-dd hh:mm:ss").toStdString() << " ";
		for (long j = 0; j < _DTdata.rows[i].data.size(); j++) {
			cout << _DTdata.rows[i].data[j] << " ";
		}
		cout << "\n";
	}
}
/* ******************************************************************************************** */
double MscsDateTimeData::getFirstValueAfter(string dt, long col, long startFrom, long* idx) {
	QDateTime qdt=QDateTime::fromString(dt.c_str(),_DTdata.dtFmt);
	return getFirstValueAfter(qdt,col,startFrom,idx);
}
/* ******************************************************************************************** */
double MscsDateTimeData::getFirstValueAfter(QDateTime qdt, long col, long startFrom, long * idx) {
	long i=startFrom;

	while (i<linesCount() and _DTdata.rows[i].dt<qdt) {
		i++;
	}
	double val=0;
	if (i<linesCount()) val=_DTdata.rows[i].data[col];
	if (idx!=NULL) *idx=i;
	
	return val;	

}
/* ******************************************************************************************** */
double MscsDateTimeData::getFirstValueAfter(double JD, long col, long startFrom, long * idx) {
	string dt=cpeds_JDToYMDhms(JD,"%li-%02li-%02li %02li:%02li:%02.0lf");
	QDateTime qdt=QDateTime::fromString(dt.c_str(),"yyyy-MM-dd hh:mm:ss");
	return getFirstValueAfter(qdt,col,startFrom,idx);	

}
/* ******************************************************************************************** */
double MscsDateTimeData::getTimeJD(long idx) {
	double hour=double(_DTdata.rows[idx].dt.time().hour())+
			double(_DTdata.rows[idx].dt.time().minute())/60+
			double(_DTdata.rows[idx].dt.time().second())/3600;
	return cpeds_julian_time(
			_DTdata.rows[idx].dt.date().year(),
			_DTdata.rows[idx].dt.date().month(),
			_DTdata.rows[idx].dt.date().day(),
			hour);
}
