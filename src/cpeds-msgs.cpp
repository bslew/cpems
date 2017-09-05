#include "cpeds-msgs.h"
#include "sstream"
#include "cpeds-math.h"

/* ************************************************************************************************************************************************************************************ */

cpedsMsgs::cpedsMsgs(string sender) {
	msgSender=sender; 
	logging = false;
	messaging=true;
	logFileName="";
	setVerbosity(CPEDS_defaultVerbosityLevel);
	setLogVerbosity(CPEDS_defaultVerbosityLevel);
	dateTimeOn=true;
	timeSinceStartOn=false;
	logFile=NULL;
	prefix="";
	
	_tmpMessage = new char[_maxialMessageLength];
	_singleMessage="";
	_saveRunWriteMode='w';
	_startTime=time(NULL);
}


cpedsMsgs::cpedsMsgs(string sender, bool log, string logfilename, cpeds_VerbosityLevel verb) {
	msgSender=sender; 
	logging = log;
	messaging=true;
	logFileName=logfilename;
//	if (logfilename!="") {
//		if (cpeds_fileExists(logfilename)) remove(logfilename.c_str());
//	}
	setVerbosity(verb);
	setLogVerbosity(verb);
	dateTimeOn=true;
	timeSinceStartOn=false;
	logFile=NULL;
	prefix="";
	/* if (logfilename!="") logFile=fopen(logFileName.c_str(),"w"); //!< TODO: make this smarter to generate an alternative name in case many threads of the same program want to write th the same logfile at the same time. Add functionality via constructor params. to do overwrights or appends */
	
	_tmpMessage = new char[_maxialMessageLength];
	_singleMessage="";
	_saveRunWriteMode='w';
	_startTime=time(NULL);
}

cpedsMsgs::cpedsMsgs(const cpedsMsgs& parent) { //TODO: make this stuff go into the operator=
	msgSender=parent.getSender();
	logging = false; //parent.getLogging(); 
	messaging=parent.getMessaging();
	logFileName=parent.getLogFileName(); // TODO: this could add .clone to the original name and go on with logging if logging is on.
	setVerbosity(parent.getVerbosity());
	setLogVerbosity(parent.getLogVerbosity());
	dateTimeOn=parent.dateTimeStatus();
	timeSinceStartOn=parent.getTimeSinceStartOn();
	logFile=NULL; //parent.logFile; - even if the parent is logging we cannot write to the same file anyway
	prefix=parent.pref();  
	/* if (logging) startLogging(); else stopLogging(); */

	_tmpMessage = new char[_maxialMessageLength];
	_singleMessage="";
	_saveRunWriteMode=parent._saveRunWriteMode;
	_startTime=time(NULL);

}
/* ************************************************************************************************************************************************************************************ */

cpedsMsgs::~cpedsMsgs() { 
	//if (logFile!=NULL) fclose(logFile);
	delete [] _tmpMessage;
}


/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::print(string msg, cpeds_messageImportance importance) {
	if (verbosity >= importance) {
		printf("%s",msg.c_str());
		sprintf(_tmpMessage,"%s",msg.c_str()); _singleMessage+=_tmpMessage;
	}
	
	if ( getLogging() ) {
		startLogging();
		if (logVerbosity >= importance) fprintf(logFile,"%s",msg.c_str());
		stopLogging();
	}
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::println(string msg, cpeds_messageImportance importance) {
	print(msg, importance);
	if (verbosity >= importance) { 
		printf("\n");
		sprintf(_tmpMessage,"\n"); _singleMessage+=_tmpMessage;
	}
	if ( getLogging() ) {
		startLogging();
		if (logVerbosity >= importance) fprintf(logFile,"\n");
		stopLogging();
	}
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::say(string from, string msg, cpeds_messageImportance importance, bool isError) {
	_singleMessage="";
	
	if (messaging) {
		if (importance == Zero && verbosity >= Top) {
			sayDateTime();
			printf("%s|%s>    - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s>    - %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
		}
		if (importance == Low && verbosity >= High) {
			sayDateTime();
			printf("%s|%s>   - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s>    - %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
		}
		if (importance == Medium && verbosity >= Medium) {
			sayDateTime();
			printf("%s|%s>  - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s>  -- %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
		}
		if ((importance == High && verbosity >= Low) or isError) { // we make visible all errors
			sayDateTime();
			printf("%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
		}
		if (importance == Top && verbosity >= Zero) {
			printf("\n***********************************************\n");
			sprintf(_tmpMessage,"\n***********************************************\n"); _singleMessage+=_tmpMessage;
			sayDateTime();
			printf("%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
			printf("***********************************************\n\n");
			sprintf(_tmpMessage,"***********************************************\n\n"); _singleMessage+=_tmpMessage;
		}
		if (importance == MustSee && verbosity >= Zero) {
			printf("\n***********************************************\n");
			printf("A MUST-SEE INFORMATION\n");
			sprintf(_tmpMessage,"\n***********************************************\n"); _singleMessage+=_tmpMessage;
			sayDateTime();
			printf("%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			sprintf(_tmpMessage,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str()); _singleMessage+=_tmpMessage;
			printf("***********************************************\n\n");
			sprintf(_tmpMessage,"***********************************************\n\n"); _singleMessage+=_tmpMessage;
		}
	}
	
	if ( getLogging() ) {
		startLogging();
		
		if (importance == Zero && logVerbosity >= Top) {
			logDateTime();
			fprintf(logFile,"%s|%s>   - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
		}
		if (importance == Low && logVerbosity >= High) {
			logDateTime();
			fprintf(logFile,"%s|%s>  - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
		}
		if (importance == Medium && logVerbosity >= Medium) {
			logDateTime();
			fprintf(logFile,"%s|%s> - %s\n",prefix.c_str() , from.c_str(), msg.c_str());
		}
		if (importance == High && logVerbosity >= Low) {
			logDateTime();
			fprintf(logFile,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
		}
		if (importance == Top && logVerbosity >= Zero) {
			fprintf(logFile,"\n***********************************************\n");
			logDateTime();
			fprintf(logFile,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			fprintf(logFile,"***********************************************\n\n");
		}
		if (importance == MustSee && logVerbosity >= Zero) {
			fprintf(logFile,"\n***********************************************\n");
			logDateTime();
			fprintf(logFile,"%s|%s> * %s\n",prefix.c_str() , from.c_str(), msg.c_str());
			fprintf(logFile,"***********************************************\n\n");
		}
		stopLogging();
	}
	
	return _singleMessage;
}

/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::say(string msg, cpeds_messageImportance importance) {
//#pragma omp critical
//	{
		say(msgSender,msg,importance);
//	}
	return _singleMessage;
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, const char* strch, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),strch);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, double d0, cpeds_messageImportance importance ) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),d0);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, long l0, cpeds_messageImportance importance ) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),l0);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, int val0, cpeds_messageImportance importance ) {
	return say(fmt,long(val0),importance);
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, double val0, double val1, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, double val0, double val1, double val2, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1,val2);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, double val0, double val1, double val2, double val3, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1,val2,val3);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, double val0, double val1, double val2, double val3, double val4, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1,val2,val3,val4);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, long val0, long val1, long val2, long val3, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1,val2,val3);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);		
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, long val0, long val1, long val2, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1,val2);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, long val0, long val1, cpeds_messageImportance importance) {
	bzero(_tmpMessage,_maxialMessageLength);
	sprintf(_tmpMessage,fmt.c_str(),val0,val1);
	_singleMessage=_tmpMessage;
	return say(_singleMessage,importance);	
	
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, int val0, int val1, cpeds_messageImportance importance) {
	return say(fmt,long(val0),long(val1),importance);
}
/***************************************************************************************/
string cpedsMsgs::say(string fmt, bool val0, cpeds_messageImportance importance) {
	_singleMessage=fmt;
	if (val0) {
		_singleMessage+="True";
	}
	else
		_singleMessage+="False";

	return say(_singleMessage,importance);	
}


/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::error(string msg, cpeds_messageImportance importance) {
	_singleMessage="";

//#pragma omp critical
	{
		sayError();
		say(msgSender,msg,importance,true);
	}
	return _singleMessage;
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::criticalError(string msg, cpeds_messageImportance importance) {
//#pragma omp critical
	{
		sayError();
		say(msgSender,msg,importance,true);
	exit(-1);
	}
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::warning(string msg, cpeds_messageImportance importance) {
	println("",importance);
	println("",importance);
	println("*** WARNING ***",importance);
	say(msg,importance);
	println("",importance);
	println("",importance);
	
	return _singleMessage;
}
/* ************************************************************************************************************************************************************************************ */
/* void cpedsMsgs::programInfo(string msg, cpeds_messageImportance importance) { */
/*   say(msgSender,msg,importance); */
/* } */



/* ************************************************************************************************************************************************************************************ */
/* void cpedsMsgs::setLogFileName(string lf) { logFileName=lf; if (lf=="") stopLogging(); else startLogging(); } */

// correction 2010/01/21 16:11:01 
void cpedsMsgs::setLogFileName(string lf) { logFileName=lf;  }

/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::getLogFileName() const { return logFileName; }

/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::setDateTimeOn(bool dton) {dateTimeOn=dton; }
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::setTimeSinceStartOn(bool tf) {timeSinceStartOn=tf; }
/* ************************************************************************************************************************************************************************************ */
bool cpedsMsgs::getDateTimgOn() const { return dateTimeOn; }
bool cpedsMsgs::getTimeSinceStartOn() const { return timeSinceStartOn; }

/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::setSender(string sender) { msgSender=sender; }

/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::getSender() const { return msgSender; }

/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::setVerbosity(cpeds_VerbosityLevel level) { verbosity=level; } //printf("********* changing verbosity to: %li, sender: %s\n",long(level),msgSender.c_str()); }

/* ************************************************************************************************************************************************************************************ */
cpeds_VerbosityLevel cpedsMsgs::getVerbosity() const { return verbosity; }

/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::setLogVerbosity(cpeds_VerbosityLevel level) { logVerbosity=level; }

/* ************************************************************************************************************************************************************************************ */
cpeds_VerbosityLevel cpedsMsgs::getLogVerbosity() const { return logVerbosity; }










/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStrf(double val, long acc) {
	string s,format="";
	filenamestr tmpch;
	format="%."+toStr(acc)+"lf";
	sprintf(tmpch,format.c_str(),val);
	s=tmpch;
	return s;
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(double val,long acc) {
	string s,format="";
	filenamestr tmpch;
	format="%."+toStr(acc)+"lE";
	/* sprintf(format,"%%.%lilE",acc); */
	sprintf(tmpch,format.c_str(),val);
	s=tmpch;
	return s;
}
/* ************************************************************************************************************************************************************************************ */
/* string cpedsMsgs::toStr(double val) { */
/*   return toStr(val,5); */
/* } */

/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(int val, string fmt) {
	return toStr(long(val),fmt);
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(long val, string fmt) {
	string s;
	filenamestr tmpch;
	sprintf(tmpch,fmt.c_str(),val);
	s=tmpch;
	return s;
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(long long val, string fmt) {
	string s;
	filenamestr tmpch;
	sprintf(tmpch,fmt.c_str(),val);
	s=tmpch;
	return s;
}

/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(cpeds_direction n, string cs) {
	string s;
	filenamestr tmpch;
	
	sprintf(tmpch," %.20lE,  %.20lE",n.l,n.b);
	s=tmpch;
	return s;
	
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(bool val) {
	string s;
	if (val) s="true"; else s="false";
	return s;
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::toStr(const char* c) {
	string s=c;
	return s;
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::sayDateTime() {
	if (dateTimeOn) {
		tm *t;
		time_t s=time(NULL);
		t=localtime(&s);
//#pragma omp critical
//		{
			printf("|%i/%02i/%02i %02i:%02i:%02i|",(int)(t->tm_year+1900), (int)(t->tm_mon)+1, (int)(t->tm_mday), (int)(t->tm_hour), (int)(t->tm_min), (int)(t->tm_sec));
			sprintf(_tmpMessage,"|%i/%02i/%02i %02i:%02i:%02i|",(int)(t->tm_year+1900), (int)(t->tm_mon)+1, (int)(t->tm_mday), (int)(t->tm_hour), (int)(t->tm_min), (int)(t->tm_sec)); 
			_singleMessage+=_tmpMessage;
//		}

	}
	if (timeSinceStartOn) {
//		time_t s=time(NULL);
//		printf("|%i|",(int)(s-_startTime));
//		sprintf(_tmpMessage,"|%i|",(int)(s-_startTime)); _singleMessage+=_tmpMessage;
//#pragma omp critical
//		{
			printf("|%i|",(int)timeElapsed());
			sprintf(_tmpMessage,"|%i|",(int)(timeElapsed())); 
			_singleMessage+=_tmpMessage;
//		}
	}
}
/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::getDateTime() {
	tm *t;
	time_t s=time(NULL);
	t=localtime(&s);
	sprintf(_tmpMessage,"|%i/%i/%i %02i:%02i:%02i|",(int)(t->tm_year+1900), (int)(t->tm_mon)+1, (int)(t->tm_mday), (int)(t->tm_hour), (int)(t->tm_min), (int)(t->tm_sec));
	return _tmpMessage;
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::logDateTime() {
	if (logging) {
		tm *t;
		time_t s=time(NULL);
		t=localtime(&s); 
		fprintf(logFile,"|%i/%i/%i %i:%i:%i|",(int)(t->tm_year+1900), (int)(t->tm_mon)+1, (int)(t->tm_mday), (int)(t->tm_hour), (int)(t->tm_min), (int)(t->tm_sec));    
	}
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::sayError() {
	perror("Error!: ");
//	printf("Error!: ");
	sprintf(_tmpMessage,"Error!: ");
	
	if (logging) fprintf(logFile,"Error!: ");
}

/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::clearLogFile() {
	bool wasLogging=getLogging();
	
	stopLogging();
	if (logFileName !="") {
		logFile = fopen(logFileName.c_str(),"w");
		fclose(logFile);
		logFile=NULL;
	}
	if (wasLogging)  startLogging();
	
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::startLogging() {
	if (logFileName !="") {
		if (logFile==NULL) { // start a new log file
			logFile = fopen(logFileName.c_str(),"a");
			logging=true;
			if (logFile == NULL) { stopLogging(); logging=false; } // if something went wrong then stop logging
		}
		else { 
			// stop logging to this file
			stopLogging();
			// start logging again
			logFile = fopen(logFileName.c_str(),"a");
			logging=true;
			if (logFile == NULL) { stopLogging(); logging=false;  } // if something went wrong then stop logging
		}
	}
	else { stopLogging(); return; }
}
/* ************************************************************************************************************************************************************************************ */
void cpedsMsgs::stopLogging() {
	if (logFile != NULL) { 
		fclose(logFile);
		logFile=NULL;
	}
	/* logging=false; */
}


/* ************************************************************************************************************************************************************************************ */
string cpedsMsgs::saveThisRun(int argc, char** argv, string fname) {
	string s="";
	stringstream ss;
	for (long i = 0; i < argc; i++) {
		ss << string(argv[i]) << " ";
	}
	s=ss.str();
	string tmp;
	if (fname.find('/',0)!=string::npos)
		tmp=fname;
	else
		tmp=getSender()+fname;
	printf("%s\n",s.c_str());
	if (tmp!="") {
		FILE *f=fopen(tmp.c_str(),&_saveRunWriteMode);
		if (f!=NULL) {
			fprintf(f,"%s> %s\n\n",getDateTime().c_str(),s.c_str());
			fclose(f);
		}
		else {
			say("Could not open run log file",Top);
		}
	}
	return s;
}
/***************************************************************************************/
long cpedsMsgs::timeElapsed() const {
	time_t s=time(NULL);
	return long(s-_startTime);
}
/***************************************************************************************/
void cpedsMsgs::restartTime() {
	_startTime=time(NULL);
}
