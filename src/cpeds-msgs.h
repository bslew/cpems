#include <iostream>
/* #include <iomanip> */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cpeds-common.h"
#include "cpeds-consts.h"

using namespace std;
#ifndef _cpedsMsgs
#define _cpedsMsgs

/*!
  \class cpedsMsgs 
  \brief A %class for managing the output messages and logging.
  \details 
  
  This class enters many other classes and is a base for user communication and logging.
  However since verbosity match for outputting given message is checked every time something is said,
  this probably isn't the very optimal way of doing that. 

  One possibility for optimizations is to add macros depending on compilation flags -D to decide the certain 
  verbosity level and the compilation stage. This will prevent though flexible verbosity manipulations at the execution 
  time.

  One of the TODOs here is the implement sending errors to the stderr rather than stdout
  There are several other pending TODOs for improvements in this class, search the code for full description.

  \date 2009/05/26 20:55:55 
  \author Bartosz Lew
 */
class cpedsMsgs {
		
	public:
		
		//! This constructor by default will not log into file at all.
		cpedsMsgs(string sender="");
		//! This constructor will by default overwrite the existing log files
		cpedsMsgs(string sender, bool log, string logfilename, cpeds_VerbosityLevel verb);
		/*   //! This constructor will do with the log file as indicated  by the lfop char */
		/*   cpedsMsgs(string sender, bool log, string logfilename, string lfop, cpeds_VerbosityLevel verb); */
		//! copy constructor;
		cpedsMsgs(const cpedsMsgs& parent);
		
		~cpedsMsgs();
		
		//! prints the msg to the screen and possibly to a log file if loogging is enabled
		/*!  This does not go into the new line */
		void print(string msg, cpeds_messageImportance importance);
		
		//! prints the msg to the screen and possibly to a log file if loogging is enabled
		/*!  This DOES go to the new line */
		void println(string msg, cpeds_messageImportance importance);
		//! Says a message msg from sender from in a style for a requested importance
		string say(string from, const char* strch, cpeds_messageImportance importance);
		//! Says a message msg  in a style for a requested importance; the sender must me previously defined; otherwise there will be no sender by the message
		string say(string msg, cpeds_messageImportance importance);
		string say(string from, string str, cpeds_messageImportance importance, bool isError=false);
		string say(string fmt, double d0, cpeds_messageImportance importance);
		string say(string fmt, long l0, cpeds_messageImportance importance);
		string say(string fmt, int val0, cpeds_messageImportance importance);
		string say(string fmt, double val0, double val1, cpeds_messageImportance importance);
		string say(string fmt, double val0, double val1, double val2, cpeds_messageImportance importance);
		string say(string fmt, double val0, double val1, double val2, double val3, cpeds_messageImportance importance);
		string say(string fmt, double val0, double val1, double val2, double val3, double val4, cpeds_messageImportance importance);
		string say(string fmt, long val0, long val1, long val2, long val3, cpeds_messageImportance importance);
		string say(string fmt, long val0, long val1, long val2, cpeds_messageImportance importance);
		string say(string fmt, long val0, long val1, cpeds_messageImportance importance);
		string say(string fmt, int val0, int val1, cpeds_messageImportance importance);
		string say(string fmt, bool val0, cpeds_messageImportance importance);
		/* //! Says a message msg  in a style for a requested importance; the sender must me previously defined; otherwise there will be no sender by the message */
		/* void say(string msg, cpeds_messageImportance importance); */
		
		//! Says a message msg  in a style for a requested importance but doesn't stop the program; the sender must me previously defined; otherwise there will be no sender by the message
		string error(string msg, cpeds_messageImportance importance);
		//! Says a message msg  in a style for a requested importance and stops the program; the sender must me previously defined; otherwise there will be no sender by the message
		void criticalError(string msg, cpeds_messageImportance importance);
		//! Issues a warning message in a style for a requested importance 
		string warning(string msg, cpeds_messageImportance importance);
		
		
		//! Sets or returns the name of the log file
		//! \details if "" then no logging is done.
		void setLogFileName(string lf);
		string getLogFileName() const;
		
		//! switches on the messaging to the screen
		void messagingOn() { messaging=true; }
		//! switches off the messaging to the screen
		void messagingOff() { messaging=false; }
		//! switches on the logging to file; does not change the status of the log file pointer
		void loggingOn() { logging=true; }
		//! switches off the logging flag; does not change the status of the log file pointer etc. -eg. the log file is not automatically closed
		void loggingOff() { logging=false; }
		bool getLogging()  const { return logging; }
		bool getMessaging()  const { return messaging; }
		
		bool dateTimeStatus()  const { return dateTimeOn; }
		string pref()  const { return prefix; }
		
		//! Sets or returns the sender of the msg.
		void setSender(string sender);
		string getSender() const;
		
		//! Sets or returns the date/time printout
		void setDateTimeOn(bool dton);
		bool getDateTimgOn() const;
		
		bool getTimeSinceStartOn() const;
		/*!
			\brief returns the number of seconds from the moment of object creation
			\details 
			@return number of seconds from the moment of object creation
		
			\date Nov 16, 2012, 1:44:50 PM
			\author Bartosz Lew
		*/
		long timeElapsed() const;
		/*!
			\brief resets the timer to current time. 
			\details 

			This will impact the returned timeElapsed value.
		
			\date Dec 20, 2013, 4:38:05 PM
			\author Bartosz Lew
		*/
		void restartTime();
		void setTimeSinceStartOn(bool tf);

		/*!
			\brief sets the default mode for saving the command line parameters used in the last run of a program
			\details 
			@param param m - can be wither 'a' for appending to the log file or, 'w' (default) for overwriting the run log file
		
			\date Jun 26, 2012, 10:25:49 PM
			\author Bartosz Lew
		*/
		void setSaveRunWriteMode(char m) { _saveRunWriteMode=m; }
		
		//! Sets or returns the verbosity level: the levels are defined in cpeds-common.h
		void setVerbosity(cpeds_VerbosityLevel level);
		cpeds_VerbosityLevel getVerbosity() const;
		
		//! Sets or returns the log verbosity level: the levels are defined in cpeds-common.h
		void setLogVerbosity(cpeds_VerbosityLevel level);
		cpeds_VerbosityLevel getLogVerbosity() const;
		
		/*!
			\brief saves to file how the program was called with all arguments.
			\details 
			@param argc - command line parameters count
			@param argv - array of char* with command line parameters
			@param fname - string with the preferred file name where the last cmd run line is to be saved. If not given then the file will not be saved.
			@return cmdline returned as string
		
			\date Jun 21, 2012, 2:08:08 PM
			\author Bartosz Lew
		*/
		string saveThisRun(int argc, char** argv, string fname=".last_run.log");
		
		
		//! Returns the string containing the value val in the %lf format
		string toStrf(double val,long acc=5);
		//! Returns the string containing the value val with acc digig precision (not implemented yet)
		string toStr(double val,long acc=5);
		//! Returns the string containing the value val
		/* string toStr(double val); */
		//! Returns the string containing the value val
		string toStr(long val, string fmt="%li");
		//! Returns the string containing the value val
		string toStr(long long val, string fmt="%lli");
		//! Returns the string containing the value val
		string toStr(int val, string fmt="%li");
		//! Returns the string containing the value val being a cpeds_direction as defined in cpeds-common.h
		string toStr(cpeds_direction n, string cs);
		//! Returns the string containing true or false as defined by val
		string toStr(bool val);
		//! Returns the string containing converted from array of chars
		string toStr(const char* c);
		
		string getDateTime();
		
		
		//! 
		void startLogging();
		void clearLogFile();
		void stopLogging();


		
		// define useful operators here
		
		
	protected:
		//! TODO: make this stuff structured for easier assignments.
		FILE* logFile;
		void sayDateTime();
		void sayError();
		void logDateTime();
		
		bool dateTimeOn, timeSinceStartOn;
		bool logging, messaging;
		
		
		string logFileName;
		string msgSender;
		string prefix;
		
		char _saveRunWriteMode;
		
		cpeds_VerbosityLevel verbosity;
		cpeds_VerbosityLevel logVerbosity;
		
	private:
		static const int _maxialMessageLength=1000;
		char *_tmpMessage;
		string _singleMessage;
		time_t _startTime;
		
};

#endif
