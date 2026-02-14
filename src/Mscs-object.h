/*!
  \file defines an empty object frame
*/

#ifndef MSCS_OBJECT
#define MSCS_OBJECT
/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */
#include <string.h>
#include "cpeds-common.h"
#include "cpeds-msgs.h"


/* **************************************************************************************************** */
/* USING NAMESPACES */
/* **************************************************************************************************** */
using namespace std;


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsObject
  \brief Encapsulates basic object functionality
  \details 
  Functionality like, object name, messaging, and logging, managing verbosity
  levels.
  \date 2009/05/28 19:25:28 
  \author Bartosz Lew
*/
class mscsObject {

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:
/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  //! an empty constructor; the default name of the object is "object"
  mscsObject() { 
    object_name="mscsObject";  // not using the handler to avoid messages on initialization
    initiate_object_variables();
  }

  //! constructor that assigns a name to the object
  mscsObject(string name, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel) {
    object_name=name;  // not using the handler to avoid messages on initialization
    initiate_object_variables();
    msgs->setVerbosity(verbosityLoc);
  }

  mscsObject(const mscsObject& parent) {
//	  initiate_object_variables();
	  msgs=NULL;
	  *this=parent;
  }


  //! destructor
  virtual ~mscsObject() { if (msgs!=NULL) { delete msgs; msgs=NULL; }  }


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  //! returns the object's name
  string getName() const { return object_name; }
  //! sets a new name to the object
  void setName(string newName) { 
//    msgs->say("changing name to: "+newName,High);
    object_name = newName; 
    msgs->setSender(newName);
  }

  //! enable logging
  void enableLogging(string filename) {     
    msgs->say("enabling logging to file: "+msgs->getLogFileName(),Medium);
    msgs->setLogFileName(filename); 
  }
  //! disenable logging
  void disableLogging() { 
    msgs->say("enabling logging",Medium);    
    msgs->setLogFileName(""); 
  }

  cpeds_VerbosityLevel getVerbosityLevel() {  return msgs->getVerbosity(); }
  //! sets the verbosity level to verb
  void setVerbosityLevel(cpeds_VerbosityLevel verb) {     if (verb>=High) msgs->say("Setting verbosity level to "+msgs->toStr(long(verb)),Low);  msgs->setVerbosity(verb); }
  //! sets the logging verbosity level to verb
  void setLogVerbosityLevel(cpeds_VerbosityLevel verb) { if (verb>=High) msgs->say("Setting log verbosity level to "+msgs->toStr(long(verb)),Low); msgs->setLogVerbosity(verb); }


  /*!
    \brief manages the object's space
    \details 
    @param whattodo - defines the action to be performed: "make" (equivalent with "load") or "kill"
    @param for_what - defines a kind of space to be created: 
    "F" - space for alms
    "T" - space for temperature map 
    "m" - space for map mask
    "N" - space for number of pixel observations, 
    "C" - space for coordinates, 
    "rot" - space for rotation map
    @param how - defines subaction - how the action should be performed;
    used to matter only for case of "C" and "T" only but to distinguish between 
    the spherical temperatures (coordinates) and those on a flatmap (how = 1 and 2 respectively)
    but this is now depreciated since flat map is to be removed to a separate object
    
    \date 2009/06/04 20:17:19 
    \author Bartosz Lew
  */
  /* virtual void  makekill_space_manager(string whattodo, string for_what, int how) {} */
  /* virtual void  makekill_space_manager() {} */

/*   virtual void loadsave_manager(string whattodo, string fileformat, string what, int how, string where, long whatmultipole) {} */
/*   virtual void loadsave_manager() {} */


  const cpedsMsgs& getMessenger() const { if (msgs==NULL) { printf("!!!!!!!!!!!!!! WARNING RETURNING NULL MSGS. You should never see this message. \n"); exit(0); }return *msgs; }
  const cpedsMsgs* getMessengerP() const { return msgs; }
  cpedsMsgs& messenger() { return *msgs; }

  mscsObject& operator=(const mscsObject& rhs) {
	  if (this == &rhs) return *this;

	  
	  const cpedsMsgs* msg=rhs.getMessengerP();
	  if (msg!=NULL) {
/*
		  if (msgs!=NULL) { delete msgs; msgs=NULL; } // WHY WHEN THIS IS UNCOMMENTED THERE'S SOMETIMES SEGMENTATION FAULT ??
*/
		  /*
		 * Comment:
		 * 
		 * This is a constructor, so msgs is not yet initialized, so comparing its value
		 * against anything in the constructor doesn't make any sense.
		 * 
		 * author: blew
		 * date: Oct 16, 2020 4:08:16 PM
		 *
		 */
		
		  
//		  else {
//			  printf("!!!!!!!!!!!!!! WARNING THe msgs is NULL. This should not be !!! Will exit now. \n"); exit(0);			  
		  // msgs can actually be null here when entering  from the copy constructor
//		  }
	  /* msgs->setVerbosity(parent.verbosity); */
	  /* logVerbosity=parent.logVerbosity; */
		  msgs=new cpedsMsgs(rhs.getMessenger());	
	  }
	  else {
		  // this also should never happen but just in case of some bug...
		  printf("!!!!!!!!!!!!!! WARNING THe parent msgs is NULL. This should not be !!! Will exit now. \n"); exit(0);
	  }
	  //    object_name=rhs.getName(); // commented out on 2010-10-06 - we don't want functions to change names after eg. f=g
	  setName(rhs.getName()); 
	  
	  return *this;
  }
  
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:

/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  //! initiates local variables
  void initiate_object_variables() {
    msgs=new cpedsMsgs(getName());
  }


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
  string object_name;
  /* cpeds_VerbosityLevel verbosity; */
  /* cpeds_VerbosityLevel logVerbosity; */

  cpedsMsgs* msgs;


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


 };

#endif
