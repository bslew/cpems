/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <string.h>
//#include <cpgplot.h>
//#include <fitsio.h>
//#include "map-scs-map.h"
//#include "cpeds-consts.h"
//#include "cpeds-math.h"
//#include "chealpix.h"



void map_scs::gui_set_nside()
{
    double input = lineEdit_nside->text().toDouble();
    lCDNumber1_nside_status.setValue(input);
}
