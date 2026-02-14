/*!
  \file qtbase.h - 
*/


#ifndef SRC_QTBASE_H_
#define SRC_QTBASE_H_

#include "QtCore/QtGlobal"
#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#warning QT_VERSION
#warning Qt>=5.0 is prefered. You are including Qt headers of version <5. 
#warning cpems can work with this version of Qt but when building other projects that link against cpems libraries
#warning qt5 dependence may be needed depending on what is available on the system. Consider defining QT_INC_DIR compiler flag.
#endif

#define QT_INC_DIR

#include QT_INC_DIR"QtCore/QString"
#include QT_INC_DIR"QtCore/QStringList"
#include QT_INC_DIR"QtCore/QtAlgorithms"
#include QT_INC_DIR"QtCore/QPointF"
#include QT_INC_DIR"QtCore/QDateTime"
#include QT_INC_DIR"QtCore/QMap"
#include QT_INC_DIR"QtCore/QFileInfo"

#endif /* SRC_QTBASE_H_ */ 



