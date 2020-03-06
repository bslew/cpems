/*!
 * \file cpeds-point3d.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: blew
 */

#include "cpeds-point3d.h"

ostream& operator<<(ostream& s, cpedsPoint3D& p) {
	return s << "(" << p.x() << "," << p.y() << "," << p.z() << ")";
}
