/*!
 * \file mscsLine.h
 *
 *  Created on: Jan 16, 2012
 *      Author: blew
 */

#ifndef MSCSLINE_H_
#define MSCSLINE_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-object.h"
#include "cpeds-point3d.h"
#include "cpeds-point_set.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates straight line in 3d defined in a parametric way by two points and parametrized by z-value
 \details 
 This class provides a contained for the analytical form of line. This class can basically be also treated as a ray of light
 or a vector since there is an implicit convention that the direction of the vector associated with this line is
 given by p2-p1 where p1 and p2 are the two points given to define the line via the eg. constructor.
 This matters for instance when calculating angles between lines.
 
 \date Jan 16, 2012, 10:58:34 AM 
 \author Bartosz Lew
 */

class mscsLine : public mscsObject {
	public:
		mscsLine(string name="mscsLine");
		mscsLine(cpedsPoint3D p1, cpedsPoint3D p2, string name="mscsLine");
		mscsLine(double x1,double y1,double z1,double x2,double y2,double z2, double zmin=0, double zmax=0, string name="mscsLine");		
//		mscsLine(double a,double b,double c,double x1,double y1,double z1, double zmin=0, double zmax=0);
		mscsLine(const mscsLine& parent);
		virtual ~mscsLine();
		void setInitialValues();
		void initialize_abc();

		
		/*!
			\brief returns the function point on the line defined by the z parameter in the z-parametrization
			\details 
			@param z - parameter
			@return point
			
			The parametrization of line in this case is:
			X= a (z-z1)/c + x1
			Y= a (z-z1)/c + y1
			Z= z
		
			\date Jan 16, 2012, 1:10:57 PM
			\author Bartosz Lew
		*/
		cpedsPoint3D getValue_zparam(double z);
		
		double x_zparam(double z) const;
		double y_zparam(double z) const;

		/*!
			\brief returns the function point on the line defined by the z parameter in the z-parametrization
			\details 
			@param z - parameter
			@return point
			
			The parametrization of line in this case is:
			X= a t + x1
			Y= b t + y1
			Z= c t + z1
		
			\date Jan 16, 2012, 1:10:57 PM
			\author Bartosz Lew
		*/
		cpedsPoint3D getValue(double t);
		
		cpedsPoint3D getP1() const { return cpedsPoint3D(_params.x1,_params.y1,_params.z1); }
		cpedsPoint3D getP2() const { return cpedsPoint3D(_params.x2,_params.y2,_params.z2); }
//		const cpedsPoint3D getP1() const { return cpedsPoint3D(_params.x1,_params.y1,_params.z1); }
//		const cpedsPoint3D getP2() const { return cpedsPoint3D(_params.x2,_params.y2,_params.z2); }
		
		double getP1x() const { return _params.x1; }
		double getP1y() const { return _params.y1; }
		double getP1z() const { return _params.z1; }
		double getP2x() const { return _params.x2; }
		double getP2y() const { return _params.y2; }
		double getP2z() const { return _params.z2; }

		double& x1() { return _params.x1; }
		double& y1() { return _params.y1; }
		double& z1() { return _params.z1; }
		double& x2() { return _params.x2; }
		double& y2() { return _params.y2; }
		double& z2() { return _params.z2; }
		double& a() { return _params.a; }
		double& b() { return _params.b; }
		double& c() { return _params.c; }
		
		/*!
			\brief calculates line between this and provided line after translating them to the origin of the coordinate system
			\details 
			@param l - reference line
			@return angle in radians
		
			\date Jan 17, 2012, 11:49:52 AM
			\author Bartosz Lew
		*/
		double angle(mscsLine& l) const;
		/*!
			\brief as angle but takes smaller angle between two lines
			\details 
			@return angle in rad.
			This returns value always within range [0, PI/2]
		
			\date Jan 18, 2012, 10:47:28 AM
			\author Bartosz Lew
		*/
		double angleNoVect(mscsLine& l) const;
		double dot(mscsLine& l) const;
		mscsLine cross(mscsLine& l) const;
		
		
		/*!
			\brief rotates line treated as vector about x axis.
			\details 
			@param ax - rotation angle in radians
			@return returns this as rotated object
			
			The convention is that for positive angles rotation is from Y axis towards Z axis.
		
			\date Jan 17, 2012, 5:22:25 PM
			\author Bartosz Lew
		*/
		mscsLine& Rx(double ax);
		/*!
			\brief rotates line treated as vector about Y axis.
			\details 
			@param ay - rotation angle in radians
			@return returns this as rotated object
			
			The convention is that for positive angles rotation is from Z axis towards X axis.
		
			\date Jan 17, 2012, 5:22:25 PM
			\author Bartosz Lew
		*/
		mscsLine& Ry(double ay);
		/*!
			\brief rotates line treated as vector about Z axis.
			\details 
			@param ax - rotation angle in radians
			@return returns this as rotated object
			
			The convention is that for positive angles rotation is from X axis towards Y axis.
		
			\date Jan 17, 2012, 5:22:25 PM
			\author Bartosz Lew
		*/
		mscsLine& Rz(double az);
		mscsLine& translate(const cpedsPoint3D& v);
		
		void printLine(string info="") const;
		
		cpedsPointSet3D tabulateZ(double zmin, double zmax, double dz);
		/*!
			\brief tabulates the line with N points between the points defining the line
			\details 
			@param N - number of points in the resulting point set.
			@return
		
			\date Jan 17, 2012, 5:49:50 PM
			\author Bartosz Lew
		*/
		cpedsPointSet3D tabulateN(long N);
		mscsLine& operator=(const mscsLine& rhs);
		
		
	protected:
		typedef struct {
			double x1,x2,y1,y2,z1,z2;
			double a,b,c;
		} _line_params_t;
		
		_line_params_t _params;

		
		
		
	private:

};

#endif /* MSCSLINE_H_ */






