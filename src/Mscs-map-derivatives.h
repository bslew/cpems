/*
 * \file Mscs-map-derivative.h - implements mscsMap class derivatives
 * UNDER CONSTRUCTION
 *
 */

#include "Mscs-map.h"

#ifndef __Mscs_map_derivatives__
#define __Mscs_map_derivatives__

/*!
\class map_derivatives
\brief Encapsulates the derivatives of a mscsMap class maps
\details
This class is dedicated to accurately calculate derivatives of the CMB shperical maps given in an
map_class object masked with a given mask - also given in the map_class object.

It is designed to take as small memory as possible. The class allows to get derivatives ring by ring
if the pixelization scheme of the map allows for that and it remembers its last state so that sequential
calls to rows subsequent row numbers were efficient

If the mask was supplied in the map_class object then the derivatives will not be calculated
in the masked pixels nor in the directly neighobouring pixles in either direction
The extended mask can also be returned

\date 2009/05/27 02:14:56
\author Bartosz Lew
*/
class mapDerivative : public mscsObject {
	public:
		mapDerivative();

		/*!
		\brief Constructor that initiate the object for subsequent calculations
		\details
		@param th - whether or not it is planned to derive the derivative in theta direction
		@param phi - whether or not it is planned to derive the derivative in phi direction
		@param th2 - whether or not it is planned to derive the second derivative in theta direction
		@param phi2 - whether or not it is planned to derive the second derivative in phi direction
		@param thphi - whether or not it is planned to derive the second derivative in theta and phi directions respectively
		@param map - map to be differenciated
		@param grid - true - the derivatives will be done on an equispaced grid, false - derivatives are calculated on
		@param Nth - number of grid pixels in latitude
		@param Nphi - number of grid pixels in longitude
		@return

		\date May 21, 2010, 10:36:01 AM
		\author Bartosz Lew
		*/
		mapDerivative(bool th, bool phi, bool th2, bool phi2, bool thphi, map_class* map, bool grid, long Nphi, long Nth);

		//! Clean up destructor
		~mapDerivative();

		/*!   This method defines how much there's interpolation between the pixels */
		/*!   @param N - defines the accuray level: 1 - single pixel accuracy */
		void set_accuracy(long N);


		/*!   Derives requested derivatives in the N'th row in the map */
		/*!   The result is returned in the matrix object where each row of the matrix corresponds to  */
		/*!   different type of the derivative in ordering as given to the constructor. */
		/*!   The number of columns in the matrix corresponds to the number of pixels in the N'th row in the map */
		/*!   Note that it is not checked the consistency between the ... */
		matrix<double>* get_allDrow(long N);

		//! This and the following routines returns a derivative of the initial map in an new object
		map_class* Dth();
		map_class* Dth2();
		map_class* Dphi();
		map_class* Dphi2();
		map_class* Dthphi();

		//! Returns the extended mask actually used in the derivative calculations
		map_class* get_extended_mask();


	private:
		string object_name;
		bool TH, PHI, TH2, PHI2, THPHI;
		map_class* m;

		long ACC;

		// map information variables
		long rownum;
		long rownumleo;
		bool masked;


		// run variables
		//
		// X,Y: the x and y coordinates of the pixels in row, note that Y is same for all X values hence it's a scalar
		// Z: the value in those pixels
		// Zint: the interpolated values onto the grid defined by the biggest row (or the bottom row in case they are of the same size)
		// M: the mask values in the row

		double *X1, Y1, *Z1, *Z1int, *M1, *M1int;
		double *X2, Y2, *Z2, *Z2int, *M2, *M2int;
		double *X3, Y3, *Z3, *Z3int, *M3, *M3int;

		// S number of pixels in the corresponding row in the above structures
		long S1, S2, S3;

		long* ringbeg;

		// state  variables
		long currow;
		long curpix;
		long curpixinrow;


		void resolve_pix_rings();
		void get_row(long fromPix, double*  X, long* pixinrow);


};
#endif
