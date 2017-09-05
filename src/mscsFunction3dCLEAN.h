/*!
 \file mscsFunction3dCLEAN.h - 
 */

#ifndef MSCSFUNCTION3DCLEAN_H_
#define MSCSFUNCTION3DCLEAN_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function3dregc.h"
#include "cpeds-direction_set.h"
/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class mscsFunction3dCLEAN
 \brief Encapsulates 
 \details 
 
 \date created: Oct 28, 2014, 4:37:48 PM 
 \author Bartosz Lew
 */
class mscsFunction3dCLEAN : public mscsFunction3dregc {
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
			mscsFunction3dregc data,psf,psfclean; // dirty map,input psf, cleaned psf- all of the same sizes
			mscsFunction3dregc b,br; // size-optimized psf and a container for rotated optimized psf
			mscsFunction3dregc bc; // size-optimized clean beam (psf)
			mscsFunction3dregc src; // clean map
			mscsFunction3dregc resid; // residual noise
			mscsFunction3dregc srcresid; // clean map with residual noise
			mscsFunction3dregc srcbcn; // cleaned map convolved with clean beam and with noise added
			long MaxIter;
			double beamOptimizationThres; // threshold used for beam optimization
			double loopGain;
			string algorithm;
			bool fitBeam;
			cpedsDirectionSet radec; // cleaned positions, value stores cleaned amplitudes
			cpedsDirectionSet srcij; // cleaned source positions, value stores cleaned amplitudes
//			cpedsList<double> alpha; // beam orientations for each cleaned source
		} mscsFunction3dCLEAN_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		mscsFunction3dCLEAN();
		~mscsFunction3dCLEAN();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		/*!
			\brief initialize the object data 
			\details 
			@param data - dirty map; the imaginary part should contain the position angles of the DD beam for that DD intensity pointing
			@param psf - point spread function. Should be centered in the field.
			@param MaxIter - maximal number of iterations 
			@param loopGain - loop gain value of the Hogbom algorithm
			@return
		
			\date Oct 28, 2014, 5:07:03 PM
			\author Bartosz Lew
		*/
		void initiate(mscsFunction3dregc& data, mscsFunction3dregc& psf, long MaxIter, double loopGain, string algorithm="CLEAN");
		void setBeamOptimizationThres(double thres) { _mscsFunction3dCLEAN_data.beamOptimizationThres=thres; }
		void setFitBeamOrientationAmplitude(bool tf) { _mscsFunction3dCLEAN_data.fitBeam=tf; }
		void clean();
//		void clean2();
		
		mscsFunction3dregc calculateCleanBeam(mscsFunction3dregc& psf);
		void rotatePSF(double angle);
		
		// cleaned map
		mscsFunction3dregc& cleanMap() { return _mscsFunction3dCLEAN_data.src; }
		// full-size clean beam
		mscsFunction3dregc& cleanBeam() { return _mscsFunction3dCLEAN_data.psfclean; }
		// size optimized clean beam
		mscsFunction3dregc& cleanBeamOptimized() { return _mscsFunction3dCLEAN_data.bc; }
		mscsFunction3dregc& cleanMapBeamNoise() { return _mscsFunction3dCLEAN_data.srcbcn; }
		mscsFunction3dregc& cleanMapNoise() { return _mscsFunction3dCLEAN_data.srcresid; }
		mscsFunction3dregc& residual() { return _mscsFunction3dCLEAN_data.resid; }
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		
		void fitBeamOrientation();
		void fitBeamAmplitude();

		
		/*!
			\brief shrinks the size of the provided psf to the minimal possible size
			\details 
			@param
			@return
		
			\date Oct 28, 2014, 4:57:23 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc optimizeDDPSFsize(mscsFunction3dregc psf, double thres);
		
		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */

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
		mscsFunction3dCLEAN_t _mscsFunction3dCLEAN_data;
		static const long Nsteps=1000;

		double vMax,vMaxB,vMinB;
		double rx0,ry0,x,y; 
		double g,dI, angle,rms;
		long iMax,jMax,kMax,id;
		long iMaxB,jMaxB,kMaxB;
		long iMinB,jMinB,kMinB;
		long iBoff,jBoff; // indexes offsets of the center of the rotated beam - pointing the place where the source should be - from the positive beam location indexes

		
};
#endif /* MSCSFUNCTION3DCLEAN_H_ */ 

