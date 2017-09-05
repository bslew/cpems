/*!
 \file pointsDensity.h - 
 */

#ifndef POINTSDENSITY_H_
#define POINTSDENSITY_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */
/* STANDALONE HEADERS */
#include "cpeds-point_set.h"
#include "subdomain.h"

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class pointsDensity
 \brief Encapsulates 
 \details 
 
 \date created: Feb 26, 2013, 10:04:19 PM 
 \author Bartosz Lew
 */
class pointsDensity : public cpedsPointSet3D {

	/* ---------------------------------------------------------------------------------------------------- */
	/* CLASS PUBLIC MEMBERS */
	/* ---------------------------------------------------------------------------------------------------- */
	public:

	/* ------------- */
	/* CLASS TYPES */
	/* ------------- */
	typedef struct {
		mscsVector<double> hsml;
		mscsVector<double> rho;
		subDomain *tree;
		double hsmlMax; //!< maximal hsml found when the density calculation was done
	}pointsDensity_t;

	/* ---------------------------- */
	/* CONSTRUCTORS AND DESTRUCTORS */
	/* ---------------------------- */
	pointsDensity();
	pointsDensity(const mscsVector<cpedsPoint3D>& points);
	void initiateVars();
	~pointsDensity();

	/* ---------------------------- */
	/* PUBLIC METHODS */
	/* ---------------------------- */
	/*!
		\brief calculate the points number density from the points distribution at the locations of those points.
		\details 
		@param NeighborsMin - minimal number of neighbors for smoothing length calculation
		@param NeighborsMax - maximal number of neighbors for smoothing length calculation
		@param is2dcase - flag that indicates whether the calculation should be done as for 2D points distribution. In this case the z coordinates
		of the points will be ignored
		@param smKernel - name of the smoothing kernel. Currently only "gadget2" name is the valid name.
		@param treeSheme - tree scheme to be used for domain decomposition and the definition of the domain that the points set occupies. 
		If not given then oct-tree and the full domain that the points set occupies is assumed. 
		@return
		
		The units of the resulting density is lu^(-dim) where
		lu is the length unit of the input points coordinates and dim is the dimension of the problem depending on the is2dcase flag. 
		So to obtain the density one needs to multiply the result by the mass particle, but it is assumed that all particles
		have the same masses. If this is not the case then in principle the method should be changed to account for that in the neighbors search so as to
		have a constant mass within smoothing radius rather then constant number of neighbors within that radius.
	
		\date Feb 26, 2013, 10:22:11 PM
		\author Bartosz Lew
	*/
	void calculateDensity(long NeighborsMin, long NeighborsMax, bool is2dcase=false, string smKernel="gadget2",  subDomain_region_t *treeScheme=NULL, double MassMin=0, double MassMax=0);
	
	double sml(long i) const { return _pointsDensity_data.hsml[i]; }
	double& sml(long i) { return _pointsDensity_data.hsml[i]; }
	double rho(long i) const { return _pointsDensity_data.rho[i]; }
	
	mscsVector<double>& density() { return _pointsDensity_data.rho; }
	mscsVector<double>& sml() { return _pointsDensity_data.hsml; }
	double getMaxHSML() { return _pointsDensity_data.hsmlMax; }
	subDomain* getTree() { return _pointsDensity_data.tree; }
	
	/* ---------------------------------------------------------------------------------------------------- */
	/* CLASS PROTECTED MEMBERS */
	/* ---------------------------------------------------------------------------------------------------- */
	protected:

	/* ---------------------------- */
	/* PROTECTED METHODS */
	/* ---------------------------- */

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
	pointsDensity_t _pointsDensity_data;

};
#endif /* POINTSDENSITY_H_ */ 

