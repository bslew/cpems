/*!
  \file implementation of various cosmological functions
*/

#ifndef CPEDS_FUNCTION_COSMO
#define CPEDS_FUNCTION_COSMO

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "Mscs-function.h"
#include "Mscs-map-window_function.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsFunctionCosmo
  \brief Implements various interesting cosmological functions
  \details 
  
  \date 2010/02/11 20:36:09 
  \author Bartosz Lew
*/
class cpedsFunctionCosmo : public mscsFunction {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


/* ------------- */
/* CLASS FRIENDS */
/* ------------- */
  

/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  cpedsFunctionCosmo() : mscsFunction() {}
  cpedsFunctionCosmo(const string name) : mscsFunction(name) {}
  ~cpedsFunctionCosmo() {}


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  /*!
    \brief generates age of the universe redshift relation for a given cosmology
    \details 
    @param N - number of points per redshift decade
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @return
    Age is given in Gyr
      
    \date 2010/02/11 20:43:21 
    \author Bartosz Lew
  */
  cpedsFunctionCosmo& make_age_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long N=0);

  
  /*!
    \brief generates scale factor time dependence  for a given cosmology
    \details 
    @param N - number of points per redshift decade
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @return
    Time is given in Gyr
      
    \date 2010/02/15 11:19:19 
    \author Bartosz Lew
  */
  cpedsFunctionCosmo& make_a_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0, long N);


  /*!
    \brief generates age of the universe redshift relation for a given cosmology
    \details 
    @param N - number of points per redshift decade
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @return
    Age is given in Gyr
      
    \date 2010/02/11 20:43:21 
    \author Bartosz Lew
  */
  /* cpedsFunctionCosmo& make_a_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double H0,long N=0); */
  

  /*!
    \brief generates the deceleration parameter - redshift relation for a given cosmology
    \details 
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @param N - number of points per redshift decade
    @return
    Age is given in Gyr
      
    \date 2010/02/11 20:43:21 
    \author Bartosz Lew
  */
  cpedsFunctionCosmo& make_q_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long N=0, cpedsFunctionCosmo* qvsz=NULL);

  /*!
    \brief generates the visibility - redshift relation for a given cosmology
    \details 
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @param N - number of points per redshift decade
    @return
    Age is given in Gyr
      
    \date 2010/02/11 20:43:21 
    \author Bartosz Lew
  */
  cpedsFunctionCosmo& make_tau_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long N=0, cpedsFunctionCosmo* visibz=NULL);


  cpedsFunctionCosmo& make_sigma0_vs_R(double Rst, double Ren, double dR, double Wr0, double Wm0, double Wl0, double w0, double h0,long N=0);
  
  /*!
	\brief calculate angular size distance vs redshift function
	\details 
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @param Nout - number of points to store the calculation (This is not the calculation precision but the number of output points)
    Use the cpeds function cpeds_set_points_number_per_logz to set the integration precision. -- THIS PARAMETER IS NOT USED AT THE MOMENT
    @param chi - comoving distance vs z; if X!=NULL then it is returned as a byproduct [Mpc]
	@return function angular size distance vs redshift

	\date Jan 19, 2011, 11:44:55 AM
	\author Bartosz Lew
  */
  cpedsFunctionCosmo& make_dA_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long Nout=200, cpedsFunctionCosmo* chi=NULL);
  
  
  /*!
	\brief calculate look-back time distance vs redshift relation
	\details 
    @param zSt - probe from redshift
    @param zEn - probe to redshift
    @param Wr0 - current density of radiation [in units of critical density]
    @param Wm0 - current density of mattter [in units of critical density]
    @param Wl0 - current density of DE [in units of critical density]
    @param w0 -  current value of EoS parameter 
    @param H0 - current value of the Hubble constant [km/s/Mpc]
    @param Nout - number of points to store the calculation (This is not the calculation precision but the number of output points)
    Use the cpeds function cpeds_set_points_number_per_logz to set the integration precision. -- THIS PARAMETER IS NOT USED AT THE MOMENT
    @param chi - comoving distance vs z; if X!=NULL then it is returned as a byproduct [Mpc]
	@return function look-back time distance vs redshift

	\date Jan 19, 2011, 11:44:55 AM
	\author Bartosz Lew
  */
  cpedsFunctionCosmo& make_lbtd_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long Nout=200);

  cpedsFunctionCosmo& make_Wr_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0);
  cpedsFunctionCosmo& make_Wm_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0);
  cpedsFunctionCosmo& make_Wl_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0);
  cpedsFunctionCosmo& make_Wtot_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0);
  cpedsFunctionCosmo& make_Wk_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0);

  /*!
	\brief calculate BBKS standard CDM model transfer function (linear prediction for today if today Omega values are given)
	\details 
	@param kmin - h/Mpc
	@param kmax - h/Mpc
	@param dlogk - log10 of the dk increment [log10(h/Mpc)]
	@param Wb 
	@return

	\date Aug 1, 2013, 2:55:46 PM
	\author Bartosz Lew
   */
  cpedsFunctionCosmo& make_BBKS_transferFunction(double kmin, double kmax, double dlogk, double Wb, double Wtot, double h);
  
  /*!
	\brief generates Tinker mass function at redshift z=0 for the requested over-density deltaRho
	\details 
	@param sigmaM - sigma(M) relation for given cosmology (unit of mass in [kg])
	@param Wm0 - Wm0
	@param deltaRho - over-density enclosing mass M
	@return returns tinker mass function multiplied by factor M^2/rho_m
	
	The mass function is given by:
	
	dn/dM = f(sigma) * rho_m/M * dln(sigma^-1)/dM
	
	where f(sigma) = A [ (sigma/b)^-a + 1 ] exp(-c/sigma^2)
	
	where A,a,b,c are constants interpolated for the requested overdensity.
	
	The overdensity should be within range [200, 3200].
	
	The mass function is calculated for redshift z=0.
	
	The created mass function is not the one that is returned. The created mass function is simply dn/dM
	whereas the returned mass function is M^2/rho_m * dn/dM

	Units of the created (not returned) mass function: [Mpc^-3 h^3 M_sol^-1 h^1 ]

	\date Aug 3, 2013, 12:19:43 AM
	\author Bartosz Lew
   */
  cpedsFunctionCosmo make_tinker_mass_function(mscsFunction& sigmaM, double Wm0, double deltaRho);

  /*!
	\brief generates Press-Schechter mass function at redshift 0 for the provided sigmaM function
	\details 
	@param sigmaM - sigma(M) relation for given cosmology (unit of mass in [kg])
	@param Wm0 - Wm0
	@param z - NOT IMPLEMENTED YET
	@return

	The mass function is given by:
	
	dn/dM = -sqrt(2/PI) * rho_m/M * delta_c/sigma(M) * 1/sigma(M) * dsigma(M)/dM
	
	where delta_c = 1.69
		
	The mass function is calculated/tested for redshift z=0 although extension to other redshifts should be easy.
	
	The created mass function is not the one that is returned. The created mass function is simply dn/dM
	whereas the returned mass function is M^2/rho_m * dn/dM (consistent with Tinker paper)
	
	Units: [Mpc^-3 h^3 M_sol h^-1 ]

	
	REVISION: Mar 4, 2014, 11:33:18 PM
	delta_c - should be model dependent 1.69 value is only for the EdS model


	\date Aug 3, 2013, 10:45:37 AM
	\author Bartosz Lew
   */
  cpedsFunctionCosmo make_PS_mass_function(mscsFunction& sigmaM, double Wm0, double z=0);

  
  /*!
	\brief calculate matter power spectrum of density contrast from transfer function and power-law curvature pert. spectrum
	\details 
	@param tf - matter transfer function (calibrated to unity at large scales) k units are assumed to be [h/Mpc]
	@param A - scalar pert. amplitude
	@param k0 - pivot scale [Mpc^-1]
	@param ns - spectral indes
	@param h - hubble (unitless) needed to convert from k/h to k 
	@return P(k_h) in units of [h^3/Mpc^3], k_h [h/Mpc]

	\date Aug 3, 2013, 2:51:18 PM
	\author Bartosz Lew
   */
  cpedsFunctionCosmo make_Pk(mscsFunction& tf, double A, double k0, double ns, double h);

  /*!
	\brief calculate mass function from the provided halo mass data
	\details 
	@param haloMass - function containing halo masses in units of 10^10 Msol_SI/h. Arguments are halo ID numbers.
	@param vol - volume containing the halo masses [Mpc^3/h^3]
	@param Npts - number of points in the mass function
	@param factors - mass function calibration modifiers; Currently not supported
	@param from - calculate mass function from this mass [10^10 Msol_SI/h]
	@param to - calculate mass function to this mass [10^10 Msol_SI/h]
	@return calculated mass function 

	The units of the output mass function will be [ Mpc^-3 h^3 Msol^-1 h^1 vs Msol/h ] 
	The output mass function is given in log10 - log10 space.
	

	\date Aug 22, 2013, 1:22:00 PM
	\author Bartosz Lew
   */
  cpedsFunctionCosmo calculateMassFunction(mscsFunction haloMass, double vol, long Npts=10, string factors="", double from=0, double to=0);

  cpedsFunctionCosmo& calculateBeamPowerPattern(double d,double ds, double lambda, double taper, double thetaMax, double dTheta);
  
  const cpedsFunctionCosmo& operator=(const mscsFunction& rhs);

  

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
	 typedef double (*fptrWx_t)(double , double , double , double , double);

/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */
	  cpedsFunctionCosmo& make_Wx_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, cpedsFunctionCosmo::fptrWx_t fWx);


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };
#endif /* CPEDS_FUNCTION_COSMO */ 

