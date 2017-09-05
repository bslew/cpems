#include "cpeds-math.h"
/*! \mainpage 
  <H1> CPEDS is a cosmological parametres evolution derivation software for the standard FLRW cosmology.</H1>

  \section Purpose
  The package consists of number of files and classes dedicated to different purposes. These are:

  \li cpeds-cosmo - number of useful cosmological parameters in Friedman-Lemaitre-Robertson-Walker cosmology 
  \li cpeds-pixelization - implements the Healpix sphere pixelization scheme
  according to formulas described in 
  <a href="http://adsabs.harvard.edu/abs/2004astro.ph..9533R>astro-ph/0409533">http://adsabs.harvard.edu/abs/2004astro.ph..9533R>astro-ph/0409533</a>
  \li cpeds-math - contains several basic numerical routines for simple and useful %math 
  \li cpeds-msgs - class for providing a user communication interface and logging functionality
  \li cpeds-consts - number of useful mathematical and physical and astronomical constants definitions
  \li cpeds-pdf - class for managing probability distribution functions
  \li cpeds-templates - re-implements some templates also present in standard libraries 
  (eg. a queue template for standard C types)
  
  \section Compilation and Installation
  To compile you need a number of external libraries available to link against. These include:
  \li cpgplog
  \li gsl

  \section Conventions
  There are some conventions that are generally followed in the code. 

  \li \b angles are generally kept in radians unless otherwise noted
  \li <b> coordinate systems</b> - The generic convention is that CPEDS
  operates in the spherical polar coordinate system \f$(\theta,\phi )\f$ with 
  \f$ \theta\in <0,\pi> \f$ and \f$ \theta \in <0,2\pi) \f$, with
  \f$ \theta=0\f$ on the north pole (a convention consistent with the Healpix
  package.
  

  \section License
  LGPL
  
  \section Author
 
  Bartosz Lew <blew at umk.pl>
  
  \date Last modification: 2009/05/27 01:20:12 
*/

#ifndef CPEDS_COSMO_DEFS
#define CPEDS_COSMO_DEFS

/*!
  \brief Returns the frequency  at which CMB energy density spectrum peaks  [Hz]
  \details 
  @param CMB temperature for which to derive the frequency
*/
double cpeds_cmbfMax(double T);

//! Returns the \f$\Omega_{tot} = \Omega_{r0} +\Omega_{m0} + \Omega_{\Lambda0}\f$
double cpeds_Wtot0(double Wr0, double Wm0, double Wl0);

//! set the number of points to probe in integrals over the logarithimc redshift space
/*! @param np - number of points that probe the log10 interval of the redshift space */
void cpeds_set_points_number_per_logz(long np);
//! This function gives step in redshift, for different integrals the be calculated; it is uniform in logarithmic space
/*! The step is controlled by the CPEDS_NP global variable */
/*! @param z - redshift for which to evaluate the apropriate step */
double cpeds_deltaz(double z);
//! returns the E(z) facror defined as: \f$E(z)=\sqrt{ \Omega_{r0}(1+z)^{-4} +\Omega_{m0}(1+z)^{-3} + \Omega_{\Lambda0}(1+z)^{3*(1+w)}+ (1-\Omega_{tot0})*(1+z)^2}\f$
double cpeds_Efactor(double Wr0, double Wm0, double Wl0, double z, double w);

//! Particle horizon at given redshift. in units [Gpc h^-1]
/*! Derives \f$d_p(z) =  \int_\infty^z \frac{dz}{E(z,\Omega_{r0},\Omega_{m0},\Omega_{\Lambda0},w,z)}\f$ */
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w);

/*!
	\brief Particle horizon at given redshift. in units [Gpc h^-1]
	\details 
	@param Z - table of redshift at which the dp is tabulated
	@param dp - table or redshifts at which dp is calculated
	@param N - size of the Z and dp arrays
	@return
	The Z and dp arrays should not be allocated. They will be allocated inside of this function and their size will be returned.

	\date Jan 19, 2011, 12:26:37 PM
	\author Bartosz Lew
*/
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w, double** Z, double** dp, long* N);

//! Derives particle horizon at given redshift. in units [Gpc]
/*! Derives \f$d_p(z) = \frac{1}{h} \int_\infty^z \frac{dz}{E(z,\Omega_{r0},\Omega_{m0},\Omega_{\Lambda0},w,z)}\f$ */
double cpeds_dp(double Wr0, double Wm0, double Wl0, double z, double w, double h);
//! Derive the age of the universe in [Gyr h^-1]. 
/*! The parameters have their usual meaning  */
/*! The age of the universe is calculated as: \f$t(\Omega_{r0},\Omega_{m0},\Omega_{l0},z,w,h) = \int_{\inf}^0\frac{dz}{(1+z)E(\Omega_{r0},\Omega_{m0},\Omega_{l0},z)}\f$ */
double cpeds_age_of_universe(double Wr0, double Wm0,double Wl0,double z, double w);
//! Derive the age of the universe in hubble time [Gyr]. 
/*! The parameters have their usual meaning  */
/*! The age of the universe is calculated as: \f$t(\Omega_{r0},\Omega_{m0},\Omega_{l0},z,w,h) =\frac{1}{h} \int_{\inf}^0\frac{dz}{(1+z)E(\Omega_{r0},\Omega_{m0},\Omega_{l0},z)}\f$ */
double cpeds_age_of_universe(double Wr0, double Wm0,double Wl0,double z, double w, double h);
//! Derives the lookback time distance (time travel distance) to object of a redshift z in [Gpc h^-1]: ie.:
/*! \f$ d_{lbtd}=c ( age(\Omega_{r0},\Omega_{m0},\Omega_{l0},w,z=0)-age(\Omega_{r0},\Omega_{m0},\Omega_{l0},w,z) ) \f$ */
double cpeds_lbtd(double Wr0, double Wm0,double Wl0,double z, double w);
//! lookback time distance to object of a redshift z in [Gpc ]
double cpeds_lbtd(double Wr0, double Wm0,double Wl0,double z, double w, double h);
//! comoving distance to object of a redshift z in [Gpc h^-1]
/*! \f$ \xi = \xi(z_0)-\xi(z) \f$ */
double cpeds_comoving_distance(double Wr0, double Wm0,double Wl0,double z, double w);
//! comoving distance (conformal distance) to object of a redshift z in [Gpc ]
/*! \f$ \xi = \xi(z_0)-\xi(z) \f$ */
double cpeds_comoving_distance(double Wr0, double Wm0,double Wl0,double z, double w, double h);
//! comoving distance to SLS controlled by f parameter in [Gpc h^-1]
/*! for f<1 distance is larger than the lbtd_SLS */
/*! for f>1 distance is smaller than the lbtd_SLS */
double cpeds_comoving_distance_rec(double Wr0, double Wm0,double Wl0, double w, double f);
//! lookback time distance to SLS controlled by f parameter in [Gpc]
/*! for f<1 distance is larger than the lbtd_SLS */
/*! for f>1 distance is smaller than the lbtd_SLS */
double cpeds_comoving_distance_rec(double Wr0, double Wm0,double Wl0, double w, double h, double f);

//! Derives the position of the accoustic peak in the CMB power spectrum \f$C_l\f$:
/*! \f$ l_p = d_A/d_c \f$ */
/*!  where  */
/*! \f$d_A\f$ - angular diameter distance to SLS */
/*! \f$ d_c \f$ - physical sound horizon size at SLS ~ 1/Sqrt(3) d_p(SLS)/(1+z) */
/*! This depends on the sound speed during recombination: here it is assumed to be 1/sqrt(3) of the light speed in vacuum */
double cpeds_acoustic_horizon_multipole(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h);
//! the accoustic angular scale in the C_l: \f$a_{acc} = \frac{180}{PI} \frac{1}{l_p}\f$ in [deg]
double cpeds_acoustic_horizon_angle(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h);

//! Derives the curvature radius of the Universe in [Gpc*h^-1] 
/*! \f$ R_c = \frac{c}{H_0}\frac{1}{\sqrt{|\Omega_{tot0}|}} \f$ */
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0);
//! Derives the curvature radius of the Universe in [Gpc] 
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0, double h);
//! Derives the curvature radius of the Universe  at redshift z [Gpc]
double cpeds_curvature_radius(double Wr0, double Wm0,double Wl0, double z, double w, double h);

//**************************************************************************************
//! Derives a comoving angular diameter distance to a given redshift [Gpc]
/*! \f$ d_A = \frac{f(X)}{(1+z)  }\f$ where:
\f$ f(X) = \Biggl\{\begin{array}{ccc}
            R_c \sin( X/R_c ) & for& W_k < 0 \\
            X &for & W_k=0 \\
            R_c \sinh( -X/R_c ) & for & W_k < 0
\end{array}\f$

 where X is the comoving distance  and 
\f$ R_c\f$ - is the curvature raduis 
*/
double cpeds_angular_diameter_distance(double Wr0, double Wm0,double Wl0, double z, double w, double h);

double cpeds_get_dA_from_ang(double Wr0, double Wb0, double Wm0, double  Wl0, double* Z, double w0, double h, double size);

//! Derives the comoving luminosity distance dL = (1+z)^2 dA to an object at redshift z [Gpc]
/*! \f$ D_L(z) = (1+z)^2 d_A \f$ where \f$ d_A \f$ is the angular size distance to that object. */
double cpeds_luminosity_distance(double Wr0, double Wm0,double Wl0, double z, double w, double h);

//! Derives the \f$\Omega_l\f$ value at given redshift for given cosmology
double cpeds_Wl(double Wr0, double Wm0, double Wl0, double z, double w);
//! Derives the \f$\Omega_b\f$ value at given redshift for given cosmology
double cpeds_Wb(double Wr0, double Wb0, double Wl0, double z, double w);
//! Derives the \f$\Omega_m\f$ value at given redshift for given cosmology
double cpeds_Wm(double Wr0, double Wm0, double Wl0, double z, double w);
//! Derives the \f$\Omega_r\f$ value at given redshift for given cosmology
double cpeds_Wr(double Wr0, double Wm0, double Wl0, double z, double w);
double cpeds_Wtot(double Wr0, double Wm0, double Wl0, double z, double w);
//! Derives the \f$\Omega_{k}\f$ value  for given cosmology at given redshift
double cpeds_Wk(double Wr0, double Wm0,double Wl0, double z, double w);
//! Derives the \f$\Omega_{k0}\f$ value  for given cosmology, defined as \f$\Omega_{k0} = 1-\Omega_{tot}\f$
double cpeds_Wk0(double Wr0, double Wm0,double Wl0);

//! Derives the sound horizon at given redshift for given cosmology (if it makes sense) [Gpc ] 
//\f$ \rho_C = d_p * 
double cpeds_sound_horizon(double Wr0, double Wb0, double Wm0,double Wl0, double z, double w, double h);

//double cpeds_sound_horizon(double Wr0, double Wm0,double Wl0, double z, double w, double h);




//! sound speed in units of c: 
/*! 
  \f$ c_s = \frac{1}{\sqrt{3(1+\frac{3}{4}\frac{ \rho_b}{\rho_r})}} \f$
  \warning these routines are not yet tested
*/
double cpeds_sound_speed(double rhor, double rhob);
//! Derives the critical density of the universe: \f$ \rho_C\f$
/*! 
  \f$ \rho_C = \frac{3H^2}{8\pi G} \f$ 
  return Result is in [kg/m^3 * h^2] 
*/
double cpeds_rhoC(double Wr0, double Wm0,double Wl0, double z, double w);

//! critical density of the universe: rhoC = 3H^2/8PIG [kg/m^3]
double cpeds_rhoC(double Wr0, double Wm0,double Wl0, double z, double w, double h);

//! current critical density of the Universe [kg/m^3]
double cpeds_rhoC0(double h=1);

//! Derives the Hubble constant at redshift z: H = h*100 * E(z,Wr0,Wm0,Wl0,w) [km/s/Mpc]
double cpeds_hubble(double Wr0, double Wm0,double Wl0, double z, double w, double h);


double fX_e(double T, double eta);
//**************************************************************************************
// funkcja zwraca energie wiazania atomu wdoropodobnego w elektronowoltach w zaleznosci od 
// podanej temperatury w kelwinach, ktorej odpowiada generalnie jakistam stan wzgudzenia 
// atomow w osrodku.  np dla temperatury powyzej 160 000 K mymy jonizacje z n=1 a dla 
// T = 50 0000 K mamny jonizacje z n=2
double cpeds_BB(double T);
//**************************************************************************************
// funkcja zwraca m - stan podstawowy dla danej temperatury. Jesli temperarura jest wieksza
// od energii jonizacji to funkcja zwraca neff. czyli zaklada sie w programie ze atom w stanie
// neff jest juz zjonizowany, chyba dlatego tez dalem sumowanie do neff-1 w funkcji r_w
double cpeds_m_base (double T);
//**************************************************************************************
// funkcja zwraca wartosc efektywnej liczny stanow wzbudzonych (liczbe mozliwych serii)
// ktora zalezy od temperatury osrodka - a dokladniej powinna zalezec jak gestosc^-1/3
double cpeds_neff(double T);

/*!
	\brief calculates mean molecular weight per particle mass for the stantard cosmology with only hydrogen and helium
	\details 
	@param Xp - hydrogen mass fraction (typically 0.76=1-Yp)
	@return

	This procedure could be improved to include cases of heavy elements both in neutral and in fully ionized medium.

	\date Apr 26, 2013, 1:48:26 PM
	\author Bartosz Lew
*/
double cpeds_calculate_meanMolecularWeight(double Xp);

/*!
	\brief convert gadget internal energy to temperature in K
	\details 
	@param u - particle internal energy in gadget default units (defined in cpeds-consts)
	@param gamma - adiabatic index (typically 5/3)
	@param Xp - hydrogen mass fraction (typically 0.76=1-Yp) used to calculate the mean molecular mass
	@return

	IMPORTANT:
		For the moment, this procedure calculates the temperature using mean molecular weight for the ionized IGM only.
		The result will be 2x too small in the limit of small temperatures where mu = ~1.22 (rather than 0.58) for the mixture of 
		neutral hydrogen and helium


	\date Apr 26, 2013, 1:50:27 PM
	\author Bartosz Lew
*/
double cpeds_convert_gadget_internalEnergy_to_temp(double u, double gamma, double Xp);

/*!
	\brief calculates the y0 constant for selected cosmology
	\details 
	@param Wb0 - baryon density today
	@param h - hubble parameter in units of 100 km/s/Mpc
	@param mueinv - number of electrons per proton (default 1.136)
	@return compton_y0 parameter [K^-1 Mpc^-1]

	\date Apr 29, 2013, 5:31:55 PM
	\author Bartosz Lew
*/
double cpeds_calculate_comptonY(double Wb0, double h, double mueinv=1.136);

#endif

