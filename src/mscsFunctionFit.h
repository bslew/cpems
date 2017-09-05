/*!
 \file mscsFunctionFit.h - 
 */

#ifndef MSCSFUNCTIONFIT_H_
#define MSCSFUNCTIONFIT_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "Mscs-function.h"

/* FORWARD DECLARATIONS */

int mscsFunctionFit_exponential_model_f(const gsl_vector * x, void *data, gsl_vector * f);
int mscsFunctionFit_exponential_model_df (const gsl_vector * x, void *data, gsl_matrix * J);
int mscsFunctionFit_exponential_model_fdf (const gsl_vector * x, void *data,	gsl_vector * f, gsl_matrix * J);


/*!
	\brief 
	\details 

	The function is:
	gauss * (1.0-gsl_sf_erf(nsig/sq2))
	where
		gauss=A * exp(-(x-m)*(x-m)/(2*s*s));
	and
		nsig=(m-x)*alpha/s;
	

	THe fitted parameters are A,m,s,alpha
	param0 - A
	param1 - m
	param2 - s
	param3 - alpha
	

	\date Nov 14, 2013, 1:34:26 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_skewgauss_4param_model_f(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_skewgauss_4param_model_df(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_skewgauss_4param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);

/*!
	\brief 
	\details 
	@param
	@return

	The function is:
		f=A * exp(-(x-m)*(x-m)/(2*s*s)) + B;

	THe fitted parameters are A,m,s,B
	param0 - A
	param1 - m
	param2 - s
	param3 - B

	\date Jun 2, 2015, 5:35:16 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_gauss_4param_model_f(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_gauss_4param_model_df(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_gauss_4param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);

/*!
	\brief 
	\details 
	@param
	@return

	The function is:
		f=A1 * exp(-(x-m1)*(x-m1)/(2*s1*s1)) + A2 * exp(-(x-m2)*(x-m2)/(2*s2*s2)) + B;

	THe fitted parameters are A1,m1,s1,A2,m2,s2,B
	param0 - A1
	param1 - m1
	param2 - s1
	param3 - A2
	param4 - m2
	param5 - s2
	param6 - B

	\date Jun 2, 2015, 5:35:16 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_2gauss_7param_model_f(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_2gauss_7param_model_df(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);
int mscsFunctionFit_2gauss_7param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);


int mscsFunctionFit_line_2param_model_f(const gsl_vector * x, void *data, gsl_vector * f);
int mscsFunctionFit_line_2param_model_df (const gsl_vector * x, void *data, gsl_matrix * J);
int mscsFunctionFit_line_2param_model_fdf (const gsl_vector * x, void *data,	gsl_vector * f, gsl_matrix * J);

int mscsFunctionFit_line_1paramB_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_line_1paramB_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_line_1paramB_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);

/*!
	\brief 
	\details 
	y = A*(x/c0)^B
	This model requires setting 1 constant mscsFunctionFit_model_const0
	THe fitted parameters are A and B

	\date Nov 14, 2013, 1:34:26 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_powerLaw_2param_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_powerLaw_2param_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_powerLaw_2param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);

/*!
	\brief 
	\details 
	y = A*(x/c0)^c1
	This model requires setting 2 constants mscsFunctionFit_model_const0 and mscsFunctionFit_model_const1
	THe fitted parameters are A

	\date Nov 14, 2013, 1:34:26 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_powerLaw_1paramA_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_powerLaw_1paramA_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_powerLaw_1paramA_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);


/*!
	\brief beta model of the FNW galaxy density profile
	\details 

	n(r) = n0 ( 1+ (r/rc)^2 )^(-3/2*beta)
	
	param1 - n0
	param2 - rc
	param3 - beta
	
	\date Nov 28, 2013, 2:59:40 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_beta_3param_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_beta_3param_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_beta_3param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);


/*!
	\brief 2D gaussian function with possibly non-zero covatiance
	\details 

	f(x,y) = A exp(- 1/(2(1-rho^2)[ (x-x0)^2/(sx^2) + (y-y0)^2/(sy^2) - 2rho (x-x0)(y-y0)/(sx * sy) ] )
	
	There are no constant values in this model.
	
	param1 - A
	param2 - x0
	param3 - sx
	param4 - y0
	param5 - sy
	param6 - rho
	
	\date Nov 28, 2013, 2:59:40 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_gauss2D_6param_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_gauss2D_6param_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_gauss2D_6param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);

/*!
	\brief power-law 3 parameter
	\details 
	y = A*( (x-c0) / c1 )^n + B
	This model requires setting 2 constants mscsFunctionFit_model_const0, and mscsFunctionFit_model_const1
	THe fitted parameters are A, B and n

	\date Nov 14, 2013, 1:34:26 PM
	\author Bartosz Lew
*/
int mscsFunctionFit_powerLaw_3param_model_f(const gsl_vector * param, void *data, gsl_vector * f);
int mscsFunctionFit_powerLaw_3param_model_df (const gsl_vector * param, void *data, gsl_matrix * J);
int mscsFunctionFit_powerLaw_3param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J);


extern double mscsFunctionFit_model_const0;
extern double mscsFunctionFit_model_const1;
extern double mscsFunctionFit_model_const2;
extern double mscsFunctionFit_model_const3;
/* USING NAMESPACES (only for inside-header files implementations) */



/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class mscsFunctionFit
 \brief Encapsulates 
 \details 
 
 \date created: Apr 10, 2013, 11:57:07 PM 
 \author Bartosz Lew
 */
class mscsFunctionFit : public mscsFunction {

	/* ---------------------------------------------------------------------------------------------------- */
	/* CLASS PUBLIC MEMBERS */
	/* ---------------------------------------------------------------------------------------------------- */
	public:

	/* ------------- */
	/* CLASS TYPES */
	/* ------------- */
	typedef struct {
		cpedsList<double> fitedParam;
		cpedsList<double> fitedParamErr;
		matrix<double> cov;
		mscsFunction bestFitModel;
		mscsFunction resid;
		double chisq; // chi^2/dof
		double residStDev;
		double DOF;
		double Pvalue;
		string modelName;
		long parametersCount;
		long constsCount;
		gsl_multifit_function_fdf fit;
		int status;
		
	}mscsFunctionFit_t;

	
	typedef struct  {
			size_t n;
			double * x;
			double * y;
			double * sigma;
	} fitData_t;

	
	/* ---------------------------- */
	/* CONSTRUCTORS AND DESTRUCTORS */
	/* ---------------------------- */
	mscsFunctionFit(string model="");
	~mscsFunctionFit();

	/* ---------------------------- */
	/* PUBLIC METHODS */
	/* ---------------------------- */
	/*!
		\brief fits loaded data with a specified model
		\details 
		@param model - name of the model. See below for the valid model names
		@param initGuessVals - initial values to start minimization from (the size should correspond to the number of parameters in the model)
		@param errors - measurement errors. vector of size of the number of function points.
		@return
	
		Valid models:\n
		"exponential_3param" - A*exp(-lambda*x) + b (params ordering: A,lambda,b) \n
		"skewgauss_4param" - A*exp(-(x-m)^2/s^2) * (1.0-erf((m-x)*alpha/(sqrt(2)*s))) (params ordering: A,m,s,alpha) \n
		"line_2param" - A*x + B (params ordering: A,B) \n
	
		\date Apr 11, 2013, 2:32:28 AM
		\author Bartosz Lew
	*/
	cpedsList<double> fitData(cpedsList<double>& initGuessVals, cpedsList<double>& errors, string model="" );
	
	/*!
		\brief generates the best fit model based on the found parameters
		\details 
		@param acc - best fit model resolution parameter. If 1, the resolution of the returned model will be same as the input data. For 0.1 it will be 10 times denser, and so on.
		@return
	
		\date Apr 11, 2013, 12:27:07 PM
		\author Bartosz Lew
	*/
	mscsFunction getBestFitModel(double acc=1,mscsFunction* args=NULL);
	mscsFunction getBestFitModel(double from, double to, double acc);
	
	
	void setModel(string model) { select_model(model); }
	string getModel() const { return _mscsFunctionFit_data.modelName; }
	int getFitStatus() const { return _mscsFunctionFit_data.status; }
	long getModelParametersCount() { return _mscsFunctionFit_data.parametersCount; }
	long getModelConstsCount() { return _mscsFunctionFit_data.constsCount; }
	void setModelConst(int constParam, double value);
	double getModelConst(int constParam);
	double& chisq() { return _mscsFunctionFit_data.chisq; }
	double chisqPerDOF() { return chisq()/DOF(); }
	double& pValue();
	double& DOF() { return _mscsFunctionFit_data.DOF; }
	
	cpedsList<double>& params() { return _mscsFunctionFit_data.fitedParam; }
	cpedsList<double>& paramsErr() { return _mscsFunctionFit_data.fitedParamErr; }
	/*!
		\brief convenience function to make fitting in a single call.
		\details 
		@param data - function with the data points
		@param yerr - list with measurement errors. if a single value is passed then the same value will be used for all data points
		@param model - name of the model to fit
		@param Pini - list with initial parameter guess values
		@param bfModelRes - best fit model resolution - defines the resolution of the output function tabulated for the best fit model
		@param Srad - search radius for selecting multiple starting points for independent runs for improve the best fit. Should 
		be a list of search radii - one for each parameter
		@param Nstates - number of random initial point draws for Levenberg-Marquardt fitting.
		@return
	
		\date Jun 24, 2013, 8:20:52 PM
		\author Bartosz Lew
	*/
	mscsFunction fit(mscsFunction& data, cpedsList<double> yerr, string model, cpedsList<double>& Pini, double bfModelRes, cpedsList<double>* Srad=NULL, long Nstates=1000);
	void print_fit_result(cpedsList<double>& p,cpedsList<double>& pErr);
	
//	void setFittedParams(cpedsList<double>& fitedParam, cpedsList<double>& fitedParamErr) { 
//		_mscsFunctionFit_data.fitedParam=fitedParam; 
//		_mscsFunctionFit_data.fitedParamErr=fitedParamErr; 
//	}
	
	/* ---------------------------------------------------------------------------------------------------- */
	/* CLASS PROTECTED MEMBERS */
	/* ---------------------------------------------------------------------------------------------------- */
	protected:

	/* ---------------------------- */
	/* PROTECTED METHODS */
	/* ---------------------------- */
	void print_state (size_t iter, gsl_multifit_fdfsolver * s);
	
	void select_model(string model);
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
	mscsFunctionFit_t _mscsFunctionFit_data;
//	long i;
//	long NP;
//	double P[10];
//	double *y;
//	double *sigma;
//	long Ndata;
//	double Yi;
	

	
};

void mscsFunctionFit_radialGridSearch(mscsFunctionFit& fit, cpedsList<double> yerr, cpedsList<double> Pini, cpedsList<double> radii, long Nstates);

#endif /* MSCSFUNCTIONFIT_H_ */ 

