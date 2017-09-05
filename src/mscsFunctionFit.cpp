/*!
 * \file mscsFunctionFit.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: blew
 */

#include "mscsFunctionFit.h"
#include "Mscs-function.h"
#include <gsl/gsl_sf_erf.h>


double mscsFunctionFit_model_const0;
double mscsFunctionFit_model_const1;
double mscsFunctionFit_model_const2;
double mscsFunctionFit_model_const3;

/***************************************************************************************/
void mscsFunctionFit_radialGridSearch(mscsFunctionFit& fit, cpedsList<double> yerr, cpedsList<double> Pini, cpedsList<double> radii, long Nstates) {
	long chisqMin_idx,chisqMax_idx;
	double chisqMin,chisqMax;
	cpedsList<double> chisq;
	vector< cpedsList<double> > params;
	vector< cpedsList<double> > paramsErr;
	cpedsList<double> fitParams;
	long Nparam=fit.getModelParametersCount();

	double tmp=0.1*fit.stdev();
	cpedsList<double> Srad=radii;
//	for (unsigned long i = 0; i < Pini.size(); i++) {
//		Srad.append(Pini[i]*radius);
//	}
		
	vector< cpedsList<double> > params_ini;
	cpedsRNG rng("uniform");
	cpedsList<double> rns;
	rns.makeLength(Nstates);

	for (long j = 0; j < Nparam; j++) {
		if (Srad[j]>0) {
			rng.setMinMax(Pini[j]-Srad[j],Pini[j]+Srad[j]);
			params_ini.push_back(rng.getRNs(Nstates));				
		}
		else {
			rns=Pini[j];
			params_ini.push_back(rns);
		}
	}
		
	for (long i = 0; i < Nstates; i++) {
		for (long j = 0; j < Nparam; j++) {
			Pini[j]=params_ini[j][i];
		}
		
		fitParams=fit.fitData(Pini,yerr);
		chisq.append(fit.chisq());
		params.push_back(fitParams);
		paramsErr.push_back(fit.paramsErr());
//		if ((i % 1000) == 0) msgs.say("done: %lf %%",double(i)/_Nstates*100, Low);
	}

	chisq.getMinMaxValues(&chisqMin,&chisqMax,&chisqMin_idx,&chisqMax_idx);
	fit.params()=params[chisqMin_idx];
	fit.paramsErr()=paramsErr[chisqMin_idx];
	fit.chisq()=chisq[chisqMin_idx];
}


/***************************************************************************************/
mscsFunctionFit::mscsFunctionFit(string model) {
	_mscsFunctionFit_data.parametersCount=0;
	_mscsFunctionFit_data.chisq=-1;
	_mscsFunctionFit_data.modelName="";
	_mscsFunctionFit_data.parametersCount=-1;
	select_model(model);
	
}

mscsFunctionFit::~mscsFunctionFit() {
	// TODO Auto-generated destructor stub
}



/***************************************************************************************/
int mscsFunctionFit_exponential_model_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *sigma = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	double b = gsl_vector_get (x, 2);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = A * exp(-lambda * i) + b */
		double t = i;
		double Yi = A * exp (-lambda * t) + b;
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_exponential_model_df (const gsl_vector * x, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *sigma = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double t = i;
		double s = sigma[i];
		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, e/s); 
		gsl_matrix_set (J, i, 1, -t * A * e/s);
		gsl_matrix_set (J, i, 2, 1/s);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_exponential_model_fdf(const gsl_vector * x, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_exponential_model_f(x, data, f);
	mscsFunctionFit_exponential_model_df(x, data, J);
	return GSL_SUCCESS;
}
/***************************************************************************************/

int mscsFunctionFit_skewgauss_4param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double m = gsl_vector_get (param, 1);
	double s = gsl_vector_get (param, 2);
	double alpha = gsl_vector_get (param, 3);
	
	mscsFunction fn("skewgauss",Zero);
//	printf("generating skewgauss from: %lE to %lE, dx:%lE, A:%lE m: %lE s:%lE alpha: %lE\n",x[0],x[n-1],x[1]-x[0],A,m,s,alpha)	;
	fn.mkSkewGauss(x[0],x[n-1],x[1]-x[0],A,m,s,alpha,(x[1]-x[0])/10);
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set (f, i, (fn.getY(i) - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_skewgauss_4param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double m = gsl_vector_get (param, 1);
	double s = gsl_vector_get (param, 2);
	double alpha = gsl_vector_get (param, 3);
	double sq2=1.41421356237310;
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double arg = x[i]-m;
		double g = A*exp(-arg*arg/(2*s*s));
		double nsig=-arg*alpha/s;
		double erfc = (1.0-gsl_sf_erf(nsig/sq2));
		
		double g1 = A*A*exp(-arg*arg*(1.0+alpha*alpha)/(2*s*s));
		double g2 = exp(arg*arg*alpha*alpha/(2*s*s));
		
		
		gsl_matrix_set (J, i, 0, 2.0*g*sqrt(2*PI)*s*erfc/yerr[i]); 
		gsl_matrix_set (J, i, 1, -g1*(2.0*alpha*s+g2*sqrt(2*PI)*(-arg)*erfc)/(yerr[i]*s));
		gsl_matrix_set (J, i, 2, (g1*(2.0*alpha*s*(-arg)+g2*sqrt(2*PI)*(arg*arg+s*s)*erfc))/(yerr[i]*s*s));
		gsl_matrix_set (J, i, 3, 2.0*g1*arg/yerr[i]);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_skewgauss_4param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_skewgauss_4param_model_f(param, data, f);
	mscsFunctionFit_skewgauss_4param_model_df(param, data, J);
	return GSL_SUCCESS;
}
/***************************************************************************************/

int mscsFunctionFit_gauss_4param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double m = gsl_vector_get (param, 1);
	double s = gsl_vector_get (param, 2);
	double B = gsl_vector_get (param, 3);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A*exp(-((x[i]-m)*(x[i]-m))/(2.0*s*s)) + B - y[i])/yerr[i]);
	}
	
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_gauss_4param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double m = gsl_vector_get (param, 1);
	double s = gsl_vector_get (param, 2);
	double B = gsl_vector_get (param, 3);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-(x-m)^2/(2 s^2))  */
		/* and the xj are the parameters (A,m,s,B) */
		double arg = x[i]-m;
		double g = A*exp(-arg*arg/(2*s*s)); //+B;
		
		gsl_matrix_set (J, i, 0, g/A/yerr[i]); 
		gsl_matrix_set (J, i, 1, g*arg/(yerr[i]*s*s));
		gsl_matrix_set (J, i, 2, g*arg*arg/(yerr[i]*s*s*s));
		gsl_matrix_set (J, i, 3, 1.0/yerr[i]);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_gauss_4param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_gauss_4param_model_f(param, data, f);
	mscsFunctionFit_gauss_4param_model_df(param, data, J);
	return GSL_SUCCESS;
}
/***************************************************************************************/
/***************************************************************************************/

int mscsFunctionFit_2gauss_7param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A1 = gsl_vector_get (param, 0);
	double m1 = gsl_vector_get (param, 1);
	double s1 = gsl_vector_get (param, 2);
	double A2 = gsl_vector_get (param, 3);
	double m2 = gsl_vector_get (param, 4);
	double s2 = gsl_vector_get (param, 5);
	double B  = gsl_vector_get (param, 6);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A1*exp(-((x[i]-m1)*(x[i]-m1))/(2.0*s1*s1)) + A2*exp(-((x[i]-m2)*(x[i]-m2))/(2.0*s2*s2)) + B - y[i])/yerr[i]);
	}
	
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_2gauss_7param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A1 = gsl_vector_get (param, 0);
	double m1 = gsl_vector_get (param, 1);
	double s1 = gsl_vector_get (param, 2);
	double A2 = gsl_vector_get (param, 3);
	double m2 = gsl_vector_get (param, 4);
	double s2 = gsl_vector_get (param, 5);
	double B  = gsl_vector_get (param, 6);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A1 * exp(-(x-m1)^2/(2 s1^2)) + A2 * exp(-(x-m2)^2/(2 s2^2))  */
		/* and the xj are the parameters (A1,m1,s1,A2,m2,s2,B) */
		double arg1 = x[i]-m1;
		double arg2 = x[i]-m2;
		double g1 = A1*exp(-arg1*arg1/(2*s1*s1));
		double g2 = A2*exp(-arg2*arg2/(2*s2*s2));
		
		gsl_matrix_set (J, i, 0, g1/A1/yerr[i]); 
		gsl_matrix_set (J, i, 1, g1*arg1/(yerr[i]*s1*s1));
		gsl_matrix_set (J, i, 2, g1*arg1*arg1/(yerr[i]*s1*s1*s1));
		gsl_matrix_set (J, i, 3, g2/A2/yerr[i]); 
		gsl_matrix_set (J, i, 4, g2*arg2/(yerr[i]*s2*s2));
		gsl_matrix_set (J, i, 5, g2*arg2*arg2/(yerr[i]*s2*s2*s2));
		gsl_matrix_set (J, i, 6, 1.0/yerr[i]);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_2gauss_7param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_gauss_4param_model_f(param, data, f);
	mscsFunctionFit_gauss_4param_model_df(param, data, J);
	return GSL_SUCCESS;
}


/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/


int mscsFunctionFit_line_1paramB_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double B = gsl_vector_get (param, 0);
	

	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (mscsFunctionFit_model_const0*(x[i]-mscsFunctionFit_model_const1)+B - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_line_1paramB_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
//	double B = gsl_vector_get (param, 0);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = const0 * (x-const1) + B */
		
		gsl_matrix_set (J, i, 0, 1.0/yerr[i]);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_line_1paramB_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_line_1paramB_model_f(param, data, f);
	mscsFunctionFit_line_1paramB_model_df(param, data, J);
	return GSL_SUCCESS;
}

/***************************************************************************************/
int mscsFunctionFit_line_2param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A*(x[i]-mscsFunctionFit_model_const0)+B - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
}

int mscsFunctionFit_line_2param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * (x-const) + B */
		
		gsl_matrix_set (J, i, 0, (x[i]-mscsFunctionFit_model_const0)/yerr[i]); 
		gsl_matrix_set (J, i, 1, 1.0/yerr[i]);
	}
	return GSL_SUCCESS;
}

int mscsFunctionFit_line_2param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_line_2param_model_f(param, data, f);
	mscsFunctionFit_line_2param_model_df(param, data, J);
	return GSL_SUCCESS;
}


/***************************************************************************************/

int mscsFunctionFit_powerLaw_2param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A*pow(x[i]/mscsFunctionFit_model_const0,B) - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
}
int mscsFunctionFit_powerLaw_2param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * (x/const0)^B */
		double arg = pow(x[i]/mscsFunctionFit_model_const0,B);
		double arg2 = A*pow(x[i]/mscsFunctionFit_model_const0,B)*log(x[i]/mscsFunctionFit_model_const0);
		
		gsl_matrix_set (J, i, 0, arg/yerr[i]); 
		gsl_matrix_set (J, i, 1, arg2/yerr[i]);
	}
	return GSL_SUCCESS;

}
int mscsFunctionFit_powerLaw_2param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_powerLaw_2param_model_f(param, data, f);
	mscsFunctionFit_powerLaw_2param_model_df(param, data, J);
	return GSL_SUCCESS;	
}
/***************************************************************************************/
int mscsFunctionFit_beta_3param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double n0 = gsl_vector_get (param, 0);
	double rc = gsl_vector_get (param, 1);
	double beta = gsl_vector_get (param, 2);
	double alpha=-1.5*beta;
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (n0*pow(1.0 + (x[i]/rc)*(x[i]/rc),alpha) - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
}
int mscsFunctionFit_beta_3param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double n0 = gsl_vector_get (param, 0);
	double rc = gsl_vector_get (param, 1);
	double beta = gsl_vector_get (param, 2);
	double alpha=-1.5*beta;
	double rrc2;
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * (x/const0)^B */
		rrc2=(x[i]/rc)*(x[i]/rc);
		double dfdn0 = pow(1.0 + rrc2,alpha);
		double dfdrc = 3.0*x[i]*x[i]*n0*pow(1.0 + rrc2,alpha-1);
		double dfdbeta = 3.0*n0*pow(1.0 + rrc2,alpha) * log10(1.0+rrc2);
		
		gsl_matrix_set (J, i, 0, dfdn0/yerr[i]); 
		gsl_matrix_set (J, i, 1, dfdrc/(rc*rc*rc*yerr[i]));
		gsl_matrix_set (J, i, 2, dfdbeta/(2.0*yerr[i]));
	}
	return GSL_SUCCESS;

}
int mscsFunctionFit_beta_3param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_beta_3param_model_f(param, data, f);
	mscsFunctionFit_beta_3param_model_df(param, data, J);
	return GSL_SUCCESS;	
}

/***************************************************************************************/

int mscsFunctionFit_powerLaw_1paramA_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);

	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A*pow(x[i]/mscsFunctionFit_model_const0,mscsFunctionFit_model_const1) - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
	
}
int mscsFunctionFit_powerLaw_1paramA_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
//	double B = gsl_vector_get (param, 0);
	
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * (x/const0)^const1 */
		double arg = pow(x[i]/mscsFunctionFit_model_const0,mscsFunctionFit_model_const1);
		
		gsl_matrix_set (J, i, 0, arg/yerr[i]);
	}
	return GSL_SUCCESS;
	
}
int mscsFunctionFit_powerLaw_1paramA_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_powerLaw_1paramA_model_f(param, data, f);
	mscsFunctionFit_powerLaw_1paramA_model_df(param, data, J);
	return GSL_SUCCESS;		
}
/***************************************************************************************/
int mscsFunctionFit_powerLaw_3param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	double nn = gsl_vector_get (param, 2);

	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (A*pow( (x[i]-mscsFunctionFit_model_const0)/mscsFunctionFit_model_const1, nn) + B - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
		
}
int mscsFunctionFit_powerLaw_3param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double A = gsl_vector_get (param, 0);
	double B = gsl_vector_get (param, 1);
	double nn = gsl_vector_get (param, 2);
	double tmp1,tmp2;
	double dfdA, dfdB, dfdn;
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * ((x-const0)/const1)^nn + B */
//		rrc2=(x[i]/rc)*(x[i]/rc);
		tmp1=(x[i]-mscsFunctionFit_model_const0)/mscsFunctionFit_model_const1;
		tmp2=pow( tmp1,nn);
		dfdA = tmp2/yerr[i];
		dfdB = 1/yerr[i];
		if (tmp1==0) dfdn=0;
		else dfdn= A*tmp2*log(tmp1)/yerr[i];
		
		gsl_matrix_set (J, i, 0, dfdA); 
		gsl_matrix_set (J, i, 1, dfdB);
		gsl_matrix_set (J, i, 2, dfdn);
	}

	return GSL_SUCCESS;
}
int mscsFunctionFit_powerLaw_3param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_powerLaw_3param_model_f(param, data, f);
	mscsFunctionFit_powerLaw_3param_model_df(param, data, J);
	return GSL_SUCCESS;			
}

/***************************************************************************************/
int mscsFunctionFit_poly_4param_model_f(const gsl_vector * param, void *data, gsl_vector * f) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double a0 = gsl_vector_get (param, 0);
	double a1 = gsl_vector_get (param, 1);
	double a2 = gsl_vector_get (param, 2);
	double a3 = gsl_vector_get (param, 3);

	size_t i;
	
	for (i = 0; i < n; i++) {
		gsl_vector_set(f, i, (
				a0*mscsFunctionFit_model_const0 + 
				a1*mscsFunctionFit_model_const1*x[i] + 
				a2*mscsFunctionFit_model_const2*x[i]*x[i] + 
				a3*mscsFunctionFit_model_const3*x[i]*x[i]*x[i] 
				 - y[i])/yerr[i]);
	}
	
	return GSL_SUCCESS;
		
}
int mscsFunctionFit_poly_4param_model_df (const gsl_vector * param, void *data, gsl_matrix * J) {
	size_t n = ((mscsFunctionFit::fitData_t*)data)->n;
	double *x = ((mscsFunctionFit::fitData_t*)data)->x;
	double *y = ((mscsFunctionFit::fitData_t*)data)->y;
	double *yerr = ((mscsFunctionFit::fitData_t*) data)->sigma;
	
	double a0 = gsl_vector_get (param, 0);
	double a1 = gsl_vector_get (param, 1);
	double a2 = gsl_vector_get (param, 2);
	double a3 = gsl_vector_get (param, 3);
//	double tmp1,tmp2;
	double dfda0, dfda1, dfda2, dfda3;
	size_t i;
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = a0c0 + a1c1 x + a2c2 x^2 +a3c3 x ^3  */
//		rrc2=(x[i]/rc)*(x[i]/rc);
//		tmp1=(x[i]-mscsFunctionFit_model_const0)/mscsFunctionFit_model_const1;
//		tmp2=pow( tmp1,nn);
		dfda0 = mscsFunctionFit_model_const1/yerr[i];
		dfda1 = mscsFunctionFit_model_const1*x[i]/yerr[i];
		dfda2 = mscsFunctionFit_model_const2*x[i]*x[i]/yerr[i];
		dfda3 = mscsFunctionFit_model_const3*x[i]*x[i]*x[i]/yerr[i];
		
		gsl_matrix_set (J, i, 0, dfda0); 
		gsl_matrix_set (J, i, 1, dfda1);
		gsl_matrix_set (J, i, 2, dfda2);
		gsl_matrix_set (J, i, 3, dfda3);
	}

	return GSL_SUCCESS;
}
int mscsFunctionFit_poly_4param_model_fdf(const gsl_vector * param, void *data,	gsl_vector * f, gsl_matrix * J) {
	mscsFunctionFit_poly_4param_model_f(param, data, f);
	mscsFunctionFit_poly_4param_model_df(param, data, J);
	return GSL_SUCCESS;			
}


/***************************************************************************************/
/******************************************************************************************/
/***************************************************************************************/
/******************************************************************************************/
void mscsFunctionFit::print_state (size_t iter, gsl_multifit_fdfsolver * s) {
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2), 
			gsl_blas_dnrm2 (s->f));
	//	for (long i = 0; i < s->f->size; i++) {
	//		printf("x:%li y:%lE\n",i,(s->f)[i]);
	//	}
}



/***************************************************************************************/
cpedsList<double> mscsFunctionFit::fitData(cpedsList<double>& initGuessVals, cpedsList<double>& errors, string model ) {
	_mscsFunctionFit_data.fitedParam.clear();
	_mscsFunctionFit_data.fitedParamErr.clear();
	
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	unsigned int i, iter = 0;
	const size_t n = pointsCount();
	const size_t paramCount = initGuessVals.size();

	gsl_matrix *covar = gsl_matrix_alloc (paramCount, paramCount);
	gsl_matrix *J;  /* BLcomment (Aug 28, 2017, 7:49:04 PM): 
		A new matrix to hold covariance matrix since the new GSL gsl_multifit_fdfsolver structure
		does not have this member anymore.
	*/
	double *x=extractArguments();
	double *y=extractValues();
	double *sigma=errors.toCarray();
	
	fitData_t d = { n, x, y, sigma};
	
	double *param_ini = initGuessVals.toCarray();
	gsl_vector_view paramVec = gsl_vector_view_array (param_ini, paramCount);
	
	
	if (model!="") select_model(model);
	_mscsFunctionFit_data.fit.n = n;
	_mscsFunctionFit_data.fit.p = paramCount;
	_mscsFunctionFit_data.fit.params = &d;
	
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, paramCount);
	gsl_multifit_fdfsolver_set (s, &_mscsFunctionFit_data.fit, &paramVec.vector);
	
//	print_state (iter, s);
	
	do {
		iter++;
		_mscsFunctionFit_data.status = gsl_multifit_fdfsolver_iterate (s);
		if (_mscsFunctionFit_data.status!=0 and msgs->getVerbosity()>Low)
			printf ("status = %s\n", gsl_strerror (_mscsFunctionFit_data.status));
		
//		print_state (iter, s);
		
		if (_mscsFunctionFit_data.status)
			break;
		
		_mscsFunctionFit_data.status = gsl_multifit_test_delta (s->dx, s->x,1e-5, 1e-5);
	} while (_mscsFunctionFit_data.status == GSL_CONTINUE && iter < 5000);
	
	 gsl_multifit_covar (s->J, 0.0, covar); // BLcomment (Aug 28, 2017, 7:50:58 PM): this does not compile with gsl>=2.3.1
//	gsl_multifit_covar(J,0.0, covar);
	
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	
	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - paramCount;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
//	printf("chi %lf, dof %lf, iter: %li, status: %i\n",chi,dof, iter, _mscsFunctionFit_data.status);
//	printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
	//	
	//	printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
	//	printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
	//	printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	
//	printf ("status = %s\n", gsl_strerror (_mscsFunctionFit_data.status));
	
	_mscsFunctionFit_data.cov.SetSize(initGuessVals.size(),initGuessVals.size());
	for (long i = 0; i < initGuessVals.size(); i++) {
		_mscsFunctionFit_data.fitedParam.append(FIT(i));
		_mscsFunctionFit_data.fitedParamErr.append(c*ERR(i));
		for (long j = 0; j < initGuessVals.size(); j++) {
			_mscsFunctionFit_data.cov(i,j)=gsl_matrix_get(covar,i,j);
		}
	}
	_mscsFunctionFit_data.chisq=pow(chi, 2.0) / dof;
	_mscsFunctionFit_data.DOF=dof;
	//	_mscsFunctionFit_data.bestFitModel.importFunction(d.x,d.y,d.n);
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	delete [] x;
	delete [] y;
	delete [] sigma;
	delete [] param_ini;
	
	if (cpeds_isnan(_mscsFunctionFit_data.chisq)) { 
		_mscsFunctionFit_data.chisq=-1;
		_mscsFunctionFit_data.status=-1; 
		printf("chisq is nan, something must be wrong\n");
	}
	else {
		_mscsFunctionFit_data.status=0;
	}
	
	return _mscsFunctionFit_data.fitedParam;
}
/***************************************************************************************/
void mscsFunctionFit::select_model(string model) {
	if (model=="exponential_3param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_exponential_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_exponential_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_exponential_model_fdf;		
		_mscsFunctionFit_data.parametersCount=3;
		_mscsFunctionFit_data.constsCount=0;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="skewgauss_4param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_skewgauss_4param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_skewgauss_4param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_skewgauss_4param_model_fdf;		
		_mscsFunctionFit_data.parametersCount=4;
		_mscsFunctionFit_data.constsCount=0;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="gauss_4param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_gauss_4param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_gauss_4param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_gauss_4param_model_fdf;		
		_mscsFunctionFit_data.parametersCount=4;
		_mscsFunctionFit_data.constsCount=0;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="2gauss_7param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_2gauss_7param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_2gauss_7param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_2gauss_7param_model_fdf;		
		_mscsFunctionFit_data.parametersCount=7;
		_mscsFunctionFit_data.constsCount=0;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="line_2param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_line_2param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_line_2param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_line_2param_model_fdf;
		_mscsFunctionFit_data.parametersCount=2;
		_mscsFunctionFit_data.constsCount=1;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="line_1paramB") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_line_1paramB_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_line_1paramB_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_line_1paramB_model_fdf;
		_mscsFunctionFit_data.parametersCount=1;
		_mscsFunctionFit_data.constsCount=2;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="powerLaw_2param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_powerLaw_2param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_powerLaw_2param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_powerLaw_2param_model_fdf;
		_mscsFunctionFit_data.parametersCount=2;
		_mscsFunctionFit_data.constsCount=1;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="powerLaw_3param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_powerLaw_3param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_powerLaw_3param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_powerLaw_3param_model_fdf;
		_mscsFunctionFit_data.parametersCount=3;
		_mscsFunctionFit_data.constsCount=2;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="powerLaw_1paramA") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_powerLaw_1paramA_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_powerLaw_1paramA_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_powerLaw_1paramA_model_fdf;
		_mscsFunctionFit_data.parametersCount=1;
		_mscsFunctionFit_data.constsCount=2;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="beta_3param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_beta_3param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_beta_3param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_beta_3param_model_fdf;		
		_mscsFunctionFit_data.parametersCount=3;
		_mscsFunctionFit_data.constsCount=0;
		_mscsFunctionFit_data.modelName=model;
	}
	if (model=="poly_4param") {
		_mscsFunctionFit_data.fit.f = &mscsFunctionFit_poly_4param_model_f;
		_mscsFunctionFit_data.fit.df = &mscsFunctionFit_poly_4param_model_df;
		_mscsFunctionFit_data.fit.fdf = &mscsFunctionFit_poly_4param_model_fdf;		
		_mscsFunctionFit_data.parametersCount=4;
		_mscsFunctionFit_data.constsCount=4;
		_mscsFunctionFit_data.modelName=model;
	}
	
}
/***************************************************************************************/
mscsFunction mscsFunctionFit::getBestFitModel(double from, double to, double acc) {
	mscsFunction m;
	m.mkConst(from,to,acc);
	return getBestFitModel(acc,&m);
}
mscsFunction mscsFunctionFit::getBestFitModel(double acc,mscsFunction* args) {
	_mscsFunctionFit_data.bestFitModel.clearFunction();
	checkRanges();
	
	if (getModel()=="exponential_3param") {
		_mscsFunctionFit_data.parametersCount=3;
	}
	if (getModel()=="skewgauss_4param") {
//		_mscsFunctionFit_data.bestFitModel.mkSkewGauss(getMinArg(),getMaxArg(),acc*(getX(pointsCount()-1)-getX(0))/pointsCount(),
		_mscsFunctionFit_data.bestFitModel.mkSkewGauss(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				_mscsFunctionFit_data.fitedParam[1],
				_mscsFunctionFit_data.fitedParam[2],
				_mscsFunctionFit_data.fitedParam[3]);
	}
	/***************************************************************************************/
	if (getModel()=="gauss_4param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
		_mscsFunctionFit_data.bestFitModel.mkGauss(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				_mscsFunctionFit_data.fitedParam[1],
				_mscsFunctionFit_data.fitedParam[2],
				_mscsFunctionFit_data.fitedParam[3]
				);
	}
	if (getModel()=="2gauss_7param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
		_mscsFunctionFit_data.bestFitModel.mkGauss(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				_mscsFunctionFit_data.fitedParam[1],
				_mscsFunctionFit_data.fitedParam[2],
				_mscsFunctionFit_data.fitedParam[6]
				);
		mscsFunction g2=_mscsFunctionFit_data.bestFitModel;
		g2.mkGauss(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[3],
				_mscsFunctionFit_data.fitedParam[4],
				_mscsFunctionFit_data.fitedParam[5],
				0.0
		);
		_mscsFunctionFit_data.bestFitModel+=g2;

	}
	if (getModel()=="line_2param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
//		shiftX(-mscsFunctionFit_model_const0);
//		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc*(getX(pointsCount()-1)-getX(0))/pointsCount(),
		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				_mscsFunctionFit_data.fitedParam[1],
				getModelConst(0),
				0);
//		shiftX(mscsFunctionFit_model_const0);
	}
	
	if (getModel()=="line_1paramB") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
//		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc*(getX(pointsCount()-1)-getX(0))/pointsCount(),
		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc,
				getModelConst(0),
				_mscsFunctionFit_data.fitedParam[0],
				getModelConst(1),
				0
				);
	}
	/***************************************************************************************/
	if (getModel()=="powerLaw_2param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
//		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc*(getX(pointsCount()-1)-getX(0))/pointsCount(),
		_mscsFunctionFit_data.bestFitModel.mkPowerLaw(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				getModelConst(0),
				_mscsFunctionFit_data.fitedParam[1]);
	}
	if (getModel()=="powerLaw_3param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
		_mscsFunctionFit_data.bestFitModel.mkPowerLaw(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				getModelConst(0), getModelConst(1),
				_mscsFunctionFit_data.fitedParam[2],
				_mscsFunctionFit_data.fitedParam[1]);
	}
	
	if (getModel()=="powerLaw_1paramA") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
//		_mscsFunctionFit_data.bestFitModel.mkLine(getMinArg(),getMaxArg(),acc*(getX(pointsCount()-1)-getX(0))/pointsCount(),
		_mscsFunctionFit_data.bestFitModel.mkPowerLaw(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				getModelConst(0),
				getModelConst(1)
				);
	}

	if (getModel()=="beta_3param") {
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
		}
		_mscsFunctionFit_data.bestFitModel.mkBetaModel(getMinArg(),getMaxArg(),acc,
				_mscsFunctionFit_data.fitedParam[0],
				_mscsFunctionFit_data.fitedParam[1],
				_mscsFunctionFit_data.fitedParam[2]
		);
	}
	if (getModel()=="poly_4param") {
		cpedsList<double> c;
		for (unsigned long i = 0; i < 4; i++) {
			c.append(_mscsFunctionFit_data.fitedParam[i]*getModelConst(i));			
		}
		if (args!=NULL) {
			_mscsFunctionFit_data.bestFitModel=(*args);
			_mscsFunctionFit_data.bestFitModel.mkPolynomial(c);
//			printf("not null %li\n",_mscsFunctionFit_data.bestFitModel.pointsCount());
			
		}
		else {			
			_mscsFunctionFit_data.bestFitModel.mkPolynomial(getMinArg(),getMaxArg(),acc,c);
//			printf("from: %lf to: %lf, dx: %lf\n",getMinArg(),getMaxArg(),acc);
//			printf("null %li\n",_mscsFunctionFit_data.bestFitModel.pointsCount());
		}
	}

	return _mscsFunctionFit_data.bestFitModel; 
}
/***************************************************************************************/
double& mscsFunctionFit::pValue() { 
	cpeds_chisq_accepted(0.05,1,_mscsFunctionFit_data.chisq,&_mscsFunctionFit_data.Pvalue);
//	cpeds_chisq_accepted(0.05,DOF(),chisq()*DOF(),&_mscsFunctionFit_data.Pvalue);
	return _mscsFunctionFit_data.Pvalue; 
}

/***************************************************************************************/
mscsFunction mscsFunctionFit::fit(mscsFunction& data, cpedsList<double> yerr, string model, cpedsList<double>& Pini, double bfModelRes, cpedsList<double>* Srad, long Nstates) {
	
	long Nparam;
	bool convOK;
	
	cpedsList<double> fitParams;

	mscsFunction bfModel;
	// load data and set model
	mscsFunction::operator=(data);
	setModel(model);
	Nparam=getModelParametersCount();
	msgs->say("The model has %li parameters",Nparam,High);
	
	// load errors
	msgs->say("Loading errors",Top);
	double err;
	if (yerr.size()==1) { 
		err=yerr[0];
		yerr.makeLength(pointsCount()); 
		yerr=err;
		msgs->say("Making same errors for all data points",High);
	}
	
	
	// define initial parameter values
	msgs->say("Setting initial parameter guess values",Top);
	
	if (Pini.size()!=Nparam) { msgs->say("The selected model requires %li parameters, but only %li are provided. Cannot continue",Nparam, long(Pini.size()), Top); exit(0); }

	
	// make fit

	msgs->say("Making fit",Top);
	fitParams=fitData(Pini,yerr,model);

	// saving fit
	
	long chisqMin_idx,chisqMax_idx;
	double chisqMin,chisqMax;
	cpedsList<double> chi_sq;
	vector< cpedsList<double> > params_val;
	vector< cpedsList<double> > params_err;

	if (_mscsFunctionFit_data.status==0) {
		params_val.push_back(fitParams);
		params_err.push_back(paramsErr());
		chi_sq.push_back(chisq());
	}

	msgs->say("Best fit parameters",High);
	if (getVerbosityLevel()>=Medium) fitParams.print();

	
	// define search radius
	if (Srad!=NULL) {
		
		if (Srad->size()!=Nparam) { msgs->say("The selected model requires %li parameters, but only %li are provided for search radius. Cannot continue",Nparam, long(Srad->size()), Top); exit(0); }
		

		vector< cpedsList<double> > params_ini;
		cpedsRNG rng("uniform");
		cpedsList<double> rns;
		rns.makeLength(Nstates);

		for (long j = 0; j < Nparam; j++) {
			if (Srad->at(j)>0) {
				rng.setMinMax(Pini[j]-Srad->at(j),Pini[j]+Srad->at(j));
				params_ini.push_back(rng.getRNs(Nstates));				
			}
			else {
				rns=Pini[j];
				params_ini.push_back(rns);
			}
		}
		
		for (long i = 0; i < Nstates; i++) {
			for (long j = 0; j < Nparam; j++) {
				Pini[j]=params_ini[j][i];
			}
			
			fitParams=fitData(Pini,yerr);
			if (_mscsFunctionFit_data.status==0) {
				chi_sq.append(chisq());
				params_val.push_back(fitParams);
				params_err.push_back(paramsErr());
			}
			if ((i % 1000) == 0) msgs->say("done: %lf %%",double(i)/Nstates*100, Low);
		}
	}

	// print results
	if (chi_sq.size() >0) {
		chi_sq.getMinMaxValues(&chisqMin,&chisqMax,&chisqMin_idx,&chisqMax_idx);
		params()=params_val[chisqMin_idx];
		paramsErr()=params_err[chisqMin_idx];
		chisq()=chi_sq[chisqMin_idx];
		msgs->say("Best fit chisq per DOF: %lE",chisqMin,High);
		msgs->say("DOF: %lf",DOF(),High);
		msgs->say("p-value [%%]: %lE ",pValue()*100,High);
		msgs->say("Best fit params: ",High);
		print_fit_result(params_val[chisqMin_idx],params_err[chisqMin_idx]);
		msgs->say("Data range used: [%lE, %lE]",getX(0),getX(pointsCount()-1),High);
		// save results
		
		_mscsFunctionFit_data.resid=getBestFitModel(1)-data;
		_mscsFunctionFit_data.residStDev=_mscsFunctionFit_data.resid.stdev();
		bfModel=getBestFitModel(bfModelRes);
	//	bestfit.save("bestFitModel");
	//	resid.save("resid");
		msgs->say("residual stdev: %lE",_mscsFunctionFit_data.residStDev,High);
	}
	else {
		msgs->say("Fitting failed",High);		
	}

	return bfModel;
}
/***************************************************************************************/
void mscsFunctionFit::print_fit_result(cpedsList<double>& p,cpedsList<double>& pErr) {
	for (long i = 0; i < p.size(); i++) {
		msgs->say("p%.0lf = %lE +- %lE",double(i),p[i],pErr[i],Low);
	}
}
/***************************************************************************************/
void mscsFunctionFit::setModelConst(int constParam, double value) {
	switch (constParam) {
		case 0:
			mscsFunctionFit_model_const0=value;
			break;
		case 1:
			mscsFunctionFit_model_const1=value;
			break;
		case 2:
			mscsFunctionFit_model_const2=value;
			break;
		case 3:
			mscsFunctionFit_model_const3=value;
			break;
		default:
			break;
	}
}
/***************************************************************************************/
double mscsFunctionFit::getModelConst(int constParam) {
	switch (constParam) {
		case 0:
			return mscsFunctionFit_model_const0;
			break;
		case 1:
			return mscsFunctionFit_model_const1;
			break;
		case 2:
			return mscsFunctionFit_model_const2;
			break;
		case 3:
			return mscsFunctionFit_model_const3;
			break;
		default:
			return 0;
			break;
	}
}
