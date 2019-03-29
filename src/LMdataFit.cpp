/*!
 * \file LMdataFit.cpp
 *
 *  Project: Mscs
 *  Created on: Apr 11, 2013 11:23:56 AM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
//#include <QtCore/QStringList>
//#include <QtCore/QString>
#include "cpeds-msgs.h"
#include "Mscs-global-defs.h"
#include "mscsFunctionFit.h"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;

string _data, _model, _Pini,_modelConst,_searchRadius,_calcPvalue, _outNamePref, _dataCols;
double _yerr, _bfModelRes,_bfModelFrom,_bfModelTo;
long _Nstates,_bfModelPts;
bool _dontSaveThisRun, _logX, _logY;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

void calculatePvalue();

void print_fit_result(cpedsList<double>& p,cpedsList<double>& pErr, cpedsMsgs& msgs);
void save_fit_result(cpedsList<double>& p,cpedsList<double>& pErr, string fname);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".LMdataFit",false,"LMdataFit.log",High);
	msgs.setSaveRunWriteMode('a');

	if (_dontSaveThisRun==false)
		msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	
	if (_calcPvalue!="") {
		calculatePvalue();
		exit(0);
	}
	
	mscsFunction data;
	vector<long> dataCols=Mscs_stringToLongList(_dataCols,",");
	mscsFunctionFit fit;
	long Nparam;
	long Nconsts;
	bool convOK;
	cpedsList<double> yerr;
	mscsFunction yerrTmp;
	cpedsList<double> fitParams;
	string dataFile=basename(_data.c_str());
	string outFilePrefix=dataFile+"--fit_"+_model;

	// load data and set model
	msgs.say("Loading data",Top);
	data.load(_data,true,dataCols[0],dataCols[1]);
	if (dataCols.size()==3)
		yerrTmp.load(_data,true,dataCols[0],dataCols[2]);
	if (_logX) {
		// check for sanity of data
		data.removeValue(0);
		data.removeSmaller(0);
		yerrTmp.removeValue(0);
		yerrTmp.removeSmaller(0);
		data.logX(10);
	}
	if (_logY) {
		data.removeValue(0);
		data.removeSmaller(0);
		yerrTmp.removeValue(0);
		yerrTmp.removeSmaller(0);
		data.logY(10);
	}
	data.save(outFilePrefix+".in");
	
//	data.sortFunctionArgAscending();
//	fit.load(_data);
	fit.concatenate(data);
//	fit.sortFunctionArgAscending();
	fit.checkRanges();
	fit.setModel(_model);
	Nparam=fit.getModelParametersCount();
	Nconsts=fit.getModelConstsCount();
	msgs.say("The model has %li parameters",Nparam,High);
	
	// load errors
	msgs.say("Loading errors",Top);
	yerr.makeLength(fit.pointsCount()); 
	if (_yerr!=-1) { 
		if (_logY)
			yerr=log10(_yerr);
		else
			yerr=_yerr;

		msgs.say("Making same errors for all data points",High);
	}
	else {
		msgs.say("Loading errors from 3rd column of data file",High);
		if (cpeds_get_cols_num_first_ln(_data.c_str(),false,1000) == 2) {
			msgs.criticalError("No yerr value set, and only two columns found in the data file, so I do not know the yerrors.", Top);
		}
		yerr=yerrTmp.toYList();
//		matrix<double> m =cpeds_matrix_load(_data);
//		for (long i = 0; i < m.RowNo(); i++) {
//			if (_logY)
//				yerr.append(log10(m(i,dataCols[2])));
//			else 
//				yerr.append(m(i,dataCols[2]));
//		}
	}
	
	// define initial parameter values
	msgs.say("Setting initial parameter guess values",Top);
	cpedsList<double> Pini;
	QString qs=_Pini.c_str();
	QStringList qsl=qs.split(",");
	for (long i = 0; i < Nparam; i++) {
		Pini.append(qsl[i].toDouble(&convOK));
		if (convOK==false) msgs.criticalError("could not convert initial parameter value to double.", Top);
		msgs.say("setting initial guess value for parameter %.0lf to %lE",double(i),Pini.last(),Medium);
	}
	if (Pini.size()!=Nparam) { msgs.say("The selected model requires %li parameters, but only %li are provided. Cannot continue",Nparam, long(qsl.size()), Top); exit(0); }

	// define model constants
	msgs.say("Setting model constants",Top);
	cpedsList<double> modelConsts;
	QString qsConsts=_modelConst.c_str();
	QStringList qslConsts;
	
	if (_modelConst!="") {
		qslConsts=qsConsts.split(",");
	}
//	printf("%li, |%s|\n",qslConsts.size(),qslConsts[0].toStdString().c_str());
//	exit(0);

	for (long i = 0; i < qslConsts.size(); i++) {
		modelConsts.append(qslConsts[i].toDouble(&convOK));
		if (convOK==false) msgs.criticalError("could not convert the model constant value to double.", Top);
	}
	if (modelConsts.size()!=Nconsts) { msgs.say("The selected model requires %li consts, but only %li are provided. Cannot continue",Nconsts, long(qslConsts.size()), Top); exit(0); }
	for (unsigned long i = 0; i < modelConsts.size(); i++) {
		msgs.say("setting model constant %.0lf to %lE",double(i),modelConsts[i],Medium);
		fit.setModelConst(i,modelConsts[i]);
	}
//	if (modelConsts.size()!=Nparam) { msgs.say("The selected model requires %li parameters, but only %li are provided. Cannot continue",Nparam, long(qsl.size()), Top); exit(0); }

	

	
	// make fit

	msgs.say("Making fit",Top);
	fitParams=fit.fitData(Pini,yerr);
	msgs.say("Best fit parameters",High);
	fitParams.print();

	long chisqMin_idx,chisqMax_idx;
	double chisqMin,chisqMax;
	cpedsList<double> chisq;
	vector< cpedsList<double> > params;
	vector< cpedsList<double> > paramsErr;
	params.push_back(fitParams);
	paramsErr.push_back(fit.paramsErr());
	chisq.push_back(fit.chisq());
	
	// define search radius
	if (_searchRadius!="") {
		msgs.say("Defining grid search radius",Top);
		cpedsList<double> Srad;
		QString qs=_searchRadius.c_str();
		QStringList qsl=qs.split(",");
		for (long i = 0; i < Nparam; i++) {
			Srad.append(qsl[i].toDouble(&convOK));
			if (convOK==false) msgs.criticalError("could not convert grid search radius parameter value to double.", Top);
		}
		if (Srad.size()!=Nparam) { msgs.say("The selected model requires %li parameters, but only %li are provided for search radius. Cannot continue",Nparam, long(qsl.size()), Top); exit(0); }
		

		vector< cpedsList<double> > params_ini;
		cpedsRNG rng("uniform");
		cpedsList<double> rns;
		rns.makeLength(_Nstates);

		for (long j = 0; j < Nparam; j++) {
			if (Srad[j]>0) {
				rng.setMinMax(Pini[j]-Srad[j],Pini[j]+Srad[j]);
				params_ini.push_back(rng.getRNs(_Nstates));				
			}
			else {
				rns=Pini[j];
				params_ini.push_back(rns);
			}
		}
		
		for (long i = 0; i < _Nstates; i++) {
			for (long j = 0; j < Nparam; j++) {
				Pini[j]=params_ini[j][i];
			}
			
			fitParams=fit.fitData(Pini,yerr);
			chisq.append(fit.chisq());
			params.push_back(fitParams);
			paramsErr.push_back(fit.paramsErr());
			if ((i % 1000) == 0) msgs.say("done: %lf %%",double(i)/_Nstates*100, Low);
		}
	}

	// print results
	if (_outNamePref!="") {
		outFilePrefix=_outNamePref+"-"+outFilePrefix;
	}
	chisq.getMinMaxValues(&chisqMin,&chisqMax,&chisqMin_idx,&chisqMax_idx);
	fit.params()=params[chisqMin_idx];
	fit.paramsErr()=paramsErr[chisqMin_idx];
	fit.chisq()=chisq[chisqMin_idx];
	msgs.setLogFileName(outFilePrefix+".params");
	msgs.loggingOn();
	msgs.say("Model name: "+fit.getModel(),High);
	msgs.say("Best fit chisq per DOF: %lE",chisqMin,Medium);
	msgs.say("DOF: %lf",fit.DOF(),Medium);
	msgs.say("p-value [%%]: %lE ",fit.pValue()*100,Medium);
	msgs.say("Best fit params: ",Medium);
	print_fit_result(params[chisqMin_idx],paramsErr[chisqMin_idx],msgs);
	msgs.loggingOff();
	save_fit_result(params[chisqMin_idx],paramsErr[chisqMin_idx],outFilePrefix+".params.val_err");
	msgs.say("Data range used: [%lE, %lE]",fit.getX(0),fit.getX(fit.pointsCount()-1),Medium);

//	exit(0);
	// save results
	data.checkRanges();
	double modelRes=(data.getMaxArg()-data.getMinArg())/data.pointsCount();
	mscsFunction resid=fit.getBestFitModel(modelRes,&data)-data;
	mscsFunction residRel=resid/data;
	if (_bfModelRes==-1) {
//		data.checkRanges();
		_bfModelRes=(data.getMaxArg()-data.getMinArg())/_bfModelPts;
		msgs.say("data.getMaxArg(): %lE",data.getMaxArg(),High);
		msgs.say("data.getMinArg(): %lE",data.getMinArg(),High);
		msgs.say("setting best fit model resolution: %lE",_bfModelRes,High);
		msgs.say("setting best fit model points: %li",_bfModelPts,High);
//		exit(0);
	}
	mscsFunction bfm;
	if (_bfModelFrom!=_bfModelTo) {
		bfm=fit.getBestFitModel(_bfModelFrom,_bfModelTo,_bfModelRes);
		if (_logX) bfm.exponentX(10);
		if (_logY) bfm.exponentY(10);
		bfm.save(outFilePrefix+"-bestFitModel");
	}
	else {
		if (_bfModelRes==-2) {
			bfm=fit.getBestFitModel(_bfModelRes,&data);
			if (_logX) bfm.exponentX(10);
			if (_logY) bfm.exponentY(10);
			bfm.save(outFilePrefix+"-bestFitModel");								
		}
		else {
			bfm=fit.getBestFitModel(_bfModelRes);
			if (_logX) bfm.exponentX(10);
			if (_logY) bfm.exponentY(10);
			bfm.save(outFilePrefix+"-bestFitModel");					
		}
	}
	if (_logX) { resid.exponentX(10); residRel.exponentX(10); }
	if (_logY) { resid.exponentY(10); residRel.exponentY(10); }
	resid.save(outFilePrefix+"_resid");
	msgs.say("residual stdev: %lE",resid.stdev(),High);
	residRel.save(outFilePrefix+"_residRel");
//	msgs.say("residual stdev: %lE",resid.stdev(),High);


	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("LMdataFit\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: LMdataFit"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
		ValueArg<double> yerr("", "yerr", "amplitude of measurement uncertainty (if same for all data points) to be used and override values from the input data file", false,-1,"double");	cmd.add( yerr );
		ValueArg<double> bfModelRes("r", "bfres", "best fit model resolution parameter. (default: 1 - same as input data for a regular spacing). "
				"Fractional values indicate the fraction of the input data"
				"separation that will be used for generating the best fit model for saving."
				"If -1 then it will be set according to the input data range and --bfPts count."
				"If -2 then it will be probed at the same arguments as the input data", false,1,"double");	cmd.add( bfModelRes );
		ValueArg<double> bfModelFrom("", "bffrom", "best fit model starts from this value. (default: 0). ", false,0,"double");	cmd.add( bfModelFrom );
		ValueArg<double> bfModelTo("", "bfto", "best fit model ends at this value. (default: 0). If zero and the same as bffrom then it will be taken from the data", false,0,"double");	cmd.add( bfModelTo );
		ValueArg<long> bfModelPts("", "bfPts", "best fit points number. (default: 0 - same as input data). If not given then --bfres is used", false,0,"long");	cmd.add( bfModelPts );
		
		SwitchArg dontSaveThisRun("","dontSaveThisRun", "do not save this run -- useful for parallel runs which could conflict by writing to the same files (default: false)", false);	cmd.add( dontSaveThisRun );
		SwitchArg logX("","logX", "take log10(X) after loading data and before doing anything else (default: false)", false);	cmd.add( logX );
		SwitchArg logY("","logY", "take log10(Y) after loading data and before doing anything else (default: false)", false);	cmd.add( logY );
		
		ValueArg<string> data("d","data","input data (function type, x,y,[yerr])",true,"","string"); cmd.add(data);
		ValueArg<string> dataCols("","dataCols","comma separated list of columns in the input data file to be used for fitting."
				"The first value should be the x-data, second - the y-data, and the third - yerror data (default: 0,1). Use 0,1,2 to use the first three columns ",false,"0,1","string"); cmd.add(dataCols);
		ValueArg<string> modelConst("","const","comma separated list of model const values. This applies only to the models that support const values."
				"The only model that supports that now is line_1paramB where the const value defines the A parameter of the fitted line",false,"","string"); cmd.add(modelConst);
		ValueArg<string> Pini("p","Pini","comma separated list of initial guess values for all parameters",true,"","string"); cmd.add(Pini);
		ValueArg<string> searchRadius("","gsr","comma separated list of grid search radii. Such defined parameter space region will be searched and chisq values recorded"
				"in order to fine tune the fit, which would otherwise stack at some non-optimal point in the parameter space. If zero, then the corresponding"
				"parameter will be set to its initial value as defined by Pini option",false,"","string"); cmd.add(searchRadius);
		ValueArg<long> Nstates("N", "Nstates", "Total number of fitting runs with initial values drawn randomly from within searchRadius. (default: 10000)", false,10000,"long");	cmd.add( Nstates );

		ValueArg<string> pValue("", "pValue", "calculates p-value for given chisq,DOF value (chisq,DOF should be comma separated)", false,"","string");	cmd.add( pValue );

		ValueArg<string> outNamePref("","pref","output file names prefix. - is appended after it if the string is different than '' (default: '') ",false,"","string"); cmd.add(outNamePref);

		std::vector<string> allowedStr;
		allowedStr.push_back("exponential_3param");
		allowedStr.push_back("skewgauss_4param");
		allowedStr.push_back("gauss_4param");
		allowedStr.push_back("2gauss_7param");
		allowedStr.push_back("line_2param");
		allowedStr.push_back("line_1paramB");
		allowedStr.push_back("powerLaw_2param");
		allowedStr.push_back("powerLaw_3param");
		allowedStr.push_back("powerLaw_1paramA");
		allowedStr.push_back("beta_3param");
		allowedStr.push_back("poly_4param");
		
		ValuesConstraint<string>* allowedStrNew;
		allowedStrNew = new ValuesConstraint<string>(allowedStr);
		ValueArg<string> model("m", "model", "theoretical model name.\n"
				"\n"
				"Model descriptions:\n"
				"\n"
				"line_2param - (A*(x-const0)+B) requires 2 parameters and 1 constant\n"
				"line_1paramB - (const0*(x-const1)+B) requires 1 parameter and 2 constants\n"
				"powerLaw_2param - (A*(x/const0)^B requires 2 parameters and 1 constant\n"
				"powerLaw_3param - (A*((x-const0)/const1)^n+B requires 3 parameters (A,B,n) and 2 constants\n"
				"powerLaw_1paramA - (A*(x/const0)^const1 requires 1 parameter and 2 constants\n"
				"beta_3param - n(r) = n0 ( 1+ (r/rc)^2 )^(-3/2*beta), where n0,rc and beta are fitted parameters\n"
				"poly_4param - a0 c0 + a1 c1 x + a2 c2 x^2 + a3 c3 x^3 requires 4 parameters and 4 constants\n"
				"skewgauss_4param - THe fitted parameters are A,m,s,alpha of the model: gauss * (1.0-gsl_sf_erf(nsig/sq2)) where gauss=A * exp(-(x-m)*(x-m)/(2*s*s))"
				"and nsig=(m-x)*alpha/s. There are no constants in this model\n"
				"gauss_4param - The fitted parameters are A,m,s,B of the model: gauss=A * exp(-(x-m)*(x-m)/(2*s*s))+B\n"
				"2gauss_7param - The fitted parameters are A1,m1,s1,A2,m2,s2,B of the model: gauss=A1 * exp(-(x-m1)*(x-m1)/(2*s1*s1)) + A2*exp(-(x-m2)*(x-m2)/(2*s2*s2))+B\n"
				" ", true,"skewgauss_4param", allowedStrNew);	cmd.add( model );
//		allowedStr.clear();

		//     UnlabeledValueArg<string> input_file("input_file","file to plot",true,"nofile","string");	cmd.add( input_file );
		//     ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		//     ValueArg<string> proj("p","proj","projection to use [mall - Mollweide/ aitoff etc] (default: mall)",false,"mall","string"); cmd.add(proj);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_data=data.getValue();
		_model=model.getValue();
		_Pini=Pini.getValue();
		_modelConst=modelConst.getValue();
		_yerr=yerr.getValue();
		_bfModelRes=bfModelRes.getValue();
		_bfModelFrom=bfModelFrom.getValue();
		_bfModelTo=bfModelTo.getValue();
		_bfModelPts=bfModelPts.getValue();
		if (bfModelPts.isSet()) _bfModelRes=-1;
		_dontSaveThisRun=dontSaveThisRun.getValue();
		_logX=logX.getValue();
		_logY=logY.getValue();
		
		_searchRadius=searchRadius.getValue();
		_Nstates=Nstates.getValue();
//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";

		_calcPvalue=pValue.getValue();
		
		_outNamePref=outNamePref.getValue();
		_dataCols=dataCols.getValue();
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}

string getProgramVersionString() {
	string rev;
#ifdef GENERATE_CVS_VERSION_STRING
	rev="$Revision: 1.1 $";
	QString qrev=rev.c_str();
	QStringList qrevl=qrev.split(" ");
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
	else rev="could not retrieve version number";
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif

	return rev;
}
/***************************************************************************************/
void print_fit_result(cpedsList<double>& p,cpedsList<double>& pErr, cpedsMsgs& msgs) {
	for (long i = 0; i < p.size(); i++) {
		msgs.say("p%.0lf = %lE +- %lE",double(i),p[i],pErr[i],Low);
	}
}
/***************************************************************************************/
void save_fit_result(cpedsList<double>& p,cpedsList<double>& pErr, string fname) {
	FILE* f=fopen(fname.c_str(),"w");
	for (unsigned long i = 0; i < p.size(); i++) {
		fprintf(f,"%.10lE %.10lE\n",p[i],pErr[i]);		
	}
	fclose(f);
}
/***************************************************************************************/
void calculatePvalue() {
	QString qs=_calcPvalue.c_str();
	QStringList qsl=qs.split(",");
	double chisq=qsl[0].toDouble();
	double dof=qsl[1].toDouble();
	double p;
	printf("dof: %lf\n",dof);
	printf("chisq: %lf\n",chisq);
	printf("chisq/DOF: %lf\n",chisq/dof);

	cpeds_chisq_accepted(0.05,dof,chisq,&p);
	printf("p-value(chisq, DOF=dof): %lE\n",p);
	cpeds_chisq_accepted(0.05,1,chisq/dof,&p);
	printf("p-value(chisq/dof, DOF=1): %lE\n",p);
//	cpeds_chisq_accepted(0.05,dof,chisq/dof,&p);
//	printf("p-value(chisq/dof, DOF=dof): %lE\n",p);
//	cpeds_chisq_accepted(0.05,1,chisq,&p);
//	printf("p-value(chisq, DOF=1): %lE\n",p);
}
