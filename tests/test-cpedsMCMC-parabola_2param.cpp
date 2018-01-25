#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"

class parabolaMCMC : public cpedsMCMC {
	public:
		parabolaMCMC(cpedsMCMC& parent) : cpedsMCMC(parent) { }
		parabolaMCMC() { msgs->setSender("parabolaMCMC"); }
		parabolaMCMC(long Ndim) : cpedsMCMC(Ndim) { msgs->setSender("parabolaMCMC"); }
//		parabolaMCMC(mscsFunction& data) { _data=data; }
		
		//		virtual double chisq(matrix<double>& data, const MClink &th, matrix<double>& cov) {
		double chisq(const MClink& th) {
			//
			// define the model data
			//
			chisqData().model.mkPowerLaw(th[0],1,2);
			chisqData().model+=th[1];
			
			//
			// calculate chisq for this model
			//
			mscsFunction X=getData()-getModel();
			X.power(2.0);
			return X.sumY()/getCov();
		}
		
		
};


void calculateFisher(mscsFunction& data, double yerr, string outdir, mscsFunction3dregc estimated, double bfa, double bfb) {
	matrix<double> F(2,2);
	long N=data.pointsCount();
	double sigma2=yerr*yerr;

	F=0.0;
	// 11 term of the Fisher matrix. Parameter a
	for (long k = 0; k < N; k++) {	F(0,0)+=pow(data.getx(k),4)/sigma2;	}
	// 12 and 21 terms of the Fisher matrix. 
	for (long k = 0; k < N; k++) {	F(0,1)+=pow(data.getx(k),2)/sigma2;	}
	F(1,0)=F(0,1);
	// 22 term of the Fisher matrix. Parameter b
	F(1,1)=double(N)/sigma2;

	cout << "Fisher information matrix\n";
	cpeds_print_matrix(F);
	cpeds_matrix_save(F,outdir+"Fisher.mat");

	matrix<double> C=F.Inv();
	cout << "Parameters correlation matrix\n";
	cpeds_print_matrix(C);
	cpeds_matrix_save(C,outdir+"Correlation.mat");
	
	printf("Fisher parameter1 sigma: %lE\n",sqrt(C(0,0)));
	printf("Fisher parameter2 sigma: %lE\n",sqrt(C(1,1)));
	
	//calculate PDF
	mscsFunction3dregc pdf=estimated; // set size and ranges as in the reconstructed PDF
	
//	pdf.setSizeRange(500,500,1,0.986920832971,0.317682574679,0,1.01512912914,1.59310328613,0);
//	pdf.allocFunctionSpace();

	pdf.mkGauss2D(bfa,bfb,C(0,0),C(1,1),C(0,1),1,2,0);
	pdf.saveHDF5(outdir+"Fisher2d.hdf5","L");
	
}

double getPDF_fwhm(MscsPDF1D pdf) {
	double fwhm;
	
	pdf-=0.5;
	mscsVector<double> r=pdf.findRoot();
	
	cpedsList<double> l(r);
	l.sort(12);
	if (l.size()<2) {
		printf("cannot calculate PDF FWHM\n"); exit(-1);
	}
	fwhm=l.last()-l[0];
	
	return fwhm;
}

void calculate_sigma_1Dposterior(MscsPDF1D p1, MscsPDF1D p2) {

	// calculate sigma 
	double s1,s2;
	s1=cpeds_fwhm2sigma(getPDF_fwhm(p1));
	s2=cpeds_fwhm2sigma(getPDF_fwhm(p2));
	
	printf("parameter1 sigma: %lE\n",s1);
	printf("parameter2 sigma: %lE\n",s2);
}

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int main() {
	string outdir="parabola_test2/";
	mkdir(outdir.c_str(),0755);
	
	//
	// PREPARE THE OBSERVATIONAL DATA
	//
	double xmin=-10;
	double xmax=10;
	double a=1, b=1;
	mscsFunction data, model;
	double yerr=10;
	
	data.mkPowerLaw(xmin,xmax,0.1,a,1,2); // generate y=a*x^2+b with 
	data+=b;
	model=data;
	model.save(outdir+"parabola-model.txt");

	cpedsRNG rns;
	rns.setRNsType(rns.gaussian);
	rns.seed(1);
	rns.setMeanVariance(0.0,yerr);
//	rns.setPDF(data.pointsCount(),data.extractArguments(),data.extractValues());
	//	rns.getRNs(100).save("parabola-data.txt");
	
	for (long i = 0; i < data.pointsCount(); i++) {		data.f(i)+=rns.getRN();	}
	data.save(outdir+"parabola-data.txt");
	

	
	
	
	// set the dimensionality of the problem
	cpedsMCMC parab_conf(2), parabMC;

	//
	// set MC parameters
	//
	parab_conf.setOutputDir(outdir);
	parab_conf.setCoolingRate(1);
//	parab_conf.setInitialStepSize(0.5);
//	parab_conf.setInitialWalkStepSize(0.05);
	parab_conf.setBurnInLength(1000);
	parab_conf.setChisqSignature("parabola fit");
	parab_conf.setMaximalRejectionsCount(50);
	parab_conf.setWalkInfoOutputFrequency(500);
	
	
	// define the parameter space
	parab_conf.addParameter("a", -20,20,2);
	parab_conf.addParameter("b", -20,20,2);

//	parabmc.saveStepPDFs("parabola-stepPDFini");

	// set the observational data
	parab_conf.setData(data);
	parab_conf.setDiagonalCovarianceMatrix(yerr*yerr);
	parab_conf.printInfo();
	parab_conf.setCalculate2Dposteriors(false);
	parabMC=parab_conf;
	
	//
	// RUN THE MCMC CHAIN
	//
	long id;
#pragma omp parallel for private(id)
	for (id = 0; id < 20; ++id) {
		parabolaMCMC mc(parab_conf);
		mc.setVerbocity(Medium);
		mc.setID(id);
		mc.msgs->setSender("parabolaMCMC."+mc.msgs->toStr(id));
		mc.startChain();
#pragma omp critical
		{
		parabMC.append(mc);
		}
	}
	parabMC.saveResults();
	
//	cpeds_fwhm2sigma()
	
	//
	// calculate Fisher information matrix
	//
	calculateFisher(data,yerr, outdir, parabMC.get2Dposterior("a","b"), parabMC.getBestFitLink().getParam(0),parabMC.getBestFitLink().getParam(1));
	
	//
	// calculate comparizon between 1D calculated posterior and the theoretical predictions
	//
	calculate_sigma_1Dposterior(parabMC.get1Dposterior(0), parabMC.get1Dposterior(1));
	
	
	return 0;
}





