#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "Mscs-function3dregc.h"
#include "cpedsMCMC.h"
#include <boost/program_options.hpp>
#include "omp.h"

mscsFunction3dregc getGaussModel(double x0, double y0) {
	mscsFunction3dregc model;
	mscsFunction f;
	model.setSizeRange(21,21,1,-10,-10,0,10,10,0);
	model.allocFunctionSpace();
	model.mkGauss2D(x0,y0,3,3,1,2,0);
	return model;
}
class gaussxyMCMC : public cpedsMCMC {
	public:
		gaussxyMCMC(cpedsMCMC& parent) : cpedsMCMC(parent) { }
		gaussxyMCMC() { msgs->setSender("gaussxy"); }
		gaussxyMCMC(long Ndim) : cpedsMCMC(Ndim) { msgs->setSender("gaussxy"); }
//		parabolaMCMC(mscsFunction& data) { _data=data; }
		
		//		virtual double chisq(matrix<double>& data, const MClink &th, matrix<double>& cov) {
		double chisq(const MClink& th) {
			//
			// define the model data
			//
			mscsFunction3dregc model=getGaussModel(th[0],th[1]);
			mscsFunction f;
			for (long i = 0; i < model.Nx(); i++) {		
				for (long j = 0; j < model.Ny(); j++) {		
					f.newPoint(0,model(i,j));
				}
			}

			chisqData().model=f;
			
			//
			// calculate chisq for this model
			//
			mscsFunction X=getData()-getModel();
			X.power(2.0);
			return X.sumY()/getCov();///X.pointsCount();
		}
		
		
};


double gauss2d(double x, double y, double x0, double y0) {
	return exp(-pow(x-x0,2)-pow(y-y0,2));
}

void calculateFisher(mscsFunction3dregc& data, double zerr, string outdir, mscsFunction3dregc estimated, double x0, double y0) {
	matrix<double> F(2,2);
	double sigma2=zerr*zerr;

	F=0.0;
	// 11 term of the Fisher matrix. Parameter x0
	for (long i = 0; i < data.Nx(); i++) {
		double x=data.getX(i);
		for (long j = 0; j < data.Ny(); j++) {
			double y=data.getY(j);
			double fij=gauss2d(x,y,x0,y0);
			F(0,0)+=1./sigma2*(
					2*fij*(data(i,j)-fij) +
					2*fij*fij*pow(x-x0,2)-
					2*fij*pow(x-x0,2)*(data(i,j)-fij)
					);
		}
	}
	// 12 and 21 terms of the Fisher matrix. 
	for (long i = 0; i < data.Nx(); i++) {
		double x=data.getX(i);
		for (long j = 0; j < data.Ny(); j++) {
			double y=data.getY(j);
			double fij=gauss2d(x,y,x0,y0);
			F(0,1)+=2*fij*(x-x0)*(y-y0)/sigma2 -
					2*fij*(x-x0)*(y-y0)*(data(i,j)-fij)/sigma2;
		}
	}

	F(1,0)=F(0,1);
	// 22 term of the Fisher matrix. Parameter y0
	for (long i = 0; i < data.Nx(); i++) {
		double x=data.getX(i);
		for (long j = 0; j < data.Ny(); j++) {
			double y=data.getY(j);
			double fij=gauss2d(x,y,x0,y0);
			F(1,1)+=1./sigma2*(
					2*fij*(data(i,j)-fij) +
					2*fij*fij*pow(y-y0,2)-
					2*fij*pow(y-y0,2)*(data(i,j)-fij)
					);
		}
	}

	cout << "Fisher information matrix\n";
	cpeds_print_matrix(F);
	cpeds_matrix_save(F,outdir+"Fisher.mat");

	matrix<double> C=F.Inv();
	cout << "Parameters correlation matrix\n";
	cpeds_print_matrix(C);
	cpeds_matrix_save(C,outdir+"Covariance.mat");
	
	printf("Fisher parameter1 sigma: %lE\n",sqrt(C(0,0)));
	printf("Fisher parameter2 sigma: %lE\n",sqrt(C(1,1)));
	
	//calculate PDF
	mscsFunction3dregc pdf=estimated; // set size and ranges as in the reconstructed PDF
	
//	pdf.setSizeRange(500,500,1,0.986920832971,0.317682574679,0,1.01512912914,1.59310328613,0);
//	pdf.allocFunctionSpace();

	pdf.mkGauss2D(x0,y0,C(0,0),C(1,1),C(0,1),1,2,0);
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
boost::program_options::variables_map parser(int argc, char** argv) {
	namespace po = boost::program_options;
	po::variables_map vm;

	
	po::options_description desc("test-cpedsMCMC-gaussXY-independent\n"
			"\n"
			"This program tests the parameter CR reconstruction for the "
			"case of two-parameter model with zero covariance.\n"
			"The reconstructed CR should have a circle like shape.\n\n"
			"The size of the 1-sigma CR along a given parameter dimension"
			"on 2-D plot is wider by a factor \n\n"
			"    The 2-D confidence region along a given parameter direction is broader"
			"by a factor\n\n "
    
			"Solve[Integrate[f[r] r, {r, 0, x}]/Integrate[f[r] r, {r, 0, \[Infinity]}] == Erf[Nsigma/Sqrt[2]], x]"
			"\n"
			"which gives x=Sqrt[-2 Log[Erfc[Nsigma/Sqrt[2]]]]\n"
			"\n"
			"which gives ~1.52 for 68\% confidence region.\n"
			"This is what is indeed observed."
			"\n\n"
			"Allowed options");
	
	int opt=1;
	long longopt=1;
	bool ompnested;
	long Bin;
	double dbl;
	string stropt;
	desc.add_options()
	    ("help", "produce help message")
	    ("Nchains,N", po::value<int>(&opt)->default_value(1), "Number of mcmc chains")
	    ("ompnested", po::value<bool>(&ompnested)->default_value(false), "Enables nested omp")
	    ("num_threads", po::value<int>(), "sets maximal number of omp threads")
	    ("Bin", po::value<long>(&Bin)->default_value(10000), "Burn-in length")
	    ("Bout", po::value<long>(&Bin)->default_value(100000), "Burn-out length")
	    ("ctol", po::value<double>(&dbl)->default_value(1.e-10), "X2 improvement convergence tolerance")
	    ("odir", po::value<string>(&stropt)->default_value("test-gauss_2param_indep"), "output directory")
	    ("zerr", po::value<double>(&dbl)->default_value(0.01), "Independent variable error. ")
	    ("x0", po::value<double>(&dbl)->default_value(1.), "model parameter 0: x0")
	    ("y0", po::value<double>(&dbl)->default_value(1.), "model parameter 1: y0")
	    ("outFreq", po::value<int>(&opt)->default_value(500), "outout frequency")
		("Nrej", po::value<long>(&longopt)->default_value(1), "number of rejections to cool down")
	;

	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
	    cout << desc << "\n";
	    exit(1);
	}

/*
	if (vm.count("compression")) {
	    cout << "Compression level was set to " 
	 << vm["compression"].as<int>() << ".\n";
	} else {
	    cout << "Compression level was not set.\n";
	}	
*/

	return vm;
}

int main(int argc, char** argv) {
	
	boost::program_options::variables_map opt=parser(argc,argv);
	
	
	string outdir=opt["odir"].as<string>()+"/";
	mkdir(outdir.c_str(),0755);
	
	//
	// PREPARE THE OBSERVATIONAL DATA
	//
	double xmin=-10;
	double xmax=10;
	double x0=opt["x0"].as<double>(), y0=opt["y0"].as<double>();
	mscsFunction model, data;
	double zerr=opt["zerr"].as<double>();
//	double yerr=1e-10;
	
	mscsFunction3dregc data3d=getGaussModel(x0,y0);

	cpedsRNG rns;
	rns.setRNsType(rns.gaussian);
	rns.seed(1);
	rns.setMeanStd(0.0,zerr);
//	rns.setPDF(data.pointsCount(),data.extractArguments(),data.extractValues());
	//	rns.getRNs(100).save("parabola-data.txt");

	cpedsList<double> noise;
	long k=0;
	if (opt["zerr"].as<double>()>0)
		for (long i = 0; i < data3d.Nx(); i++) {		
			for (long j = 0; j < data3d.Ny(); j++) {		
				double n=rns.getRN();
				noise.append(n);
				data3d(i,j)+=n;	
				data.newPoint(k++,data3d(i,j));
			}
		}
//	data=data3d.e
	data3d.savetxtlin(outdir+"gaussxy-data.txt",true);
	cout << "noise variance: "<< noise.variance() << endl;
	cout << "noise std: "<< sqrt(noise.variance())  << endl;
//	exit(0);
	
	
	
	
	// set the dimensionality of the problem
	cpedsMCMC confMC(2), gaussMC;

	//
	// set MC parameters
	//
	confMC.setOutputDir(outdir);
	confMC.setUpHillClimbing(true);
	confMC.setUpHillGradient(true);
	confMC.setConvergenceThres(opt["ctol"].as<double>());
	confMC.setAcceptWorsePvalues(0,0);
//	parab_conf.setCoolingRate(1);
	confMC.setInitialStepSize(1);
	confMC.setInitialWalkStepSize(1);
	confMC.setMaximalRejectionsCount(opt["Nrej"].as<long>());
	confMC.setBurnInLength(opt["Bin"].as<long>());
	confMC.setBurnOutLength(opt["Bout"].as<long>());
	confMC.setTemperatures(1e10,1e-30);
//	confMC.setChisqSignature("parabola fit");
	confMC.setWalkInfoOutputFrequency(opt["outFreq"].as<int>());
	
	
	// define the parameter space
	confMC.addParameter("x0", -10,10,0.001,"$x_0$");
	confMC.addParameter("y0", -10,10,0.001,"$y_0$");

//	parabmc.saveStepPDFs("parabola-stepPDFini");

	// set the observational data
	confMC.setData(data);
	confMC.setDiagonalCovarianceMatrix(zerr*zerr);
	confMC.printInfo();
	confMC.setCalculate2Dposteriors(true);
	gaussMC=confMC;
	
	//
	// RUN THE MCMC CHAIN
	//
	long id;
	omp_set_nested(int(opt["ompnested"].as<bool>()));
	if (opt.count("num_threads"))
		omp_set_num_threads(opt["num_threads"].as<int>());

#pragma omp parallel for private(id)
	for (id = 0; id < opt["Nchains"].as<int>(); ++id) {
		gaussxyMCMC mc(confMC);
		mc.setVerbocity(Medium);
		mc.setID(id);
		mc.msgs->setSender("gaussMCMC."+mc.msgs->toStr(id));
		mc.startChain();
#pragma omp critical
		{
		gaussMC.append(mc);
		}
	}
	gaussMC.saveResults();
	
//	cpeds_fwhm2sigma()
	
	//
	// calculate Fisher information matrix
	//
	calculateFisher(data3d,zerr, outdir, gaussMC.get2Dposterior("x0","y0"), gaussMC.getBestFitLink().getParam(0),gaussMC.getBestFitLink().getParam(1));
	
	//
	// calculate comparizon between 1D calculated posterior and the theoretical predictions
	//
	calculate_sigma_1Dposterior(gaussMC.get1Dposterior(0), gaussMC.get1Dposterior(1));
	
	
	return 0;
}





