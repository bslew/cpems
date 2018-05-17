/*!
 * \file test-mscsfn.cpp
 *
 *  Created on: Sep 13, 2012
 *      Author: blew
 */

#include "stdio.h"
#include "stdlib.h"
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map-window_function.h"
#include "Mscs-function3dregc.h"
#include "mscsFunction3dCLEAN.h"
#include "mscsFunctionFit.h"
#include "cpeds-point_set.h"
#include "cpeds-msgs.h"
#include "cpeds-list.h"
#include "mscsVector.h"
#include "MscsPDF2D.h"
#include "MscsPDF1D.h"
#include "MscsDateTimeData.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#ifdef HAVE_OPENMP
omp_lock_t  _omp_lock;
#endif

void makeLMFitTest();
void goOMP();
void gofxy();

double test51_fn(double x, double y) {
	return x*x+y*y;
}

int main(int argc, char** argv) {
#ifdef HAVE_OPENMP
	omp_init_lock(&_omp_lock);
#endif

	if (argc==1) {
		printf("USAGE: test-mscsfn testNo\n");
		printf("TESTS\n");
		printf("0 - air mass elevation relation test\n");
		printf("1 - mscsfn3dregc cutaway block test\n");
		printf("2 - mscsfn cut from to test\n");
		printf("3 - mscsfn interpolation with different types depending on the arguments separation threshold test\n");
		printf("4 - mscsFunction: averaging within same arguments test\n");
		printf("5 - mscsFunction: make atm. pressure profile\n");
		printf("6 - mscsFunction3dregc: matrix/slice saving convensions test\n");
		printf("7 - mscsFunction: skewgauss test\n");
		printf("8 - makeLMFitTest\n");
		printf("9 - dNdS fitting test on real data from north6cm catalog\n");
		printf("10 - 3d gauss convolution test\n");
		printf("11 - test reaching function arguments in python style by negative indexes\n");
		printf("12 - test integrateX_Xmax\n");
		printf("13 - 3d real-space convolution test \n");
		printf("15 - test 3d function intersections\n");
		printf("16 - test function hdf5 functionality\n");
		printf("17 - test function3dregc concatenate capability\n");
		printf("18 - mscsFunction: finter test for n<5\n");
		printf("19 - test function3dregc-hdf5 IO with openMP\n");
		printf("20 - mscsFunction3dregc: fxy test\n");
		printf("21 - mscsFunction3dregc: gauss fwhm test\n");
		printf("22 - mscsFunction3dregc: matrix orientation for plotting with plot_matrix.py\n");
		printf("23 - mscsFunction: test root finder\n");
		printf("24 - mscsFunction: test log space generation\n");
		printf("25 - mscsFunction3dregc: test RL deconvolution\n");
		printf("26 - mscsFunction3dCLEAN: test CLEAN deconvolution\n");
		printf("27 - mscsFunction: test mkPowerLaw\n");
		printf("28 - mscsFunction: test LMdatafit poly\n");
		printf("29 - mscsFunction: test interpolation (periodic)\n");
		printf("30 - test arithmetic\n");
		printf("31 - mscsFunctionFit: test gauss data fit with grid search\n");
		printf("32 - mscsFunction: use case - average step-like function\n");
		printf("33 - mscsFunction3dregc: test power spectrum calculation (white noise)\n");
		printf("34 - mscsFunction3dregc: test power spectrum calculation (top hat/gauss)\n");
		printf("35 - mscsFunction3dregc: test power spectrum calculation (Pk)\n");
		printf("36 - mscsFunction3dregc: test correlation function calculation\n");
		printf("37 - mscsFunction3dregc: test multi-scale power spectrum calculation\n");
		printf("38 - mscsFunction3dregc: test multi-scale power spectrum calculation2\n");
		printf("39 - mscsFunction: binning function with calculation of statistics in bins\n");
		printf("40 - mscsFunction: testing hdf5 non-closed files issue\n");
		printf("41 - mscsFunction3dregc: test integration accuracy\n");
		printf("42 - mscsFunction3dregc: mk 100 x 50 hdf5 dataset\n");
		printf("43 - mscsFunction: Fourier coefficients\n");
		printf("44 - mscsFunction: Fourier series test\n");
		printf("45 - LINE 1\n");
		printf("46 - test cpeds_point3d rotations\n");
		printf("47 - test cpeds_sort_data - new qsort implementation\n");
		printf("48 - test mscsFunction distributions\n");
		printf("49 - test cpeds list median implementation\n");
		printf("50 - test mscsFunction binFunction fixed interval\n");
		printf("51 - test mscsFunction3dregc triangulation interpolation\n");
		printf("52 - test mscsPDF2D contours\n");
		printf("53 - test mscsPDF1D CR\n");
		printf("54 - test mscsFunction derivative\n");
		printf("55 - test cpeds shift array\n");
		printf("56 - test cpeds bilinear interpolation\n");
		printf("57 - test mscsFunction3dregc hole filling by degrade/prograde\n");
		printf("58 - test MscsDataTimeData\n");
		printf("59 - mscsFunction3dregc: test hdf5 string attribute\n");
		printf("60 - shiftYwrtX\n");
		exit(0);
	}
	long testNo=strtol(argv[1],NULL,10);
	mscsWindowFunction wfn;
	mscsFunction3dregc f3d,g3d,h3d;
	mscsFunction3dCLEAN cleanfn;
	mscsFunction f,f2,f3;
	cpedsPointSet3D ps;
	subDomain_region_t sd;
	long i,j;
	long n,Nx,Ny, Nz;
	double Lx,Ly,Lz;
	mscsFunctionFit fit;
	cpedsList<double> ini,fitparam;
	cpedsList<double> err;
	cpedsList<double> alpha,weight;
	cpedsList<long> binsize;
	vector< cpedsList<double> > stats;
	double fwhm,sigma,nu,T0;
	cpedsMsgs msgs;
	mscsVector<double> roots;

	
	switch (testNo) {
		case 0:
			wfn.mkAirMassElevationRelation();
			wfn.save("airMassElevationRelation.txt");
			wfn.clearFunction();
			wfn.newPoint(1,cpeds_air_mass(1));
			wfn.newPoint(2,cpeds_air_mass(2));
			wfn.save("airMassElevationRelation-cpeds-test.txt");
			
			break;
			
		case 1:
			f3d.setSize(2,3,4,1,2,3,0,0,0);
			f3d.allocFunctionSpace();
			f3d.generateRandomGaussianField(0,1,true,true,1,0);
			f3d.printInfo();
			f3d.savetxtlin("in",true);
			f3d=f3d.cutAwayBlock(0,0,0,1,1,2);
			f3d.printInfo();
			f3d.savetxtlin("cut",true);
			break;
			
		case 2:
			f.mkSin(0,1,0.1,1,0,1);
			f.save("cutfromto-test.txt");
			f2=f.cut(0.4,0.6);
			f2.save("cutfromto-test-check.txt");
			break;
		case 3:
			f.mkSin(0,1,0.1,1,0,1);
			f.cut(0.4,0.6);
			f.save("int-test.txt");
			f2=f.interpolate(0.01,"cspline",0.15,"linear");
			f2.save("int-test-check.txt");
			break;
		case 4:
			f.mkConst(0,100,1,1);
			f2.mkConst(0,100,1,2);
			f3.mkConst(0,100,1,4);
			f.concatenate(f2);
			f.concatenate(f3);
			f.sortFunctionArgAscending();
			f.save("average");
			f.average_sameArgs();
			f.save("average.done");
			break;
		case 5:
			f.mkPressureProfile(0,100,0.1,1013,300,0);
			f.save("atmPressureProfile");
			break;

		case 6:
			f3d.setSize(2,3,1);
			f3d.allocFunctionSpace();
			f3d=double(0);
			f3d.fRe(0,0,0)=0;
			f3d.fRe(1,0,0)=1;
			f3d.fRe(0,1,0)=2;
			f3d.fRe(1,1,0)=3;
			f3d.fRe(0,2,0)=4;
			f3d.fRe(1,2,0)=5;
			f3d.printInfo();
//			f3d.savetxtlin("in",true);
//			f3d=f3d.cutAwayBlock(0,0,0,1,1,2);
//			f3d.printInfo();
			f3d.saveSlice(2,0,"f3d-test.slice",0);
//			f3d.savetxtlin("cut",true);
			break;
		case 7:
			f.mkSkewGauss(-5,5,0.1,1,0,1,0);
			f.save("skewGauss--A1_m0_s1_a0");
			f.clearFunction();
			f.mkSkewGauss(-5,5,0.1,1,0,1,2);
			f.save("skewGauss--A1_m0_s1_a2");
			f.clearFunction();
			f.mkSkewGauss(-5,5,0.1,1,0,1,-2);
			f.save("skewGauss--A1_m0_s1_a-2");
			f.clearFunction();
			f.mkSkewGauss(-5,5,0.1,1,-2,1,-2);
			f.save("skewGauss--A1_m-2_s1_a-2");

			break;
			
		case 8:
			makeLMFitTest();
			
			fit.load("data");
			ini.append(1);
			ini.append(0);
			ini.append(1);
			for (long i = 0; i < fit.pointsCount(); i++) {
				err.append(0.1);
			}
			fit.print();
			printf("now fitting data");
			fit.fitData(ini,err,"exponential_3param");
			break;

		case 9:
			fit.load("/home/blew/cosmo-data/radio-catalogs/rtv3c15ghz/rtv9c15ghz_cross_north6cm--alpha.h");
			fit.checkRanges();
			ini.append(30);
			ini.append(0.22);
			ini.append(0.22);
			ini.append(0);
			for (long i = 0; i < fit.pointsCount(); i++) {
				err.append(5);
			}
			fit.print();
			printf("now fitting data");
			fitparam=fit.fitData(ini,err,"skewgauss_4param");
			fitparam.print();
//			fit.getBestFitModel().save("bestFitModel");
			f.mkSkewGauss(fit.getMinArg(),fit.getMaxArg(),0.01,fitparam[0],fitparam[1],fitparam[2],fitparam[3]);
			f.save("bestFitModel");
			break;
		case 10:
			printf("Begin smoothing test\n");
			printf("Signal is 1 pixel in the center with value 1. All other pixels are 0.\n");
			printf("Smoothing 501x501 pixels field. Field size is 501x501 units\n\n");
			printf("Gaussian smoothing length: fwhm=20 units\n\n");
			f3d.setSizeRange(501,501,1,0,0,0,501,501,0); // 2D function
			f3d.allocFunctionSpace();
			f3d=double(0);
			f3d.fRe(250,250,0)=1;
			f3d.saveSlice(2,0,"signal.L501",0);
			f3d.printInfo();
			printf("integral before convolution: %lE\n",f3d.integrateRe());
			fwhm=20.0/(2.0*sqrt(2.0*log(2.0)));
			f3d.smooth3DGauss(fwhm,fwhm,fwhm);
			printf("integral after convolution: %lE\n",f3d.integrateRe());
			f3d.saveSlice(2,0,"signalConvolved_fwhm20.L501re",0);
			f3d/=f3d.getMaxValue();
			f3d.saveSlice(2,0,"signalConvolved_fwhm20.L501renorm",0);
			f3d.saveSlice(2,0,"signalConvolved_fwhm20.L501imnorm",1);
			printf("Test END\n\n");

			printf("--------------------\n");
			printf("Begin smoothing test\n");
			printf("Signal is 1 pixel in the center with value 1. All other pixels are 0.\n");
			printf("Smoothing 501x501 pixels field. Field size is 1x1 units\n\n");
			printf("Gaussian smoothing length: fwhm=0.1 units\n\n");
			//			exit(0);
			f3d.setSizeRange(501,501,1,0,0,0,1,1,0); // 2D function
			f3d.allocFunctionSpace();
			f3d=double(0);
			f3d.fRe(250,250,0)=1;
			printf("maximal field amplitude before smoothing: %lE\n",f3d.getMaxValue());
			f3d.saveSlice(2,0,"signal.L1",0);
			fwhm=0.1/(2.0*sqrt(2.0*log(2.0)));
			f3d.smooth3DGauss(fwhm,fwhm,fwhm);
			f3d.saveSlice(2,0,"signal.L1.convolved_fwhm0.1.re",0);
			printf("maximal field amplitude after smoothing: %lE\n",f3d.getMaxValue());
			printf("Test END\n\n");

			printf("--------------------\n");
			printf("Begin smoothing test\n");
			printf("Signal is 1 pixel in the center with value 1. All other pixels are 0.\n");
			printf("Smoothing 401x401 pixels field. Field size is 1x1 units\n\n");
			printf("Gaussian smoothing length: fwhm=0.1 units\n\n");
			//			exit(0);
			f3d.setSizeRange(401,401,1,0,0,0,1,1,0); // 2D function
			f3d.allocFunctionSpace();
			f3d=double(0);
			f3d.fRe(200,200,0)=1;
			printf("maximal field amplitude before smoothing: %lE\n",f3d.getMaxValue());
			printf("integral before convolution: %lE\n",f3d.integrateRe());
			f3d.saveSlice(2,0,"signal.L1",0);
			fwhm=0.1/(2.0*sqrt(2.0*log(2.0)));
			f3d.smooth3DGauss(fwhm,fwhm,fwhm);
			f3d.saveSlice(2,0,"signal.L1.convolved_fwhm0.1.re",0);
			printf("maximal field amplitude after smoothing: %lE\n",f3d.getMaxValue());
			printf("integral after convolution: %lE\n",f3d.integrateRe());
			printf("Test END\n\n");
			
			exit(0);
			f3d.setSizeRange(501,501,1,0,0,0,3,5,0);
			f3d.allocFunctionSpace();
			f3d=double(0);
			f3d.fRe(250,250,0)=1;
			f3d.saveSlice(2,0,"signal.L1",0);
			f3d.printInfo();
			f3d.fft(true);
			f3d.fft(false);
			f3d.saveSlice(2,0,"signal.L1.fftcheck",0);
			printf("fftw check fwd-bkw: %lE\n",f3d.getMaxValue());
			printf("integral before convolution: %lE\n",f3d.integrateRe());
			fwhm=0.1/(2.0*sqrt(2.0*log(2.0)));
			f3d.smooth3DGauss(fwhm,fwhm,fwhm);
			printf("integral after convolution: %lE\n",f3d.integrateRe());
			f3d.saveSlice(2,0,"signalConvolved_fwhm20.L1re",0);

			break;

		case 11:
			f.newPoint(1,1);
			f.newPoint(2,2);
			f.newPoint(3,3);
			f.newPoint(4,4);
			f.newPoint(5,5);
			f.print(true);
			printf("f(-1) = %lE\n",f.fpsty(-1));
			printf("f(-2) = %lE\n",f.fpsty(-2));
			printf("f(-5) = %lE\n",f.fpsty(-5));
			printf("f(0) = %lE\n",f.fpsty(0));
			printf("f(-6) = %lE\n",f.fpsty(-6));
			printf("f(6) = %lE\n",f.fpsty(6));

			break;

		case 12:
			f.newPoint(1,1);
			f.newPoint(2,2);
			f.newPoint(3,3);
			f.newPoint(4,4);
			f.newPoint(5,5);
			f.print(true);
			f2=f.integrateX_Xmax();
			f2.print(true);
			f.save("test_XXmax");
			f2.save("test_intXXmax");

			break;

		case 13:
			n=30;
			f3d.setSizeRange(n,n,n,0,0,0,2,2,2);
			f3d.allocFunctionSpace();
			f3d.mkGauss3D(0.1,0.1,0.1);
			f3d.saveHDF5("convtest.hdf5","in");
			g3d=f3d;
			g3d.normalizeRe();
//			exit(0);
			f3d=f3d.convolve(g3d);
			f3d.saveHDF5("convtest.hdf5","out");

			break;
//		case 14:
//			ps.append(cpedsPoint3D(0.1,0.1,0.1));
//			n=30;
//			f3d.setSizeRange(n,n,n,0,0,0,2,2,2);
//			f3d.allocFunctionSpace();
//			f3d.mkGauss3D(0.1,0.1,0.1);
//			f3d.saveHDF5("convtest.hdf5","in");
//			g3d=f3d;
//			g3d.normalizeRe();
////			exit(0);
//			f3d=f3d.convolve(g3d);
//			f3d.saveHDF5("convtest.hdf5","out");

			break;
			
		case 15:
			n=5;
			f3d.setSizeRange(n,n,1,0,0,0,5,5,1);
			f3d.allocFunctionSpace();

			g3d.setSizeRange(n,n,1,2,2,0,7,7,1);
			g3d.allocFunctionSpace();

			f3d.mkGauss3D(1,1,5);
			g3d.mkGauss3D(1,1,5);

			f3d.saveHDF5("intersectiong-test.hdf5","f");
			g3d.saveHDF5("intersectiong-test.hdf5","g");
			
			if (f3d.intersects(g3d)) printf("intersects: true\n"); else printf("intersects: false\n");

			f3d.intersection_set(g3d,-1);
			f3d.saveHDF5("intersectiong-test.hdf5","fout");
			g3d.saveHDF5("intersectiong-test.hdf5","gout");

			break;
			
		case 16:
			f.mkGauss(-5,5,1,1,0,1);
			f.print();
			f.saveHDF5("fn.hdf5","fout","[m]","[cm]","testx", "testy");
			
			f2.loadHDF5("fn.hdf5","fout");
			f2.print();

			break;

		case 17:
			n=5;
			f3d.setSize(2,n,1,1,1,1,0,0,0);
			f3d.allocFunctionSpace();

			n-=2;
			g3d.setSize(2,n,1,1,1,1,0,0,0);
			g3d.allocFunctionSpace();

			f3d.mkGauss3D(1,1,1);
			f3d.printInfo();
			g3d=1.0;
			for (i = 0; i < n; ++i) {
				g3d.setf(0,i,0,i,0);
				g3d.setf(1,i,0,i*2,0);
			}

			h3d=f3d.concatenate(g3d,1);
			
			h3d.saveHDF5("concatenate-test.hdf5","joint");
			

			break;
			
		case 18:
			f.newPoint(1,1);
			f.newPoint(4,4);
			f.print(true);
			printf("x: %lf, finter: %lf\n",3.0,f.finter(3.0));
			printf("x: %lf, finter: %lf\n",0.5,f.finter(0.5));
			printf("x: %lf, finter: %lf\n",4.5,f.finter(4.5));
			
			f.clearFunction();
			
			f.newPoint(1,1);
			printf("single point interpolation\n");
			f.print(true);
			printf("x: %lf, finter: %lf\n",4.5,f.finter(4.5));


		case 19:
			goOMP();
			break;
		
		case 20:
			gofxy();
			break;

		case 21:
			n=100;
			f3d.setSizeRange(n,n,1,0,0,0,n,n,0);
			f3d.allocFunctionSpace();

			fwhm=10; //deg
			sigma=fwhm/(2.0*sqrt(2.0*log(2.0)));
//			printf("")
			f3d.mkGauss2D(n/2,n/2,sigma,sigma,1.0,2,0);

			f3d.saveHDF5("gauss_fwhm_test.hdf5","fwhm_"+msgs.toStrf(fwhm,5));

			break;

		case 22:
			Nx=20;
			Ny=40;
			f3d.setSizeRange(Nx,Ny,1,0,0,0,Nx,Ny,0);
			f3d.allocFunctionSpace();

			for (i = 0; i < Nx; i++) {
				for (j = 0; j < Ny; j++) {
					f3d.fRe(i,j,0)=i*Ny+j;
				}
			}

			f3d.saveHDF5("matrix.hdf5","mat");

			break;

		case 23:
			f.mkGauss(-5,5,0.01,1,0,2);
			f.save("gauss");
			f-=0.5;
			roots=f.findRoot();
			roots.printVector();
			break;

		case 24:
			f.mkLogSpace(1,1000,10,10);
			f.save("logSpace");
			break;

		case 25:
			Nx=100;
			Ny=100;
			f3d.setSizeRange(Nx,Ny,1,0,0,0,Nx,Ny,0);
			f3d.allocFunctionSpace();
			f3d=0.0;
			f3d.setf(Nx/2,Ny/2,0,1,0);
			f3d.smooth3DGauss(2.5,5,1);
			g3d=f3d;
			g3d=0.0;
			g3d.mkGauss2D(Nx/2,Ny/2,2.5,5,1,2,0);
			
			f3d=f3d.RLdeconvolution(f3d,g3d,20);

			f3d.saveHDF5("RLestimate.hdf5","estimate");

			break;
			
		case 26:
			Nx=100;
			Ny=100;
			f3d.setSizeRange(Nx,Ny,1,0,0,0,Nx,Ny,0);
			f3d.allocFunctionSpace();
			f3d.mkGauss2D(Nx/2,Ny/2,5,5,1,2,0);
//			f3d=0.0;
//			f3d.setf(Nx/2,Ny/2,0,1,0);
//			f3d.smooth3DGauss(5,5,1);
			g3d=f3d;
			g3d*=-1.0;
			f3d.shift(-10,0,0);
			g3d.shift(10,0,0);
			g3d+=f3d;
			
			g3d.saveSlice(2,0,"beam",0);
			g3d=g3d.rotateSlice(2,0,PIsnd/2,Nx/2,Ny/2,0);
			g3d.saveSlice(2,0,"rot",0);
			
			cleanfn.setVerbosityLevel(Top);
			cleanfn.initiate(g3d,g3d,5000,0.2);			
			cleanfn.clean();

		case 27:
			f.mkConst(0,36,1,0);
			f.mkPowerLaw(0,36,1,0.3,0,1,2,1);
			f.save("powerLaw");
			break;

		case 28:
			f.mkConst(-2,2,0.1,1);
			f.save("const.dat");
			fitparam.append(1); 
			fitparam.append(-1);
			fitparam.append(0.5); 
			fitparam.append(2.0);
			f.mkPolynomial(fitparam);
			f.save("poly.dat");
			f.clearFunction();
			f.mkPolynomial(-2,2,1,fitparam);
			f.save("poly2.dat");

			break;

		case 29:
			f.mkLine(0,1,0.1,1,1);
			f.checkRanges();
//			f.cut(0.4,0.6);
			f.newPoint(f.getMaxArg()+0.1,f.f(0));
			f.save("interpolation.in");
			f2=f.interpolate(0.001,false,"akima_periodic");
			f2.save("interpolation.out");

			break;
			
		case 30:
			printf("Specific intensity of 2.726 K black body at 30 GHz should be: 5.72156423-19 J/s/m^2/Hz/sr\n");
			nu=30.0e9;
			T0=2.726;
			printf("We have: %lE\n",2.0*CPEDS_h*nu*nu*nu/CPEDS_c/CPEDS_c/(exp(CPEDS_h*nu/CPEDS_kB/T0)-1.0));
//			printf("We have: %lE\n",2.0*CPEDS_h*pow(nu,3)/pow(CPEDS_c,2)/(exp(CPEDS_h*nu/CPEDS_kB/T0)-1.0));
//			printf("We have (2.0*CPEDS_h*nu*nu*nu/CPEDS_c): %lE\n",2.0*CPEDS_h*nu*nu*nu/CPEDS_c);
//			printf("We have (2.0*CPEDS_h*nu*nu*nu/300 000 000): %lE\n",2.0*CPEDS_h*nu*nu*nu/300000000.0);
			
			break;

//		case 31:
//			fit.load("/home/blew/tmp/1/data");
//			fit.checkRanges();
//			ini.append(1);
//			ini.append(1);
//			ini.append(1);
//			ini.append(1);
//			for (long i = 0; i < fit.pointsCount(); i++) {
//				err.append(0.1);
//			}
//			fit.print();
//			printf("now fitting data\n");
//			fitparam=fit.fitData(ini,err,"gauss_4param");
//			fitparam.print();
//			printf("performing grid search\n");
//			mscsFunctionFit_radialGridSearch(fit,err,ini,20,1000);
////			fit.getBestFitModel().save("bestFitModel");
//			f.mkGauss(fit.getMinArg(),fit.getMaxArg(),0.01,fit.params()[0],fit.params()[1],fit.params()[2],fit.params()[3]);
//			fit.params().print();
//			f.save("/home/blew/tmp/1/bestFitModel");
//
//
//			break;

		case 32:
			f.load(argv[2],true,13,9);
			f.sortFunctionArgAscending();
			f.unique();
			f.add(-f.finter(0,"linear"));
			f.save("test-mscsfn-usecase32");
			f.invert();
			f.average_sameArgs();
			f.invert();
			f.newPoint(0,0);
			f.sortFunctionArgAscending();
//			f2=f.interpolate(1.0,0,"linear");
			f2=f.extrapolate(0,85,1,false,"linear"); //interpolate(1.0,0,"linear");
			
			f.save("test-mscsfn-usecase32.avg");
			f2.save("test-mscsfn-usecase32.int");
			break;

		case 33:
			Nx=100;
			Ny=100;
			Nz=100;
			Lx=1000; // 1000 Mpc
			Ly=1000; // dx=1000/100=10 Mpc
			Lz=1000;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.generateRandomGaussianField(0,1,true,false,1,0);
			g3d=f3d*f3d;
			f=f3d.powerSpectrum(1);
			printf("sum fi^2: %lf\n",g3d.sumRe());
//			printf("V Sum Pk: %lf\n",g3d.sumRe());
			
			printf("kmin: %lf\n",1.0/Lx);
			printf("kmax: %lf\n",1.0/(2.0 * Lx/Nx));
			
//			g3d.saveSlice(2,0,"beam",0);
			f.save("Pk");
			break;

		case 34:
			Nx=500;
			Ny=500;
			Nz=1;
			Lx=1000; // 1000 Mpc
			Ly=1000; // dx=1000/100=10 Mpc
			Lz=0;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
//			f3d.mkBall3D(100);
			f3d.mkGauss2D(Lx/2,Ly/2,10,10,1,2,0);
			f3d.shift(Nx/2,Ny/2,0);
			f3d.saveSlice(2,0,"real");
			f3d.fft(true);
			f3d.saveSlice(2,0,"fourierRe",0);
			f3d.saveSlice(2,0,"fourierIm",1);
			f3d.shift(Nx/2,Ny/2,0);
			f3d.saveSlice(2,0,"fourierReShift",0);
			f3d.saveSlice(2,0,"fourierImShift",1);
			f3d.shift(Nx/2,Ny/2,0);
			f3d.absoluteValue();
			f3d*=f3d;
			f3d.saveSlice(2,0,"power",0);
			f3d.shift(Nx/2,Ny/2,0);
			f3d.saveSlice(2,0,"powerShift",0);
			break;
		case 35:
			Nx=400;
			Ny=400;
			Nz=400;
			Lx=700; // 100 Mpc
			Ly=700; // dx=1000/100=10 Mpc
			Lz=700;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
			f3d.setRe(0);
//			f.mkPowerLaw(0.001,10,0.001,2.41e-9,0.002,0.961-1.0);
//			f.mkPowerLaw(0.001,10,0.001,1,0.01,0);
			f.mkPowerLaw(0.001,10,0.001,1,0.01,-1);
			f.save("inPk");
			f3d.generateRandomGaussianField(f,0);
//			f3d.generateRandomGaussianField(0,1,true,false,0,0);
			f3d.saveSlice(2,0,"grf",0);
			f3d.saveSlice(2,0,"grfim",1);
			g3d=f3d;
			g3d*=f3d;
			printf("sum fi^2: %lf\n",g3d.sumRe());

			f=f3d.powerSpectrum(1);
			
			f.save("outPk");
			break;

		case 36:
			Nx=50;
			Ny=50;
			Nz=50;
			Lx=100; // 100 Mpc
			Ly=100; // dx=100/100=1 Mpc
			Lz=100;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
//			f.mkPowerLaw(0.001,10,0.001,2.41e-9,0.002,0.961-1.0);
//			f.mkPowerLaw(0.001,10,0.001,1,0.01,-2);
//			f.save("inPk");
//			f3d.generateRandomGaussianField(f);
			f3d.generateRandomGaussianField(0,1,true,false,0,0);
//			g3d=f3d*f3d;
			f=f3d.correlationFunctionR();
			
			f.save("corr");
			break;

		case 37:
			Nx=400;
			Ny=400;
			Nz=400;
			Lx=3000; // 100 Mpc
			Ly=3000; // dx=1000/100=10 Mpc
			Lz=3000;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
			f3d.setRe(0);
//			f.mkPowerLaw(0.001,10,0.001,2.41e-9,0.002,0.961-1.0);
//			f.mkPowerLaw(0.001,10,0.001,1,0.01,0);
			f.mkPowerLaw(0.001,10,0.001,1,0.01,-1);
			f.save("inPk");
//			f3d.generateRandomGaussianField(f,0);
			f3d.generateRandomGaussianField(0,1,true,false,0,0);
//			f3d.saveSlice(2,0,"grf",0);
//			f3d.saveSlice(2,0,"grfim",1);
//			g3d=f3d;
//			g3d*=f3d;
//			printf("sum fi^2: %lf\n",g3d.sumRe());
			printf("rms fi: %lf\n",f3d.RMSre());

//			mscsVector< mscsFunction > msPk;
			f=f3d.powerSpectrum(1);
			g3d=f3d.cutAwayBlock(0,0,0,99,99,99);
			f2=g3d.powerSpectrum(1);
			
			f.save("outPk");
			f2.save("outPksub");
			break;
		case 38:
			sd.xmin=0;
			sd.ymin=0;
			sd.zmin=0;
			sd.xmax=3000;
			sd.ymax=3000;
			sd.zmax=3000;
			ps.generateRandomPoints(10000000,sd);

			Nx=200;
			Ny=200;
			Nz=200;
			Lx=3000; // 100 Mpc
			Ly=3000; // dx=1000/100=10 Mpc
			Lz=3000;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();

			// calculate density field 
			f3d.populateField(ps,false,0,0,false,true);
			for (i = 0; i < f3d.getN(); i++) {
				f3d.fRe(i)*=f3d.fIm(i)/f3d.getPixelVolume(0,0,0,false);
				f3d.fRe(i)*=f3d.fIm(i);
			}
			printf("dV: %lf\n",f3d.getPixelVolume(0,0,0,false));
			f=f3d.powerSpectrum(1);			
			f.save("outPk2");

			f3d.freeFunctionSpace();
			Nx=100;
			Ny=100;
			Nz=100;
			Lx=300; // 100 Mpc
			Ly=300; // dx=1000/100=10 Mpc
			Lz=300;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();

			// calculate density field 
			f3d.populateField(ps,false,0,0,false,true);
			for (i = 0; i < f3d.getN(); i++) {
				f3d.fRe(i)*=f3d.fIm(i)/f3d.getPixelVolume(0,0,0,false);
				f3d.fRe(i)*=f3d.fIm(i);
			}
			printf("dV: %lf\n",f3d.getPixelVolume(0,0,0,false));
			f=f3d.powerSpectrum(1);			
			f.save("outPksub");

			break;
		case 39:
			f.mkGaussianNoise(10000,0,1,0);
			f.save("noise");
			binsize << 1000 << 1000 << 1000 << 1000 << 1000 << 1000 << 1000 << 1000 << 1000 << 1000 ;
//			weight.makeLength(10000);
//			weight=double(1.0);
//			alpha << 0.0026998 << 0.0455003 << 0.317311 << 0.5 << 0.682689 << 0.9545 << 0.9973;
			alpha << 0.0455003/2 << 0.317311/2 << 0.5 << 1.0-0.317311/2 << 1.0-0.0455003/2 ;
//			stats.push_back(alpha);
			f.binFunctionLin2(0,binsize,weight,&err,&alpha, &stats);
			f.save("noise.bin");
			stats[0].save("2sl");
			stats[1].save("1sl");
			stats[2].save("center");
			stats[3].save("1su");
			stats[4].save("2su");
			stats[5].save("mx");
			stats[6].save("my");
			stats[7].save("sx");
			stats[8].save("sy");
			stats[9].save("minx");
			stats[10].save("maxx");
			stats[11].save("miny");
			stats[12].save("maxy");
			break;
			
		case 40:
			f.mkGaussianNoise(10000,0,1,0);
			f.saveHDF5("noise.hdf5","/dupa/test");
			while (1) {
				sleep(1);
			}

		case 41:
		{
			Nx=100;
			Ny=100;
			Nz=100;
			Lx=1000; // 1000 Mpc
			Ly=1000; // dx=1000/100=10 Mpc
			Lz=1000;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
			f3d.setRe(0);
			double sigma=cpeds_fwhm2sigma(100);
			f3d.mkGauss3D(sigma,sigma,sigma);
			g3d=f3d.integrateRe(2);
			printf("integral at max: %lf\n",g3d.fRe(Nx/2,Ny/2,0));
			printf("theoretical at max: %lf\n",sqrt(twoPI)*sigma);

			f3d.saveHDF5("gauss3d","in");
			g3d.saveSlice(2,0,"gauss3d.int");
			f3d=g3d;
			f3d.mkGauss3D(sigma,sigma,0);
			f3d*=sqrt(2.0*PI)*sigma;
			f3d.saveSlice(2,0,"model");
			g3d-=f3d;
			g3d.saveSlice(2,0,"resid");
		}	
			break;

		case 42:
		{
			Nx=100;
			Ny=50;
			Nz=1;
			Lx=8; // 1000 Mpc
			Ly=4; // dx=1000/100=10 Mpc
			Lz=0;
			f3d.setSizeRange(Nx,Ny,Nz,0,0,0,Lx,Ly,Lz);
			f3d.allocFunctionSpace();
			f3d.setIm(0);
			f3d.setRe(0);
			double sigma=cpeds_fwhm2sigma(3);
			f3d.mkGauss3D(sigma*2,sigma,0);

			f3d.saveHDF5("data.hdf5","gauss");
//			g3d.saveSlice(2,0,"gauss3d.int");
//			f3d=g3d;
//			f3d.mkGauss3D(sigma,sigma,0);
//			f3d*=sqrt(2.0*PI)*sigma;
//			f3d.saveSlice(2,0,"model");
//			g3d-=f3d;
//			g3d.saveSlice(2,0,"resid");
		}	
			break;

			
		case 43:
		{
			f.mkSin(0,4,0.001,1,0.6,1);
//			f2.mkGaussianNoise(f.pointsCount(),0,0.1);
//			f+=f2;
			f.save("sig");
			mscsFunction re,im;
			f2=f.powerSpectrum(&re,&im);
			f2.save("power");
			re.save("re");
			im.save("im");
			double RE=re.f(4);
			double IM=im.f(4);
			double phi=atan2(RE,IM);
			printf("phi= %lf f0 = %lf sig(0)=%lf\n",phi, 2.0*sqrt(RE*RE+IM*IM)*sin(phi),f.f(0));
			
			break;
		}

		case 44:
		{
			f.mkConst(0,1,0.01,0);
			f.mkConst(1,2,0.01,1);
			f.mkConst(2,3,0.01,0);
			f.mkConst(3,4,0.01,1);
			f.mkConst(4,5,0.01,0);
			f.mkConst(5,6,0.01,1);
			f.shiftYwrtX(20);
			f.shiftX(0.3);
			f+=2.0;
			f.save("sig");
			mscsFunction A,phi;
			f2=f.FourierSeries(500,0, &A,&phi);
			f2.save("series");
			A.save("A");
			phi.save("phi");
			
			break;
		}

		case 45:
		{
			long double a0=0;
			long double a1=1;
			long double an;
			long n=8181;
			for (long i = 0; i <=n; i++) {
				an=a0+a1;
				printf("a%li: %.15LE\n",i,an);
				a0=a1;
				a1=an;
			}
				
			break;
		}
		case 46:
		{
			cpedsPoint3D p(15,-15,0);
			p.print_point();
			p.Rz(11.8*PI180);
			p.print_point();
				
			break;
		}

		case 47:
		{
			f.mkGaussianNoise(1000,0,1,0);
			cpedsList<double> l=f.toYList();
			l.save("unsorted");
			l.sort(12);
			l.save("sorted12");
			l=f.toYList();
			l.sort(21);
			l.save("sorted21");
			
				
			break;
		}

		case 48:
		{
			f.mkGumbelDistr(-5,20,0.1,2,0.5);
			f.save("Gumbel-0.5-2.0.txt");
			
			break;
		}
		case 49:
		{
			cpedsList<double> l;
			l.append(1);
			l.append(2);
			l.append(3);
			l.print();
			cout << "Median: " << l.median() << "\n";

			l.clear();
			l.append(1);
			l.append(2);
			l.print();
			cout << "Median: " << l.median() << "\n";
			
				
			break;
		}

		case 50:
		{
			f.newPoint(1,1);
			f.newPoint(2,2);
			f.newPoint(3,3);
			f.newPoint(15,5);
			f.newPoint(15.1,6);
			f.save("data");
			cpedsList<long> bins;
			f2=f.binFunction(2,bins,"mean","mean");
			f2.save("data-meanmean");
			bins.save("meanmean.bins");
			bins.clear();
			f2=f.binFunction(2,bins,"mean","max");
			f2.save("data-meanmax");
			f2=f.binFunction(2,bins,"bin_center","max");
			f2.save("data-bincentermax");
			break;
		}
		case 51:
		{
			mscsVector<cpedsPoint3D> positions;
			cpedsPoint3D p;
			cpedsPointSet3D ps;
			for (long x = 0; x < 4; x++) {
				for (long y = 0; y < 4; y++) {
					positions.push_back(cpedsPoint3D(x,y,test51_fn(x,y)));
					ps.append(cpedsPoint3D(x,y,test51_fn(x,y)));
					
				}
			}
//			positions.push_back(cpedsPoint3D(x,y,test51_fn(x,y)));
			ps.append(cpedsPoint3D(1,1,test51_fn(2,2)));
			ps.append(cpedsPoint3D(0.0001,0.00001,test51_fn(1,1)));
			ps.save("points");

			f3d.setSizeRange(100,100,1,0,0,0,4,4,0);
			f3d.allocFunctionSpace();
			f3d.setRe(0);
//			f3d.mkInterpolatedFieldTriangSibsonGrad2D(positions);
//			f3d.mkInterpolatedFieldTriangLinear2D(positions);
			f3d.mkInterpolatedFieldTriangLinear2D(ps);
			f3d.saveHDF5("triang-linear-interpolation.hdf5","1");
//			f3d.saveHDF5("triang-sibson-interpolation.hdf5","2");
			break;
		}
		case 52:
		{
			n=500;
			f3d.setSizeRange(n,n,1,-5,-5,0,5,5,0);
			f3d.allocFunctionSpace();

			fwhm=2; //deg
			sigma=fwhm/(2.0*sqrt(2.0*log(2.0)));
			printf("fwhm: %lf\n",fwhm);
			printf("sigma: %lf\n",sigma);
			f3d.mkGauss2D(0,0,sigma,sigma,1.0,2,0);
//			f3d+=1.0;
			MscsPDF2D pdf=f3d;
			pdf.saveSlice(2,0,"gauss");
			pdf.getContour(0.68).save("CR68");

			break;
		}
		case 53:
		{
			n=500;
			fwhm=2; //deg
			sigma=fwhm/(2.0*sqrt(2.0*log(2.0)));
			printf("fwhm: %lf\n",fwhm);
			printf("sigma: %lf\n",sigma);
			f.mkGauss(-5,5,0.01,1,0,sigma);
//			f3d+=1.0;
			MscsPDF1D pdf=f;
			pdf.save("pdf_good");
			cout << "TESTING WELL DEFINED PDF\n";
			cout << "calcuating 1-sigma CR\n";
			cpedsList<double> v=pdf.getCR(0.682689);
			v.print();
			cout << "\n\n";
			cout << "calcuating 2-sigma CR\n";
			v=pdf.getCR(0.9545);
			v.print();
			cout << "\n\n";
			
			cout << "TESTING ILL DEFINED PDF 1\n";
			f.clearFunction();
			f.mkGauss(-1,5,0.01,1,0,sigma);
			pdf=f;
			pdf.save("pdf_bad1");
			v=pdf.getCR(0.9545);
			v.print();
			cout << "\n\n";

			cout << "TESTING ILL DEFINED PDF 2\n";
			f.clearFunction();
			f.mkGauss(-5,1,0.01,1,0,sigma);
			pdf=f;
			pdf.save("pdf_bad2");
			v=pdf.getCR(0.9545);
			v.print();
			cout << "\n\n";

			break;
		}
		case 54:
		{
			double T=1;
			f.mkSin(0,0.999,0.01,T,0,1);
			f.save("sin");
			long i=0;
//			printf("%li\n",(i-1) % 10);
//			exit(0);
			f.derivative(true,&T);
			f.save("sin-deriv");
			break;
		}
		case 55:
		{
			f.mkSin(0,100,1,1,0,1);
			double* d=f.extractArguments();
			long n=f.pointsCount();
			f.save("in");
			cpeds_save_matrix(d,n,1,"in0");
			cpeds_shift_array(d,n,1,true);
			cpeds_save_matrix(d,n,1,"in1");
			cpeds_shift_array(d,n,11,true);
			cpeds_save_matrix(d,n,1,"in12");
			break;
		}
		case 56:
		{
			double x1=1;
			double x2=2;
			double y1=1;
			double y2=2;
			double f11=1;
			double f21=2;
			double f12=3;
			double f22=4;
			cout << "(x,y)_1 = (" << x1 << ", " << y1 << "), f = "<< f11 << "\n";
			cout << "(x,y)_2 = (" << x2 << ", " << y1 << "), f = "<< f21 << "\n";
			cout << "(x,y)_3 = (" << x1 << ", " << y2 << "), f = "<< f12 << "\n";
			cout << "(x,y)_4 = (" << x2 << ", " << y2 << "), f = "<< f22 << "\n";
			double dx=x2-x1;
			double dy=y2-y1;
			double x=1.5;
			double y=1.5;
			cout << "interpolating at: (x,y) = (" << x <<", " << y << ")\n";				
			cout << "f(x,y) = " << cpeds_bilinear_interpolation(x1,x2,y1,y2,f11,f12,f21,f22,x,y) <<"\n";
			
			x=1;
			y=1;
			cout << "interpolating at: (x,y) = (" << x <<", " << y << ")\n";				
			cout << "f(x,y) = " << cpeds_bilinear_interpolation(x1,x2,y1,y2,f11,f12,f21,f22,x,y) <<"\n";

			x=2;
			y=2;
			cout << "interpolating at: (x,y) = (" << x <<", " << y << ")\n";				
			cout << "f(x,y) = " << cpeds_bilinear_interpolation(x1,x2,y1,y2,f11,f12,f21,f22,x,y) <<"\n";

			cout << "\n\n";
			
			mscsFunction3dregc F;
			mscsFunction3dregc Fint;
			long N=10;
			long Nint=200;
			F.setSizeRange(N,N,1,-5,-5,0,5,5,0);
			Fint.setSizeRange(Nint,Nint,1,-5,-5,0,5,5,0);
			F.allocFunctionSpace();
			Fint.allocFunctionSpace();
			
			F.mkGauss2D(0,0,1,1,0,1,2,0);
//			F.shift(N/2,N/2,0);
			
//			long i=5;
//			long j=5;
//			x=Fint.getX(i);
//			y=Fint.getY(j);
//			printf("x: %lE, y: %lE\n",x,y);
//			printf("i: %li j:%li  F: %lE,  Fint: %lE\n",i,j,F.fRe(i,j,0),F.fxyLin(x,y,0,0));
//			printf("\n\n");
//			exit(0);
			
			for (long i = 0; i < Nint; i++) {
				x=Fint.getX(i);
				for (long j = 0; j < Nint; j++) {
					y=Fint.getY(j);
					printf("x: %lE, y: %lE\n",x,y);
					Fint.setf(i,j,0,F.fxyLin(x,y,0,0),0.0);
					printf("\n\n");
				}
			}
			
	
			F.saveSlice(2,0,"F");
			Fint.saveSlice(2,0,"Fint");
			
			
			
			break;
		}
		case 57:
		{
			
			mscsFunction3dregc F;
			mscsFunction3dregc Fint;
			long N=20;
			long Nint=100;
			F.setSizeRange(N,N,1,-5,-5,0,5,5,0);
			Fint.setSizeRange(Nint,Nint,1,-5,-5,0,5,5,0);
			F.allocFunctionSpace();
			Fint.allocFunctionSpace();

			// make signal
			F.mkGauss2D(0,0,1,1,0,1,2,0);
			// make hole (this is not needed really because the mask does the thing)
			F.fRe(4,4,0)=0.0;
			
			// make mask
			F.setIm(1);
			F.fIm(4,4,0)=0.0;
			
			
			
			
			
			F.saveSlice(2,0,"F");
			Fint.saveSlice(2,0,"Fint");

			break;
		}			
		case 58: {
			cout << "line params count: " << argc << "\n";
			cout << "param 0: " << argv[0] << "\n";
			if (argc<3) { cout << "provide file name as argument\n"; return -1; }
			cout << "param 2: " << argv[2] << "\n";
			MscsDateTimeData dtd;
//			dtd.load(argv[2],0,"%Y-%m-%d %H:%M:%S");
			dtd.load(argv[2],0,"yyyy-MM-dd HH:mm:ss");
			dtd.printData();
			
			long idx;
			double val=dtd.getFirstValueAfter("2016-11-27 19:05:48",0,2,&idx);
			cout << "found at idx: " << idx << " value: " << val << "\n";

			val=dtd.getFirstValueAfter("2016-11-28 19:05:48",0,2,&idx);
			cout << "found at idx: " << idx << " value: " << val << "\n";

			double JD=cpeds_julian_time(2016,11,27,19.0+6.0/60+22.1/3600);
			val=dtd.getFirstValueAfter(JD,0,2,&idx);
			cout << "\ntesting from JD.\nFound at idx: " << idx << " value: " << val;
			
			break;
		}
		case 59: {
			mscsFunction3dregc F;
			F.setSize(10,10);
			F.allocFunctionSpace();
			F.saveHDF5("test.hdf5","test");
			F.setHDF5_scalarStringAttribute("test.hdf5","test","param","$\\alpha$");
			
			break;
		}
		case 60:
			f.newPoint(1,1);
			f.newPoint(2,2);
			f.newPoint(3,3);
			f.newPoint(4,4);
			f.save("inf");
			f.shiftYwrtX(1);
			f.save("inf.shift1");
			f.shiftYwrtX(-1);
			f.save("inf.shift1.shift-1");

			break;
		default:
			break;
	}
	
#ifdef HAVE_OPENMP
		omp_destroy_lock(&_omp_lock);
#endif
	return 0;
}

/***************************************************************************************/
struct data {
		size_t n;
		double * y;
		double * sigma;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data *)data)->n;
	double *y = ((struct data *)data)->y;
	double *sigma = ((struct data *) data)->sigma;
	
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

int expb_df (const gsl_vector * x, void *data, 	gsl_matrix * J) {
	size_t n = ((struct data *)data)->n;
	double *sigma = ((struct data *) data)->sigma;
	
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

int
expb_fdf (const gsl_vector * x, void *data,
		gsl_vector * f, gsl_matrix * J)
{
	expb_f (x, data, f);
	expb_df (x, data, J);
	
	return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multifit_fdfsolver * s) {
	printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2), 
			gsl_blas_dnrm2 (s->f));
}

/***************************************************************************************/
void makeLMFitTest() {
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t N = 40;
	const size_t n = N;
	const size_t p = 3;
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double y[N], sigma[N];
	struct data d = { n, y, sigma};
	gsl_multifit_function_fdf f;
	double x_init[3] = { 1.0, 0.0, 0.0 };
	gsl_vector_view x = gsl_vector_view_array (x_init, p);
	const gsl_rng_type * type;
	gsl_rng * r;
	
    gsl_matrix*  J = gsl_matrix_alloc(n, p);

	gsl_rng_env_setup();
	
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	
	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;
	
	/* This is the data to be fitted */
	
	for (i = 0; i < n; i++)
	{
		double t = i;
		y[i] = 1.0 + 5 * exp (-0.1 * t) 
		+ gsl_ran_gaussian (r, 0.1);
		sigma[i] = 0.1;
		printf ("data: %u %g %g\n", i, y[i], sigma[i]);
	};
	
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	
	print_state (iter, s);
	
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		
		printf ("status = %s\n", gsl_strerror (status));
		
		print_state (iter, s);
		
		if (status)
			break;
		
		status = gsl_multifit_test_delta (s->dx, s->x,
				1e-4, 1e-4);
	}
	while (status == GSL_CONTINUE && iter < 500);
	
//	gsl_multifit_covar (s->J, 0.0, covar);
    gsl_multifit_fdfsolver_jac (s, J);
	gsl_multifit_covar(J,0.0, covar);
	
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	
	{ 
		double chi = gsl_blas_dnrm2(s->f);
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
		
		printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
		
		printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
		printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
	}
	
	printf ("status = %s\n", gsl_strerror (status));
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	gsl_rng_free (r);

}
/***************************************************************************************/
void goOMP() {
	mscsFunction3dregc f3d,h;
	cpedsMsgs msgs;
	long n=100;
	f3d.setSize(n,n,n,1,1,1,0,0,0);
	f3d.allocFunctionSpace();
	long Ndone=0;
	long M=1000;
	long i;
	long Nx=100;
	h.setSize(Nx,M,1);
	h.allocFunctionSpace();
	h=0.0;
	
	omp_set_nested(0);
	omp_set_dynamic(0);
	omp_set_num_threads(100);
	msgs.say("entering parallel for",High);
#pragma omp parallel for schedule(guided) private(i) firstprivate(Nx,f3d,msgs) shared(h,Ndone,M)
	for (i = 0; i < M; i++) {
		mscsFunction3dregc f;
		mscsFunction g;
		string name;
		name="halo_"+msgs.toStr(i);
		for (long j = 0; j < 10; j++) {
			f3d=0.0;
			f3d.mkGauss3D(1,1,1);
			f=f3d;
			g.mkSin(0,1,0.01,1,1,1);
			
		}
		subDomain_region_t r,t;
		r.xmin=0;
		r.xmax=1;
		r.ymin=0;
		r.ymax=1;
		r.zmin=0;
		r.zmax=1;
		r.subx=50;
		r.suby=50;
		r.subz=50;
		t=r;
		t.subx=2;
		t.suby=2;
		t.subz=2;
		
		mscsVector<cpedsPoint3D> pos;
		mscsVector<double> val;
		pos.clear();
		val.clear();
		cpedsPoint3D p;
		cpedsRNG rns("uniform");
		rns.setMinMax(0,1);
		for (unsigned long k = 0; k < 20; k++) {
			p.set(rns.getRN(),rns.getRN(),rns.getRN());
			pos.push_back(p);
			val.push_back(rns.getRN());
		}
		
		f.mkInterpolatedFieldScatter(r,t,pos,val,"gadget2",10,12);
//		omp_set_lock(&_omp_lock);
#pragma omp critical (A)
		{
			f.saveHDF5("halos.hdf",name);
			g.saveHDF5("profile.hdf",name);
			Ndone++;
			msgs.say("saved halo %i/%i",Ndone,M,High);
		}

		if (Ndone % 10 == 0) {
#pragma omp critical (B)
			{
			msgs.say("saving slice",Medium);
			h.saveSlice(2,0,"slice",0);
			}
		}

//		omp_unset_lock(&_omp_lock);

		for (unsigned long k = 0; k < Nx; k++) {
//			printf("setting h(i=%li)\n",i);
			h.setf(k,i,0,double(i)*Nx+k,0);			
		}

		
		
	}
	

}
/***************************************************************************************/
void gofxy() {
	mscsFunction3dregc f,g;
//	long i,j;
//	long Nx=90,Ny=90;
//	long Nx=850,Ny=600;
	long Nx=170,Ny=120;
	double x,y;
	
	
//	f.setSizeRange(3,3,1,1,1,1,2,2,2);
//	f.setSizeRange(3,3,1,0,0,0,Nx,Ny,0);
//	f.allocFunctionSpace();
//	f.fRe(0,0,0)=0;
//	f.fRe(0,1,0)=0;
//	f.fRe(0,2,0)=0;
//	f.fRe(1,0,0)=1;
//	f.fRe(1,1,0)=2;
//	f.fRe(1,2,0)=3;
//	f.fRe(2,0,0)=1;
//	f.fRe(2,1,0)=1;
//	f.fRe(2,2,0)=1;
	
	f.loadHDF5("comptonY-2D-res_50.0-scale_1.0.hdf","halo_609");
	f.setVerbosityLevel(Top);
	f.printInfo();
//	f.get
	g.setSizeRange(Nx,Ny,1,f.getMinX(),f.getMinY(),f.getMinZ(),f.getMaxX(),f.getMaxY(),f.getMaxZ());
	g.allocFunctionSpace();
	g.setVerbosityLevel(Top);
	g.printInfo();
	for (unsigned long i = 0; i < Nx; i++) {
		x=f.getMinX()+double(i)/Nx*f.lengthX()+g.getDxo2();
		for (unsigned long j = 0; j < Ny; j++) {
			y=f.getMinY()+double(j)/Ny*f.lengthY()+g.getDyo2();
			printf("i,j: %li, %li,   x,y: %lf, %lf\n",i,j,x,y);
			g.fRe(i,j,0)=f.fxy(x,y,0,0);
		}
	}
	
	f.saveHDF5("fxy_test.hdf","f");
	g.saveHDF5("fxy_test.hdf","g");
	
	
	
}
