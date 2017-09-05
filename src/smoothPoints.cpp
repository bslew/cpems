/*!
 * \file smoothPoints.cpp - scattered points smoother and regular grid smooth field generator
 *
 *  Created on: Feb 2, 2012
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <string.h>
#include <tclap/CmdLine.h>
#include "cpeds-msgs.h"
#include "cpeds-point_set.h"
#include "Mscs-function3dregc.h"
#include "lssKitFileHandler.h"
#include "pointsDensity.h"
#include <omp.h>

#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


#define HALO_INFO_FILE_COL_HaloID 27
#define HALO_INFO_FILE_COL_DomainID 25
#define HALO_INFO_FILE_COL_isInsideFlag 22


//bool _testHyperboloid, _testParaboloid, _testRT32model, _testDerivative, _testLine;
bool _calcDensity, _interpolate;
bool _gather,_scatter, _cube, _dims2,_inFormat3NN, _inFormatHDF5, _mkTestData,_mkTestData3d, _save3NNvalues, _shiftPointsToGridCells, _calcMass, _unitVolume, _useProvidedHSML, _twoTrees;
bool _processAllHalos;
long _haloInterpolationStartWithID;

long _Nx,_Ny, _Nz, _NneighborsMin,_NneighborsMax, _maxDepth, _ptsno;
double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax; // defines the region of interpolation in the deep field coordinates
double _haloSizeScaleFactor;
double _haloInterpolationComovingResolution; // [comoving Mpc]
double _maxHSML;
string _outfilePrefix, _infile, _outDirPrefix, _outDatasetName, _infileFormat, _haloFileName, _hdf5dset, _haloParticlesFileName;
subDomain_region_t _treeSubDomain, _treeSubDomain3;
long _from, _to;
long _currentHaloIdx=-1;
long _currentHaloID=-1;
long _currentDomainID=-1;
long _loadedDomainID=-1;
long _requestedDomainID=-1;
mscsFunction3dregc _halos;

lssPyramidSliceFileHeader_t _h;

void parseOptions(int argc, char** argv);

mscsFunction3dregc calcDensity(cpedsMsgs& msgs, bool calcMass=false, bool is2D=false, bool isSpherical=false, bool isUnitVolume=false);
mscsFunction3dregc makeInterpolation(cpedsMsgs& msgs);
void lss3NNloadPointsSet(cpedsPointSet3D& ps, cpedsMsgs& msgs, mscsVector<double>* hsml=NULL, mscsVector<double>* dens=NULL);
void lssHDF5loadPointsSet(cpedsPointSet3D& ps, cpedsMsgs& msgs, mscsVector<double>* hsml=NULL, mscsVector<double>* dens=NULL, long haloID=-1);
void setupGlobalVariables(cpedsMsgs& msgs);
void save3NNvalues(cpedsMsgs& msgs);


/*
 * Comment:
 *  These functions were designed for SZproject.
 *  
 * author: blew
 * date: Oct 10, 2013 3:09:21 PM
 *
 */
void performInSparseGrid(cpedsMsgs& msgs);
void saveToHDF(mscsFunction3dregc& f);
void makeSparseInterpolation(mscsFunction3dregc& f, cpedsMsgs& msgs);
//-----------------------------------

int main(int argc, char **argv) {
	cpedsMsgs msgs("smoothPoints");
	mscsFunction3dregc f("field");
	
	parseOptions(argc,argv);
	msgs.saveThisRun(argc,argv);
	
	if (_mkTestData) {
//		long pdone=0;
//		long i;
//#pragma omp parallel for  default(none) private(i) reduction(+:pdone)
//		for (i = 0; i < 10000; i++) {
//			pdone=pdone+1;
//			if (omp_get_thread_num()==0) printf("pdone: %li\n",pdone);
//			
//		}
//		exit(0);

		mscsVector<double> vals;
		pointsDensity dens;
		dens.append(cpedsPoint3D(0,0,0)); 		vals.push_back(1);
		dens.append(cpedsPoint3D(0.1,0.1,0)); 	vals.push_back(1);
		dens.append(cpedsPoint3D(0.2,0.2,0)); 		vals.push_back(2);
		dens.append(cpedsPoint3D(0.2,0.15,0)); 	vals.push_back(3);
//		dens.append(cpedsPoint3D(1,0,0)); 		vals.push_back(1);
//		dens.append(cpedsPoint3D(1,1,0)); 		vals.push_back(1);
//		dens.calculateDensity(2,2,true);		
		subDomain_region_t r,t;
		r.xmin=0; 		r.xmax=1;
		r.ymin=0; 		r.ymax=1;
		r.zmin=0; 		r.zmax=0;
		r.subx=_Nx;
		r.suby=_Nx;
		r.subz=1;
		t=r;
		t.ymin=-1;
		t.subx=2;
		t.suby=2;
		t.subz=1;
		if (_scatter) 
			f.mkInterpolatedFieldScatter(r,t,dens,vals,"gadget2",2,2);
		if (_cube) {
			f.mkDensityFieldCube(r,t,dens,&vals,_maxDepth,false,_calcMass,_unitVolume);				
		}
		f.saveSlice(2,0,"inter",0);

		exit(0);
		
		long nmax=_Nx;
		long pts=_ptsno;
		f.setSizeRange(nmax,nmax,1,0,0,0,1,1,0);
		f.allocFunctionSpace();
		f=double(0);
		
		mscsVector<cpedsPoint3D> ps;
		cpedsPointSet3D ps3d;
		for (long i = 0; i < pts; i++) {
//			ps.push_back(cpedsPoint3D(cpeds_random_uniform_number(0,1),cpeds_random_uniform_number(0,1),cpeds_random_uniform_number(0,1)));
			ps.push_back(cpedsPoint3D(cpeds_random_uniform_number(0,1),cpeds_random_uniform_number(0,1),0));
		}
		ps3d=ps;
		ps3d.save("testData.pts");

		subDomain_region_t reg;
		nmax++;
		reg.xmin=0;	reg.ymin=0; 	reg.zmin=0;
		reg.xmax=1;	reg.ymax=1; 	reg.zmax=0;
		reg.subx=nmax; reg.suby=nmax; reg.subz=1;
		_treeSubDomain.subz=1;
		f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,NULL,"gadget2",_NneighborsMin,_NneighborsMax);
		f.saveSlice(2,0,"density2d",0);
		exit(0);
		
	}

	if (_mkTestData3d) {

		mscsVector<double> vals;
		pointsDensity dens;
		dens.append(cpedsPoint3D(0,0,0)); 		vals.push_back(1);
		dens.append(cpedsPoint3D(0.1,0.1,0)); 	vals.push_back(2);
		dens.append(cpedsPoint3D(0.2,0.2,0)); 		vals.push_back(2);
		dens.append(cpedsPoint3D(0.2,0.15,0)); 	vals.push_back(3);
		subDomain_region_t r,t;
		r.xmin=0; 		r.xmax=1;
		r.ymin=0; 		r.ymax=1;
		r.zmin=0; 		r.zmax=1;
		r.subx=_Nx;
		r.suby=_Ny;
		r.subz=_Nz;
		t=r;
		t.subx=2;
		t.suby=2;
		t.subz=2;
		if (_scatter) 
			f.mkInterpolatedFieldScatter(r,t,dens,vals,"gadget2",2,2);
		if (_cube) {
			f.mkDensityFieldCube(r,t,dens,&vals,_maxDepth,false,_calcMass,_unitVolume);				
		}
		f.saveSlice(2,0,"density-slice0",0);
		f.saveHDF5("density.hdf5","cube");

		exit(0);
		
	}

	if (_inFormat3NN and _save3NNvalues) {
		save3NNvalues(msgs);
		exit(0);
	}
	
	if (_inFormat3NN or _inFormatHDF5) {
		if (_haloFileName!="") {
			performInSparseGrid(msgs);
			exit(0);			
		}
		else
			setupGlobalVariables(msgs);
	}

	if (_calcDensity) f=calcDensity(msgs); 
	if (_interpolate) f=makeInterpolation(msgs); 
	
	
	msgs.say("Saving to file", High);
	if (_inFormat3NN) {
		
						/*
						 * Comment: This format is meant to serve for deriving compton-y maps for the SZ effect, and so it is a dedicated
						 * format. It is assumed that if the input data are given in this format then the output format will be an
						 * hdf5 file.
						 * 
						 * Any other usages are handled in the "else" block
						 * 
						 * author: blew
						 * date: Mar 22, 2012, 12:27:18 PM
						 *
						 */
//		f.saveHDF5(_outfilePrefix,_outDatasetName,(_h.DNx-_h.SNx)/2,(_h.DNy-_h.SNy)/2,(_h.DNz-_h.SNz)/2,  _h.SNx,_h.SNy,_h.SNz,0);
//		f.saveHDF5(_outfilePrefix,_outDatasetName+"_neigh",(_h.DNx-_h.SNx)/2,(_h.DNy-_h.SNy)/2,(_h.DNz-_h.SNz)/2,  _h.SNx,_h.SNy,_h.SNz,1);
//		f.saveAllSlices("pyramid","slice",2,0);
		f.saveHDF5(_outfilePrefix,_outDatasetName,0,0,0,  _h.SNx,_h.SNy,_h.SNz,0);
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"x0_comment","function domain lower range in x direction [comoving Mpc]");
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"y0_comment","function domain lower range in y direction [comoving Mpc]");
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"z0_comment","function domain lower range in z direction [comoving Mpc]");
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"dx_comment","function cell size in x direction [comoving Mpc]");
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"dy_comment","function cell size in y direction [comoving Mpc]");
//		f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"dz_comment","function cell size in z direction [comoving Mpc]");
//		f.saveHDF5(_outfilePrefix,_outDatasetName+"_neigh",0,0,0,  _h.SNx,_h.SNy,_h.SNz,1);
		// slice position in the deep field
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SX",_h.SX,"slice position in the deep field in x direction [comoving Mpc]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SY",_h.SY,"slice position in the deep field in y direction [comoving Mpc]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SZ",_h.SZ,"slice position in the deep field in z direction [comoving Mpc]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdx",_h.Sdx,"slice size in x direction [comoving Mpc]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdy",_h.Sdy,"slice size in y direction [comoving Mpc]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdz",_h.Sdz,"slice size in z direction [comoving Mpc]");
		// cosmology parameters
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"redshift",_h.z);
//		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"X",cpeds_comoving_distance(_h.Wb0,_h.Wm0,_h.Wl0,_h.z,_h.w0,_h.h),"comoving distance to slice redshift";
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"hubble",_h.h,"unitless Hubble parameter");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Wb0",_h.Wb0,"baryon density today in units of critical density");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"compton_y0",_h.y0,"integral constant: y0 = sigma_T kB rhoC0(h=1) Wb0_WMAP / ( mu_e_BBN m_p * m_e c^2 ) * CPEDS_MPC  = 1.50582*10^-16 [ K^-1 Mpc^-1 ]");
		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"mu_e",_h.mu_e, "electron number per proton mass");
		
	}
	else {
		if (_dims2) {
			f.saveSlice(2,0,_outfilePrefix,0);
		}
		else {
			f.saveAllSlices(_outDirPrefix,_outfilePrefix,2,0);
		}		
	}
	msgs.say("done", Low);
	
	
	return 0;
}




void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd("density calculator. Calculates adaptive density estimate from the distribution of points, optionally with masses (weights).\n"
				"It can also perform adaptive interpolations on a set of points. "
				"\n"
				"For example, in order to smooth points stored in the file pointSet defined over the range from 0 to 10 in every dimension using the octTree scheme"
				"and 3 neighbors for defining the smoothing length, and with the output density field consisting of 30x30x30 pixels you can use:\n\n"
				"smoothPoints -i pointSet --sdZ 2 --Nneighbors 3 --zmax 10 --zmin 0 --ymax 10 --ymin 0 --xmax 10 --xmin 0 --rho --Nx 30 --Ny 30 --Nz 30 -o density\n\n"
				""
				"The output file names will be named density and will be placed in the default out subdirectory created in the current directory.\n\n"
				""
				"To smooth points in 2D set --sdZ 1 and Nz 1."
				"smoothPoints -i pointSet --sdZ 1 --Nneighbors 1 --zmax 10 --zmin 0 --ymax 10 --ymin 0 --xmax 10 --xmin 0 --rho --Nx 10 --Ny 10 --Nz 1 -o tmp\n\n"
				""
				"If you only smooth 1 particle then the hsml is not well defined in the scatter scheme. Currently it is defined arbitrarily to be 10xsqrt(dx^2+dy^2) for the case"
				"of 2d smoothing and 10xsqrt(dx^2+dy^2+ dz^2) for the case of 3d smoothing. There is no way to control it for the moment.\n"
				"\n"
				"This program has a special usage (--inFormat3NN) for interpolating deep fields for the deepSky project. It is started from mkPyramidGrid.py script"
				"and and in this mode the interpolation parameters are read from the header of the input 3NN format file and the outputs are stored"
				"in hdf5 file. The data stored in that file is density weighted temperature which can be used for integrating along the line of sight.\n"
				"Since the geometry of the preparred input 3NN files is such that there is the same number of pixels in x,y directions for every slice"
				"the LOS integration in the output hdf5 cube is simply an integration along z direction in the cube.", ' ', "" );
		
		// 
		// Define arguments
		//
		
		// tests
//		SwitchArg testHyperboloid("","testHyperboloid", "tests the hyperboloid generation", false);	cmd.add( testHyperboloid );
//		SwitchArg testParaboloid("","testParaboloid", "tests the paraboloid generation", false);	cmd.add( testParaboloid );
//		SwitchArg testRT32model("","testRT32model", "tests the RT32 model generation with basic reflections off primary and secondary mirrors and taper shaping.", false);	cmd.add( testRT32model );
//		SwitchArg testDerivative("","testDerivative", "tests derivatives calculation on the interpolated surfaces.", false);	cmd.add( testDerivative );
//		SwitchArg testLine("","testLine", "tests mscsLine class.", false);	cmd.add( testLine);
		
		// switches
		SwitchArg mkTestData("","mkTestData", "Generate test data - do not load input data file", false);	cmd.add(mkTestData);
		SwitchArg mkTestData3d("","mkTestData3d", "Generate 3d test data - do not load input data file", false);	cmd.add(mkTestData3d);
		ValueArg<long> ptsno("", "ptsno", "number of generated points for test data (1000)", false,1000,"long");	cmd.add( ptsno );

		SwitchArg interpolate("","int", "Interpolate on a regular grid a 2D field represented by a distribution of points in space.", false);	cmd.add(interpolate);
		SwitchArg calcDensity("","rho", "Derive the density distribution of points on a regular grid.", false);	cmd.add(calcDensity);
		SwitchArg gather("","gather", "Use gather approach for density calculation.", false);	cmd.add(gather);
		SwitchArg scatter("","scatter", "Use scatter approach for density calculation.", false);	cmd.add(scatter);
		SwitchArg cube("","cube", "Use cube approach for density calculation. This is way faster, as this only does domain decomposition into requested"
				"number of subdomains, and then the returned density function gives simply the number of points in each subdomain divided"
				"by the subdomain volume.", false);	cmd.add(cube);
		SwitchArg inFormat3NN("","inFormat3NN", "Indicates the 3NN input file format. All the required information will be taken from the header of the input file. "
				"All domain related options are ignored.", false);	cmd.add(inFormat3NN);
		SwitchArg inFormatHDF5("","inFormatHDF5", "Indicates the hdf5 input file format. All the required information will be taken from the header of the input file. "
				"All domain related options are ignored.", false);	cmd.add(inFormatHDF5);		
		ValueArg<string> hdf5dset("", "hdf5dset", "name of the hdf5 dataset to use ", false,"","string");	cmd.add( hdf5dset );
		

		SwitchArg save3NNvalues("","save3NNvals", "to be used only with --inFormat3NN. Dumps all values data to txt file.", false);	cmd.add(save3NNvalues);
		SwitchArg shiftPointsToGridCells("","shiftToGrid", "Shifts particles to grid cell centers.", false);	cmd.add(shiftPointsToGridCells);
		SwitchArg calcMass("","calcMass", "Derive the mass distribution on a regular grid - this a --rho and --int options modifier only", false);	cmd.add(calcMass);
		SwitchArg unitVolume("","unitVolume", "assume a unit volume in density calculations or in interpolations", false);	cmd.add(unitVolume);
		SwitchArg useProvidedHSML("","useProvidedHSML", "triggers using pre-calculated hsml for SPH (scatter) interpolations. Eg if data are loaded from 3NNN files then "
				"the last columns is taken to be hsml and it can be used to speed up the density calculations needed for interpolations.", false);	cmd.add(useProvidedHSML);
		
		
		// options
//		ValueArg<long> Nside("N", "Nside", "number of points on the regular grid on the mirrors along side (1000)", false,1000,"long");	cmd.add( Nside );
		ValueArg<long> Nx("", "Nx", "number of points on the regular grid along X direction (1000)", false,1000,"long");	cmd.add( Nx );
		ValueArg<long> Ny("", "Ny", "number of points on the regular grid along Y direction (1000)", false,1000,"long");	cmd.add( Ny );
		ValueArg<long> Nz("", "Nz", "number of points on the regular grid along Z direction (1000)", false,1000,"long");	cmd.add( Nz );
		ValueArg<double> xmin("", "xmin", "minimal x coordinate on the grid (0)", false,0,"double");	cmd.add( xmin );
		ValueArg<double> xmax("", "xmax", "maximal x coordinate on the grid (1)", false,1,"double");	cmd.add( xmax );
		ValueArg<double> ymin("", "ymin", "minimal y coordinate on the grid (0)", false,0,"double");	cmd.add( ymin );
		ValueArg<double> ymax("", "ymax", "maximal y coordinate on the grid (1)", false,1,"double");	cmd.add( ymax );
		ValueArg<double> zmin("", "zmin", "minimal z coordinate on the grid (0)", false,0,"double");	cmd.add( zmin );
		ValueArg<double> zmax("", "zmax", "maximal z coordinate on the grid (1)", false,1,"double");	cmd.add( zmax );
		ValueArg<long> NneighborsMin("", "NneighborsMin", "number of neighbors to use in calculation (30)", false,30,"long");	cmd.add( NneighborsMin );
		ValueArg<long> NneighborsMax("", "NneighborsMax", "number of neighbors to use in calculation (30)", false,30,"long");	cmd.add( NneighborsMax );
		ValueArg<long> subDomainX("", "sdX", "number of subdomain divisions along X direction. 2 for octTree (2)", false,2,"long");	cmd.add( subDomainX );
		ValueArg<long> subDomainY("", "sdY", "number of subdomain divisions along Y direction. 2 for octTree (2)", false,2,"long");	cmd.add( subDomainY );
		ValueArg<long> subDomainZ("", "sdZ", "number of subdomain divisions along Z direction. 2 for octTree. "
				"For 2D case use sdZ=1. (2)", false,2,"long");	cmd.add( subDomainZ );

		ValueArg<long> maxDepth("", "maxDepth", "maximal depth of the tree. Only with --cube option. (default: -1 = no end)", false,-1,"long");	cmd.add( maxDepth );

		ValueArg<long> from("", "from", "number to start dumping from (for save3NNvalues only)", false,0,"long");	cmd.add( from );
		ValueArg<long> to("", "to", "number to end dumping at (for save3NNvalues only)", false,100,"long");	cmd.add( to );			
		ValueArg<string> infileFormat("","infmt","input file name format, to be used with options --from and --to. To be used with --save3NNvals only.",false,"","string");	cmd.add( infileFormat );
		
		ValueArg<string> infile("i","infie","file with list of points and their masses (weights) - 4 column format for 3D case (x,y,z,w) and 3 column formar for 2D case (x,y,w)",false,"infile","string");	cmd.add( infile );
		ValueArg<string> outfilePrefix("o","outfile","output file name prefix (out)",false,"out","string");	cmd.add( outfilePrefix );
		ValueArg<string> outDirPrefix("","outDir","output directory name prefix for 3d smoothing (out)",false,"out","string");	cmd.add( outDirPrefix );
		ValueArg<string> outDatasetName("","outDatasetName","name for the dataset to be saved in HDF5 file",false,"slice","string");	cmd.add( outDatasetName);

		SwitchArg processAllHalos("","processAllHalos", "special use for the SZproject. This flag makes the program to process all halos in the haloFileName and not "
				"only those inside of the field of view (default), i.e. having flag 1 set in the proper column in the haloFileName.", false);	cmd.add(processAllHalos);
		ValueArg<string> haloFileName("","haloFileName","special use for the SZproject. With this option you can provide the information on the locations"
				"of the CDM halo inside of the grid in order to speed up the calculations by calculating the interpolated values only in the grid cells around"
				"the halos as defined by the halo position and size - both provided in the haloFileName (the format of the file should be compatible"
				"with the outputs from the calculate-halo-properties program of the lssKit package). (default: '')",false,"","string");	cmd.add( haloFileName);
		ValueArg<string> haloParticlesFileName("","haloParticlesFileName","special use for the SZproject. With this option you can provide the information on the location"
				"of the hdf5 file containing indexes of the particles inside of the input hdf5 file that build up a halo. (default: '')",false,"","string");	cmd.add( haloParticlesFileName);

		ValueArg<double> haloSizeScaleFactor("", "haloSizeScaleFactor", "used only if haloFileName is provided. The linear sizes of the parallelpipe in which the interpolation is done are"
				"increased by this factor. (Default: 1)", false,1,"double");	cmd.add( haloSizeScaleFactor );
		
		ValueArg<double> haloInterpolationComovingResolution("", "haloInterpolationComovingResolution", "used only if haloFileName is provided. The spatial resolution"
				"at which the halo is to be interpolated. This effectively defines the number of grid cells to hold each halo. Given variable halo sizes the number of cells"
				"will vary from halo to halo. [kpc]"
				"(Default: 50 kpc)", false,50,"double");	cmd.add( haloInterpolationComovingResolution );
		ValueArg<long> requestedDomainID("", "requestedDomainID", "Used only if haloFileName is provided. Indicates the domain ID corresponding to the input file. "
				"This could be extracted from parsing input file name but this is to make things general. This value will be set by the "
				"mkPyramidGrid.py program automatically when interpolating deep field. "
				"This is needed to calculate interpolated values only for the halos that lie inside of the currently loaded slice."
				"(default: -1)", false,-1,"long");	cmd.add( requestedDomainID );
		ValueArg<double> maxHSML("", "maxHSML", "used only if haloFileName is provided. THe maximal size of the hsml used for particle selection"
				"for interpolation. If the hsml from the provided data is smaller than this value then it is used, but if it is larger then"
				"maxHSML is used instead. This effectively decreases the size of the cube when searching for particles that should participate"
				"in interpolation on the grid cells. [Mpc]"
				"(Default: 2 Mpc)", false,2,"double");	cmd.add( maxHSML );

		ValueArg<long> haloInterpolationStartWithID("", "haloInterpolationStartWithID", "Used only if haloFileName is provided. You can restart a previous run with this option "
				"from such provided halo ID."
				"(default: -1 - not used)", false,-1,"long");	cmd.add( haloInterpolationStartWithID );
		
		SwitchArg twoTrees("","twoTrees", "use extra tree with a different domain splits to optimize neighbors finding in case "
				"of grid cells laying on the tree domain boundaries", false);	cmd.add(twoTrees);
		
		
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_infile=infile.getValue();
		_outfilePrefix=outfilePrefix.getValue();
		_outDirPrefix=outDirPrefix.getValue();
		_outDatasetName=outDatasetName.getValue();
		
		_inFormat3NN=inFormat3NN.getValue();
		_inFormatHDF5=inFormatHDF5.getValue();
		_hdf5dset=hdf5dset.getValue();
		_gather=gather.getValue();
		_scatter=scatter.getValue();
		_cube=cube.getValue();
		_calcDensity=calcDensity.getValue();
		_calcMass=calcMass.getValue();
		_unitVolume=unitVolume.getValue();
		_useProvidedHSML=useProvidedHSML.getValue();
		_interpolate=interpolate.getValue();
		_Nx=Nx.getValue();
		_Ny=Ny.getValue();
		_Nz=Nz.getValue();
		_xmin=xmin.getValue();
		_xmax=xmax.getValue();
		_ymin=ymin.getValue();
		_ymax=ymax.getValue();
		_zmin=zmin.getValue();
		_zmax=zmax.getValue();
		_dims2=false;
		if (_zmin==_zmax) _dims2=true;

		_maxDepth=maxDepth.getValue();
		
		_treeSubDomain.subx=subDomainX.getValue();
		_treeSubDomain.suby=subDomainY.getValue();
		_treeSubDomain.subz=subDomainZ.getValue();
		_treeSubDomain.xmin=_xmin;
		_treeSubDomain.ymin=_ymin;
		_treeSubDomain.zmin=_zmin;
		_treeSubDomain.xmax=_xmax;
		_treeSubDomain.ymax=_ymax;
		_treeSubDomain.zmax=_zmax;
		_treeSubDomain3=_treeSubDomain;
		_treeSubDomain3.subx=_treeSubDomain.subx+1;
		_treeSubDomain3.suby=_treeSubDomain.suby+1;
		_treeSubDomain3.subz=_treeSubDomain.subz+1;
		
		_NneighborsMin=NneighborsMin.getValue();
		_NneighborsMax=NneighborsMax.getValue();
		_mkTestData=mkTestData.getValue();
		_mkTestData3d=mkTestData3d.getValue();
		_ptsno=ptsno.getValue();
		_save3NNvalues=save3NNvalues.getValue();
		_from=from.getValue();
		_to=to.getValue();
		_infileFormat=infileFormat.getValue();
		_shiftPointsToGridCells=shiftPointsToGridCells.getValue();
		_haloFileName=haloFileName.getValue();
		_haloParticlesFileName=haloParticlesFileName.getValue();
		_haloSizeScaleFactor=haloSizeScaleFactor.getValue();
		_haloInterpolationComovingResolution=haloInterpolationComovingResolution.getValue()/1000; // convert to Mpc
		_requestedDomainID=requestedDomainID.getValue();
		_haloInterpolationStartWithID=haloInterpolationStartWithID.getValue();
		_maxHSML=maxHSML.getValue();
		_twoTrees=twoTrees.getValue();
		_processAllHalos=processAllHalos.getValue();
		
		/* 	_save_overplot_as = save_overplot_as.getValue(); */
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}
/***************************************************************************************/
void lss3NNloadPointsSet(cpedsPointSet3D& ps, cpedsMsgs& msgs, mscsVector<double>* hsml, mscsVector<double>* dens) {
	msgs.say("loading points from 3NN file", High);
	if (_inFormat3NN) {
		float* data=NULL;
		lssPyramidSliceFileHeader_t h;

//		lssPyramidSliceFileQuadruplet_t* data;
		long Npart, Ncol;
//		lssKitLoad3NNData(_infile,&data,&Npart);
		data=lssKitLoad3NNNData(_infile,&h);
		Npart=h.Npart;
		Ncol=h.scalarColNum;
		ps.setSize(Npart);
		//		for (long i = 0; i < Npart; i++) {
		//			ps.set(i,data[i].x,data[i].y,data[i].z,data[i].T);
		//		}
		for (long i = 0; i < Npart; i++) {
			ps.set(i,data[Ncol*i],data[Ncol*i+1],data[Ncol*i+2],data[Ncol*i+3]);
		}
		if (Ncol>=5 and hsml!=NULL) {
			hsml->setSize(Npart);
			for (long i = 0; i < Npart; i++) hsml->at(i)=data[Ncol*i+4];
		}
		if (Ncol>=6 and dens!=NULL) {
			dens->setSize(Npart);
			for (long i = 0; i < Npart; i++) dens->at(i)=data[Ncol*i+5];
		}
		
		double xmin,xmax,ymin,ymax,zmin,zmax;
		ps.getRanges(xmin,xmax,ymin,ymax,zmin,zmax);
		printf("particles lie within:\n");
		printf("[xmin,xmax] = [%lE,%lE]\n",xmin,xmax);
		printf("[ymin,ymax] = [%lE,%lE]\n",ymin,ymax);
		printf("[zmin,zmax] = [%lE,%lE]\n",zmin,zmax);
//		ps.print();
//		exit(0);
		if (data!=NULL) delete [] data;
	}
	else {
		ps.load(_infile);		
	}
	msgs.say("done", Low);
	

}
/***************************************************************************************/
void lssHDF5loadPointsSet(cpedsPointSet3D& ps, cpedsMsgs& msgs, mscsVector<double>* hsml, mscsVector<double>* dens, long haloID) {
	msgs.say("loading points from HDF5 file", High);
	static mscsFunction3dregc data,tmpdata;
	mscsFunction3dregc partIdx;
	vector<string> dSetNames;
	double shiftZ0=0;
	double simID;
	string ptype;
//	double *idx;
	if (_loadedDomainID!=_requestedDomainID ) { //or _currentDomainID==-1
		msgs.say("loading data from: "+_infile,High);
		dSetNames=data.getHDF5dataSets(_infile);
//		data.loadHDF5(_infile,_hdf5dset);
		simID=data.getHDF5_scalarDoubleAttribute(_infile,_hdf5dset,"simID",NULL);
//		printf("dataset count: %li\n",);
		
		/*
		 * Comment: load all data from the input file that are relevant to the requested simulation ID
		 * and concatenate the results into a one big array
		 * This is needed because fof typically runs over all particles from the full simulation box
		 * so in order to match particle indexes  (and store halo particle indexes) we need to load all data
		 * for the requested simulation ID iterating through domains/slabs in the correct order.
		 * 
		 * We do not need to load CDM data because all baryons are stored first in the tipsy file and then they are
		 * followed by all CDM particles.
		 * 
		 * The data function is static to avoid loading the same data for every halo located in the same
		 * simulation box. So when we return to this routine to process the next halo from the same domain/slab
		 * then this part is skipped.
		 * 
		 * author: blew
		 * date: Oct 26, 2013 11:21:02 AM
		 *
		 */
		for (unsigned long i = 0; i < dSetNames.size(); i++) {
			ptype=data.getHDF5_stringAttribute(_infile,dSetNames[i],"particle_type",NULL);
			if (ptype=="baryon") {
				printf("looking for simID=%.0lf in %s\n",simID,dSetNames[i].c_str());
				if (data.getHDF5_scalarDoubleAttribute(_infile,dSetNames[i],"simID",NULL)==simID) {
					printf(" -found simID=%.0lf in %s\n",simID,dSetNames[i].c_str());
					tmpdata.loadHDF5(_infile,dSetNames[i],0);
					printf("data size is: %li\n",data.Ny());
					data=data.concatenate(tmpdata,1);
				}				
			}
		}
		printf("data size is: %li\n",data.Ny());
//		exit(0);
//		data.saveSlice(2,0,"sphdata.tmp",0);
//		shiftZ0=data.getHDF5_scalarDoubleAttribute(_infile,_hdf5dset,"simBoxZ0",NULL);
//		_loadedDomainID=data.getHDF5_scalarDoubleAttribute(_infile,_hdf5dset,)
	}
	

	long l,jj,Npart=0;
	
	if (_haloParticlesFileName!="") {
		partIdx.loadHDF5(_haloParticlesFileName,"halo_"+msgs.toStr(haloID));
//		partIdx.printFunction();
//		idx=partIdx.getSlice1D(1,0,0,0);
//		Npart=partIdx.Ny();
		Npart=0;
		for (long j = 0; j < partIdx.Ny(); j++) {
			if (partIdx.fRe(1,j,0)==0) Npart++; // count only gas particles
		}

		ps.setSize(Npart);
//		ps.clear();
		jj=0;
		printf("loading gas particles\n");
		for (long j = 0; j < partIdx.Ny(); j++) {
			if (partIdx.fRe(1,j,0)==0) { // count only gas particles
				l=partIdx.fRe(0,j,0);
				ps.set(jj,data.fRe(0,l,0),data.fRe(1,l,0),data.fRe(2,l,0),data.fRe(10,l,0));
				jj++;
				
			}
//			l=partIdx.fRe(0,j,0);
////			ps.set(j,data.fRe(0,l,0),data.fRe(1,l,0),data.fRe(2,l,0)-shiftZ0,data.fRe(10,l,0));
//			ps.set(j,data.fRe(0,l,0),data.fRe(1,l,0),data.fRe(2,l,0),data.fRe(10,l,0));
		}
		if (hsml!=NULL) {
			hsml->setSize(Npart);
			jj=0;
			for (long j = 0; j < partIdx.Ny(); j++) {
				if (partIdx.fRe(1,j,0)==0) { // count only gas particles
					l=partIdx.fRe(0,j,0);
					hsml->at(jj)=data.fRe(9,l,0);
					jj++;
				}
			}
		}
		if (dens!=NULL) {
			dens->setSize(Npart);
			jj=0;
			for (long j = 0; j < partIdx.Ny(); j++) {
				if (partIdx.fRe(1,j,0)==0) { // count only gas particles
					l=partIdx.fRe(0,j,0);
					dens->at(jj)=data.fRe(8,l,0);
					jj++;
				}
			}
		}
//		if (hsml!=NULL) {
//			hsml->setSize(Npart);
//			for (long j = 0; j < Npart; j++) {
//				l=partIdx.fRe(0,j,0);
//				hsml->at(j)=data.fRe(9,l,0);
//			}
//		}
//		if (dens!=NULL) {
//			dens->setSize(Npart);
//			for (long j = 0; j < Npart; j++) {
//				l=partIdx.fRe(0,j,0);
//				dens->at(j)=data.fRe(8,l,0);
//			}
//		}

	}
	else {		
		Npart=data.Ny();
		ps.setSize(Npart);
		for (long i = 0; i < Npart; i++) {
//			ps.set(i,data.fRe(0,i,0),data.fRe(1,i,0),data.fRe(2,i,0)-shiftZ0,data.fRe(10,i,0));
			ps.set(i,data.fRe(0,i,0),data.fRe(1,i,0),data.fRe(2,i,0),data.fRe(10,i,0));
		}
		if (hsml!=NULL) {
			hsml->setSize(Npart);
			for (long i = 0; i < Npart; i++) hsml->at(i)=data.fRe(9,i,0);
		}
		if (dens!=NULL) {
			dens->setSize(Npart);
			for (long i = 0; i < Npart; i++) dens->at(i)=data.fRe(8,i,0);
		}
	}
	
	
	double xmin,xmax,ymin,ymax,zmin,zmax;
	ps.getRanges(xmin,xmax,ymin,ymax,zmin,zmax);
//	ps.print();
	printf("loaded %li particles\n",ps.size());
	printf("particles lie within:\n");
	printf("[xmin,xmax] = [%lE,%lE]\n",xmin,xmax);
	printf("[ymin,ymax] = [%lE,%lE]\n",ymin,ymax);
	printf("[zmin,zmax] = [%lE,%lE]\n",zmin,zmax);
//	printf("[zmin,zmax] = [%lE,%lE]\n",ps.vzmin,zmax);
//	ps.values().printVector();
//	dens->printVector();
//	hsml->printVector();
	msgs.say("done", Low);	
}
/***************************************************************************************/
mscsFunction3dregc calcDensity(cpedsMsgs& msgs, bool calcMass, bool is2D, bool isSpherical, bool isUnitVolume) {
	msgs.say("allocating function space", High);
	double dx,dy,dz;
	dx=(_xmax-_xmin)/_Nx;
	dy=(_ymax-_ymin)/_Ny;
	dz=(_zmax-_zmin)/_Nz;
	mscsFunction3dregc f(_Nx,_Ny,_Nz,dx,dy,dz,_xmin,_ymin,_zmin);
	msgs.say("done", Low);
	subDomain_region_t reg;
	reg.xmin=_xmin;	reg.ymin=_ymin; 	reg.zmin=_zmin;
	reg.xmax=_xmax;	reg.ymax=_ymax; 	reg.zmax=_zmax;
	reg.subx=_Nx; reg.suby=_Ny; reg.subz=_Nz;
	
	
	cpedsPointSet3D ps;
	lss3NNloadPointsSet(ps,msgs);
	
	
	if (_shiftPointsToGridCells) {
		ps.save("smooth-points.ps_before_shift");
		for (long i = 0; i < ps.size(); i++) {
			ps[i].setX(long(ps[i].x() / dx) + dx/2);
			ps[i].setY(long(ps[i].y() / dy) + dy/2);
			ps[i].setZ(long(ps[i].z() / dz) + dz/2);
		}		
		ps.save("smooth-points.ps_after_shift");
	}
	
	if (_gather) {		
		msgs.say("smoothing with gather method, using "+msgs.toStr(_NneighborsMin)+" neighbors.", High);
		if (ps.values().size()==0)	f.mkDensityFieldGather(reg,_treeSubDomain,ps,NULL,0,"gadget2",_NneighborsMin);
		else
			f.mkDensityFieldGather(reg,_treeSubDomain,ps,&(ps.values()),0,"gadget2",_NneighborsMin);
	}
	if (_scatter) {
		msgs.say("smoothing with scatter method, using [%li,%li] neighbors.",_NneighborsMin,_NneighborsMax, High);
		if (ps.values().size()==0)	f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,NULL,"gadget2",_NneighborsMin,_NneighborsMax);
		else
			f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,&(ps.values()),"gadget2",_NneighborsMin,_NneighborsMax);			
	}
	if (_cube) {
		msgs.say("smoothing with cube method.", High);
		if (ps.values().size()==0)	{
//			f.mkDensityFieldCube(reg,_treeSubDomain,ps,NULL,0,_maxDepth);
			f.setVerbosityLevel(Low);
			f.mkDensityFieldCube(_Nx,_Ny,_Nz,ps,_maxDepth,is2D,isSpherical,calcMass,isUnitVolume);
		}
		else {
			f.setVerbosityLevel(Low);
			f.mkDensityFieldCube(reg,_treeSubDomain,ps,&(ps.values()),_maxDepth,isSpherical,calcMass,isUnitVolume);	
		}
	}

	msgs.say("done", Low);
	
	return f;
}

/***************************************************************************************/
void setupGlobalVariables(cpedsMsgs& msgs) {
	msgs.say("resetting global variables",High);
	lssPyramidSliceFileHeader_t h;
	if (_inFormat3NN) {
		h=lssKitLoad3NNData(_infile);
		lssKitPrint3NNFileHeader(h);

	}
	if (_inFormatHDF5) {
		mscsFunction3dregc f;
		string dsetName="gas-";
		h.DNx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"DNx",NULL);
		h.DNy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"DNy",NULL);
		h.DNz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"DNz",NULL);
		h.Ddx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Ddx",NULL);
		h.Ddy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Ddy",NULL);
		h.Ddz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Ddz",NULL);
		h.Dx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Dx",NULL);
		h.Dy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Dy",NULL);
		h.Dz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Dz",NULL);
		h.SNx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SNx",NULL);
		h.SNy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SNy",NULL);
		h.SNz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SNz",NULL);
		h.SX=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SX",NULL);
		h.SY=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SY",NULL);
		h.SZ=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"SZ",NULL);
		h.Sdx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sdx",NULL);
		h.Sdy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sdy",NULL);
		h.Sdz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sdz",NULL);
		h.Sx=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sx",NULL);
		h.Sy=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sy",NULL);
		h.Sz=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Sz",NULL);
		h.Wb0=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"Wb0",NULL);
		h.h=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"hubble",NULL);
		h.mu_e=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"mu_e",NULL);
		h.simBoxZ0=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"simBoxZ0",NULL);
		h.simID=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"simID",NULL);
		h.y0=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"compton-y0",NULL);
		h.z=f.getHDF5_scalarDoubleAttribute(_infile,dsetName+msgs.toStr(_requestedDomainID,"%03i"),"redshift",NULL);

		
	}
	_h=h;

	_xmin=h.Sx;
	_ymin=h.Sy;
	_zmin=h.Sz;
	_xmax=(h.Sx+h.Sdx);
	_ymax=(h.Sy+h.Sdy);
	_zmax=(h.Sz+h.Sdz);

	_Nx=h.SNx;
	_Ny=h.SNy;
	_Nz=h.SNz;
	
	_treeSubDomain.xmin=_h.Dx; // the unit is comoving Mpc
	_treeSubDomain.ymin=_h.Dy;
	_treeSubDomain.zmin=_h.Dz;
	_treeSubDomain.xmax=(_h.Dx+_h.Ddx);
	_treeSubDomain.ymax=(_h.Dy+_h.Ddy);
	_treeSubDomain.zmax=(_h.Dz+_h.Ddz);

	
	if (_currentHaloIdx!=-1) {
		_xmin=_halos.fRe(0,_currentHaloIdx,0)-_halos.fRe(3,_currentHaloIdx,0)*_haloSizeScaleFactor/2;
		_xmax=_halos.fRe(0,_currentHaloIdx,0)+_halos.fRe(3,_currentHaloIdx,0)*_haloSizeScaleFactor/2;
		_ymin=_halos.fRe(1,_currentHaloIdx,0)-_halos.fRe(4,_currentHaloIdx,0)*_haloSizeScaleFactor/2;
		_ymax=_halos.fRe(1,_currentHaloIdx,0)+_halos.fRe(4,_currentHaloIdx,0)*_haloSizeScaleFactor/2;
		_zmin=_halos.fRe(2,_currentHaloIdx,0)-_halos.fRe(5,_currentHaloIdx,0)*_haloSizeScaleFactor/2; // this is in the deep field coordinates
		_zmax=_halos.fRe(2,_currentHaloIdx,0)+_halos.fRe(5,_currentHaloIdx,0)*_haloSizeScaleFactor/2; // this is in the deep field coordinates
		
		_Nx=int((_xmax-_xmin)/_haloInterpolationComovingResolution);
		_Ny=int((_ymax-_ymin)/_haloInterpolationComovingResolution);
		_Nz=int((_zmax-_zmin)/_haloInterpolationComovingResolution);

//		// translation along Z to convert from deep field coordinates to sim box cooridnates
//		_zmin-=h.simBoxZ0;
//		_zmax-=h.simBoxZ0;
	
		
		_treeSubDomain.xmin=_xmin; // the unit is comoving Mpc
		_treeSubDomain.ymin=_ymin;
		_treeSubDomain.zmin=_zmin;
		_treeSubDomain.xmax=_xmax;
		_treeSubDomain.ymax=_ymax;
		_treeSubDomain.zmax=_zmax;

	}
//	if (_Nz==1) {
//		_treeSubDomain.subz=1;
//	}


	_treeSubDomain3=_treeSubDomain;
	_treeSubDomain3.subx=_treeSubDomain.subx+1;
	_treeSubDomain3.suby=_treeSubDomain.suby+1;
	_treeSubDomain3.subz=_treeSubDomain.subz+1;

}
/***************************************************************************************/
mscsFunction3dregc makeInterpolation(cpedsMsgs& msgs) {
	msgs.say("allocating function space for interpolation on regular grid", High);
	double dx,dy,dz;
	dx=(_xmax-_xmin)/_Nx;
	dy=(_ymax-_ymin)/_Ny;
	dz=(_zmax-_zmin)/_Nz;
	mscsFunction3dregc f(_Nx,_Ny,_Nz,dx,dy,dz,_xmin,_ymin,_zmin);
	msgs.say("done", Low);
	subDomain_region_t reg;
	reg.xmin=_xmin;	reg.ymin=_ymin; 	reg.zmin=_zmin;
	reg.xmax=_xmax;	reg.ymax=_ymax; 	reg.zmax=_zmax;
	reg.subx=_Nx; reg.suby=_Ny; reg.subz=_Nz;
	
	
	cpedsPointSet3D ps;
	
	if (_scatter) {
		msgs.say("interpolating with scatter method, using [%li,%li] neighbors.",_NneighborsMin,_NneighborsMax, High);
		mscsVector<double> hsml;
		if (_useProvidedHSML) {
			lss3NNloadPointsSet(ps,msgs,&hsml);
			f.mkInterpolatedFieldScatter(reg,_treeSubDomain,ps,ps.values(),"gadget2",_NneighborsMin,_NneighborsMax,&hsml);
		}
		else {
			lss3NNloadPointsSet(ps,msgs,NULL);
			f.mkInterpolatedFieldScatter(reg,_treeSubDomain,ps,ps.values(),"gadget2",_NneighborsMin,_NneighborsMax,NULL);			
		}
//		if (ps.values().size()==0)	f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,NULL,"gadget2",_NneighborsMin,_NneighborsMax);
//		else
//			f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,&(ps.values()),"gadget2",_NneighborsMin,_NneighborsMax);			
	}
	if (_cube) {
		msgs.say("interpolating with cube method.", High);
//		f=calcDensity(msgs,true,true,false,false);
		if (ps.values().size()==0)	{
//			f.mkDensityFieldCube(_Nx,_Ny,_Nz,ps,_maxDepth,is2D,isSpherical,calcMass,isUnitVolume);
			f.setVerbosityLevel(Low);
			f.mkDensityFieldCube(_Nx,_Ny,_Nz,ps,_maxDepth,_dims2,false,_calcMass,_unitVolume);
		}
		else {
			f.setVerbosityLevel(Low);
			f.mkDensityFieldCube(reg,_treeSubDomain,ps,&(ps.values()),_maxDepth,false,_calcMass,_unitVolume);	
		}

//		if (ps.values().size()==0)	f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,NULL,"gadget2",_NneighborsMin,_NneighborsMax);
//		else
//			f.mkDensityFieldScatter2(reg,_treeSubDomain,ps,&(ps.values()),"gadget2",_NneighborsMin,_NneighborsMax);			
	}
//	f.saveAllSlices("yinter-slices","yinter",2,0);

	msgs.say("done", Low);
	return f;
}
/***************************************************************************************/
void save3NNvalues(cpedsMsgs& msgs) {
	msgs.say(string("saving data to txt file")+_outfilePrefix,Medium);
	FILE* fout=fopen(_outfilePrefix.c_str(), "w");
	string fname;
	char fnametmp[100];

	lssPyramidSliceFileQuadruplet_t* data;
	long Npart;
	
	for (long i =_from; i < _to; i++) {
		sprintf(fnametmp,_infileFormat.c_str(),i);
		printf("processing file: %s\n", fnametmp);
		fname=fnametmp;
		lssKitLoad3NNData(fname,&data,&Npart);
		for (long i = 0; i < Npart; i++) {
			fprintf(fout,"%lE\n",data[i].T);
		}
		delete [] data;
	}
	fclose(fout);
	msgs.say("done",Low);
}
/***************************************************************************************/
void performInSparseGrid(cpedsMsgs& msgs) {
	_halos.loadMatrix(_haloFileName);
	mscsFunction3dregc f("field");
	f.setVerbosityLevel(High);
	string datasetName;
	double isInside;
	double t0=msgs.timeElapsed();
	double speed;
	double halosDone=0;
	long jStart=0;
	if (_haloInterpolationStartWithID!=-1) {
		for (long j = 0; j < _halos.Ny(); j++) {
			_currentHaloID=_halos.fRe(HALO_INFO_FILE_COL_HaloID,j,0);
			if (_currentHaloID!=_haloInterpolationStartWithID) {
				jStart++;	
				halosDone++;
			}
			else break;
		}
	}

	for (long j = jStart; j < _halos.Ny(); j++) {
		
		_currentHaloIdx=j;
		_currentHaloID=_halos.fRe(HALO_INFO_FILE_COL_HaloID,j,0);
		_currentDomainID=_halos.fRe(HALO_INFO_FILE_COL_DomainID,j,0);
		isInside=_halos.fRe(HALO_INFO_FILE_COL_isInsideFlag,j,0);

//		printf("current domain ID: %li\n",_currentDomainID);
//		printf("requested domain ID: %li\n",_requestedDomainID);
//		printf("_currentHaloIdx: %li\n",_currentHaloIdx);
//		printf("_currentHaloID: %li\n",_currentHaloID);
//		printf("isInside: %lE\n",isInside);
//		printf("\n");
		/*
		 * Comment: isInside should be replaced with a flag isPartiallyInside which should be implemented in the calculate-halo-properties
		 * program to also calculate halos that only partially overlap with the FOV.
		 * 
		 * author: blew
		 * date: Oct 20, 2013 2:42:26 AM
		 *
		 */
		
		if (_processAllHalos) isInside=1;
		
		if (_currentDomainID==_requestedDomainID and isInside==1) {
			msgs.say("Performing interpolation for halo: %li of %li",long(j),_halos.Ny(),Top);
			msgs.say("halo ID: %li",long(_currentHaloID),High);
			msgs.say("current domain: %li",long(_currentDomainID),High);
			setupGlobalVariables(msgs);
			
			if (_interpolate) makeSparseInterpolation(f,msgs); 
			else {
				msgs.criticalError("performing something else than 'interpolation' calculation in a sparse grid mode - this is probably not what you wanted. "
						"It is not implemented anyway. Please check the code",Top);
			}
			msgs.say("Saving to file", High);
			datasetName="halo_"+msgs.toStr(_currentHaloID);
			_outDatasetName=datasetName;
			saveToHDF(f);
			
			halosDone++;
			speed=halosDone/(msgs.timeElapsed()-t0);
			msgs.say("Interpolation speed [halos/min]: %lf",speed*60, High);
//			exit(0);
		}
		else {
			msgs.say("Skipping interpolation for halo: %li of %li",long(j),_halos.Ny(),Medium);
			msgs.say("current domain: %li, requested domain: %li",long(_currentDomainID),long(_requestedDomainID),High);
		}
	}
	
	// clean up
	_currentHaloIdx=-2;
	makeSparseInterpolation(f,msgs);
}
/***************************************************************************************/
void saveToHDF(mscsFunction3dregc& f) {
	f.saveHDF5(_outfilePrefix,_outDatasetName);
	// slice position in the deep field
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SX",_h.SX,"slice position in the deep field in x direction [comoving Mpc]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SY",_h.SY,"slice position in the deep field in y direction [comoving Mpc]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"SZ",_h.SZ,"slice position in the deep field in z direction [comoving Mpc]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdx",_h.Sdx,"slice size in x direction [comoving Mpc]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdy",_h.Sdy,"slice size in y direction [comoving Mpc]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Sdz",_h.Sdz,"slice size in z direction [comoving Mpc]");
	// cosmology parameters
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"slice_redshift",_h.z);
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"halo_redshift",_halos.fRe(17,_currentHaloIdx,0));
//		f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"X",cpeds_comoving_distance(_h.Wb0,_h.Wm0,_h.Wl0,_h.z,_h.w0,_h.h),"comoving distance to slice redshift";
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"hubble",_h.h,"unitless Hubble parameter");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"Wb0",_h.Wb0,"baryon density today in units of critical density");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"compton_y0",_h.y0,"integral constant: y0 = sigma_T kB rhoC0(h=1) Wb0_WMAP / ( mu_e_BBN m_p * m_e c^2 ) * CPEDS_MPC  = 1.50582*10^-16 [ K^-1 Mpc^-1 ]");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"mu_e",_h.mu_e, "electron number per proton mass");
	f.setHDF5_scalarDoubleAttribute(_outfilePrefix,_outDatasetName,"haloSizeScaleFactor",_haloSizeScaleFactor, "halo linear comoving size scale factor used to define function space for interpolation");
	f.setHDF5_scalarStringAttribute(_outfilePrefix,_outDatasetName,"haloFileName",_haloFileName);

}
/***************************************************************************************/
void makeSparseInterpolation(mscsFunction3dregc& f, cpedsMsgs& msgs) {
	static mscsVector<double> hsml;
	static mscsVector<double> dens;
	static cpedsPointSet3D ps;
	static subDomain *D=NULL;
	static subDomain *D3=NULL;
		
	if (_currentHaloIdx==-2) {
		msgs.say("cleaning up interpolation data",Low);
		delete D;
	}
	else {
		msgs.say("region: xmin: %lf, xmax: %lf",_xmin,_xmax,Medium);
		msgs.say("region: ymin: %lf, ymax: %lf",_ymin,_ymax,Medium);
		msgs.say("region: zmin: %lf, zmax: %lf",_zmin,_zmax,Medium);
		msgs.say("halo will be resolved with: Nx: %li, Ny: %li, Nz: %li grid cells",_Nx,_Ny,_Nz,Medium);
		
		if (_Nx==0) { msgs.warning("Nx was detected to be zero. Changing to 1", Top);	_Nx=1;	}
		if (_Ny==0) { msgs.warning("Ny was detected to be zero. Changing to 1", Top);	_Ny=1;	}
		if (_Nz==0) { msgs.warning("Nz was detected to be zero. Changing to 1", Top);	_Ny=1;	}
		if (_xmin==_xmax) { msgs.criticalError("halo linear size along X was detected to be zero. Stopping", Top);	}
		if (_ymin==_ymax) { msgs.criticalError("halo linear size along Y was detected to be zero. Stopping", Top);	}
		if (_zmin==_zmax) { msgs.criticalError("halo linear size along Z was detected to be zero. Stopping", Top);	}

		subDomain_region_t reg;
		reg.xmin=_xmin;	reg.ymin=_ymin; 	reg.zmin=_zmin;
		reg.xmax=_xmax;	reg.ymax=_ymax; 	reg.zmax=_zmax;
		reg.subx=_Nx; reg.suby=_Ny; reg.subz=_Nz;


		if (_loadedDomainID!=_requestedDomainID or _currentHaloIdx==0 or _inFormatHDF5) {
			// load domain data
			if (_inFormat3NN)
				lss3NNloadPointsSet(ps,msgs,&hsml,&dens);
			if (_inFormatHDF5)
				lssHDF5loadPointsSet(ps,msgs,&hsml,&dens,_currentHaloID);
			_loadedDomainID=_requestedDomainID;
//			exit(0);
			//
			// prepare tree
			//
			long minPointsInDomain=2;
			if (D!=NULL) delete D;
			D=new subDomain(&ps,&_treeSubDomain,minPointsInDomain); 
			msgs.say("building tree",Medium);
			D->tree();

			if (_twoTrees) {
				if (D3!=NULL) delete D3;
				D3=new subDomain(&ps,&_treeSubDomain,minPointsInDomain); 
				msgs.say("building tree 3",Medium);
				D3->tree();				
			}
		}
		//	f.mkInterpolatedFieldScatter(reg,_treeSubDomain,ps,ps.values(),hsml,dens,"gadget2");
//		ps.values().printVector();
//		f.mkInterpolatedFieldScatter(reg,D,ps,ps.values(),hsml,dens,_maxHSML, "gadget2",D3);		
		f.mkInterpolatedFieldScatter(reg,_treeSubDomain,ps,ps.values(),"gadget2",30,33,&hsml);
//		double interpErr=0;
//		for (unsigned long i = 0; i < ps.size(); i++) {
//			interpErr=abs(ps.val(i)-f.fReCoordPeriodic(ps[i].x(),ps[i].y(),ps[i].z()));
//			printf("ps val: %lE, fRe val: %lE,  interpErr: %lE\n",ps.val(i), f.fReCoordPeriodic(ps[i].x(),ps[i].y(),ps[i].z()), interpErr);
//		}
	}
	
}
