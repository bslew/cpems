/*!
 * \file mscsFunction3dCLEAN.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: blew
 */

#include "mscsFunction3dCLEAN.h"
#include "cpeds-direction.h"

mscsFunction3dCLEAN::mscsFunction3dCLEAN() {
	_mscsFunction3dCLEAN_data.MaxIter=5000;
	_mscsFunction3dCLEAN_data.loopGain=0.001;
	_mscsFunction3dCLEAN_data.beamOptimizationThres=0.01;
}

mscsFunction3dCLEAN::~mscsFunction3dCLEAN() {
}
/***************************************************************************************/
void mscsFunction3dCLEAN::initiate(mscsFunction3dregc& data, mscsFunction3dregc& psf, long MaxIter, double loopGain, string algorithm) {
	_mscsFunction3dCLEAN_data.MaxIter=MaxIter;
	_mscsFunction3dCLEAN_data.loopGain=loopGain;
	_mscsFunction3dCLEAN_data.data=data;
	_mscsFunction3dCLEAN_data.psf=psf;	
	_mscsFunction3dCLEAN_data.psf/=psf.getMaxValue();
	
	_mscsFunction3dCLEAN_data.algorithm=algorithm;
//	if (algorithm=="hogbom") _mscsFunction3dCLEAN_data.fitBeam=false;
//	if (algorithm=="fit") _mscsFunction3dCLEAN_data.fitBeam=true;
	
}
/***************************************************************************************/
void mscsFunction3dCLEAN::clean() {
	msgs->say("START CLEANING",High);
	bool stop=false;
	long ii,jj,kk;
	long iter=0;
	mscsFunction3dregc& data=_mscsFunction3dCLEAN_data.data;
	mscsFunction conv;

	
	// optimize the size of PSF
	_mscsFunction3dCLEAN_data.b=optimizeDDPSFsize(_mscsFunction3dCLEAN_data.psf,_mscsFunction3dCLEAN_data.beamOptimizationThres);
	rx0=(_mscsFunction3dCLEAN_data.b.getMinX()+_mscsFunction3dCLEAN_data.b.getMaxX())/2;
	ry0=(_mscsFunction3dCLEAN_data.b.getMinY()+_mscsFunction3dCLEAN_data.b.getMaxY())/2;
	_mscsFunction3dCLEAN_data.b.saveSlice(2,0,"b");
	_mscsFunction3dCLEAN_data.psf.saveSlice(2,0,"input_b");
	
	
	// estimate the initial stop condition threshold
	double stopThres=sqrt(data.varianceRe());
	msgs->say("stopThres: %lE",stopThres,Medium);
	
	data.saveSlice(2,0,"input_map.mat");

	
	while (iter<_mscsFunction3dCLEAN_data.MaxIter and stop==false) {
		msgs->say("iteration: %li",iter,Medium);
		
		// find maximal value in the data
		vMax=data.getMaxValue(&id);
		data.idx2ijk(id,iMax,jMax,kMax);
		x=data.getX(iMax);
		y=data.getY(jMax);
		msgs->say("found maximal map value: %lE",vMax,Medium);
		msgs->say("at: %lE, %lE",x,y,Low);
		msgs->say("idx: %li, %li",iMax,jMax,Low);
		
		// check if the source should be cleaned
		if (vMax > stopThres) {
			//
			// generate DD beam
			//
			if (_mscsFunction3dCLEAN_data.fitBeam) {
				fitBeamOrientation();
				fitBeamAmplitude();
				_mscsFunction3dCLEAN_data.br.saveSlice(2,0,"brfit");
				g=1;
				stop=true;
			}
			else {
				angle=data.fIm(id);
	//			angle=0*PI180;
				msgs->say("Forming beam",Low);
				_mscsFunction3dCLEAN_data.br=_mscsFunction3dCLEAN_data.b.rotateSlice(2,0,angle,rx0,ry0,0);
				vMaxB=_mscsFunction3dCLEAN_data.br.getMaxValue(&id);
				_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMaxB,jMaxB,kMaxB);
				vMinB=_mscsFunction3dCLEAN_data.br.getMinValue(&id);
				_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMinB,jMinB,kMinB);
				iBoff=(iMinB-iMaxB)/2;
				jBoff=(jMinB-jMaxB)/2;

				g=_mscsFunction3dCLEAN_data.loopGain;
				_mscsFunction3dCLEAN_data.br*=g;
				dI=vMax*g;
				msgs->say("Loop gain is: %lE",g,Low);
				msgs->say("dI=vMax * loog gain: %lE",dI,Low);
				msgs->say("iBoff: %li, jBoff: %li",iBoff,jBoff,Low);
			}

			msgs->say("Subtracting source",Low);
			// position the beam at the source and subtract its fraction according to the loop gain
//			_mscsFunction3dCLEAN_data.br.saveSlice(2,0,"br");
//			exit(0);
			for (long i = 0; i < _mscsFunction3dCLEAN_data.br.Nx(); i++) {
				ii=iMax-iMaxB+i;
				if (ii>=0 and ii<data.Nx()) {
					for (long j = 0; j < _mscsFunction3dCLEAN_data.br.Ny(); j++) {
						jj=jMax-jMaxB+j;
						if (jj>=0 and jj<data.Ny()) {
							data.fRe(ii,jj,0)-=_mscsFunction3dCLEAN_data.br.fRe(i,j,0)*dI;
						}
					}
				}
			}
			
			
			// store cleaned source data
			_mscsFunction3dCLEAN_data.radec.append(cpedsDirection(x+iBoff*_mscsFunction3dCLEAN_data.br.getDx(),y+jBoff*_mscsFunction3dCLEAN_data.br.getDy(),dI));
			_mscsFunction3dCLEAN_data.srcij.append(cpedsDirection(iMax+iBoff,jMax+jBoff,dI));
			
			
			// check the convergence
			rms=sqrt(data.varianceRe());
			conv.newPoint(double(iter),rms);
			msgs->say("rms: %lE",rms,Low);
//			if (vMax < rms) stop=true;
			stopThres=1.5*rms;
		}
		else {
			stop=true;
		}
		
		iter++;
		
	}
	_mscsFunction3dCLEAN_data.br.saveSlice(2,0,"br");
	
	// save residual map
	_mscsFunction3dCLEAN_data.resid=data;
	_mscsFunction3dCLEAN_data.resid.saveSlice(2,0,"residual_noise.mat");
	//save convergence track
	conv.save("convergence");
	//save source positions
	_mscsFunction3dCLEAN_data.radec.save("positions","",true);
	_mscsFunction3dCLEAN_data.srcij.save("positionsij","",true);
	
	
	// calculate clean beam
	_mscsFunction3dCLEAN_data.bc=calculateCleanBeam(_mscsFunction3dCLEAN_data.b);
	_mscsFunction3dCLEAN_data.bc.saveSlice(2,0,"clean_beam.mat");
	
	// construct the clean map
	_mscsFunction3dCLEAN_data.src=_mscsFunction3dCLEAN_data.data;
	_mscsFunction3dCLEAN_data.src=0.0;
	for (long i = 0; i < _mscsFunction3dCLEAN_data.radec.size(); ++i) {
		_mscsFunction3dCLEAN_data.src.fRe(_mscsFunction3dCLEAN_data.srcij[i].lon(),_mscsFunction3dCLEAN_data.srcij[i].lat(),0)+=_mscsFunction3dCLEAN_data.srcij[i].val()*g;
	}
	_mscsFunction3dCLEAN_data.src.saveSlice(2,0,"src");
	
	// convolve with beam
	_mscsFunction3dCLEAN_data.psfclean=_mscsFunction3dCLEAN_data.data;
	_mscsFunction3dCLEAN_data.psfclean=0.0;
	_mscsFunction3dCLEAN_data.psfclean.pasteAdd(_mscsFunction3dCLEAN_data.bc,0,0,0);
	_mscsFunction3dCLEAN_data.psfclean.centerMax();
//	_mscsFunction3dCLEAN_data.psfclean.shift(_mscsFunction3dCLEAN_data.psfclean.Nx(),_mscsFunction3dCLEAN_data.psfclean.Ny(),0);
//	_mscsFunction3dCLEAN_data.psfclean/=_mscsFunction3dCLEAN_data.psfclean.integrateRe();
//	_mscsFunction3dCLEAN_data.srcbcn=_mscsFunction3dCLEAN_data.src.convolve(_mscsFunction3dCLEAN_data.psfclean);
	mscsFunction3dregc ker=_mscsFunction3dCLEAN_data.psfclean;
//	ker/=ker.integrateRe();
	ker.shift(ker.Nx()/2,ker.Ny()/2,0);
	ker*=ker.size();
	ker.fft(true);
	_mscsFunction3dCLEAN_data.srcbcn=_mscsFunction3dCLEAN_data.src.convolve_fft(ker);
	_mscsFunction3dCLEAN_data.srcbcn.saveSlice(2,0,"src_beam.mat");
	_mscsFunction3dCLEAN_data.srcbcn+=_mscsFunction3dCLEAN_data.resid;
	_mscsFunction3dCLEAN_data.srcbcn.saveSlice(2,0,"src_beam_noise.mat");
	
	_mscsFunction3dCLEAN_data.srcresid=_mscsFunction3dCLEAN_data.src+_mscsFunction3dCLEAN_data.resid;
	_mscsFunction3dCLEAN_data.srcresid.saveSlice(2,0,"src_noise.mat");
	
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dCLEAN::optimizeDDPSFsize(mscsFunction3dregc psf, double thres) {

	psf/=psf.getMaxValue();
//	psf.saveSlice(2,0,"psf",0);
	long iMax,jMax,kMax,id;
	long iMin,jMin,kMin;
	long ii,jj;
	psf.getMaxValue(&id);
	psf.idx2ijk(id,iMax,jMax,kMax);
	psf.getMinValue(&id);
	psf.idx2ijk(id,iMin,jMin,kMin);
//	cpedsPoint3D ijMax=psf.getMaxValueCell(true);
//	cpedsPoint3D ijMin=psf.getMinValueCell(true);
//	cpedsPoint3D ctr=(ijMax+ijMin);
//	ctr/=2.0;
//	ijMax.x()/=psf.getDx();
//	ijMax.y()/=psf.getDy();
//	ctr.x()/=psf.getDx();
//	ctr.y()/=psf.getDy();
//	printf("%lE\n",psf.getDx());
	
	ii=(iMax+iMin)/2;
	jj=(jMax+jMin)/2;
	
	long dir=1;
	if (iMax>ii) dir=1; else dir=-1;
//	ijMax.print_point("ijMax");
//	ijMin.print_point("ijMin");
	
	long i=iMax;
	long j=jMax;
	long k=0;
	bool cont=true;
	while (i<psf.Nx()-1 and i>0 and cont) {
		if (psf.fRe(i,j,k)>thres) i+=dir;
		else cont=false;
	}
	long iSt,iEn,jSt,jEn;
	if (dir==1) { iSt=2*ii-i; iEn=i; } else { iSt=i; iEn=2*ii-i; }
	iSt=cpeds_get_max(iSt,0);
	iEn=cpeds_get_min(iEn,psf.Nx()-1);
	jSt=iSt;
	jEn=iEn;

	if (iSt==iEn or jSt==jEn) {
		printf("dir: %li, iSt: %li, iEn: %li\n",dir,iSt,iEn);
		msgs->criticalError("optimizeDDPSFsize: ist==iEn or jst==jEn, this is weird - check this out",Top);
		
	}
	if (iSt>iEn) {
		printf("dir: %li, iSt: %li, iEn: %li\n",dir,iSt,iEn);
		msgs->criticalError("optimizeDDPSFsize: ist>iEn this is weird - check this out",Top);
	}
	return psf.cutAwayBlock(iSt,jSt,0,iEn,jEn,0);
//	_mscsFunction3dCLEAN_data.b.saveSlice(2,0,"psf-optimized",0);
}
/***************************************************************************************/
void mscsFunction3dCLEAN::rotatePSF(double angle) { 
	rotateSlice(2,0,angle,0,0,0);
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dCLEAN::calculateCleanBeam(mscsFunction3dregc& psf)  {
	msgs->say("Calculating clean beam",High);
	mscsFunction3dregc c=psf;
	mscsFunction3dregc g;
	
	for (long i = 0; i < c.Nx(); i++) {
		for (long j = 0; j < c.Ny(); j++) {
			if (c.fRe(i,j,0)<0) c.fRe(i,j,0)=0;
		}
	}

//	c.saveSlice(2,0,"nulled");
	c.centerMax();
//	c.saveSlice(2,0,"center");
	// now the beam is in the center of the field
//	c.setVerbosityLevel(Top);
//	c.printInfo();
//	c=optimizeDDPSFsize(c,_mscsFunction3dCLEAN_data.beamOptimizationThres);
	g=c;
//	c.printInfo();
//	c.saveSlice(2,0,"clean-optimized");
//	exit(0);
	// find best fitting gaussian beam
	
	mscsFunction beamSize;
	double bs,dbs,bsMin,bsMax,s;
	double chisq;
	
	dbs=c.lengthX()/Nsteps;
	bsMin=dbs;
	bsMax=c.lengthX();
	
//	printf("dbs: %lE\n",dbs);
//	printf("bsMax: %lE\n",bsMax);

	
	bs=bsMin;
	while (bs<bsMax) {
//		s=bs;
		s=cpeds_fwhm2sigma(bs);
		g.mkGauss2D(g.getMinX()+g.lengthX()/2,g.getMinY()+g.lengthY()/2,s,s,1,2,0);
//		g=0.0;
//		g.setf(g.Nx()/2,g.Ny()/2,0,1,0);
//		g.smooth3DGauss(s,s,1);
//		g.saveSlice(2,0,"iter1");
//		if (g.getMaxValue()!=0)	g/=g.getMaxValue();
		
		chisq=0;
		for (long i = 0; i < g.Nx(); i++) {
			for (long j = 0; j < g.Ny(); j++) {
				chisq+=pow(g.fRe(i,j,0)-c.fRe(i,j,0),2);
			}
		}
//		printf("chisq for bs=%lE: %lE\n",bs,chisq);
		beamSize.newPoint(bs,chisq);
		bs+=dbs;
	}
	
	beamSize.save("beamSizes");
	beamSize.checkRanges();
	bs=beamSize.getX(beamSize.getMinValueIdx());
	msgs->say("best fit beam size: %lE",bs,Low);
	
	// form clean beam
	s=cpeds_fwhm2sigma(bs);
	g.mkGauss2D(g.getMinX()+g.lengthX()/2,g.getMinY()+g.lengthY()/2,s,s,1,2,0);
	return g;
}
/***************************************************************************************/
//void mscsFunction3dCLEAN::clean2() {
//	msgs->say("START CLEANING",High);
//	bool stop=false;
//	double vMax,vMaxB,vMinB;
//	double rx0,ry0,x,y; 
//	double g,dI, angle,rms;
//	long iMax,jMax,kMax,id;
//	long iMaxB,jMaxB,kMaxB;
//	long iMinB,jMinB,kMinB;
//	long iBoff,jBoff; // indexes offsets of the center of the rotated beam - pointing the place where the source should be - from the positive beam location indexes
//	long ii,jj,kk;
//	long iter=0;
//	mscsFunction3dregc& data=_mscsFunction3dCLEAN_data.data;
//	mscsFunction conv;
//
//	
//	// optimize the size of PSF
//	_mscsFunction3dCLEAN_data.b=optimizeDDPSFsize(_mscsFunction3dCLEAN_data.psf,_mscsFunction3dCLEAN_data.beamOptimizationThres);
//	rx0=(_mscsFunction3dCLEAN_data.b.getMinX()+_mscsFunction3dCLEAN_data.b.getMaxX())/2;
//	ry0=(_mscsFunction3dCLEAN_data.b.getMinY()+_mscsFunction3dCLEAN_data.b.getMaxY())/2;
//	_mscsFunction3dCLEAN_data.b.saveSlice(2,0,"b");
//	_mscsFunction3dCLEAN_data.psf.saveSlice(2,0,"input_b");
//	
//	
//	// estimate the initial stop condition threshold
//	double stopThres=sqrt(data.varianceRe());
//	msgs->say("stopThres: %lE",stopThres,Medium);
//	
//	data.saveSlice(2,0,"input_map.mat");
//
//	
//	while (iter<_mscsFunction3dCLEAN_data.MaxIter and stop==false) {
//		msgs->say("iteration: %li",iter,Medium);
//		
//		// find maximal value in the data
//		vMax=data.getMaxValue(&id);
//		data.idx2ijk(id,iMax,jMax,kMax);
//		x=data.getX(iMax);
//		y=data.getY(jMax);
//		msgs->say("found maximal map value: %lE",vMax,Medium);
//		msgs->say("at: %lE, %lE",x,y,Low);
//		msgs->say("idx: %li, %li",iMax,jMax,Low);
//		
//		// check if the source should be cleaned
//		if (vMax > stopThres) {
//			fitBeamOrientation();
//			fitBeamAmplitude();
//			
//			g=_mscsFunction3dCLEAN_data.loopGain;
//			_mscsFunction3dCLEAN_data.br*=g;
//			dI=vMax*g;
//			msgs->say("Loop gain is: %lE",g,Low);
//			msgs->say("dI=vMax * loog gain: %lE",dI,Low);
//			msgs->say("iBoff: %li, jBoff: %li",iBoff,jBoff,Low);
//
//			msgs->say("Subtracting source",Low);
//			// position the beam at the source and subtract its fraction according to the loop gain
////			_mscsFunction3dCLEAN_data.br.saveSlice(2,0,"br");
////			exit(0);
//			for (long i = 0; i < _mscsFunction3dCLEAN_data.br.Nx(); i++) {
//				ii=iMax-iMaxB+i;
//				if (ii>=0 and ii<data.Nx()) {
//					for (long j = 0; j < _mscsFunction3dCLEAN_data.br.Ny(); j++) {
//						jj=jMax-jMaxB+j;
//						if (jj>=0 and jj<data.Ny()) {
//							data.fRe(ii,jj,0)-=_mscsFunction3dCLEAN_data.br.fRe(i,j,0)*dI;
//						}
//					}
//				}
//			}
//			
//			
//			// store cleaned source data
//			_mscsFunction3dCLEAN_data.radec.append(cpedsDirection(x+iBoff*_mscsFunction3dCLEAN_data.br.getDx(),y+jBoff*_mscsFunction3dCLEAN_data.br.getDy(),dI));
//			_mscsFunction3dCLEAN_data.srcij.append(cpedsDirection(iMax+iBoff,jMax+jBoff,dI));
//			
//			
//			// check the convergence
//			rms=sqrt(data.varianceRe());
//			conv.newPoint(double(iter),rms);
//			msgs->say("rms: %lE",rms,Low);
////			if (vMax < rms) stop=true;
//			stopThres=1.5*rms;
//		}
//		else {
//			stop=true;
//		}
//		
//		iter++;
//		
//	}
//	_mscsFunction3dCLEAN_data.br.saveSlice(2,0,"br");
//	
//	// save residual map
//	_mscsFunction3dCLEAN_data.resid=data;
//	_mscsFunction3dCLEAN_data.resid.saveSlice(2,0,"residual_noise.mat");
//	//save convergence track
//	conv.save("convergence");
//	//save source positions
//	_mscsFunction3dCLEAN_data.radec.save("positions","",true);
//	_mscsFunction3dCLEAN_data.srcij.save("positionsij","",true);
//	
//	
//	// calculate clean beam
//	_mscsFunction3dCLEAN_data.bc=calculateCleanBeam(_mscsFunction3dCLEAN_data.b);
//	_mscsFunction3dCLEAN_data.bc.saveSlice(2,0,"clean_beam.mat");
//	
//	// construct the clean map
//	_mscsFunction3dCLEAN_data.src=_mscsFunction3dCLEAN_data.data;
//	_mscsFunction3dCLEAN_data.src=0.0;
//	for (long i = 0; i < _mscsFunction3dCLEAN_data.radec.size(); ++i) {
//		_mscsFunction3dCLEAN_data.src.fRe(_mscsFunction3dCLEAN_data.srcij[i].lon(),_mscsFunction3dCLEAN_data.srcij[i].lat(),0)+=_mscsFunction3dCLEAN_data.srcij[i].val()*_mscsFunction3dCLEAN_data.loopGain;
//	}
//	_mscsFunction3dCLEAN_data.src.saveSlice(2,0,"src");
//	
//	// convolve with beam
//	_mscsFunction3dCLEAN_data.psfclean=_mscsFunction3dCLEAN_data.data;
//	_mscsFunction3dCLEAN_data.psfclean=0.0;
//	_mscsFunction3dCLEAN_data.psfclean.pasteAdd(_mscsFunction3dCLEAN_data.bc,0,0,0);
//	_mscsFunction3dCLEAN_data.psfclean.centerMax();
////	_mscsFunction3dCLEAN_data.psfclean.shift(_mscsFunction3dCLEAN_data.psfclean.Nx(),_mscsFunction3dCLEAN_data.psfclean.Ny(),0);
////	_mscsFunction3dCLEAN_data.psfclean/=_mscsFunction3dCLEAN_data.psfclean.integrateRe();
////	_mscsFunction3dCLEAN_data.srcbcn=_mscsFunction3dCLEAN_data.src.convolve(_mscsFunction3dCLEAN_data.psfclean);
//	mscsFunction3dregc ker=_mscsFunction3dCLEAN_data.psfclean;
//	ker.shift(ker.Nx()/2,ker.Ny()/2,0);
//	ker*=ker.size();
//	ker.fft(true);
//	_mscsFunction3dCLEAN_data.srcbcn=_mscsFunction3dCLEAN_data.src.convolve_fft(ker);
//	_mscsFunction3dCLEAN_data.srcbcn.saveSlice(2,0,"src_beam.mat");
//	_mscsFunction3dCLEAN_data.srcbcn+=_mscsFunction3dCLEAN_data.resid;
//	_mscsFunction3dCLEAN_data.srcbcn.saveSlice(2,0,"src_beam_noise.mat");
//	
//	_mscsFunction3dCLEAN_data.srcresid=_mscsFunction3dCLEAN_data.src+_mscsFunction3dCLEAN_data.resid;
//	_mscsFunction3dCLEAN_data.srcresid.saveSlice(2,0,"src_noise.mat");
//	
//}
/***************************************************************************************/
void mscsFunction3dCLEAN::fitBeamOrientation() {
	mscsFunction3dregc& data=_mscsFunction3dCLEAN_data.data;

	//
	// fit DD beam by rotation
	//
	msgs->say("Fitting beam orientation",Low);
	mscsFunction angCorrelation;
	for (long ang = 0; ang < 359; ang++) {
		// generate rotated beam
		angle=ang*PI180;
		_mscsFunction3dCLEAN_data.br=_mscsFunction3dCLEAN_data.b.rotateSlice(2,0,angle,rx0,ry0,0);
		vMaxB=_mscsFunction3dCLEAN_data.br.getMaxValue(&id);
		_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMaxB,jMaxB,kMaxB);
		vMinB=_mscsFunction3dCLEAN_data.br.getMinValue(&id);
		_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMinB,jMinB,kMinB);
//		iBoff=(iMinB-iMaxB)/2;
//		jBoff=(jMinB-jMaxB)/2;
		
		mscsFunction3dregc f=data.cutAwayBlock(iMax-iMaxB,jMax-jMaxB,0,iMax-iMaxB+_mscsFunction3dCLEAN_data.br.Nx()-1,jMax-jMaxB+_mscsFunction3dCLEAN_data.br.Ny()-1,0);
		if (f.Nx()!=_mscsFunction3dCLEAN_data.b.Nx()) { msgs->criticalError("fitBeamOrientation>> cutting away too close the X boundary, please increase the field margin",Top); }
		if (f.Ny()!=_mscsFunction3dCLEAN_data.b.Ny()) { msgs->criticalError("fitBeamOrientation>> cutting away too close the Y boundary, please increase the field margin",Top); }
		f.saveSlice(2,0,"fcut");

		f*=_mscsFunction3dCLEAN_data.br;
		
		
		angCorrelation.newPoint(angle,f.sumRe());
		
	}
	
	// find the optimal angle
	angCorrelation.checkRanges();
	angCorrelation.save("angCorrelation");
	angle=angCorrelation.getX(angCorrelation.getMaxValueIdx());
	_mscsFunction3dCLEAN_data.br=_mscsFunction3dCLEAN_data.b.rotateSlice(2,0,angle,rx0,ry0,0);
	vMaxB=_mscsFunction3dCLEAN_data.br.getMaxValue(&id);
	_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMaxB,jMaxB,kMaxB);
	vMinB=_mscsFunction3dCLEAN_data.br.getMinValue(&id);
	_mscsFunction3dCLEAN_data.br.idx2ijk(id,iMinB,jMinB,kMinB);
	iBoff=(iMinB-iMaxB)/2;
	jBoff=(jMinB-jMaxB)/2;
	
}

/***************************************************************************************/
void mscsFunction3dCLEAN::fitBeamAmplitude() {
	mscsFunction3dregc& data=_mscsFunction3dCLEAN_data.data;

	//
	// fit DD beam amplitude
	//
	msgs->say("Fitting beam amplitude",Low);
	_mscsFunction3dCLEAN_data.br/=_mscsFunction3dCLEAN_data.br.getMaxValue();
	double Imax=1.5*vMax;
	double I=_mscsFunction3dCLEAN_data.loopGain*vMax;
	double delta=I;
	mscsFunction3dregc br;
	mscsFunction3dregc f0=data.cutAwayBlock(iMax-iMaxB,jMax-jMaxB,0,iMax-iMaxB+_mscsFunction3dCLEAN_data.br.Nx()-1,jMax-jMaxB+_mscsFunction3dCLEAN_data.br.Ny()-1,0);
	f0.saveSlice(2,0,"ampdata");

	mscsFunction chisq;
	while (I<Imax)  {
		// generate rotated beam
		br=_mscsFunction3dCLEAN_data.br;
		br*=double(I);
		br-=f0;
		br*=br;
		chisq.newPoint(I,br.sumRe());
		I+=delta;
	}
	chisq.save("Ichisq");
	
	// find the optimal amplitude
	chisq.checkRanges();
	I=chisq.getX(chisq.getMinValueIdx());
	dI=I;

	br=_mscsFunction3dCLEAN_data.br;
	br*=double(I);
	br.saveSlice(2,0,"brfitamp");

}
