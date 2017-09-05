/*
 * mscsFunction3d.cpp
 *
 *  Created on: Dec 21, 2010, 3:31:32 PM
 *      Author: blew
 */
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <fftw3.h>
#include "Mscs-function2dreg.h"

mscsFunction2dreg::mscsFunction2dreg() : mscsFunction() {

}
mscsFunction2dreg::mscsFunction2dreg(string name) : mscsFunction(name) { }

mscsFunction2dreg::~mscsFunction2dreg() {
}

mscsFunction2dreg::mscsFunction2dreg(long Nx, long Ny, double dx,double dy) : mscsFunction() {
	_sizeX=Nx*dx;
	_sizeY=Ny*dy;
	_dx=dx;
	_dy=dy;
	_dxo2=_dx/2.0;
	_dyo2=_dy/2.0;
	_x0=0;
	_y0=0;
	_Nx=Nx;
	_Ny=Ny;
	
	// initiate the function domain
	double tmp=_dxo2;
	for (int i = 0; i < Nx; i++) {
		_X.append(tmp);
		tmp+=dx;
	}
	
	tmp=_dyo2;
	for (int i = 0; i < Ny; i++) {
		_Y.append(tmp);
		tmp+=dy;
	}
	

	
	// allocate space for the function
	for (int i = 0; i < _sizeX; i++) {
		for (int j = 0; j < _sizeY; j++) {
			newPoint();
		}
	}
	
}
/***************************************************************************************/
void mscsFunction2dreg::setf(long i, long j, double x, double y) {
	long idx=ij2idx(i,j);
	mscsFunction::setf(idx,x,y);
}

/***************************************************************************************/
long mscsFunction2dreg::ij2idx(long i, long j) const {
//	return k*_sizeXY+j*_sizeY+i;
	return j+i*_Ny;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction2dreg::load(string filename, bool commentedFile) {
	FILE* f;
	string format;
	double x,y;
	int res;
	msgs->say("Loading function from file "+filename,High);
	long cols=cpeds_get_cols_num_first_ln(filename.c_str(),commentedFile);
	msgs->say("Number of columns in the input file: "+msgs->toStr(cols),Medium);
	if (cols==1) { msgs->say("Assuming real-valued function.", Medium);	}
	if (cols==2) { msgs->say("Assuming complex-valued function.", Medium);	}

	switch (cols) {
		case 1:
			format="%lE";
			break;
		case 2:
			format="%lE %lE";

			break;
		default:
			break;
	}
//	for (int i = 2; i < cols; ++i) {	format+=" %*E";  }
	if (commentedFile) {
		ifstream ifs(filename.c_str());
		if (ifs.is_open()) {
			string s;
			switch (cols) {
				case 1:
					while (getline(ifs,s)) {
						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x); newPoint(x,0); }
					}
					break;
				case 2:
					while (getline(ifs,s)) {
						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x,&y); newPoint(x,y); }
					}

					break;
			}
		}
		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	}
	else {
		f= fopen(filename.c_str(),"r");
		if (f!=NULL) {
			switch (cols) {
				case 1:
					while (!feof(f)) { res=fscanf(f,format.c_str(),&x); if (res!=EOF) newPoint(x,0); }
					break;
				case 2:
					while (!feof(f)) { res=fscanf(f,format.c_str(),&x,&y); if (res!=EOF) newPoint(x,y); }

					break;
			}
			fclose(f);
			msgs->say("Read "+msgs->toStr((long)_f.count())+" lines",Medium);
		}
		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	}
	
		
	long n=pointsCount();
	msgs->say("Loaded "+msgs->toStr(pointsCount())+" points.",Medium);
	if (n!=getN()) {
		msgs->warning("The expected number of loaded points is different from the declared grid size: "+msgs->toStr(getN()),High);
	}
	return cpedsSuccess;
	
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::generateRandomGaussianField(double m, double s, long seed) {
	cpedsRNG *rns = new cpedsRNG("gaussian_circle");
	rns->seed(seed);
		
	double *a= rns->getRN(getN());
	delete rns;
	importFunction(a,NULL,getN());
	delete [] a;
	return *this;
}


/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::generateRandomGaussianFField(const mscsFunction& Pk) {		
	return *this;
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::generateRandomGaussianField(const mscsFunction& Pk) {		
	return *this;
}

/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::fft_r2c(fftw_complex** dk) {
	double* in=extractArguments();

	long M=Nx()*(long(Ny()/2)+1);
	fftw_complex* out=(fftw_complex*)fftw_malloc((M)*sizeof(fftw_complex)); 
	
	fftw_plan p;
	p=fftw_plan_dft_r2c_2d(Nx(),Ny(), in,out, FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_free(in);
	*dk=&out[0];
	return *this;
	
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::fft_c2c(fftw_complex** dk,bool fwd) {
//	fftw_complex* in= new fftw_complex[getN()];
	fftw_complex* in=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	long n=getN();
	printf("%li %li\n",n,pointsCount());
	for (long i = 0; i < n; i++) {
		in[i][0]=mscsFunction::getx(i)/n; // normalization to correct the fundamental term to be average
		in[i][1]=mscsFunction::Y(i)/n;
	}
//	exit(0);
	fftw_complex* out=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	
	fftw_plan p;
	if (fwd)
		p=fftw_plan_dft_2d(Nx(),Ny(), in,out, FFTW_FORWARD,FFTW_ESTIMATE); 
	else
		p=fftw_plan_dft_2d(Nx(),Ny(), in,out, FFTW_BACKWARD,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_free(in);
	*dk=&out[0];
	return *this;
	
}
/***************************************************************************************/
mscsFunction2dreg mscsFunction2dreg::slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky) {
//	fftw_complex* in=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	long nx=Nx();
	long ny=Ny();
	long nkx=kx.size();
	long nky=ky.size();
	mscsFunction2dreg g("FT");
	double re,im;
	
	for (long ki = 0; ki < nkx; ki++) {
		for (long kj = 0; kj < nky; kj++) {
				re=im=0;
				
				for (long i = 0; i < nx; i++) {
					for (long j = 0; j < ny; j++) {
							re+= fRe(i,j)*cos(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j)  )   );
							im+=-fIm(i,j)*sin(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j)  )   );
					}
				}
				g.newPoint(re,im);
		}
	}
	
	g._X=kx;
	g._Y=ky;
	g.setSize(nx,ny);
	
	return g;
}

/***************************************************************************************/
mscsFunction mscsFunction2dreg::powerSpectrum() {
	fftw_complex* dk;
	mscsFunction Pk("power spectrum");
	cpedsList<double> Pveck;

	// derive the fourier modes
//	fft_r2c(&dk);
	fft_c2c(&dk);
	
	// derive the power spectrum
	
	long nx=Nx();
	long ny=Ny();
//	long nz=long(Nz()/2)+1;

//	long nx=int(Nx()/2)+1;
//	long ny=int(Ny()/2)+1;

	long n=getN();

	printf("nz %li ny %li\n",nx,ny);
	double kx,ky, modk;
	
	
	//
	// test
	//
	double* tmpa=new double[getN()];
	for (long i = 0; i < getN(); i++) {	tmpa[i]=dk[i][0];	}
	cpeds_save_matrix(tmpa,Nx(),Ny(),"fftre");
	for (long i = 0; i < getN(); i++) {	tmpa[i]=dk[i][1];	}
	cpeds_save_matrix(tmpa,Nx(),Ny(),"fftim");
	for (long i = 0; i < getN(); i++) {	tmpa[i]=dk[i][0]*dk[i][0]+dk[i][1]*dk[i][1];	}
	cpeds_save_matrix(tmpa,Nx(),Ny(),"fftpow");
	
	
	
	//
	// calculate power P(vec k) 
	//
	long l=0;
	double tmp;
//	
//	for (long i = 0; i < nx; i++) {
//		for (long j = 0; j < ny; j++) {
//				tmp=(dk[l][0]*dk[l][0] + dk[l][1]*dk[l][1]); tmp*=tmp;
//				Pveck.append(tmp);
//				l++;
//		}
//	}
//	fftw_free(dk);
//	Pveck.save("fft");
	
	//
	// calculate power P(k) = integrate sqrt(kx^2+ky^2) f(vec k) dkx dky
	//
	l=0;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
				// derive k=sqrt(kx^2+ky^2+kz^2)
				kx=double(i)/double(_sizeX);
				ky=double(j)/double(_sizeY);
				modk=sqrt(kx*kx + ky*ky);
				tmp=(dk[l][0]*dk[l][0] + dk[l][1]*dk[l][1]);
				// derive the power
				Pk.newPoint(modk,tmp);
				l++;
		}
	}
//	Pk/=double(_sizeX*_sizeY);
	
	Pk.sortFunctionArgAscending();
	return Pk;	
}

/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::mkSin2Drad(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y,z;
	long i,j,k;

	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			newPoint(sin(twoPI/T*sqrt(x*x + y*y )),0);
			//				newPoint(x*x + y*y + z*z,0);
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=from;
	_dx=_dy=dr;
	_dxo2=_dyo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
		
	return *this;	
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::mkSin2D(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y;
	long i,j;

	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			newPoint(sin(twoPI/T*(x + y)),0);
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=from;
	_dx=_dy=dr;
	_dxo2=_dyo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
	
	return *this;	
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::mkSin2Dx(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y;
	long i,j;

	x=from; i=0;
	while (x<to) {
		y=from; j=0;
		while (y<to) {
			newPoint(sin(twoPI/T*(x)),0);
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=from;
	_dx=_dy=dr;
	_dxo2=_dyo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
	
	return *this;	
}
/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::mkSin2Dy(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y;
	long i,j;

	x=from; i=0;
	while (x<to) {
		y=from; j=0;
		while (y<to) {
			newPoint(sin(twoPI/T*(y)),0);
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=from;
	_dx=_dy=dr;
	_dxo2=_dyo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
	
	return *this;	
}

/***************************************************************************************/
void mscsFunction2dreg::setSize(long Nx, long Ny, double dx, double dy, double x0, double y0) {
	setSize(Nx,Ny);
	_dx=dx; _dxo2=dx/2.0;
	_dy=dy; _dyo2=dy/2.0;
	_x0=x0;
	_y0=y0;
}
/***************************************************************************************/
void mscsFunction2dreg::setSize(long Nx, long Ny) {
	_sizeX=Nx;
	_sizeY=Ny;
}


/***************************************************************************************/
void mscsFunction2dreg::printInfo() const {
	printf("sizeX sizeY: %lE %lE \n",_sizeX,_sizeY);
	printf("Nx Ny: %li %li\n",_Nx,_Ny);	
	printf("dx dy: %lE %lE\n",_dx,_dy);	
	printf("x0 y0: %lE %lE\n",_x0,_y0);	
}


/***************************************************************************************/
cpedsStatusCodes mscsFunction2dreg::saveAsMatrix(string filename) {
	double *tmp=extractArguments();
	cpeds_save_matrix(tmp,Nx(),Ny(),filename+".re");
	delete tmp;
	tmp=extractValues();
	cpeds_save_matrix(tmp,Nx(),Ny(),filename+".im");
	delete tmp;	
	return cpedsSuccess;
}


/***************************************************************************************/
void mscsFunction2dreg::antysymmetrize() {
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j <= Ny()/2; j++) {
			setf((Nx()-i) % Nx(),(Ny()-j) % Ny(),fRe(i,j),-fIm(i,j));
			if ((Nx()-i) % Nx()==i and (Ny()-j) % Ny()==j) { fIm(i,j)=0.0; }
		}
	}
//	mscsFunction::f(0)=0.0;
}

/***************************************************************************************/
mscsFunction2dreg& mscsFunction2dreg::fft(bool fwd) {
	fftw_complex* dk;
	fft_c2c(&dk,fwd);

	long n=getN();
	for (long i = 0; i < n; i++) {
		mscsFunction::setf(i,dk[i][0],dk[i][1]);
	}
//	if (fwd==false) {
//		divide(n);
//		scaleX(1.0/n);
//	}
	fftw_free(dk);
	return *this;
}










