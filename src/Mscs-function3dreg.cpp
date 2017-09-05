/*
 * mscsFunction3d.cpp
 *
 *  Created on: Sep 27, 2010
 *      Author: blew
 */
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <fftw3.h>
#include "Mscs-function3dreg.h"

mscsFunction3dreg::mscsFunction3dreg() : mscsFunction() {

}
mscsFunction3dreg::mscsFunction3dreg(string name) : mscsFunction(name) { }

mscsFunction3dreg::~mscsFunction3dreg() {
}

mscsFunction3dreg::mscsFunction3dreg(long Nx, long Ny, long Nz, double dx,double dy, double dz) : mscsFunction() {
	_sizeX=Nx*dx;
	_sizeY=Ny*dy;
	_sizeZ=Nz*dz;
	_dx=dx;
	_dy=dy;
	_dz=dz;
	_dxo2=_dx/2.0;
	_dyo2=_dy/2.0;
	_dzo2=_dz/2.0;
	_x0=0;
	_y0=0;
	_z0=0;
	_Nx=Nx;
	_Ny=Ny;
	_Nz=Nz;
	_NXY=Nx*Ny;
	_NYZ=Ny*Nz;
	
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
	
	tmp=_dzo2;
	for (int i = 0; i < Nz; i++) {
		_Z.append(tmp);
		tmp+=dz;
	}

	
	// allocate space for the function
	for (int i = 0; i < _sizeX; i++) {
		for (int j = 0; j < _sizeY; j++) {
			for (int k = 0; k < _sizeZ; k++) {
				newPoint();
			}
		}
	}
	
}
/***************************************************************************************/
void mscsFunction3dreg::setf(long i, long j, long k, double x, double y) {
	long idx=ijk2idx(i,j,k);
	mscsFunction::setf(idx,x,y);
}

/***************************************************************************************/
long mscsFunction3dreg::ijk2idx(long i, long j, long k) const {
//	return k*_sizeXY+j*_sizeY+i;
	return k+j*_Nz+i*_NYZ;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dreg::load(string filename, bool commentedFile) {
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
mscsFunction3dreg& mscsFunction3dreg::generateRandomGaussianField(double m, double s, long seed) {
	cpedsRNG *rns = new cpedsRNG("gaussian_circle");
	rns->seed(seed);
		
	double *a= rns->getRN(getN());
	delete rns;
	importFunction(a,NULL,getN());
	delete [] a;
	return *this;
}


/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::generateRandomGaussianFField(const mscsFunction& Pk) {
		
	return *this;
}
/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::generateRandomGaussianField(const mscsFunction& Pk) {
		
	return *this;
}

/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::fft_r2c(fftw_complex** dk) {
	double* in=extractArguments();
	long n=pointsCount();
	for (long i = 0; i < n; i++) { // normalizing here by n so that the fundamental mode is the average
		in[i]/=n;
	}

	long M=Nx()*Ny()*(long(Nz()/2)+1);
	fftw_complex* out=(fftw_complex*)fftw_malloc((M)*sizeof(fftw_complex)); 
	
	fftw_plan p;
	p=fftw_plan_dft_r2c_3d(Nx(),Ny(),Nz(), in,out, FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_free(in);
	*dk=&out[0];
	return *this;
	
}
/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::fft_c2c(fftw_complex** dk) {
//	fftw_complex* in= new fftw_complex[getN()];
	fftw_complex* in=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	long n=getN();
	printf("%li %li\n",n,pointsCount());
	for (long i = 0; i < n; i++) {
		in[i][0]=mscsFunction::getx(i)/n; // normalizing here by n so that the fundamental mode is the average
		in[i][1]=mscsFunction::Y(i)/n;
	}
//	exit(0);
	fftw_complex* out=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	
	fftw_plan p;
	p=fftw_plan_dft_3d(Nx(),Ny(),Nz(), in,out, FFTW_FORWARD,FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_free(in);
	*dk=&out[0];
	return *this;
	
}
/***************************************************************************************/
mscsFunction3dreg mscsFunction3dreg::slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky, const cpedsList<double>& kz) {
//	fftw_complex* in=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	msgs->warning("slowft_r2c: THIS IS NOT TESTED, PLEASE CHECK THE CODE",Top);
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	long nkx=kx.size();
	long nky=ky.size();
	long nkz=kz.size();
	mscsFunction3dreg g("FT");
	double re,im;
	
	for (long ki = 0; ki < nkx; ki++) {
		for (long kj = 0; kj < nky; kj++) {
			for (long kk = 0; kk < nkz; kk++) {
				re=im=0;
				
				for (long i = 0; i < nx; i++) {
					for (long j = 0; j < ny; j++) {
						for (long k = 0; k < nz; k++) {
							re+= fRe(i,j,k)*cos(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j) + kz[kk]*getZ(k) )   );
							im+=-fIm(i,j,k)*sin(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j) + kz[kk]*getZ(k) )   );
						}
					}
				}
				g.newPoint(re,im);
			}
		}
	}
	
	g._X=kx;
	g._Y=ky;
	g._Z=kz;
	g.setSize(nx,ny,nz);
	
	return g;
}

/***************************************************************************************/
mscsFunction mscsFunction3dreg::powerSpectrum() {
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
	long nz=Nz();

//	long nx=int(Nx()/2)+1;
//	long ny=int(Ny()/2)+1;

	long n=getN();

	printf("nz %li ny %li nz %li\n",nx,ny,nz);
	double kx,ky,kz, modk;
	
	
	//
	// calculate power P(vec k) 
	//
	long l=0;
	double tmp;
	
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				tmp=(dk[l][0]*dk[l][0] + dk[l][1]*dk[l][1]); // tmp*=tmp; // BLcomment (Nov 11, 2011, 9:58:39 AM):  commented out. Don't understand why squared ?
				Pveck.append(tmp);
				l++;
			}
		}
	}
	fftw_free(dk);
	Pveck.save("fft");
	
	//
	// calculate power P(k) = integrate sqrt(kx^2+ky^2+kz^2) f(vec k) dkx dky dkz
	//
	l=0;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				// derive k=sqrt(kx^2+ky^2+kz^2)
				kx=double(i)/double(_sizeX);
				ky=double(j)/double(_sizeY);
				kz=double(k)/double(_sizeZ);
				modk=sqrt(kx*kx + ky*ky + kz*kz);
				// derive the power
				Pk.newPoint(modk,Pveck[l]*modk);
				l++;
			}
		}
	}
	Pk/=double(_sizeX*_sizeY*_sizeZ);
	
	Pk.sortFunctionArgAscending();
	return Pk;	
}

/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::mkSin3Drad(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y,z;
	long i,j,k;

	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			z=from; k=0;
			while (z<=to) {
				newPoint(sin(twoPI/T*sqrt(x*x + y*y + z*z)),0);
//				newPoint(x*x + y*y + z*z,0);
				z+=dr;
				k++;
			}
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=_z0=from;
	_dx=_dy=_dz=dr;
	_dxo2=_dyo2=_dzo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	_Nz=k;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
	_sizeZ=_Nz*_dz;
	
	_NYZ=_Ny*_Nz;
	_NXY=_Nx*_Ny;
	
	return *this;	
}
/***************************************************************************************/
mscsFunction3dreg& mscsFunction3dreg::mkSin3D(double from, double to, double dr, double T, double phi) {
	clearFunction();
	double x,y,z;
	long i,j,k;

	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			z=from; k=0;
			while (z<=to) {
				newPoint(sin(twoPI/T*(x + y + z)),0);
//				newPoint(x*x + y*y + z*z,0);
				z+=dr;
				k++;
			}
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	_x0=_y0=_z0=from;
	_dx=_dy=_dz=dr;
	_dxo2=_dyo2=_dzo2=dr/2.0;
	
	_Nx=i;
	_Ny=j;
	_Nz=k;
	
	_sizeX=_Nx*_dx;
	_sizeY=_Ny*_dy;
	_sizeZ=_Nz*_dz;
	
	_NYZ=_Ny*_Nz;
	_NXY=_Nx*_Ny;
	
	return *this;	
}

/***************************************************************************************/
void mscsFunction3dreg::setSize(long Nx, long Ny, long Nz, double dx, double dy, double dz, double x0, double y0, double z0) {
	setSize(Nx,Ny,Nz);
	_dx=dx; _dxo2=dx/2.0;
	_dy=dy; _dyo2=dy/2.0;
	_dz=dz; _dzo2=dz/2.0;
	_x0=x0;
	_y0=y0;
	_z0=z0;
}
/***************************************************************************************/
void mscsFunction3dreg::setSize(long Nx, long Ny, long Nz) {
	_sizeX=Nx;
	_sizeY=Ny;
	_sizeZ=Nz;
}


/***************************************************************************************/
matrix<double> mscsFunction3dreg::getSlice(long plane, long coord, long part) const {
	
	matrix<double> m;
	long n1,n2;
	
	switch (plane) {
		case 0:
			n1=Ny();
			n2=Nz();
			break;
		case 1:
			n1=Nx();
			n2=Nz();
			break;
		case 2:
			n1=Nx();
			n2=Ny();
			break;
	}

	m.SetSize(n1,n2);
	for (long i1 = 0; i1 < n1; i1++) {
		for (long i2 = 0; i2 < n2; i2++) {
//			printf("idx: %li, i1 %li, i2 %li, n1: %li, n2: %li, _NYZ: %li\n",ijk2idx(coord,i1,i2),i1,i2,n1,n2,_NYZ);
			
			if (plane==0) {
				if (part==0) m(i1,i2)=getx(ijk2idx(coord,i1,i2)); else  m(i1,i2)=Y(ijk2idx(coord,i1,i2)); }
			if (plane==1) {
				if (part==0) m(i1,i2)=getx(ijk2idx(i1,coord,i2)); else  m(i1,i2)=Y(ijk2idx(i1,coord,i2)); }
			if (plane==2) {
				if (part==0) m(i1,i2)=getx(ijk2idx(i1,i2,coord)); else  m(i1,i2)=Y(ijk2idx(i1,i2,coord)); }

		}
	}
	return m;
}
/***************************************************************************************/
void mscsFunction3dreg::printInfo() const {
	printf("sizeX sizeY sizeZ: %lE %lE %lE\n",_sizeX,_sizeY,_sizeZ);
	printf("Nx Ny Nz: %li %li %li\n",_Nx,_Ny,_Nz);	
	printf("dx dy dz: %lE %lE %lE\n",_dx,_dy,_dz);	
	printf("x0 y0 z0: %lE %lE %lE\n",_x0,_y0,_z0);	
}









