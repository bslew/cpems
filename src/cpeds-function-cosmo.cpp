#include "cpeds-function-cosmo.h"
#include "cpeds-cosmo.h"
#include "gsl/gsl_sf_bessel.h"

const cpedsFunctionCosmo& cpedsFunctionCosmo::operator=(const mscsFunction& rhs) {
	if (this != &rhs) {
		clearFunction();
		long N=rhs.pointsCount();
		for (long i=0;i<N;i++) {
			//      newPoint(rhs.getPoint(i)); ---- why this doesn't work ???
			//      newPoint(rhs.getX(i),rhs.f(i)); -- this works fine
			newPoint(rhs.getQPoint(i)); // this works fine - but faster
		}
		return *this; 
	}
	return *this;
}
// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_age_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0, long N) {
	double z=zSt;
	clearFunction();
	if (N>0) cpeds_set_points_number_per_logz(N);
	if (zSt<0) {
		while (z<=-CPEDS_ZMIN) {
			newPoint(z,cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w0,h0));
			z+=cpeds_deltaz(z);
			msgs->say("zSt: "+msgs->toStr(zSt)+" zEn: "+msgs->toStr(zEn)+" now at z: "+msgs->toStr(z),Medium);
		}
		z=CPEDS_ZMIN;
	}
	while (z<=zEn) {
		newPoint(z,cpeds_age_of_universe(Wr0,Wm0,Wl0,z,w0,h0));
		z+=cpeds_deltaz(z);
		msgs->say("zSt: "+msgs->toStr(zSt)+" zEn: "+msgs->toStr(zEn)+" now at z: "+msgs->toStr(z),Medium);
	}
	//save("t_vs_z.txt");
	
	return *this;
}
// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_q_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long N, cpedsFunctionCosmo* qvsz) {
	
	// // q = (-2 z'(t)^2 + (1+z(t)) z''(t))/ z'(t)^2
	// // this implementation is probably incorrect. why ?
	//   mscsFunction Z,Zp,Zpp;
	
	//   make_age_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,H0,N);
	//   invert();
	//   sortFunctionArgAscending();
	//   Z=unique();
	//   Zp=Z.derivative();
	//   Zp.unique();
	//   Zpp=Zp.derivative();
	//   Zpp.unique();
	
	//   // Z.save("z_vs_t.txt");
	//   // Zp.save("zp_vs_t.txt");
	//   // Zpp.save("zpp_vs_t.txt");
	
	//   *this= (Zp*Zp*2.0 + (Z+1.0)*Zpp)/ (Zp * Zp);
	
	// if (qvsz !=NULL) {
	//   qvsz->clearFunction();
	//   qvsz->importFunction(Z.extractValues(), extractValues(), pointsCount());
	// }
	
	//-------------------------
	
	cpedsFunctionCosmo a,ap,app,q;
	a.make_a_vs_t(zSt,zEn,Wr0,Wm0,Wl0,w0,h0,N);
	ap=a.derivative();
	app=ap.derivative();
	
	// a.save("a.txt");
	// ap.save("ap.txt");
	// app.save("app.txt");
	
	*this=app*a*double(-1)/(ap*ap);
	//  print();
	//  exit(0);
	//  app=ap;
	//  app*=double(-1.0);
	//  app/=ap;
	//  app/=ap;
	//  *this=app;
	
	if (qvsz !=NULL) {
		long n = pointsCount();
		*qvsz=*this;
		for (long i=0;i<n;i++) {    qvsz->setarg(i,1/a.f(i)-1.0);  }
	}
	return *this;
}

// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_a_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0, long N) {
	make_age_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,h0,N);
	long n = pointsCount();
	for (long i=0;i<n;i++) {    setarg(i,1/(1+X(i)));  }
	invert();
	
	return *this;
}
// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_sigma0_vs_R(double Rst, double Ren, double dR, double Wr0, double Wm0, double Wl0, double w0, double h0,long N) {
	
	return *this;
	
}
// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_tau_vs_t(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long N, cpedsFunctionCosmo* visibz) {
	//	  double zz,d,dz;
	//
	//	  zz=zSt;  d=0.0;
	//	  if (zSt==0.0) zz=CPEDS_ZMIN;
	//
	//	  do {
	//	    dz=cpeds_deltaz(zz);
	//	    zz=zz-dz;
	//	    d+=dz / cpeds_Efactor(Wr0,Wm0,Wl0,zz+0.5*dz,w0) * sigma_T/c/(h0*100000/MPC) * cpeds_fX_e(T,eta) * eta *Wr0*(1+z)*(1+z)*(1+z);
	//	  } while (zz-0.1*dz >= z);
	//
	//	  return d*c/1.0e8;  
	
	return *this;
}

// ****************************************************************************************************
cpedsFunctionCosmo& cpedsFunctionCosmo::make_dA_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long Nout, cpedsFunctionCosmo* chi) {
	double z=zSt;
	clearFunction();
	
	double *dp=NULL;
	double *Z=NULL;
	long n;
	cpeds_dp(Wr0,Wm0,Wl0,0,w0,&Z,&dp,&n); // particle horizon [Gpc]
	//	  for (long i = 0; i < n; i++) {
	//		printf("%li z %lE dp %lE\n",i,Z[i],dp[i]);
	//	  }
	importFunction(Z,dp,n);
	// remove the data outside of the requested range.
	deleteOutsideOf(zSt,zEn);
	n=pointsCount();	  
	delete [] Z;
	delete [] dp;
	
	divide(-h0);
	double dp0=cpeds_dp(Wr0,Wm0,Wl0,0,w0,h0);
	add(dp0);
	// convert to Mpc
	multiply(1000);
	// now we have calculated the comoving distance [Mpc] vs z
	
	if (chi!=NULL) *chi=*this;
	
	// derive the angular size distance
	double Wk0=cpeds_Wk0(Wr0,Wm0,Wl0);
	if (Wk0 == 0.0) { // flat case
		return *this; 
	}
	else {
		double Rc=cpeds_curvature_radius(Wr0,Wm0,Wl0,h0)* 1000; // [Mpc]
		
		if (Rc < 0.0) { // negative curvature
			// dA = -Rc * sinh( -X/Rc ) ; 
			for (long i = 0; i < n; i++) {	f(i)=-Rc*sinh( -f(i)/Rc )/(1+X(i));	}
		}
		else { // positive curvature
			//dA = Rc * sin( X/Rc ); 
			for (long i = 0; i < n; i++) {	f(i)=Rc*sinh( f(i)/Rc )/(1+X(i));	}
		}
	}
	
	
	
	//	  double bin=pointsCount()/Nout;
	//	  binFunction;
	
	
	//	  // this is very slow implementation
	//	  while (z<=zEn) {
	//		  newPoint(z,cpeds_angular_diameter_distance(Wr0,Wm0,Wl0,z,w0,h0));
	//		  z+=cpeds_deltaz(z);
	//		  msgs->say("zSt: "+msgs->toStr(zSt)+" zEn: "+msgs->toStr(zEn)+" now at z: "+msgs->toStr(z),Medium);
	//	  }
	
	
	
	return *this;
	
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_lbtd_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, double h0,long Nout) {
	clearFunction();
	make_age_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,h0,Nout);
	double age=cpeds_age_of_universe(Wr0,Wm0,Wl0,0,w0,h0);
	
	for (long i = 0; i < pointsCount(); i++) {
		f(i)=CPEDS_c*(age-f(i))*GYR/CPEDS_MPC;
	}
	
	return *this;
	
}

/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wx_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0, cpedsFunctionCosmo::fptrWx_t fWx) {
	double z=zSt;
	clearFunction();
	if (zSt<0) {
		while (z<=-CPEDS_ZMIN) {
			newPoint(z,(*fWx)(Wr0,Wm0,Wl0,z,w0));
			z+=cpeds_deltaz(z);
//			msgs->say("zSt: "+msgs->toStr(zSt)+" zEn: "+msgs->toStr(zEn)+" now at z: "+msgs->toStr(z),Medium);
		}
		z=CPEDS_ZMIN;
	}
	while (z<=zEn) {
		newPoint(z,(*fWx)(Wr0,Wm0,Wl0,z,w0));
		z+=cpeds_deltaz(z);
//		msgs->say("zSt: "+msgs->toStr(zSt)+" zEn: "+msgs->toStr(zEn)+" now at z: "+msgs->toStr(z),Medium);
	}
	
	return *this;
	
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wr_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0) {
	fptrWx_t Wx=&cpeds_Wr;
	return make_Wx_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,Wx);
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wm_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0) {
	fptrWx_t Wx=&cpeds_Wm;	
	return make_Wx_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,Wx);
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wl_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0) {
	fptrWx_t Wx=&cpeds_Wl;	
	return make_Wx_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,Wx);
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wtot_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0) {
	fptrWx_t Wx=&cpeds_Wtot;	
	return make_Wx_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,Wx);	
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_Wk_vs_z(double zSt, double zEn, double Wr0, double Wm0, double Wl0, double w0) {
	fptrWx_t Wx=&cpeds_Wk;	
	return make_Wx_vs_z(zSt,zEn,Wr0,Wm0,Wl0,w0,Wx);	
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::make_BBKS_transferFunction(double kmin, double kmax, double dlogk, double Wb, double Wtot, double h) {
	double Gamma=Wtot*h*exp(-Wb-sqrt(2.0*h)*Wb/Wtot);
	double logkmin,logkmax,dq,logk,q,k,T,tmp;
	printf("Gamma: %lE\n",Gamma);
	printf("Wtot: %lE, Wb: %lE, h: %lE\n",Wtot,Wb, h);
	logkmin=log10(kmin);
	logkmax=log10(kmax);
	printf("logkmin: %lE, logkmax: %lE, dlogk: %lE\n",logkmin,logkmax,dlogk);
	logk=logkmin;
	while (logk<logkmax) {
		k=exp10(logk);
		q=k/Gamma;
		tmp=sqrt((1.0+3.89*q+pow(16.1*q,2)+pow(5.46*q,3)+pow(6.71*q,4)));
		T=log(1.0+2.34*q)/(2.34*q)/sqrt(tmp);
		newPoint(k,T);
		logk+=dlogk;
	}
	return *this;
}

/***************************************************************************************/
cpedsFunctionCosmo cpedsFunctionCosmo::make_tinker_mass_function(mscsFunction& sigmaM, double Wm0, double deltaRho) {
	mscsFunction A, a, b, c;
	double Ai,ai,bi,ci;
	
	A.newPoint(200,0.186);
	A.newPoint(300,0.200);
	A.newPoint(400,0.212);
	A.newPoint(600,0.218);
	A.newPoint(800,0.248);
	A.newPoint(1200,0.255);
	A.newPoint(1600,0.260);
	A.newPoint(2400,0.260);
	A.newPoint(3200,0.260);

	a.newPoint(200, 1.47);
	a.newPoint(300, 1.52);
	a.newPoint(400, 1.56);
	a.newPoint(600, 1.61);
	a.newPoint(800, 1.87);
	a.newPoint(1200,2.13);
	a.newPoint(1600,2.30);
	a.newPoint(2400,2.53);
	a.newPoint(3200,2.66);

	b.newPoint(200, 2.57);
	b.newPoint(300, 2.25);
	b.newPoint(400, 2.05);
	b.newPoint(600, 1.87);
	b.newPoint(800, 1.59);
	b.newPoint(1200,1.51);
	b.newPoint(1600,1.46);
	b.newPoint(2400,1.44);
	b.newPoint(3200,1.41);

	c.newPoint(200, 1.19);
	c.newPoint(300, 1.27);
	c.newPoint(400, 1.34);
	c.newPoint(600, 1.45);
	c.newPoint(800, 1.58);
	c.newPoint(1200,1.80);
	c.newPoint(1600,1.97);
	c.newPoint(2400,2.24);
	c.newPoint(3200,2.44);
	
	Ai=A.finter(deltaRho,"linear");
	ai=a.finter(deltaRho,"linear");
	bi=b.finter(deltaRho,"linear");
	ci=c.finter(deltaRho,"linear");
	
	printf("A: %lf, a: %lf, b: %lf, c: %lf\n",Ai,ai,bi,ci);

	cpedsFunctionCosmo tinker;
	mscsFunction sM=sigmaM;
	sM.scaleX(1.0/(CPEDS_SOLAR_MASS)); // convert mass unit from kg to Msol
	mscsFunction tmp=sM;
	double rhom=cpeds_rhoC0(1.0)*Wm0 / (CPEDS_SOLAR_MASS/pow(CPEDS_MPC,3) ); // Msol/Mpc^3 h^2
//	double rhom=cpeds_rhoC0(1.0)*Wm0; // kg/m^3 h^2
	printf("rhom: %lE\n",rhom);
	// calculate: dln sigmaM^-1/dM = -1/sigmaM * d sigmaM/dM
//	tmp.save("lnsigmaMinv");
//	tmp.sortFunctionArgAscending();
//	sM.sortFunctionArgAscending();
	
	tmp.derivative(true);
	
	
	double fsigma;
	double s,M;
	for (unsigned long i = 0; i < sM.pointsCount(); i++) {
		s=sM.f(i);
		M=sM.getX(i);
		fsigma=Ai*( pow(s/bi,-ai) + 1)*exp(-ci/(s*s));
		printf("sigmaM: %lE, M: %lE, log(1/sigmaM): %lE, fsigma: %lE, log10(fsigma): %lE, rhom/M: %lE\n",s,M,log10(1.0/s),fsigma, log10(fsigma),rhom/M);
		tmp.f(i)=- fsigma * rhom/M * tmp.f(i) / s;
		tinker.newPoint(M,M*M/rhom*tmp.f(i));
	}

	(*this)=tmp;
	
	return tinker;

}

cpedsFunctionCosmo cpedsFunctionCosmo::make_PS_mass_function(mscsFunction& sigmaM, double Wm0, double z) {
	cpedsFunctionCosmo PS;
	mscsFunction sM=sigmaM;
	sM.scaleX(1.0/(CPEDS_SOLAR_MASS)); // convert mass unit from kg to Msol
	mscsFunction tmp=sM;
	double rhom=cpeds_rhoC0(1.0)*Wm0 / (CPEDS_SOLAR_MASS/pow(CPEDS_MPC,3) ); // Msol/Mpc^3 h^2
	printf("rhom: %lE\n",rhom);
	// calculate: dln sigmaM^-1/dM = -1/sigmaM * d sigmaM/dM
	
	tmp.derivative(true);
	
	
	double fact;
	double s,M;
	double delta_c = 1.686; // this value is model dependent. The 1.686 value is for the EdS model and there should be some correction for the LCDM model
	for (unsigned long i = 0; i < sM.pointsCount(); i++) {
		s=sM.f(i);
		M=sM.getX(i);
		fact=sqrt(2.0/PI)*delta_c/s * exp(-delta_c*delta_c/(2.0*s*s));
		printf("sigmaM: %lE, M: %lE, log(1/sigmaM): %lE, fsigma: %lE, log10(fsigma): %lE, rhom/M: %lE\n",s,M,log10(1.0/s),fact, log10(fact),rhom/M);
		tmp.f(i)=- fact * rhom/M * tmp.f(i) / s;
		PS.newPoint(M,M*M/rhom*tmp.f(i));
	}

	(*this)=tmp;
	
	return PS;
	
}

/***************************************************************************************/
cpedsFunctionCosmo cpedsFunctionCosmo::make_Pk(mscsFunction& tf, double A, double k0, double ns, double h) {
//cpedsFunctionCosmo cpedsFunctionCosmo::make_Pk(mscsFunction& tf, double A, double k0, double ns, double kmin, double kmax, double dlogk) {
//	double logkmin,logkmax,logk,k,tmp;
//	printf("Wtot: %lE, Wb: %lE, h: %lE\n",Wtot,Wb, h);
//	logkmin=log10(kmin);
//	logkmax=log10(kmax);
//	printf("logkmin: %lE, logkmax: %lE, dlogk: %lE\n",logkmin,logkmax,dlogk);
//	logk=logkmin;
//	while (logk<logkmax) {
//		k=exp10(logk);
//		tmp=4.0/25*pow(CPEDS_c*k/100000,4)*
//		T=log(1.0+2.34*q)/(2.34*q)/sqrt(tmp);
//		newPoint(k,T);
//		logk+=dlogk;
//	}
	
	double ktf,pk_h,k,pk;
	for (unsigned long i = 0; i < tf.pointsCount(); i++) {		
		// since provided (via tf ) k_tf= k/h  and is in h/Mpc we calculate k in [Mpc^-1]
		ktf=tf.getx(i);
		k=ktf*h;
		pk=4.0/25*pow(CPEDS_c*k/(100000*h),4)*tf.f(i)*tf.f(i) * A*pow(k/k0,ns-1.0); 
		pk_h=pk_h*h*h*h; // convert to h^3/Mpc^3
		newPoint(k,pk_h);
	}

	return *this;
	
}



/***************************************************************************************/
cpedsFunctionCosmo cpedsFunctionCosmo::calculateMassFunction(mscsFunction haloMass, double vol, long Npts, string factors, double from, double to) {
	double* masses=haloMass.logY(10).extractValues();
	long Nbin=Npts;
	double* mbin = new double[Nbin];
	double *counts;
	if (from==0 and to==0)
		counts=cpeds_bin_data(haloMass.pointsCount(),masses,Nbin,mbin,0);
	else {
		from=log10(from);
		to=log10(to);
		counts=cpeds_bin_data(haloMass.pointsCount(),masses,Nbin,mbin,0,from,to);
	}
	clearFunction();
	importFunction(mbin,counts,Nbin);
	double dM2=(mbin[1]-mbin[0])/2;
//	printf("dM2: %lE\n",dM2);
//		massFunction.logX(10);
	// normalize by volume
//		massFunction.exponent();
	divide(vol); // change units from counts to counts * Mpc^-3 h^3 / (10^10 Msol /h)
	for (unsigned long i = 0; i < pointsCount(); i++) { // divide by mass range in the bin: we  calculate dN/dM so the units must be counts/(Mpc/h)^3/(Msol/h)
		f(i)/=exp10(mbin[i]+dM2)-exp10(mbin[i]-dM2); // change units to counts * Mpc^-3 h^3 Msol^-1 h^1			
	}
	divide(1e10); // change mass units from 10^10Msol/h to Msol/h
	logY(10); // calculate log10 (dN/dM)
	shiftX(10); // add log 10^10=10 to change units from gadget units (10^10 Msol/h) to Msol/h. We do this in log10 space hence the shift.
	
	delete [] mbin;
	delete [] counts;
	delete [] masses;

	return *this;
}
/***************************************************************************************/
cpedsFunctionCosmo& cpedsFunctionCosmo::calculateBeamPowerPattern(double d,double ds, double lambda, double taper, double thetaMax, double dTheta) {
	
	double th=-thetaMax;
	double z,q,u,J1z,J2z, J1kz,J2kz;
	double kappa=ds/d;
	double k2=kappa*kappa;
	double U0betaz,U0betakkzk,U;
	double beta=1.0-pow(10,taper/20);
	printf("beta: %lf\n",beta);
	u=(k2 * (2.0-k2*beta))/(2.0-beta);
	
//	cpedsFunctionCosmo Uth;
	
	while (th<thetaMax) {
		q=sin(th*PI180)/lambda;
		z=PI*d*q;
		J1z=gsl_sf_bessel_J1(z);
		J2z=gsl_sf_bessel_Jn(2,z);
		J1kz=gsl_sf_bessel_J1(kappa*z);
		J2kz =gsl_sf_bessel_Jn(2,kappa*z);
		U0betaz=4.0/(2.0-beta)/z * ( (1.0-beta)*J1z + 2.0*beta*J2z/z );
		U0betakkzk=4.0/(2.0-beta*k2)/(z*kappa) * ( (1.0-beta*k2)*J1kz + 2.0*beta*k2*J2kz/(z*kappa) );
		
		U=(U0betaz - u*U0betakkzk)/(1.0-u);
		
//		Uth.newPoint(th,U*U);
		newPoint(th,U*U);
		
		th+=dTheta;
	}
//	return Uth;
	return *this;
}
