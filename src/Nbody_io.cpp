/*!
* \file Nbody_io.cpp
*
*  Created on: Jan 3, 2011
*      Author: blew
*/

//#include <iostream>
//#include <fstream>
#include <stdio.h>
#include "cpeds-cosmo.h"
#include "Nbody_io.h"


Nbody_io::Nbody_io() {
	msgs = new cpedsMsgs("Nbody IO");
	P=NULL;
	Id=NULL;
	gadgetMassUnit=1e10*CPEDS_SOLAR_MASS; // in SI
	
	tipsyPdark=NULL;
	tipsyPgas=NULL;
//	gadgetMassUnit=GADGET2_UnitMass_in_g;

}
/***************************************************************************************/
Nbody_io::~Nbody_io() {
	delete msgs;
	free_gadget_memory();
	free_tipsy_memory();
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::loadGadgetSnapshot(string fname, long files, int fileType) {
	msgs->say("Reading the GADGET-2 snapshot file",Top);
	FILE *fd;
	char   buf[200];
	int    i,j,k,dummy,ntot_withmasses;
	int    t,n,off,pc,pc_new,pc_sph;
	
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
	
	for(i=0, pc=1; i<files; i++, pc=pc_new)	{
		if(files>1)	sprintf(buf,"%s.%d",fname.c_str(),i);
		else sprintf(buf,"%s",fname.c_str());
		
		if(!(fd=fopen(buf,"r"))) { msgs->criticalError("cannot open the file: "+string(buf),High); }
		
		msgs->say("Reading data from file: "+string(buf),High); fflush(stdout);
		
		//
		// read in the header information from the file
		// 
		SKIP;
		fread(&gadget_header, sizeof(gadget_header), 1, fd);
		SKIP;
		
		if(files==1) {
			for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) NumPart+= gadget_header.npart[k];
			Ngas= gadget_header.npart[0];
		}
		else {
			for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) NumPart+= gadget_header.npartTotal[k];
			Ngas= gadget_header.npartTotal[0];
		}
		
		msgs->say("total number of particles in the file: "+msgs->toStr(NumPart),High);
		
		for(k=0, ntot_withmasses=0; k<5; k++) {
			if(gadget_header.mass[k]==0)
				ntot_withmasses+= gadget_header.npart[k];
		}
		
		if(i==0) allocate_gadget_memory();
		
		msgs->say("reading positions",Medium);
		SKIP;
		for (k=0,pc_new=pc;k<6;k++) {
			for (n=0;n<gadget_header.npart[k];n++) {
				fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
				pc_new++;
				//				printf("%li\n",n);
			}
		}
		SKIP;
		
		msgs->say("reading velocities",Medium);
		SKIP;
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<gadget_header.npart[k];n++) {
				fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
				pc_new++;
			}
		}
		SKIP;
		
		
		msgs->say("reading IDs",Medium);
		SKIP;
		for(k=0,pc_new=pc;k<6;k++) {
			for(n=0;n<gadget_header.npart[k];n++)	{
				fread(&Id[pc_new], sizeof(int), 1, fd);
				pc_new++;
			}
		}
		SKIP;
		
		
		if(ntot_withmasses>0) SKIP;
		msgs->say("reading or setting the masses",Medium);
		for(k=0, pc_new=pc; k<6; k++) {
			for(n=0;n<gadget_header.npart[k];n++)	{
				P[pc_new].Type=k;
				
				if(gadget_header.mass[k]==0)
					fread(&P[pc_new].Mass, sizeof(float), 1, fd);
				else
					P[pc_new].Mass= gadget_header.mass[k];
				pc_new++;
			}
		}		
		if(ntot_withmasses>0) SKIP;
		
		if(gadget_header.npart[0]>0) {
			SKIP;
			msgs->say("reading internal energy",Medium);
			for(n=0, pc_sph=pc; n<gadget_header.npart[0];n++)	{
				fread(&P[pc_sph].U, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			SKIP;
			msgs->say("reading density",Medium);
			for(n=0, pc_sph=pc; n<gadget_header.npart[0];n++) {
				fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			if(gadget_header.flag_cooling) {
				SKIP;
				msgs->say("reading Ne",Medium);
				for(n=0, pc_sph=pc; n<gadget_header.npart[0];n++)
				{
					fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
					pc_sph++;
				}
				SKIP;
			}
			else {
				msgs->say("setting Ne to 1.0",Medium);
				for(n=0, pc_sph=pc; n<gadget_header.npart[0];n++)	{
					P[pc_sph].Ne= 1.0;
					pc_sph++;
				}
			}
		}
		
		fclose(fd);
	}
	
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::loadGraficICsDensity(string dirname) {
	//	ofstream inf;
	int32_t w;
	string fname;
	
	FILE* f;
	
	
	fname=dirname+"/"+"ic_deltab";
	msgs->say("reading file:"+fname,High);
	
	grafic_header=read_grafic_header(fname);
	print_grafic_header_info(grafic_header);
	
//	f=fopen(fname.c_str(),"rb");
//	fread(&w,sizeof(w),1,f);
//	fread(&header,sizeof(grafic_level0_header_t),1,f);
//	fread(&w,sizeof(w),1,f);	
	
//	printf("%li %li %li\n", (long)header.n1, (long)header.n2, (long)header.n3);
//	printf("%f %f %f\n", header.o1, header.o2, header.o3);
//	printf("dx %f\n", header.dx);
//	printf("Om %f\n", header.omegam);
//	printf("Ov %f\n", header.omegav);
//	printf("ast %f\n", header.astart);
//	printf("H0 %f\n", header.H0);
//	printf("\n\n");
	
	
	ICdata.setSize(grafic_header.n3,grafic_header.n2,grafic_header.n1, 
			grafic_header.dx, grafic_header.dx, grafic_header.dx, 
			grafic_header.o3,grafic_header.o2,grafic_header.o1);
	ICdata.allocFunctionSpace();
	
	f=fopen(fname.c_str(),"rb");
	fseek(f,get_grafic_header_size(),SEEK_SET);
	
	double* slice=NULL; 
	long cellsInSlice;
	
	for (long k = 0; k < grafic_header.n3; k++) {
		slice=read_grafic_slice(&f,&cellsInSlice);
		
		ICdata.importSlice(slice,cellsInSlice,k);

//		fread(&sliceSize,sizeof(uint32_t),1,f);
//		printf("slice size [B]: %li\n",sliceSize);
//		
//		bool iDouble=false;
//		if (sliceSize/(h.n1*h.n2)==4) { iDouble=false; cellsInSlice=sliceSize/4; printf("data type is 4-byte floating point\n"); } 
//		else {
//			if (sliceSize/(header.n1*header.n2)==8) { cellsInSlice=sliceSize/8; iDouble=true; printf("data type is 8-byte floating point\n"); }
//			else msgs->criticalError("wrong data type size estimation.", High);
//		}
//		
//		double *slice=new double[cellsInSlice];
//		if (iDouble) {
//			fread(slice,sizeof(double),cellsInSlice,f);			
//		}
//		else {
//			float *slicef=new float[cellsInSlice];
//			fread(slicef,sizeof(float),cellsInSlice,f);
//			for (long i = 0; i < cellsInSlice; i++) {	slice[i]=slicef[i];	}
//			delete [] slicef;
//		}
//		
//		for (long i = 0; i < cellsInSlice; i++) {
//			//		printf("%li %lE ",i,slice[i]);
//			//			printf("%lE ",slice[i]);
//		}
//		
//		int count=fread(&sliceSize,sizeof(uint32_t),1,f);
//		//		printf("\nslice size [B]: %li\n",sliceSize);
		
	}
	
	
	//	inf.open(f.c_str(), std::ios::in | std::ios::binary );
	//	assert( inf.is_open() );	
	//	inf.read(reinterpret_cast<int32_t *>(&w),sizeof(w));
	//	w = sizeof(grafic_level0_header_t);
	//	inf.read(reinterpret_cast<grafic_level0_header_t *>(&header),w);
	
	
	fclose(f);
	
	return cpedsSuccess;	
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::convertGraficToGadgetIC(string dirname, string outfile, int gadgetICformat, int fileNum) {
	string fname;
	gadget_snapshot_file_header_t gadget_header;
	double velocity_factor;	
	grafic_program_parameters_t grafic_params;
	
	if (gadgetICformat==2) { msgs->criticalError("gadget2 format is not supported here yet",High); }
	//
	// read in the grafic header for baryons
	//
	fname=dirname+"/"+"ic_deltab";
	msgs->say("reading file:"+fname,High);
	grafic_header=read_grafic_header(fname);
	
	//
	// read in other parameters
	//
	velocity_factor=getVelocityFactor(dirname);
	grafic_program_parameters=read_grafic_parameters(dirname);
	grafic_cosmo_parameters=read_grafic_cosmo_parameters(dirname);
	calculateParticleMasses();
	
	// output some info
	print_grafic_header_info(grafic_header);
	msgs->say("mass per baryon particle in Gadget units: "+msgs->toStr(baryonParticleMass),Medium);
	msgs->say("mass per CDM    particle in Gadget units: "+msgs->toStr(CDMparticleMass),Medium);	
	msgs->say("velocity factor [km/s/Mpc]: "+msgs->toStr(velocity_factor),Medium);
//	msgs->say("mass per baryon particle in Gadget units: "+msgs->toStr(baryonParticleMass),Medium);

	//
	// number of particles per file
	//
	long PperF=grafic_header.n1*grafic_header.n2*grafic_header.n3/fileNum;
	if (fileNum*PperF!=grafic_header.n1*grafic_header.n2*grafic_header.n3) {
		msgs->criticalError("this method does not yet support unequal number of particles per file. The total number of particles must be the number of particles in one file times number of files. In every file there must be the same number of particles",High);
	}
	long slicesPerFile=grafic_header.n3/fileNum; // number of grafic file slices per one Gadget-format snapshot file 
	if (fileNum*slicesPerFile!=grafic_header.n3) {
		msgs->criticalError("The size of the block saved in one file times number of files must be equal the total number of slices in the grafic file.",High);
	}
	
	//
	// construct the gadget header block for baryons
	//
//	gadget_header.BoxSize=grafic_header.dx*grafic_header.n1*1000; // this assumes that the grafic box has the same sizes in all dimentions; the unit is converted from Mpc as in grafic to the gadget defauld internal units: kpc
//	gadget_header.HubbleParam=grafic_header.H0/100.0;
	gadget_header.HubbleParam=grafic_header.H0/100.0;
	gadget_header.BoxSize=grafic_header.dx*grafic_header.n1*1000*gadget_header.HubbleParam; // this assumes that the grafic box has the same sizes in all dimentions; the unit is converted from Mpc as in grafic to the gadget defauld internal units: kpc/h
	gadget_header.Omega0=grafic_header.omegam;
	gadget_header.OmegaLambda=grafic_header.omegav;
	gadget_header.flag_cooling=0; // unused
	gadget_header.flag_feedback=0; // unused
	gadget_header.flag_sfr=0; // unused
	gadget_header.flag_stellarage=0; // unused
	gadget_header.flag_metals=0; // unused
	gadget_header.flag_entropy_instead_u=0; // unused
	gadget_header.mass[0]=baryonParticleMass;
	gadget_header.mass[1]=CDMparticleMass;
	gadget_header.mass[2]=0;
	gadget_header.mass[3]=0;
	gadget_header.mass[4]=0;
	gadget_header.mass[5]=0;
	gadget_header.npart[0]=PperF;
	gadget_header.npart[1]=0;
	gadget_header.npart[2]=0;
	gadget_header.npart[3]=0;
	gadget_header.npart[4]=0;
	gadget_header.npart[5]=0;
	gadget_header.num_files=fileNum;
	gadget_header.redshift=grafic_program_parameters.zstart;
	gadget_header.time=grafic_program_parameters.astart; // this is not used as it will be defined as a parameter TimeBegin in the parameter input file for Gadget-2
	char fill[60]; for (long i = 0; i < 60; i++) { fill[i]='\0'; }
	memcpy(gadget_header.fill,fill,sizeof(char)*60); // this is to make identical IC files really identical
	// this stuff is needed for gadget to correctly read the particles number from the header
	for (long i = 0; i < 6; i++) { gadget_header.npartTotalHighWord[i]=(unsigned int) ((long long)(gadget_header.npart[i]) >> 32); }
	
	//
	// construct the gadget header block for CDM
	//

	gadget_header.npart[1]=PperF;
	
	for (long i = 0; i < 6; i++) { gadget_header.npartTotal[i]=grafic_header.n1*grafic_header.n2*grafic_header.n3; }
	string outputFileName;
	long currentOffset=0;

	for (long fidx = 0; fidx < fileNum; fidx++) { // loop over parial files with ICs
		msgs->say("writting Gadget IC file number: "+msgs->toStr(fidx),High);
		
		//
		// save gadget header block
		//
		msgs->say("writting Gadget header block",Medium);
		if (fileNum==1) outputFileName=outfile;
		else outputFileName=outfile+"."+msgs->toStr(fidx);
		FILE* gf=fopen(outputFileName.c_str(),"wb");
		int tmp=sizeof(gadget_snapshot_file_header_t);
		fwrite(&tmp,sizeof(int),1,gf);
		fwrite(&gadget_header,sizeof(gadget_snapshot_file_header_t),1,gf);
		fwrite(&tmp,sizeof(int),1,gf);
		
		
		//
		// write remaining blocks
		//
		double dx=grafic_header.dx; // In Mpc
		double dxo2=dx/2;
		double depthPerFile=slicesPerFile*dx; // defines the depth of the pack of slices that will be saved in one gadget-format snapshot file

		double L=grafic_header.n1*dx; // Box size in Mpc (assumes that it must be a cube)
		long n; // number of particles of a given type
		int blockSize;
		int usedParticleTypes=0;
		double *vx,*vy,*vz;
		FILE *fx,*fy,*fz;
		long nsx,nsy,nsz; // number of cells in slices in velocity files x,y,z
		long ii;
		int cellsInSlice;
		float* data=NULL;
		unsigned int *dataInt;
		
		for (long i = 0; i < 6; i++) {	if (gadget_header.npart[i]>0) usedParticleTypes++;	}
		gadget_IC_blocks block;
		for (block = Coordinates; block <= SmoothingLength; block=gadget_IC_blocks(block+1)) {
						
			switch (block) {
				case Coordinates:
				case Velocities: 
					if (block==Coordinates) 	msgs->say("Writing Gadget Coordinates block",Medium);
					if (block==Velocities) 		msgs->say("Writing Gadget Velocities block",Medium);
					//
					// save gadget positions and velocities
					//
					
					blockSize=usedParticleTypes*3*grafic_header.n1*grafic_header.n2*slicesPerFile*sizeof(float);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Low);
					fwrite(&blockSize,sizeof(int),1,gf); // initiate the block
					fname=dirname+"/"+"ic_velbx";
					if (grafic_data_is_double(fname)) {
						currentOffset=fidx*(PperF*sizeof(double)+2*sizeof(int32_t)*slicesPerFile);
					}
					else {
						currentOffset=fidx*(PperF*sizeof(float)+2*sizeof(int32_t)*slicesPerFile);
					}
//					printf("current slice: %li\n", currentOffset);

					for (long particles = 0; particles < 2; particles++) { // loop over particle types
						msgs->say("Processing particles of type: "+msgs->toStr(particles),Low);
						// get the right input file, open and seek to the 
						switch (particles) {
							case 0: 
								fname=dirname+"/"+"ic_velbx"; fx=fopen(fname.c_str(),"rb");	fseek(fx,get_grafic_header_size(),SEEK_SET); fseek(fx,currentOffset,SEEK_CUR);
								fname=dirname+"/"+"ic_velby"; fy=fopen(fname.c_str(),"rb");	fseek(fy,get_grafic_header_size(),SEEK_SET); fseek(fy,currentOffset,SEEK_CUR);
								fname=dirname+"/"+"ic_velbz"; fz=fopen(fname.c_str(),"rb");	fseek(fz,get_grafic_header_size(),SEEK_SET); fseek(fz,currentOffset,SEEK_CUR);
								break;
							case 1:
								fname=dirname+"/"+"ic_velcx"; fx=fopen(fname.c_str(),"rb");	fseek(fx,get_grafic_header_size(),SEEK_SET); fseek(fx,currentOffset,SEEK_CUR);
								fname=dirname+"/"+"ic_velcy"; fy=fopen(fname.c_str(),"rb");	fseek(fy,get_grafic_header_size(),SEEK_SET); fseek(fy,currentOffset,SEEK_CUR);
								fname=dirname+"/"+"ic_velcz"; fz=fopen(fname.c_str(),"rb");	fseek(fz,get_grafic_header_size(),SEEK_SET); fseek(fz,currentOffset,SEEK_CUR);
								break;				
						}
						
//						for (long k = 0; k < grafic_header.n3; k++) { // loop over slices
						for (long k = 0; k < slicesPerFile; k++) { // loop over slices
							
							vx=read_grafic_slice(&fx,&nsx);
							vy=read_grafic_slice(&fy,&nsy);
							vz=read_grafic_slice(&fz,&nsz);
							
//							cpeds_save_matrix_lin(vx,nsx,1,fname+"-x.txt",false,true);
//							cpeds_save_matrix_lin(vy,nsy,1,fname+"-y.txt",false,true);
//							cpeds_save_matrix_lin(vz,nsz,1,fname+"-z.txt",false,true);
							
							
							if (block==Coordinates) {
								// derive offsets from velocities by dividing by velocity_factor
								cpeds_divide_value(velocity_factor,vx,nsx);
								cpeds_divide_value(velocity_factor,vy,nsy);	
								cpeds_divide_value(velocity_factor,vz,nsz);	
								
								// move particles to their cell centers
								ii=0;
								for (long j = 0; j < grafic_header.n2; j++) {
									for (long i = 0; i < grafic_header.n1; i++) {
										vx[ii]+=i*dx+dxo2;
										vy[ii]+=j*dx+dxo2;
										vz[ii]+=k*dx+dxo2 + depthPerFile*fidx;
										ii++;
									}
								}
								
								// wrap periodically if needed
								ii=0;
								for (long j = 0; j < grafic_header.n2; j++) {
									for (long i = 0; i < grafic_header.n1; i++) {
										if (vx[ii]<0) vx[ii]=L+vx[ii]; if (vx[ii]>=L) vx[ii]=vx[ii]-L;
										if (vy[ii]<0) vy[ii]=L+vy[ii]; if (vy[ii]>=L) vy[ii]=vy[ii]-L;
										if (vz[ii]<0) vz[ii]=L+vz[ii]; if (vz[ii]>=L) vz[ii]=vz[ii]-L;
										ii++;
									}
								}
								
								// convert Grafic units (comoving Mpc) to Gadget units (kpc/h)
//								cpeds_mul_value(1000.0,vx,nsx);
//								cpeds_mul_value(1000.0,vy,nsy);	
//								cpeds_mul_value(1000.0,vz,nsz);	
								cpeds_mul_value(1000.0*gadget_header.HubbleParam,vx,nsx);
								cpeds_mul_value(1000.0*gadget_header.HubbleParam,vy,nsy);	
								cpeds_mul_value(1000.0*gadget_header.HubbleParam,vz,nsz);	
								// TODO should we also multiply here by h ? 
								// yes, Mar 30, 2012, 12:07:19 PM
							}
							
							if (block==Velocities) {
								// convert units here if necessary
								// we convert here the grafic velocities to the gadget-2 velocities required in gadget IC file
								// i.e. 1/sqrt(a_start)
								// this is the case for the default internal velocity units of gadget (i.e. km/s)
								double vconv=getGadgetVelocityFactor();
								cpeds_mul_value(vconv,vx,nsx);
								cpeds_mul_value(vconv,vy,nsy);	
								cpeds_mul_value(vconv,vz,nsz);	
							}
							
							
							// form the gadget-consistent linear array of data
							cellsInSlice=3*grafic_header.n1*grafic_header.n2; // cells in slice
							data= new float[cellsInSlice];
							ii=0;
							for (long j = 0; j < grafic_header.n2; j++) {
								for (long i = 0; i < grafic_header.n1; i++) {
									data[3*ii]  =vx[ii];	
									data[3*ii+1]=vy[ii];	
									data[3*ii+2]=vz[ii];	
									ii++;
								}
							}
							
							// save to output file
							fwrite(data,sizeof(float),cellsInSlice,gf);
							
							delete [] data;
							
							delete [] vx;			delete [] vy;			delete [] vz;
							
						} // end of loop over slices
						fclose(fx);
						fclose(fy);
						fclose(fz);
						
					} // end of loop over particle types
					fwrite(&blockSize,sizeof(int),1,gf); // finish the block
					
					break;
					
				case ParticleIDs:
					msgs->say("writting Gadget ParticleIDs block",Medium);
					blockSize=usedParticleTypes*grafic_header.n1*grafic_header.n2*slicesPerFile*sizeof(unsigned int);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Low);
					fwrite(&blockSize,sizeof(int),1,gf); // initiate the bloc
					for (long i = 0; i < 6; i++) {					
						n=gadget_header.npart[i];
						if (n>0) {
							dataInt= new unsigned int[n];
							for (long j = 0; j < n; j++) {
								dataInt[j]=n*usedParticleTypes*fidx+i*n+j;
							}
							fwrite(dataInt,sizeof(unsigned int),n,gf);
							delete [] dataInt;
						}
					}
					fwrite(&blockSize,sizeof(int),1,gf); // finish the block
					break;
					
				case Masses: // do not write the masses block. The masses will be set from the header information. Same for all particles of a given type
					msgs->say("skipping Gadget Masses block",Medium);
					break;
					
				case InternalEnergy:
					//			case Rho:
					//			case SmoothingLength:
					if (gadget_header.npart[0]>0) {
						if (block==InternalEnergy) 	msgs->say("writting Gadget InternalEnergy block",Medium);
						if (block==Rho) 	msgs->say("writting Gadget Density block",Medium);
						if (block==SmoothingLength) 	msgs->say("writting Gadget Smoothing Length block",Medium);					
						
						blockSize=grafic_header.n1*grafic_header.n2*slicesPerFile*sizeof(float);
						msgs->say("saving the block of size: "+msgs->toStr(blockSize),Low);
						
						fwrite(&blockSize,sizeof(int),1,gf);
						
						for (long i = 0; i <= 0; i++) {	// save only SPH (gas) particles
							n=gadget_header.npart[i];
							if (n>0) {
								data= new float[n];
								for (long j = 0; j < n; j++) {	data[j]=0;	}
								msgs->say("saving "+msgs->toStr(n)+" particles",Low);
								fwrite(data,sizeof(float),n,gf);
								delete [] data;
							}
						}
						
						fwrite(&blockSize,sizeof(int),1,gf);
					}
					break;
					
			} // end of switch
		} // end of loop over blocks
		fclose(gf);
	} // end of loop over files
	
	msgs->say("done",High);
	
	
	return cpedsSuccess;
}

/***************************************************************************************/
cpedsStatusCodes Nbody_io::convertGadgetToGadget2(string fname, string outfile, long Nfiles, bool ic) {
//	gadget_snapshot_file_header_t gadget_header;
	grafic_program_parameters_t grafic_params;
	FILE* fin;
	FILE* fout;
	string fnameIn, fnameOut;
	gadget_snapshot_type2_block_t infoBlk;		
	int usedParticleTypes=0;

	msgs->warning("This is not fully tested...",High);
	
	for (long fi = 0; fi < Nfiles; fi++) {
		usedParticleTypes=0;
		if (Nfiles>1) {
			fnameIn=fname+"."+msgs->toStr(fi);
			fnameOut=outfile+"."+msgs->toStr(fi);
		}
		else {
			fnameIn=fname;
			fnameOut=outfile;
		}
		msgs->say("Input gadget-1 file: "+fnameIn,Medium);
		msgs->say("Output gadget-2 file: "+fnameOut,Medium);
		
		msgs->say("writting Gadget header info block",Medium);
		msgs->say("header size is: %li",long(sizeof(gadget_snapshot_file_header_t)),Medium);
		fout=fopen(fnameOut.c_str(),"wb");
		infoBlk.size=16;
		infoBlk.ID="HEAD"; infoBlk.ID=infoBlk.ID.substr(0,4);
		infoBlk.next=256+8;
		infoBlk.size_check=16;
		writeGadgetInfoBlock(&infoBlk,&fout);
		
		fin=fopen(fnameIn.c_str(),"rb");
		readGadgetHeader(&gadget_header,&fin);
		//
		// save gadget header block
		//
		int tmp=sizeof(gadget_snapshot_file_header_t);
		fwrite(&tmp,sizeof(int),1,fout);
//		printf("tmp: %i\n",tmp);
		fwrite(&gadget_header,sizeof(gadget_snapshot_file_header_t),1,fout);
		fwrite(&tmp,sizeof(int),1,fout);
//		printf("tmp: %i\n",tmp);
		
//		exit(0);
		//
		// write remaining blocks
		//
		
		
		for (long i = 0; i < 6; i++) {	if (gadget_header.npart[i]>0) usedParticleTypes++;	}
		
		long blockSize,blockSize2,Npart;
		gadget_IC_blocks block;
		for (block = Coordinates; block <= InternalEnergy; block=gadget_IC_blocks(block+1)) {
			
			switch (block) {
				case Coordinates:
				case Velocities: 
					if (block==Coordinates) 	{ 
						msgs->say("writting Gadget Coordinates block",High);
						infoBlk.ID="POS "; infoBlk.ID=infoBlk.ID.substr(0,4);
					}
					if (block==Velocities) 		{
						msgs->say("writting Gadget Velocities block",High);
						infoBlk.ID="VEL "; infoBlk.ID=infoBlk.ID.substr(0,4);
					}
					//
					// save gadget positions and velocities
					//
					blockSize=getLocalNumberOfParticles(gadget_header)*3*sizeof(float);
					infoBlk.next=blockSize+8;
					writeGadgetInfoBlock(&infoBlk,&fout);					
					
					msgs->say("saving the block of size: %li bytes (%li particles)",blockSize,getLocalNumberOfParticles(gadget_header),Medium);
					fwrite(&blockSize,sizeof(int),1,fout); // initiate the block
					fread(&blockSize2,sizeof(int),1,fin); // initiate the block
					float vec[3];
					for (long pt = 0; pt < usedParticleTypes; pt++) {
						Npart=gadget_header.npart[pt];
						for (long i = 0; i < Npart; i++) {
							fread(&vec[0],sizeof(float),3,fin);
							fwrite(&vec[0],sizeof(float),3,fout);
						}
					}
					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					fread(&blockSize2,sizeof(int),1,fin); // finish the block
					
					break;
					
				case ParticleIDs:
					msgs->say("writting Gadget ParticleIDs block",High);
					infoBlk.ID="ID  "; infoBlk.ID=infoBlk.ID.substr(0,4);
					blockSize=getLocalNumberOfParticles(gadget_header)*sizeof(int);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
					infoBlk.next=blockSize+8;
					writeGadgetInfoBlock(&infoBlk,&fout);					
					fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
					fread(&blockSize2,sizeof(int),1,fin); // initiate the block
					
					int id;
					for (long pt = 0; pt < usedParticleTypes; pt++) {
						Npart=gadget_header.npart[pt];
						for (long i = 0; i < Npart; i++) {							
							fread(&id,sizeof(int),1,fin);
							fwrite(&id,sizeof(int),1,fout);
						}
					}					
					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					fread(&blockSize2,sizeof(int),1,fin); // finish the block
					break;
					
				case Masses: 
					//				if (ic==false) {
					//					msgs->say("writting Gadget Particle masses block",High);
					//					infoBlk.ID="MASS"; infoBlk.ID=infoBlk.ID.substr(0,4);
					//					blockSize=Npart*sizeof(float);
					//					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
					//					infoBlk.next=blockSize+8;
					//					writeGadgetInfoBlock(&infoBlk,&fout);					
					//					fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
					//
					//					float mass;
					//					for (long i = 0; i < Npart; i++) {							
					//						fread(&mass,sizeof(float),1,fin);
					//						fwrite(&mass,sizeof(float),1,fout);
					//					}	
					//					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					//				}
					break;
				case InternalEnergy:
					//				case Rho:
					//			case SmoothingLength:
					if (gadget_header.npart[0]>0) {
						if (block==InternalEnergy) 	{ 
							msgs->say("writting Gadget InternalEnergy block",High);
							infoBlk.ID="U   "; infoBlk.ID=infoBlk.ID.substr(0,4);
						}
						if (block==Rho) {
							msgs->say("writting Gadget Density block",High);
							infoBlk.ID="RHO "; infoBlk.ID=infoBlk.ID.substr(0,4);
						}
						if (block==SmoothingLength) {
							msgs->say("writting Gadget Smoothing Length block",High);					
							infoBlk.ID="HSML"; infoBlk.ID=infoBlk.ID.substr(0,4);
						}
						
						blockSize=Npart*sizeof(float);
						msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
						infoBlk.next=blockSize+8;
						writeGadgetInfoBlock(&infoBlk,&fout);					
						fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
						fread(&blockSize2,sizeof(int),1,fin); // initiate the block
						
						float f;
						if (block==InternalEnergy) 	{ 
//						for (long pt = 0; pt < usedParticleTypes; pt++) {
							Npart=gadget_header.npart[0];
							for (long i = 0; i < Npart; i++) {							
								fread(&f,sizeof(float),1,fin);
								fwrite(&f,sizeof(float),1,fout);
							}
//						}
						}
						fwrite(&blockSize,sizeof(int),1,fout); // finish the block
						fread(&blockSize2,sizeof(int),1,fin); // finish the block
						
					}
					break;
					
			} // end of switch
		} // end of loop over blocks
		
		
		fclose(fin);
		fclose(fout);
		msgs->say("done",High);
		
	}	
	
	return cpedsSuccess;
	
}
/***************************************************************************************/
void Nbody_io::allocate_gadget_memory() {
	
	P = new gadget_particle_data_t[NumPart];
	//	if(!(P=(gadget_particle_data_t)malloc(NumPart*sizeof(gadget_particle_data_t)))) {
	//		fprintf(stderr,"failed to allocate memory.\n");
	//		exit(0);
	//	}
	
	P--;   /* start with offset 1 */
	
	
	Id=new int[NumPart];
	
	//	if(!(Id=malloc(NumPart*sizeof(int)))) {
	//		fprintf(stderr,"failed to allocate memory.\n");
	//		exit(0);
	//	}
	
	Id--;   /* start with offset 1 */	
}
/***************************************************************************************/
void Nbody_io::free_gadget_memory() {
	if (P!=NULL) { delete [] P; P=NULL; }
	if (Id!=NULL) { delete [] Id; Id=NULL; }
	
}
/***************************************************************************************/
void Nbody_io::free_tipsy_memory() {
	if (tipsyPdark!=NULL) { delete [] tipsyPdark; tipsyPdark=NULL; }
	if (tipsyPgas!=NULL) { delete [] tipsyPgas; tipsyPgas=NULL; }
}
/***************************************************************************************/
Nbody_io::grafic_level0_header_t Nbody_io::read_grafic_header(string fname) {
	
	int32_t w;
	grafic_level0_header_t header;
	
	FILE* f;
	
	//	msgs->say("reading file:"+fname,High);
	f=fopen(fname.c_str(),"rb");
	fread(&w,sizeof(w),1,f);
	fread(&header,sizeof(grafic_level0_header_t),1,f);
	fread(&w,sizeof(w),1,f);
	fclose(f);
	
	return header;
	
}
/***************************************************************************************/
void Nbody_io::print_grafic_header_info(grafic_level0_header_t header) {
	msgs->say("Grafic header information", Top);
	printf("grid size: %li x %li x %li\n", (long)header.n1, (long)header.n2, (long)header.n3);
	printf("grid offsets: dx %f,  dy %f, dz %f\n", header.o1, header.o2, header.o3);
	printf("grid spacing: dx %f\n", header.dx);
	printf("cosmology:\n");
	printf("Om %f\n", header.omegam);
	printf("Ov %f\n", header.omegav);
	printf("a_start %f\n", header.astart);
	printf("z_start %f\n", 1.0/header.astart-1);
	printf("H0 %f\n", header.H0);
	printf("\n\n");
	
}
/***************************************************************************************/
double* Nbody_io::read_grafic_slice(FILE** f, long *cellsInSlice) {
	uint32_t sliceSize, sliceSize2;
	bool dataTypeDouble;
	double *slice=NULL;
	
	fread(&sliceSize,sizeof(uint32_t),1,*f);
	dataTypeDouble=grafic_data_is_double(grafic_header,sliceSize);
	if (dataTypeDouble) { *cellsInSlice=sliceSize/8; } else { *cellsInSlice=sliceSize/4; }
	
	slice = new double[*cellsInSlice];
	if (dataTypeDouble) {
		fread(slice,sizeof(double),*cellsInSlice,*f);
	}
	else {
		float *slicef=new float[*cellsInSlice];
		fread(slicef,sizeof(float),*cellsInSlice,*f);
		for (long i = 0; i < *cellsInSlice; i++) {	slice[i]=slicef[i];	}
		delete [] slicef;
	}
	
	fread(&sliceSize2,sizeof(uint32_t),1,*f);
	if (sliceSize!=sliceSize2) { msgs->criticalError("Internal file error. The slice sizes at the beginning and at the end of the data slice don't match.",High); }

	return slice;
}
/***************************************************************************************/
bool Nbody_io::grafic_data_is_double(grafic_level0_header_t h, uint32_t sliceSize) {
	bool iDouble=false;
//	printf("sliceSize: %li\n", sliceSize);
//	printf("n1 n2: %li %li\n", h.n1,h.n2);

	if (sliceSize/(h.n1*h.n2)==4) { iDouble=false; } 
	else {
		if (sliceSize/(h.n1*h.n2)==8) { iDouble=true; }
		else msgs->criticalError("wrong data type size estimation.", High);
	}
	return iDouble;
}
/***************************************************************************************/
bool Nbody_io::grafic_data_is_double(string fname) {
	uint32_t sliceSize;
	grafic_level0_header_t h=read_grafic_header(fname);
	FILE* f=fopen(fname.c_str(),"rb");	
	fseek(f,get_grafic_header_size(),SEEK_SET); 
	fread(&sliceSize,sizeof(uint32_t),1,f);
	fclose(f);
	return grafic_data_is_double(h,sliceSize);
}
/***************************************************************************************/
double Nbody_io::getVelocityFactor(string dirname) {
	string fname=dirname+"/velocity_factor";
	FILE* f=fopen(fname.c_str(),"r");
	if (f==NULL) { msgs->criticalError("no such file: "+fname,High); }
	double vf=0;
	fscanf(f,"%lE",&vf);
	fclose(f);
	return vf;
}

/***************************************************************************************/
Nbody_io::grafic_program_parameters_t Nbody_io::read_grafic_parameters(string dirname) {
	grafic_program_parameters_t params;
	string fname=dirname+"/program_parameters";
	FILE* f=fopen(fname.c_str(),"rb");
	if (f==NULL) { msgs->criticalError("no such file: "+fname,High); }
	fread(&params,sizeof(grafic_program_parameters_t),1,f);
	fclose(f);
	return params;
}

/***************************************************************************************/
void Nbody_io::calculateParticleMasses() {
	double volumePerParticle;
//	double z=grafic_program_parameters.zstart;
//	double opz3=z+1.0;
//	opz3=opz3*opz3*opz3;
	double h=grafic_cosmo_parameters.H0/100.0;
	double dx=grafic_header.dx * CPEDS_MPC; // cell length in meters
	volumePerParticle=dx*dx*dx *h*h*h;
	double rhoCrit0=cpeds_rhoC0(); // km/m^3 h^-2
	printf("rho critical today is: %lE [kg/m^3 h^2]\n",rhoCrit0);
	printf("volume of a cell today is [m^3 h^-3]: %lE\n",volumePerParticle);
	printf("mass within that volume is [kg h^-1]: %lE\n",volumePerParticle*grafic_cosmo_parameters.Wb*cpeds_rhoC0());
	
	baryonParticleMass=volumePerParticle*grafic_cosmo_parameters.Wb*rhoCrit0 / gadgetMassUnit;
	CDMparticleMass=volumePerParticle*grafic_cosmo_parameters.Wcdm*rhoCrit0 / gadgetMassUnit ;
	printf("baryon particle mass in gadget units: %lf\n",baryonParticleMass);
	printf("CDM particle mass in gadget units: %lf\n",CDMparticleMass);
//	baryonParticleMass=volumePerParticle*grafic_cosmo_parameters.Wb*opz3*cpeds_rhoC0(h) / gadgetMassUnit;	
//	CDMparticleMass=volumePerParticle*grafic_cosmo_parameters.Wcdm*opz3*cpeds_rhoC0(h) / gadgetMassUnit;

}
/***************************************************************************************/
Nbody_io::grafic_cosmoData_t Nbody_io::read_grafic_cosmo_parameters(string dirname) {
	grafic_cosmoData_t params;
	string fname=dirname+"/cosmological_parameters";
	FILE* f=fopen(fname.c_str(),"rb");
	if (f==NULL) { msgs->criticalError("no such file: "+fname,High); }
	fread(&params,sizeof(grafic_cosmoData_t),1,f);
	fclose(f);
	return params;	
}

/***************************************************************************************/
void Nbody_io::printGadgetInfoBlock(gadget_snapshot_type2_block_t* b) {
	msgs->say("Gadget-type-2 infomation block",High);
	msgs->say("Block ID: "+b->ID,Medium);
	msgs->say("Next block in "+msgs->toStr(b->next)+" bytes",Medium);
}
/***************************************************************************************/
Nbody_io::gadget_snapshot_file_header_t Nbody_io::readGadgetHeader(string fname, Nbody_io::gadget_snapshot_type2_block_t &infoBlk) {
	FILE* f;

	f=fopen(fname.c_str(),"rb");			if (f==NULL) { msgs->criticalError("cannot open file: "+fname,High); }
	
	gadget_snapshot_file_header_t gadget_header;
	infoBlk=readGadgetSnapshotType2InfoBlock(&f);						
	readGadgetHeader(&gadget_header,&f);
	fclose(f);
	return gadget_header;
	
}
/***************************************************************************************/
float* Nbody_io::readGadgetBlock(string fname, string blockName, int particleType, 	long* npart, Nbody_io::gadget_snapshot_type2_block_t &infoBlk) {
	FILE* f;
	int status;
	float* data;
	long endPos;
	int length;

	f=fopen(fname.c_str(),"rb");			if (f==NULL) { msgs->criticalError("cannot open file: "+fname,High); }
	fseek(f,0,SEEK_END);
	endPos=ftell(f);
//	printf("end of file position: %li\n",endPos);
	fseek(f,0,SEEK_SET);
//	printf("pozycja: %li\n",ftell(f));
	
	// read individual blocks
//	if (blockName=="HEAD") {
		infoBlk=readGadgetSnapshotType2InfoBlock(&f);						
		readGadgetHeader(&gadget_header,&f);
//		fclose(f);
//		return NULL;
//	}
	
	status=0;
	while (status==0) {
		infoBlk=readGadgetSnapshotType2InfoBlock(&f);
//		printGadgetInfoBlock(&infoBlk);
		if (infoBlk.ID==blockName) {
//			msgs->say("Found the requested block: "+blockName,High);
			fread(&length, sizeof(length), 1, f);

			status=1;
			if (particleType==1) { // seek to next type of particles
				fseek(f,gadget_header.npart[0]*infoBlk.sizeFactor*sizeof(float),SEEK_CUR);
			}
			
			*npart=gadget_header.npart[particleType]*infoBlk.sizeFactor;
			data=new float[*npart];
			fread(data, sizeof(float), *npart, f);
			fclose(f);
			return data;				
		}
		else {
//			msgs->say("The block ID doesn't match: "+infoBlk.ID,Medium);
			// seek to the next block
			status=fseek(f,infoBlk.next,SEEK_CUR);
//			if (status!=0) { msgs->criticalError("readGadgetBlock: block not found", High); }
//			if (status!=0) { return NULL; }
//			printf("pozycja: %li\n",ftell(f));
			*npart=0;
			if (ftell(f)>=endPos) return NULL;
		}
	}
	
	fclose(f);
	return NULL;
}
/***************************************************************************************/
mscsVector<float> Nbody_io::readGadgetBlockVec(string fname, long Nfiles, string blockName, int particleType, gadget_snapshot_type2_block_t &infoBlk) {
	float* data=NULL;
	mscsVector<float> dataAll;
	string inputFile;
	long npart;
	for (unsigned long l = 0; l < Nfiles; l++) {
		if (Nfiles==1)
			inputFile=fname;
		else
			inputFile=fname+"."+msgs->toStr(long(l));
		msgs->say("reading data from file: "+inputFile,Low);
		data=readGadgetBlock(inputFile,blockName,particleType,&npart,infoBlk);		

		if (blockName=="MASS") {
			if (data==NULL) {
				npart=gadget_header.npart[particleType];
				data=new float[npart];
				double m=0;
				m=gadget_header.mass[particleType];
				for (unsigned long i = 0; i < npart; i++) {
					data[i]=m;
				}
			}
		}
		
	
		for (unsigned long i = 0; i < npart; i++) {
			dataAll.push_back(data[i]);
		}
		if (data!=NULL)	delete [] data;
	}
	return dataAll;
}


/***************************************************************************************/
cpedsStatusCodes Nbody_io::dumpGadgetSnapshotData(string fname, int files, int fileType, gadget_IC_blocks data, string outputDir) {
//	fname=dirname+"/cosmological_parameters";
	FILE* f;
	gadget_snapshot_type2_block_t infoBlk;
	string datas=gadgetBlock2BlockIDStr(data);
	int status;
	string fnametmp;
	msgs->say("Will search for block with ID name:",High);
	printf(">%4.4s<\n",datas.c_str());
	for (long i = 0; i < files; i++) {
		msgs->say("PROCESSING FILE: "+msgs->toStr(i),Top);
		if(files>1)	{ fnametmp=fname+"."; fnametmp+=msgs->toStr(i); }
		else fnametmp=fname;
		
		f=fopen(fnametmp.c_str(),"rb");
		if (f==NULL) { msgs->criticalError("cannot open file: "+fnametmp,High); }

		if (fileType==2) {
			infoBlk=readGadgetSnapshotType2InfoBlock(&f);
//			printf("ID string length is :%i\n",infoBlk.ID.length());
			printGadgetInfoBlock(&infoBlk);
		}
	
		readGadgetHeader(&gadget_header,&f);
		printGadgetHeaderInfo(&gadget_header);

		// seek to the requested block in the file
		if (fileType==2) {
			status=0;
			while (status==0) {
				infoBlk=readGadgetSnapshotType2InfoBlock(&f);
				printGadgetInfoBlock(&infoBlk);
				if (infoBlk.ID==datas) {
					msgs->say("Found the requested block",High);
					status=1;
					msgs->say("Dumping data from this block to txt file",High);
					dumpGadgetSnapshotBlock(&f,outputDir+fnametmp,data);
				}
				else {
					msgs->say("The block ID doesn't match",Medium);
					// seek to the next block
					status=fseek(f,infoBlk.next,SEEK_CUR);
					if (status!=0) { msgs->criticalError("block not found", High); }
				}
			}
		}
		else {
			msgs->criticalError("Sorry, this type is not yet supported in this method",High);
		}
	}
	return cpedsSuccess;
}

/***************************************************************************************/
cpedsStatusCodes Nbody_io::dumpGadgetSnapshot2Tipsy(string fname, int files) {
	//	fname=dirname+"/cosmological_parameters";
	FILE* tf; // tipsy file descriptor
	gadget_snapshot_type2_block_t infoBlk;
	string datas;
	int status;
	string fnametmp;
	float* pos;
	float* vel;
	float* mass;
	float* rho;
	float* temp;
	long Npart;
	
	tipsy_dark_particle_t DMpart;
	tipsy_gas_particle_t GASpart;
	tipsy_star_particle_t STARpart;
	
	if(files>1)	{ fnametmp=fname+"."; fnametmp+=msgs->toStr(0); }
	else fnametmp=fname;
	readGadgetBlock(fnametmp,"HEAD",-1,&Npart,infoBlk); // this initializes gadget_header variable
	//
	// generate tipsy header block
	//
	tipsyHeader_t tipsyHeader;
	tipsyHeader.nbodies=gadget_header.npartTotal[0]+gadget_header.npartTotal[1];
	tipsyHeader.ndark=gadget_header.npartTotal[1];
	tipsyHeader.nsph=gadget_header.npartTotal[0];
	tipsyHeader.nstar=0;
	tipsyHeader.ndim=3;
	tipsyHeader.time=gadget_header.time;

	string fnametipsy=fname+".tipsy";
	tf=fopen(fnametipsy.c_str(),"wb");
	fwrite(&tipsyHeader,sizeof(tipsyHeader_t),1,tf);
	
	//
	// generate the three blocks of the tipsy native file format: sph block, dark particles bloc and stars block
	//
	
	
	for (long ptype = 0; ptype < 2; ptype++) { // loop over particle types; 0-gas, 1-halo
		 
		for (long i = 0; i < files; i++) { // loop over snapshot files

			msgs->say("PROCESSING FILE: "+msgs->toStr(i),Top);
			if(files>1)	{ fnametmp=fname+"."; fnametmp+=msgs->toStr(i); }
			else fnametmp=fname;
			
			// read particles position vector
			pos=readGadgetBlock(fnametmp,"POS ",ptype,&Npart,infoBlk);
			// read particles velocities vector
			vel=readGadgetBlock(fnametmp,"VEL ",ptype,&Npart,infoBlk);
//			// read particles mass vector
//			mass=readGadgetBlock(fname,"MASS",ptype,&Npart,infoBlk);
			// read particles density vector
			rho=readGadgetBlock(fnametmp,"RHO ",ptype,&Npart,infoBlk);

			// generate tipsy block
			if (ptype==0) {
				for (long i = 0; i < Npart; i++) {
					GASpart.hsmooth=0;
					GASpart.mass=gadget_header.mass[ptype];
					GASpart.phi=0;
					GASpart.pos[0]=pos[3*i];
					GASpart.pos[1]=pos[3*i+1];
					GASpart.pos[2]=pos[3*i+2];
					GASpart.vel[0]=vel[3*i];
					GASpart.vel[1]=vel[3*i+1];
					GASpart.vel[2]=vel[3*i+2];
					GASpart.metals=0;
					GASpart.rho=rho[i];
					GASpart.temp=0;
					fwrite(&GASpart,sizeof(tipsy_gas_particle_t),1,tf);	
				}			
			}
			if (ptype==1) {
				for (long i = 0; i < Npart; i++) {
					DMpart.eps=0;
					DMpart.mass=gadget_header.mass[ptype];
					DMpart.phi=0;
					DMpart.pos[0]=pos[3*i];
					DMpart.pos[1]=pos[3*i+1];
					DMpart.pos[2]=pos[3*i+2];
					DMpart.vel[0]=vel[3*i];
					DMpart.vel[1]=vel[3*i+1];
					DMpart.vel[2]=vel[3*i+2];
					fwrite(&DMpart,sizeof(tipsy_dark_particle_t),1,tf);	
				}
			}
			delete [] pos;
			delete [] vel;
			delete [] rho;
		}
	}
	
	fclose(tf);

/*
			f=fopen(fnametmp.c_str(),"rb");			if (f==NULL) { msgs->criticalError("cannot open file: "+fnametmp,High); }
			
			infoBlk=readGadgetSnapshotType2InfoBlock(&f);				printGadgetInfoBlock(&infoBlk);
			
			readGadgetHeader(&gadget_header,&f);			printGadgetHeaderInfo(&gadget_header);
			
			for (gadget_IC_blocks gblock = MIN; gblock < MAX; var++) {				
				datas=gadgetBlock2BlockIDStr(gblock);
				// seek to the requested block in the file
				status=0;
				while (status==0) {
					infoBlk=readGadgetSnapshotType2InfoBlock(&f);
					printGadgetInfoBlock(&infoBlk);
					if (infoBlk.ID==datas) {
						msgs->say("Found the requested block",High);
						status=1;
						msgs->say("Dumping data from this block to txt file",High);
						
						
//						if (ptype==1) { // seek to next type of particles
//							fseek(f,gadget_header.npart[0]*infoBlk.sizeFactor*sizeof(float),SEEK_CUR);
////							fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
//						}
//						
//
//						if (infoBlk.sizeFactor==1) { // get vector data 
//						
//						}
//						if (infoBlk.sizeFactor==3) { // get 3 x vector data 
//						
//						}
					}
					else {
						msgs->say("The block ID doesn't match",Medium);
						// seek to the next block
						status=fseek(f,infoBlk.next,SEEK_CUR);
						if (status!=0) { msgs->criticalError("block not found", High); }
					}
				}
				
				
				
			}

*/
	return cpedsSuccess;
	
}

/***************************************************************************************/
cpedsStatusCodes Nbody_io::loadGadgetSnapshotOrd(string fname, long files) {
	//	fname=dirname+"/cosmological_parameters";
	gadget_snapshot_type2_block_t infoBlk;
	string datas;
	int status;
	string fnametmp;
	float* pos;
	float* vel;
	float* mass;
	float* rho;
	float* temp;
	long Npart;
	
	
	if(files>1)	{ fnametmp=fname+"."; fnametmp+=msgs->toStr(0); }
	else fnametmp=fname;
	readGadgetBlock(fnametmp,"HEAD",-1,&Npart,infoBlk); // this initializes gadget_header variable

	if (P!=NULL) delete [] P;
	P = new gadget_particle_data_t[gadget_header.npartTotal[0]+gadget_header.npartTotal[1]];
	long pid=0;
	
	for (long ptype = 0; ptype < 2; ptype++) { // loop over particle types; 0-gas, 1-halo
		 
		for (long i = 0; i < files; i++) { // loop over snapshot files

			msgs->say("PROCESSING FILE: "+msgs->toStr(i),Top);
			if(files>1)	{ fnametmp=fname+"."; fnametmp+=msgs->toStr(i); }
			else fnametmp=fname;
			
			// read particles position vector
			pos=readGadgetBlock(fnametmp,"POS ",ptype,&Npart,infoBlk);
			// read particles velocities vector
			vel=readGadgetBlock(fnametmp,"VEL ",ptype,&Npart,infoBlk);
//			// read particles mass vector
//			mass=readGadgetBlock(fname,"MASS",ptype,&Npart,infoBlk);
			// read particles density vector
			rho=readGadgetBlock(fnametmp,"RHO ",ptype,&Npart,infoBlk);

			for (long i = 0; i < Npart; i++) {
				P[pid].Mass=gadget_header.mass[ptype];
				P[pid].Ne=0;
				P[pid].Pos[0]=pos[3*i];
				P[pid].Pos[1]=pos[3*i+1];
				P[pid].Pos[2]=pos[3*i+2];
				P[pid].Vel[0]=vel[3*i];
				P[pid].Vel[1]=vel[3*i+1];
				P[pid].Vel[2]=vel[3*i+2];
				P[pid].Rho=rho[i];
				P[pid].Temp=0;
				P[pid].U=0;
				pid++;
			}
			delete [] pos;
			delete [] vel;
			delete [] rho;
		}
	}

	return cpedsSuccess;

}

/***************************************************************************************/
int Nbody_io::readGadgetHeader(gadget_snapshot_file_header_t* gadget_header, FILE** f) {
	//
	// read in the header information from the file
	// 
	int length,length2;
	fread(&length, sizeof(length), 1, *f);
	fread(gadget_header, sizeof(gadget_snapshot_file_header_t), 1, *f);
	fread(&length2, sizeof(length2), 1, *f);
	if (length!=length2) {
		msgs->error("the sizes of the block opening flag and closing flag do not match",High);
		return 0;
	}
	else { return 1; }
	
}


/***************************************************************************************/
void Nbody_io::printGadgetHeaderInfo(gadget_snapshot_file_header_t* gadget_header) {
	int i;
	msgs->say("Gadget Snapshot file header information", High);
	msgs->say("Box size: "+msgs->toStr(gadget_header->BoxSize),Medium);
	msgs->say("H0: "+msgs->toStr(gadget_header->HubbleParam),Medium);
	msgs->say("Omega0: "+msgs->toStr(gadget_header->Omega0),Medium);
	msgs->say("Omega Lambda: "+msgs->toStr(gadget_header->OmegaLambda),Medium);
	
	for (i = 0; i < 6; ++i) {
//		printf("\n");
		msgs->say("  npart for type "+msgs->toStr(i)+": "+msgs->toStr(gadget_header->npart[i]),Low);
		msgs->say("  npart total for type "+msgs->toStr(i)+": "+msgs->toStr(long(gadget_header->npartTotal[i])),Low);
		msgs->say("  mass for type "+msgs->toStr(i)+": "+msgs->toStr(gadget_header->mass[i]),Low);
	}
//	printf("\n");
	msgs->say("num files: "+msgs->toStr(gadget_header->num_files),Medium);
	msgs->say("redshift: "+msgs->toStr(gadget_header->redshift),Medium);
	msgs->say("time: "+msgs->toStr(gadget_header->time),Medium);
	msgs->say("cooling flag: "+msgs->toStr(gadget_header->flag_cooling),Medium);
	msgs->say("feedback flag: "+msgs->toStr(gadget_header->flag_feedback),Medium);
	msgs->say("sfr flag: "+msgs->toStr(gadget_header->flag_sfr),Medium);
	
	
	
//	printf("\n\n");
}
/***************************************************************************************/
Nbody_io::gadget_snapshot_type2_block_t Nbody_io::readGadgetSnapshotType2InfoBlock(FILE** f) {
	gadget_snapshot_type2_block_t b;
	char id[4];
	fread(&b.size,sizeof(int),1,*f);
	fread(&id[0],sizeof(char),4,*f);
	fread(&b.next,sizeof(int),1,*f);
	fread(&b.size_check,sizeof(int),1,*f);
	b.ID=id;
	b.ID=b.ID.substr(0,4);
	b.sizeFactor=1;
	if (b.ID=="POS " or b.ID=="VEL ") b.sizeFactor=3;
//	printf("BEBUG check: If the string between >< contains no space then correct the code in Nbody_io::readGadgetSnapshotType2InfoBlock(FILE** f)\n");
//	printf("b.ID: >%s<\n",b.ID.c_str());
//	exit(0);
	return b;
}
/***************************************************************************************/
string Nbody_io::gadgetBlock2BlockIDStr(gadget_IC_blocks data) {
	string s="";
	switch (data) {
		case Header:
			s="HEAD";
			break;
		case Coordinates:
			s="POS ";
			break;
		case Velocities:
			s="VEL ";
			break;
		case ParticleIDs:
			s="ID  ";
			break;
		case Masses:
			s="MASS";
			break;
		case InternalEnergy:
		case Temperature:
			s="U   ";
			break;
		case Rho:
			s="RHO ";
			break;
		case SmoothingLength:
			s="HSML";
			break;
		case Acceleration:
			s="ACCE";
			break;
		default:
			s="";
			break;
	}
	return s;
}
/***************************************************************************************/
void Nbody_io::dumpGadgetSnapshotBlock(FILE** fd, string outputFile, gadget_IC_blocks block) {
	string fname;
	FILE* f;
	int blockSize;
	int id;
	float vecd[3];
	int n;

//			double MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * CPEDS_m_p;
	double MeanWeight= 4.0/(3*CPEDS_hydrogen_mass_fraction+1             ) * CPEDS_m_p *1e3; // 1e3 - conversion to cgs
//			double mu=MeanWeight/BOLTZMANN * (gamma-1) * u;
	double gamma= 5.0/3;
	double u,T;

	
	
	
	
	
	
	fread(&blockSize,sizeof(int),1,*fd); // initiate the bloc read
	
	switch (block) {
		case Coordinates:
		case Velocities: 
			if (block==Coordinates) 	msgs->say("dumping Gadget Coordinates block",High);
			if (block==Velocities) 		msgs->say("dumping Gadget Velocities block",High);
			msgs->say("the block size is: "+msgs->toStr(blockSize),Medium);
			//
			// save gadget positions and velocities
			//
				
			for (long i = 0; i < 6; i++) {		
				if (block==Coordinates) fname=outputFile+".Coordinates"+".Ptype_"+msgs->toStr(i);
				if (block==Velocities) fname=outputFile+".Velocities"+".Ptype_"+msgs->toStr(i);
				f=fopen(fname.c_str(), "w");
				msgs->say("processing particle type: "+msgs->toStr(i),Medium);
				n=gadget_header.npart[i];
				msgs->say("number of particles in this type: "+msgs->toStr(n),Medium);
				if (n>0) {
					for (long j = 0; j < n; j++) {
						fread(&vecd[0],sizeof(float),3,*fd);
						fprintf(f,"%E %E %E\n",vecd[0], vecd[1], vecd[2]);
					}
				}
				fclose(f);
			}
		
	
			break;
			
		case ParticleIDs:
			msgs->say("dumping Gadget ParticleIDs block to file: "+outputFile,High);
			msgs->say("the block size is: "+msgs->toStr(blockSize),Medium);
			for (long i = 0; i < 6; i++) {					
				fname=outputFile+".ParticleIDs"+".Ptype_"+msgs->toStr(i);
				f=fopen(fname.c_str(), "w");
				msgs->say("processing particle type: "+msgs->toStr(i),Medium);
				n=gadget_header.npart[i];
				msgs->say("number of particles in this type: "+msgs->toStr(n),Medium);
				if (n>0) {
					for (long j = 0; j < n; j++) {
						fread(&id,sizeof(int),1,*fd);
						fprintf(f,"%i\n",id);
					}
				}
			}
			break;
		case Rho:
			msgs->say("dumping Gadget densities block",High);
			msgs->say("the block size is: "+msgs->toStr(blockSize),Medium);
			//
			// save gadget SPH densities
			//
				
			for (long i = 0; i < 6; i++) {		
				fname=outputFile+".Rho"+".Ptype_"+msgs->toStr(i);
				f=fopen(fname.c_str(), "w");
				msgs->say("processing particle type: "+msgs->toStr(i),Medium);
				n=gadget_header.npart[i];
				msgs->say("number of particles in this type: "+msgs->toStr(n),Medium);
				if (n>0) {
					for (long j = 0; j < n; j++) {
						fread(&vecd[0],sizeof(float),1,*fd);
						fprintf(f,"%E\n",vecd[0]);
					}
				}
				fclose(f);
			}
			
			break;
		case InternalEnergy:
		case Temperature:
			if (block==InternalEnergy) msgs->say("dumping Gadget internal energy block",High);
			if (block==Temperature) msgs->say("dumping Gadget temperatures derived from internal energy ",High);
			msgs->say("the block size is: "+msgs->toStr(blockSize),Medium);
			//
			// save gadget SPH internal energy
			//
			
			for (long i = 0; i < 6; i++) {		
				if (block==InternalEnergy) fname=outputFile+".U"+".Ptype_"+msgs->toStr(i);
				if (block==Temperature) fname=outputFile+".T"+".Ptype_"+msgs->toStr(i);
				f=fopen(fname.c_str(), "w");
				msgs->say("processing particle type: "+msgs->toStr(i),Medium);
				n=gadget_header.npart[i];
				msgs->say("number of particles in this type: "+msgs->toStr(n),Medium);
				if (n>0) {
					msgs->say("This is NOT DEBUGED: PLEASE CHECK THE CODE",Top);
					if (block==InternalEnergy)
						for (long j = 0; j < n; j++) {
							fread(&vecd[0],sizeof(float),1,*fd);
							/* convert internal energy to keV units */
							u  = vecd[0] * GADGET2_UnitEnergy_in_cgs * 1e-4 * CPEDS_e; // * 1e-7 * CPEDS_e * 1e3 (conversion to Joules then to eV and then to keV) = 
							fprintf(f,"%E\n",vecd[0]);
						}
					if (block==Temperature) {
						for (long j = 0; j < n; j++) {
							/* convert internal energy to keV units */
							u  = vecd[0] * GADGET2_UnitEnergy_in_cgs / GADGET2_UnitMass_in_g; // ergs
//							MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * CPEDS_m_p;
							T= MeanWeight/(CPEDS_kB * 1e7) * (gamma-1) * u; // 1e7 - conversion from J/K to erg/K
							fread(&vecd[0],sizeof(float),1,*fd);
							fprintf(f,"%E\n",T);
						}
					}
				}
				fclose(f);
			}
			
			break;
			
		case SmoothingLength:
			msgs->say("dumping Gadget hsml block",High);
			msgs->say("the block size is: "+msgs->toStr(blockSize),Medium);
			//
			// save gadget SPH hsml
			//
				
			for (long i = 0; i < 6; i++) {		
				fname=outputFile+".hsml"+".Ptype_"+msgs->toStr(i);
				f=fopen(fname.c_str(), "w");
				msgs->say("processing particle type: "+msgs->toStr(i),Medium);
				n=gadget_header.npart[i];
				msgs->say("number of particles in this type: "+msgs->toStr(n),Medium);
				if (n>0) {
					for (long j = 0; j < n; j++) {
						fread(&vecd[0],sizeof(float),1,*fd);
						fprintf(f,"%E\n",vecd[0]);
					}
				}
				fclose(f);
			}
			
			break;
//		case Coordinates:
//			
//			break;
//		case Coordinates:
//			
//			break;
//		case Coordinates:
//			
//			break;
//		case Coordinates:
//			
//			break;
//		case Coordinates:
//			
//			break;
		default:
			msgs->error("Sorry, you didn't implement saving this block yet",High);
			break;
	}
	
	fread(&blockSize,sizeof(int),1,*fd); // finish the bloc read
	
//	for (block = Coordinates; block <= SmoothingLength; block=gadget_IC_blocks(block+1)) {
//		
//
//		switch (block) {
//
//			case Masses: // do not write the masses block. The masses will be set from the header information. Same for all particles of a given type
//					msgs->say("skipping Gadget Masses block",High);
//				break;
//				
//			case InternalEnergy:
////			case Rho:
////			case SmoothingLength:
//				if (gadget_header.npart[0]>0) {
//					if (block==InternalEnergy) 	msgs->say("writting Gadget InternalEnergy block",High);
//					if (block==Rho) 	msgs->say("writting Gadget Density block",High);
//					if (block==SmoothingLength) 	msgs->say("writting Gadget Smoothing Length block",High);					
//					
//					blockSize=grafic_header.n1*grafic_header.n2*grafic_header.n3*sizeof(float);
//					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
//					
//					fwrite(&blockSize,sizeof(int),1,gf);
//
//					for (long i = 0; i <= 0; i++) {	// save only SPH (gas) particles
//						n=gadget_header.npart[i];
//						if (n>0) {
//							data= new float[n];
//							for (long j = 0; j < n; j++) {	data[j]=0;	}
//							msgs->say("saving "+msgs->toStr(n)+" particles",Medium);
//							fwrite(data,sizeof(float),n,gf);
//							delete [] data;
//						}
//					}
//					
//					fwrite(&blockSize,sizeof(int),1,gf);
//				}
//				break;
//	
//		} // end of switch
//	} // end of loop over blocks


}
/***************************************************************************************/
void Nbody_io::writeGadgetInfoBlock(gadget_snapshot_type2_block_t* b, FILE** f) {
	char id[4];
	strcpy(id,b->ID.c_str());
	fwrite(&b->size,sizeof(int),1,*f);
	fwrite(&id[0],sizeof(char),4,*f);
	fwrite(&b->next,sizeof(int),1,*f);
	fwrite(&b->size_check,sizeof(int),1,*f);	
}
/***************************************************************************************/
long Nbody_io::getTotalNumberOfParticles(const gadget_snapshot_file_header_t& header) {
	long n=0;
	for (long i = 0; i < 6; i++) {
		n+=header.npartTotal[i];
	}
	return n;
}
/***************************************************************************************/
long Nbody_io::getLocalNumberOfParticles(const gadget_snapshot_file_header_t& header) {
	long n=0;
	for (long i = 0; i < 6; i++) {
		n+=header.npart[i];
	}
	return n;
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::writeGadget2data(matrix<double>& data, string fname, double boxSize, int particleType) {
	gadget_snapshot_file_header_t gadget_header;
	FILE* fout;
	gadget_snapshot_type2_block_t infoBlk;		
	int usedParticleTypes=0;
	
	msgs->say("writting Gadget header info block",High);
	fout=fopen(fname.c_str(),"wb");
	infoBlk.size=16;
	infoBlk.ID="HEAD"; infoBlk.ID=infoBlk.ID.substr(0,4);
	infoBlk.next=256+8;
	infoBlk.size_check=16;
	writeGadgetInfoBlock(&infoBlk,&fout);
	
	gadget_header.BoxSize=boxSize;
	gadget_header.HubbleParam=0;
	gadget_header.Omega0=0;
	gadget_header.OmegaLambda=0;
	gadget_header.flag_cooling=0;
	gadget_header.flag_entropy_instead_u=0;
	gadget_header.flag_feedback=0;
	gadget_header.flag_metals=0;
	gadget_header.flag_sfr=0;
	gadget_header.flag_stellarage=0;
	gadget_header.num_files=1;
	gadget_header.redshift=0;
	gadget_header.time=0;
	gadget_header.mass[particleType]=0;
	gadget_header.npart[particleType]=data.RowNo();
	gadget_header.npartTotal[particleType]=data.RowNo();
	for (long i = 0; i < particleType; i++) { 	
		gadget_header.mass[i]=0;	
		gadget_header.npart[i]=0;
		gadget_header.npartTotal[i]=0;
	}
	for (long i = particleType+1; i < 6; i++) { 
		gadget_header.mass[i]=0;	
		gadget_header.npart[i]=0;
		gadget_header.npartTotal[i]=0;
	}
	for (long i = 0; i < 6; i++) { gadget_header.npartTotalHighWord[i]=(unsigned int) ((long long)(gadget_header.npart[i]) >> 32); }

	//
	// save gadget header block
	//
	int tmp=sizeof(gadget_snapshot_file_header_t);
	fwrite(&tmp,sizeof(int),1,fout);
	fwrite(&gadget_header,sizeof(gadget_snapshot_file_header_t),1,fout);
	fwrite(&tmp,sizeof(int),1,fout);
	
	
	//
	// write remaining blocks
	//
	
		
	long blockSize,Npart;
	gadget_IC_blocks block;
	Npart=data.RowNo();
	int dataSpecifier=0;
	for (block = Coordinates; block <= InternalEnergy; block=gadget_IC_blocks(block+1)) {
		
		switch (block) {
			case Coordinates:
			case Velocities: 
				if (block==Coordinates) 	{ 
					msgs->say("writting Gadget Coordinates block",High);
					infoBlk.ID="POS "; infoBlk.ID=infoBlk.ID.substr(0,4);
					dataSpecifier=0;
				}
				if (block==Velocities) 		{
					msgs->say("writting Gadget Velocities block",High);
					infoBlk.ID="VEL "; infoBlk.ID=infoBlk.ID.substr(0,4);
					dataSpecifier=3;
				}
				//
				// save gadget positions and velocities
				//
				blockSize=Npart*3*sizeof(float);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the block
				float vec[3];
				for (long i = 0; i < Npart; i++) {							
					vec[0]=data(i,0+dataSpecifier);
					vec[1]=data(i,1+dataSpecifier);
					vec[2]=data(i,2+dataSpecifier);
					fwrite(&vec[0],sizeof(float),3,fout);
				}
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				
				break;
				
			case ParticleIDs:
				msgs->say("writting Gadget ParticleIDs block",High);
				infoBlk.ID="ID  "; infoBlk.ID=infoBlk.ID.substr(0,4);
				blockSize=Npart*sizeof(int);
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
				
				int id;
				for (long i = 0; i < Npart; i++) {							
					id=i;
					fwrite(&id,sizeof(int),1,fout);
				}
				
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				break;
				
			case Masses: 
				msgs->say("writting Gadget Masses block",High);
				infoBlk.ID="MASS"; infoBlk.ID=infoBlk.ID.substr(0,4);
				blockSize=Npart*sizeof(float);
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
				
				float m;
				for (long i = 0; i < Npart; i++) {							
					m=data(i,6);
					fwrite(&m,sizeof(float),1,fout);
				}
				
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				break;
			case InternalEnergy: 
				if (gadget_header.npart[0]>0) {
					if (block==InternalEnergy) {
						msgs->say("writting Gadget InternalEnergy block",High);					
						infoBlk.ID="U   "; infoBlk.ID=infoBlk.ID.substr(0,4);
					}
					
					blockSize=Npart*sizeof(float);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
					infoBlk.next=blockSize+8;
					writeGadgetInfoBlock(&infoBlk,&fout);					
					fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
					
					float f=0;
					for (long i = 0; i < Npart; i++) {							
						fwrite(&f,sizeof(float),1,fout);
					}
					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					
				}
				break;
				
		} // end of switch
	} // end of loop over blocks
	
	
	fclose(fout);
	msgs->say("done",High);
	
	
	return cpedsSuccess;
	
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::writeGadget2data(string fname, Nbody_io::gadget_snapshot_file_header_t& gadget_header, int particleType, long particleNum, float* POS, float* VEL, float* MASS,float* U, float* RHO ) {
//	gadget_snapshot_file_header_t gadget_header;
	FILE* fout;
	gadget_snapshot_type2_block_t infoBlk;		
	int usedParticleTypes=0;
	
	
//	for (long k = 0; k < particleNum; k++) {
//		printf("%li) saving particle: xyz: %f  %f  %f\n",k,POS[3*k],POS[3*k+1],POS[3*k+2]);
//	}

	fname+=".gadget2";
	msgs->say("Writing Gadget header info block in file: "+fname,High);
	fout=fopen(fname.c_str(),"wb");
	infoBlk.size=16;
	infoBlk.ID="HEAD"; infoBlk.ID=infoBlk.ID.substr(0,4);
	infoBlk.next=256+8;
	infoBlk.size_check=16;
	writeGadgetInfoBlock(&infoBlk,&fout);
	
//	gadget_header.BoxSize=boxSize;
//	gadget_header.HubbleParam=0;
//	gadget_header.Omega0=0;
//	gadget_header.OmegaLambda=0;
//	gadget_header.flag_cooling=0;
//	gadget_header.flag_entropy_instead_u=0;
//	gadget_header.flag_feedback=0;
//	gadget_header.flag_metals=0;
//	gadget_header.flag_sfr=0;
//	gadget_header.flag_stellarage=0;
	gadget_header.num_files=1;
//	gadget_header.redshift=0;
	gadget_header.time=0;
//	gadget_header.mass[particleType]=0;
	gadget_header.npart[particleType]=particleNum;
	gadget_header.npartTotal[particleType]=particleNum;
	for (long i = 0; i < particleType; i++) {
		gadget_header.mass[i]=0;	
		gadget_header.npart[i]=0;
		gadget_header.npartTotal[i]=0;
	}
	for (long i = particleType+1; i < 6; i++) { 
		gadget_header.mass[i]=0;	
		gadget_header.npart[i]=0;
		gadget_header.npartTotal[i]=0;
	}
	for (long i = 0; i < 6; i++) { gadget_header.npartTotalHighWord[i]=(unsigned int) ((long long)(gadget_header.npart[i]) >> 32); }

	//
	// save gadget header block
	//
	int tmp=sizeof(gadget_snapshot_file_header_t);
	fwrite(&tmp,sizeof(int),1,fout);
	fwrite(&gadget_header,sizeof(gadget_snapshot_file_header_t),1,fout);
	fwrite(&tmp,sizeof(int),1,fout);
	

	//
	// write remaining blocks
	//
	
		
	long blockSize,Npart;
	gadget_IC_blocks block;
	Npart=particleNum;
	float* dataSpecifier=NULL;
	
//	gadget_IC_blocks endBlock;
	
	for (block = Coordinates; block <= Rho; block=gadget_IC_blocks(block+1)) {
//	for (block = Coordinates; block <= Coordinates; block=gadget_IC_blocks(block+1)) {
		
		switch (block) {
			case Coordinates:
			case Velocities: 
				if (block==Coordinates) 	{ 
					msgs->say("Writing Gadget Coordinates block",High);
					infoBlk.ID="POS "; infoBlk.ID=infoBlk.ID.substr(0,4);
					dataSpecifier=POS;
				}
				if (block==Velocities) 		{
					msgs->say("Writing Gadget Velocities block",High);
					infoBlk.ID="VEL "; infoBlk.ID=infoBlk.ID.substr(0,4);
					dataSpecifier=VEL;
				}
				//
				// save gadget positions and velocities
				//
				blockSize=Npart*3*sizeof(float);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the block
				float vec[3];
				long i3;
				for (long i = 0; i < Npart; i++) {							
					i3=i*3;
					vec[0]=dataSpecifier[i3];
					vec[1]=dataSpecifier[i3+1];
					vec[2]=dataSpecifier[i3+2];
//					printf("saving dataSpec: %li %f %f %f\n",i,vec[0],vec[1],vec[2]);
					fwrite(&vec[0],sizeof(float),3,fout);

				}
//				exit(0);			
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				
				break;
				
			case ParticleIDs:
				msgs->say("Writing Gadget ParticleIDs block",High);
				infoBlk.ID="ID  "; infoBlk.ID=infoBlk.ID.substr(0,4);
				blockSize=Npart*sizeof(int);
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
				
				int id;
				for (long i = 0; i < Npart; i++) {							
					id=i;
					fwrite(&id,sizeof(int),1,fout);
				}
				
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				break;
				
			case Masses: 
				msgs->say("writting Gadget Masses block",High);
				infoBlk.ID="MASS"; infoBlk.ID=infoBlk.ID.substr(0,4);
				blockSize=Npart*sizeof(float);
				msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
				infoBlk.next=blockSize+8;
				writeGadgetInfoBlock(&infoBlk,&fout);					
				fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
				
				float m;
				for (long i = 0; i < Npart; i++) {							
					m=MASS[i];
					fwrite(&m,sizeof(float),1,fout);
				}
				
				fwrite(&blockSize,sizeof(int),1,fout); // finish the block
				break;
			case InternalEnergy: 
				if (gadget_header.npart[0]>0 and U!=NULL) {
					if (block==InternalEnergy) {
						msgs->say("writting Gadget InternalEnergy block",High);					
						infoBlk.ID="U   "; infoBlk.ID=infoBlk.ID.substr(0,4);
					}
					
					blockSize=Npart*sizeof(float);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
					infoBlk.next=blockSize+8;
					writeGadgetInfoBlock(&infoBlk,&fout);					
					fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
					
					float f=0;
					for (long i = 0; i < Npart; i++) {							
						f=U[i];
						fwrite(&f,sizeof(float),1,fout);
					}
					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					
				}
				break;

			case Rho: 
				if (gadget_header.npart[0]>0 and RHO!=NULL) {
					if (block==Rho) {
						msgs->say("writting Gadget Rho block",High);					
						infoBlk.ID="RHO   "; infoBlk.ID=infoBlk.ID.substr(0,4);
					}
					
					blockSize=Npart*sizeof(float);
					msgs->say("saving the block of size: "+msgs->toStr(blockSize),Medium);
					infoBlk.next=blockSize+8;
					writeGadgetInfoBlock(&infoBlk,&fout);					
					fwrite(&blockSize,sizeof(int),1,fout); // initiate the bloc
					
					float f=0;
					for (long i = 0; i < Npart; i++) {							
						f=RHO[i];
						fwrite(&f,sizeof(float),1,fout);
					}
					fwrite(&blockSize,sizeof(int),1,fout); // finish the block
					
				}
				break;
				
		} // end of switch
	} // end of loop over blocks
	
	
	fclose(fout);
	msgs->say("done",High);
	
	
	return cpedsSuccess;

}
/***************************************************************************************/
void Nbody_io::dumpGraficProperVelocities(string dirname, string outfname) {
	
	string fname;
	double velocity_factor;	
	grafic_program_parameters_t grafic_program_parameters;
	grafic_cosmoData_t grafic_cosmo_parameters;
	FILE* f;
	FILE *fx,*fy,*fz;
	double *vx,*vy,*vz;
	long nsx,nsy,nsz; // number of cells in slices in velocity files x,y,z
	long ii;

	//
	// read in the grafic header for baryons
	//
	fname=dirname+"/"+"ic_deltab";
	msgs->say("reading file:"+fname,High);
	grafic_header=read_grafic_header(fname);

	//
	// read in other parameters
	//
	velocity_factor=getVelocityFactor(dirname);
	grafic_program_parameters=read_grafic_parameters(dirname);
	grafic_cosmo_parameters=read_grafic_cosmo_parameters(dirname);
	
	// output some info
	print_grafic_header_info(grafic_header);
	msgs->say("mass per baryon particle in Gadget units: "+msgs->toStr(baryonParticleMass),Medium);
	msgs->say("mass per CDM    particle in Gadget units: "+msgs->toStr(CDMparticleMass),Medium);	
	msgs->say("velocity factor [km/s/Mpc]: "+msgs->toStr(velocity_factor),Medium);
//	msgs->say("mass per baryon particle in Gadget units: "+msgs->toStr(baryonParticleMass),Medium);




	string tmps;

	for (long particles = 0; particles < 2; particles++) { // loop over particle types
		msgs->say("Processing particles of type: "+msgs->toStr(particles),Low);
		
		switch (particles) {
			case 0: 
				tmps=outfname+".b";
				f=fopen(tmps.c_str(),"w");
				fname=dirname+"/"+"ic_velbx"; fx=fopen(fname.c_str(),"rb");	fseek(fx,get_grafic_header_size(),SEEK_SET); 
				fname=dirname+"/"+"ic_velby"; fy=fopen(fname.c_str(),"rb");	fseek(fy,get_grafic_header_size(),SEEK_SET); 
				fname=dirname+"/"+"ic_velbz"; fz=fopen(fname.c_str(),"rb");	fseek(fz,get_grafic_header_size(),SEEK_SET); 
				break;
			case 1:
				tmps=outfname+".cdm";
				f=fopen(tmps.c_str(),"w");
				fname=dirname+"/"+"ic_velcx"; fx=fopen(fname.c_str(),"rb");	fseek(fx,get_grafic_header_size(),SEEK_SET); 
				fname=dirname+"/"+"ic_velcy"; fy=fopen(fname.c_str(),"rb");	fseek(fy,get_grafic_header_size(),SEEK_SET); 
				fname=dirname+"/"+"ic_velcz"; fz=fopen(fname.c_str(),"rb");	fseek(fz,get_grafic_header_size(),SEEK_SET); 
				break;				
		}
		
		for (long k = 0; k < grafic_header.n3; k++) { // loop over slices
			
			vx=read_grafic_slice(&fx,&nsx);
			vy=read_grafic_slice(&fy,&nsy);
			vz=read_grafic_slice(&fz,&nsz);
			
//			if (block==Coordinates) {
//				// derive offsets from velocities by dividing by velocity_factor
//				cpeds_divide_value(velocity_factor,vx,nsx);
//				cpeds_divide_value(velocity_factor,vy,nsy);	
//				cpeds_divide_value(velocity_factor,vz,nsz);	
//				
//				// move particles to their cell centres
//				ii=0;
//				for (long j = 0; j < grafic_header.n2; j++) {
//					for (long i = 0; i < grafic_header.n1; i++) {
//						vx[ii]+=i*dx+dxo2;
//						vy[ii]+=j*dx+dxo2;
//						vz[ii]+=k*dx+dxo2;
//						ii++;
//					}
//				}
//				
//				// wrap periodically if needed
//				ii=0;
//				for (long j = 0; j < grafic_header.n2; j++) {
//					for (long i = 0; i < grafic_header.n1; i++) {
//						if (vx[ii]<0) vx[ii]=L+vx[ii]; if (vx[ii]>=L) vx[ii]=vx[ii]-L;
//						if (vy[ii]<0) vy[ii]=L+vy[ii]; if (vy[ii]>=L) vy[ii]=vy[ii]-L;
//						if (vz[ii]<0) vz[ii]=L+vz[ii]; if (vz[ii]>=L) vz[ii]=vz[ii]-L;
//						ii++;
//					}
//				}
//				
//				// convert Grafic units (comoving Mpc) to Gadget units (kpc/h)
//				cpeds_mul_value(1000.0,vx,nsx);
//				cpeds_mul_value(1000.0,vy,nsy);	
//				cpeds_mul_value(1000.0,vz,nsz);	
//				// TODO should we also multiply here by h ?
//			}
			
//			if (block==Velocities) {
//				// convert units here if necessary
//				// we convert here the grafic velocities to the gadget-2 velocities required in gadget IC file
//				// i.e. 1/sqrt(a_start)/h
//				// this is the case for the default internal velocity units of gadget (i.e. km/s)
//				double vconv=getGadgetVelocityFactor();
//				cpeds_mul_value(vconv,vx,nsx);
//				cpeds_mul_value(vconv,vy,nsy);	
//				cpeds_mul_value(vconv,vz,nsz);	
//			}
			
			
			// save to output file

			ii=0;
			for (long j = 0; j < grafic_header.n2; j++) {
				for (long i = 0; i < grafic_header.n1; i++) {
					fprintf(f,"%lE %lE %lE %lE\n",vx[ii],vy[ii],vz[ii],sqrt(vx[ii]*vx[ii]+vy[ii]*vy[ii]+vz[ii]*vz[ii]));
					ii++;
				}
			}
			
			
			
			delete [] vx;			delete [] vy;			delete [] vz;
			
		} // end of loop over slices
		fclose(fx);
		fclose(fy);
		fclose(fz);
		fclose(f);
		
	} // end of loop over particle types
		
}

/***************************************************************************************/
void Nbody_io::dumpDensityWeightedTemperature(string fname, int files, int fileType, double Wb0, int saveOpts, string outputFile) {
	string fname2;
	Nbody_io::gadget_snapshot_file_header_t head;
	Nbody_io::gadget_snapshot_type2_block_t infoBlk;
	long Npart=NumPart;
	float gamma= 5.0/3;
	double Xp=CPEDS_hydrogen_mass_fraction;

	
	for (long file = 0; file < files; file++) {
		//
		//	READ IN THE DATA FROM SNAPSHOT FILE
		// 
		if (files>1) {	fname2=fname+"."+msgs->toStr(file);	} else { fname2=fname; }
		head=readGadgetHeader(fname2,infoBlk);
		msgs->say("PROCESSING FILE: "+fname2,High);
		
		float* U=   readGadgetBlock(fname2,"U   ",0,&Npart,infoBlk);
		float* rho=NULL;
		if (saveOpts==0) {
			rho= readGadgetBlock(fname2,"RHO ",0,&Npart,infoBlk);
		}

	

		//
		// convert internal energy into temperature of the gas in SI units
		// 
		double *T = new double[Npart];
		float u;
	
		for (long i = 0; i < Npart; i++) {
			T[i]= cpeds_convert_gadget_internalEnergy_to_temp(U[i],gamma,Xp);
		}

		double rhoC0=cpeds_rhoC0(head.HubbleParam); // critical density today. Units: [kg/m^3]
		double rhoAv=rhoC0*Wb0*pow(head.redshift+1.0,3); 

		// convert to [kg/m^3 h^2]
		rhoAv/=(head.HubbleParam*head.HubbleParam);
		// convert to cgs
		rhoAv*=double(1.0e-3);

		// convert to gadget units
		rhoAv/=GADGET2_UnitDensity_in_cgs; /* 13^10 Msol(gram) / kpc (cm)^3 h^2
											 * Comment: conversion to the default gadget density units. 
											 * This is needed to reduce with rho which is also in gadget units.
											 * Important thing is that if the units will change in gadget - the change should also be reflected here.
											 * 
											 * author: blew
											 * date: Feb 28, 2013 1:27:43 PM
											 *
											 */

		if (saveOpts==0) {
			for (long i = 0; i < Npart; i++) {
				T[i]*= rho[i] / rhoAv; // the unit here is Kelvin
			}
			delete [] rho;
		}
		
		delete [] U;
		cpeds_save_matrix(T,Npart,1,outputFile,true,true);
	}
}
/***************************************************************************************/
cpedsStatusCodes Nbody_io::loadTipsy(string fname) {
	Nbody_io::tipsy_dark_particle_t DMpart;
	Nbody_io::tipsy_gas_particle_t GASpart;
	Nbody_io::tipsy_star_particle_t STARpart;

	FILE* tf;
//	Nbody_io::tipsyHeader_t tipsyHeader;

	if (cpeds_fileExists(fname)) { 
		// read tipsy header
		tf=fopen(fname.c_str(),"rb");
		fread(&tipsyHeader,sizeof(Nbody_io::tipsyHeader_t),1,tf);

//		printf("h.nbodies: %li\n", tipsyHeader.nbodies);
//		printf("h.ndark: %li\n", tipsyHeader.ndark);	
//		printf("h.ngas: %li\n", tipsyHeader.nsph);
//		printf("h.stars: %li\n", tipsyHeader.nstar);
//		printf("h.ndim: %li\n", tipsyHeader.ndim);
//		printf("h.time: %lE\n", tipsyHeader.time);

		long Npart;
		Npart=tipsyHeader.nsph;
		if (Npart>0) {
			tipsyPgas=new tipsy_gas_particle_t[Npart];
			for (long j = 0; j < Npart; j++) {
				fread(&GASpart,sizeof(tipsy_gas_particle_t),1,tf);	
				tipsyPgas[j]=GASpart;
//				printf("z: %lE\n",GASpart.pos[2]);
			}
		}
//		exit(0);
		Npart=tipsyHeader.ndark;
		if (Npart>0) {
			tipsyPdark=new tipsy_dark_particle_t[Npart];
			for (long j = 0; j < Npart; j++) {
				fread(&DMpart,sizeof(tipsy_dark_particle_t),1,tf);	
				tipsyPdark[j]=DMpart;
			}
		}
		

		fclose(tf);
	}
	else {
		return cpedsError;
	}

	
}
