/*!
 * \file test-fits.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: blew
 */

#include "stdio.h"
#include "stdlib.h"
#include "cpeds-math.h"
//#include "cpeds-list.h"
#include "Mscs-function.h"
#include "string.h"


int main(int argc, char** argv) {
	
	if (argc==1) {
		printf("USAGE: test-mscsfn testNo\n");
		printf("TESTS\n");
		printf("0 - air mass elevation relation test\n");
		exit(0);
	}
	long testNo=strtol(argv[1],NULL,10);

	
	mscsFunction f,g;
	double zd=0;
	double R0=6374;
	double npatm;
	double h=10;
	double fact;

	
	
	switch (testNo) {
		case 1:
			while (zd<91) {
				npatm=cpeds_nonplanar_atmospheric_layer_path(zd,R0,h,R0);
				printf("%lf %lE\n",zd,npatm);
				f.newPoint(zd,npatm/h);
				g.newPoint(zd,npatm/(h/cos(zd*PI180)));
				zd+=5;
			}
			f.save("nonPlanarAtm.zd");
			g.save("nonPlanarAtm_wrt_geometrical.zd");	
			
			break;
//		case 2:
//			while (zd<91) {
//				fact=(1.0-exp(-))
//				printf("%lf %lE\n",zd,npatm);
//				f.newPoint(zd,npatm/h);
//				g.newPoint(zd,npatm/(h/cos(zd*PI180)));
//				zd+=5;
//			}
//			f.save("nonPlanarAtm.zd");
//			
//			
//			break;
		default:
			break;
	}
	return 0;
}
