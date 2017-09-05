#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "chealpix.h"


/* // pixel numeration is from 0 to pix_num-1 */
/* // region numeration is from 0 to 11 */
/* long int cpeds_get_healpix_region(long int nside, long int pix) { */
/*   long int region; */

/*   region = (long int)floor((double)pix/(double)(nside*nside)); */
/*   //printf("*** region = %li pix = %li, **=%lf\n ",region,pix, floor((double)pix/(double)(nside*nside)));   */
/*   return region; */
/* } */

/* // returns the number of the pixel num in ordering in the region where it belongs */
/* long int cpeds_get_healpix_pixnum_in_region(long int nside, long int pix) { */
/*   long int region; */

/*   region = cpeds_get_healpix_region(nside,pix); */
/*   return pix-(region)*nside*nside; */
/* } */

/* // returns the number of the column in ordering in the region where the piexel num  belongs */
/* // the column ordering is from 0 to nside-1 */
/* long int cpeds_get_healpix_colnum_in_region(long int nside, long int pix) { */
/*   long int pix_num_in_reg; */
/*   long int pix_colnum; */

/*   pix_num_in_reg = cpeds_get_healpix_pixnum_in_region(nside,pix); */
/*   pix_colnum = (long int)floor((double)pix_num_in_reg/(double)nside); */

/*   return pix_colnum; */
/* } */

/* // returns the number of the row in ordering in the region where the piexel num belongs */
/* // the column ordering is from 0 to nside-1 */
/* long int cpeds_get_healpix_rownum_in_region(long int nside, long int pix) { */
/*   long int pix_num_in_reg; */
/*   long int pix_rownum, pix_colnum; */

/*   pix_num_in_reg = cpeds_get_healpix_pixnum_in_region(nside,pix); */
/*   pix_colnum = cpeds_get_healpix_colnum_in_region(nside,pix); */

/*   pix_rownum = pix_num_in_reg - pix_colnum*nside; */
/*   return pix_rownum; */
/* } */

/* int cpeds_pix2xy(long int nside, long int pix, double *x, double *y, int ordering) { */
/*   double half_pix_size = 1/(double)nside; */
/*   long int pix_rownum, pix_colnum, region; */
/*   int byte_num,i,sum,z,X,Y,nsteps; */
/*   double dbnum,NS; */

/*   region = (long int)floor((double)pix/(double)(nside*nside)); */

/*   if (ordering == 3) { // this is my own ordering approach allowing for nside to be any integer value */
/*     pix_colnum = cpeds_get_healpix_colnum_in_region(nside,pix); */
/*     pix_rownum = cpeds_get_healpix_rownum_in_region(nside,pix); */
/*   } */

/*   if (ordering == 1) { // this is standard nested HEALPIX ordering - the stupid one :) */
/*     dbnum = 2*cpeds_log2((double)nside); */
/*     byte_num = (int)dbnum; */

/*     /\* reducing the pixel number to the number in a partucular region *\/ */

/*     pix = pix- nside*nside*region; */

/*     /\*     pix = 10; *\/ */

/*     /\* setting up the binary representation for finding the pixel x,y position *\/ */
/*     int *tab_num = malloc(byte_num*sizeof(int)); // new int[byte_num]; */
/*     sum = 0;  */
/*     for (i=byte_num-1;i>=0;i--) { */
/*       z = (int)pow((double)2,(double)i); */
/*       if (sum+z > pix) { tab_num[byte_num-1-i] = 0; } else { tab_num[byte_num-1-i] = 1; sum = sum+z; } */
/* /\*       printf("i = %i z = %i, sum = %i, tab=%i\n",i,z, sum, tab_num[i]);  *\/ */
/*     } */
/* /\*     for (i=0;i<byte_num;i++) { *\/ */
/* /\*       printf("\n i = %i tab=%i\n",i,tab_num[i]);  *\/ */
/* /\*     } *\/ */

/* /\*     printf("\ndupa %i, %li %lf \n",byte_num,nside, dbnum); *\/ */
/* /\*     exit(0); *\/ */


/*     X = Y = 0; */
/*     NS = 1; */
/*     pix_colnum = pix_rownum = 0; */
/*     nsteps = byte_num/2; */
/*     for (i = 0; i< nsteps; i++) { */
/*       NS = pow(2,(double)(nsteps-i-1)); */
/*       if (tab_num[2*i] == 1) { pix_colnum = pix_colnum + (long int)NS; } */
/*       if (tab_num[2*i+1] == 1) { pix_rownum = pix_rownum + (long int)NS; } */
/*     } */
/*     free(tab_num); */
/*   } */
/*   printf("*** colnum = %li rownum = %li ",pix_colnum,pix_rownum);   */

/* /\* this is the simplified form of the comment below *\/ */
/*   *x = half_pix_size*( pix_rownum - pix_colnum ); */
/*   *y = half_pix_size*( pix_rownum + pix_colnum + 1 ); */

/* /\*   *x = 0.5*( (2*half_pix_size*pix_rownum + half_pix_size) - (2*half_pix_size*pix_colnum + half_pix_size) ); *\/ */
/* /\*   *y = 0.5*( (2*half_pix_size*pix_rownum + half_pix_size) + (2*half_pix_size*pix_colnum + half_pix_size) ); *\/ */

/* } */
/* /\********************************************************************************** */
/* /\* ordering -- 1 - when nested ordering *\/ */
/* /\* ordering -- 0 - when ring ordering *\/ */
/* int cpeds_pix2ang_healpix(long int nside, long int pix, double *th, double *phi, int ordering) { */
/*   double x,y,X,Y, quarter; */
/*   long int region; */

/*   cpeds_pix2xy(nside,pix,&x,&y,ordering); */
/*   region = (long int)floor((double)pix/(double)(nside*nside)); */
/*   printf("*** reg. = %li ",region); */
/*   printf("x = %lf, y = %lf, pix = %li ", x , y, pix); */

/*   if ( (region == 0) || (region == 1) || (region == 2) || (region == 3) ) { // we're in the north polar region */
/*     quarter = region; */
/*     if (y >= 1) { */
/*       *phi = quarter*PI/2 + PI*(2-y+x)/(4*(2-y)); */
/*       *th  = asin(1-(2-y)*(2-y)/3);  */
/* /\*       printf("|||th=%lf y=%lf arg=%lf",*th,y,1-(2-y)*(2-y)/3 ); *\/ */
/*     } */
/*     if (y < 1) { // we're still in the equatorial region */
/*       *phi = quarter*PI/2 + PI*(1+x)/4; */
/*       *th = asin(2*y/3); */
/*     } */
/*   } */

/*   if ((region == 4) || (region == 5) || (region == 6) || (region == 7)) { // we're in the equatorial region */
/*     quarter = region-4; */
/*     *phi = quarter*PI/2 + PI*x/4; if (*phi < 0) { *phi = 2*PI+*phi; } */
/*     *th = asin(2*(y-1)/3); */
/*   } */

/*   if ((region == 8) || (region == 9) || (region == 10) || (region == 11)) { // we're in the south polar region */
/*     quarter = region-8; */
/*     y = 2-y; */
/*     if (y >= 1) { */
/*       *phi = quarter*PI/2 + PI*(2-y+x)/(4*(2-y)); */
/*       *th = -asin((1-(2-y)*(2-y)/3)); */
/*     } */
/*     if (y < 1) { // we're still in the equatorial region */
/*       *phi = quarter*PI/2 + PI*(1+x)/4; */
/*       *th = -asin(2*y/3); */
/*     } */
/*   } */
/*   *th = PI/2-*th; */
/*   printf("th_cpeds = %.15lf, phi_cpeds = %lf\n ", 180/PI*(*th) , 180/PI*(*phi)); */
/* } */

/* /\**********************************************************************************\/ */
/* /\* ordering -- 1 - when orderinged ordering *\/ */
/* /\* ordering -- 0 - when ring ordering *\/ */
/* /\* the invers formulas cannot be used beyond the area of the initial region since they are no longer valid: exactly the formulas for the polar regions y>=1 *\/ */
/* /\* usage is ok only in some special cases in polar regions and it's ok in all equatorial regions including y<1 parts of polar regions *\/ */

/* int cpeds_ang2pix_healpix(long int nside, long int *pix, double th, double phi, int ordering) { */
/*   double TH_TROPIC = PI/2-48.189685109680696*PI/180; // boundary of the polar-equator region in corrdinates where pole has theta=0 */
/*   double half_pix_size = 1/(double)nside; */
/*   double x,y,x_eq,y_eq, x_np, y_np, x_sp, y_sp; */
/*   double phi_np,phi_sp,phi_eq,region,region_np,region_sp,region_eq, pix_rownumd, pix_colnumd,sum,dbnum,divcol,divrow; */
/*   long int pix_rownum, pix_colnum,pix_num,i,nsteps,byte_num; */
/*   th = PI/2-th; // transformation from Healpix coordinates to the coordinates in which Boud's transformations were calculated */

/*   if (th > TH_TROPIC) { // we're in the north polar region */
/*     region = floor(phi*2/PI); phi_np = phi - region* PI/2; // this is to keep x inside the inital region to avoid mistakes. */
/*     x = (-1 + 4*phi_np/PI) * sqrt(3-3*sin(th));// + half_pix_size; */
/* /\*     x = x + region*2; *\/ */
/*     y = 2 - sqrt(3-3*sin(th)); //- half_pix_size; */
/*     printf("(1) x = %lf, y=%lf \n ", x,y); */
/*     printf("-------------> x = %lf, y = %lf region = %lf\n ",x,y,region); */
/*   } */

/*   if ((th < TH_TROPIC) && (th >= -TH_TROPIC)) { // we're basically in the equatorial region but still some of the poles region are possible - so we calculate for the two cases */
/*     region_np = floor(phi*2/PI); phi_np = phi - region_np* PI/2; // this is to keep x inside the inital region to avoid mistakes. */
/*     x_np = 4*phi_np/PI - 1; y_np = 3*sin(th)/2; // this is for the polar regions in equatorial part */

/*     region_eq = floor((phi+PI/4)*2/PI); phi_eq = phi - region_eq* PI/2; if (region_eq == 4) { region_eq = 0; printf("################# %lf ###############",region_eq); }  // this is to keep x inside the inital region to avoid mistakes. */
/*     x_eq = 4*phi_eq/PI; y_eq = 0.5 * (3*sin(th) + 2); // this is a guess for equatorial region */

/*     th = -th; */
/*     region_sp = floor(phi*2/PI); phi_sp = phi - region_sp* PI/2; // this is to keep x inside the inital region to avoid mistakes. */
/*     x_sp = 4*phi_sp/PI - 1; y_sp = 2-3*sin(th)/2; // this is for the polar regions in equatorial part */

/*     printf("(2) x_eq = %lf, y_eq = %lf, x_np = %lf y_np = %lf x_sp = %lf, y_sp = %lf region_np = %lf region_eq = %lf region_sp = %lf\n ",x_eq,y_eq, x_np, y_np,x_sp, y_sp,region_np,region_eq,region_sp); */
/*     th = -th; */

/*     if ( (y_eq <= -fabs(x_eq)+2) && (y_eq >= fabs(x_eq)) ) { region = region_eq+4; x = x_eq; y = y_eq; printf("A");} else { */
/*       if (y_eq > -fabs(x_eq)+2) { region = region_np; x = x_np; y = y_np; printf("B");} else { region = region_sp+8; x = x_sp; y = y_sp; printf("C");} } */
/*     printf("-------------> x = %lf, y = %lf region = %lf\n ",x,y,region); */
/*   } */
  
/*   if (th < -TH_TROPIC) { // we're in the south polar region */
/*     th = -th; */
/*     region = floor(phi*2/PI)+8; phi_sp = phi - (region-8)* PI/2; // this is to keep x inside the inital region to avoid mistakes. */
/*     x = (-1 + 4*phi_sp/PI) * sqrt(3-3*sin(th));// + half_pix_size; */
/* /\*     x = x + region*2; *\/ */
/*     y = sqrt(3-3*sin(th));// + half_pix_size; */
/*     printf("(3) x = %lf, y=%lf \n ", x,y); */
/*     th = -th; */
/*     printf("-------------> x = %lf, y = %lf region = %lf\n ",x,y,region); */
/*   } */

/*   pix_rownumd = rint(0.5 * (y+x-half_pix_size)/half_pix_size); */
/*   pix_colnumd = rint(0.5 * (y-x-half_pix_size)/half_pix_size); */
/*   pix_rownum = (int)pix_rownumd; */
/*   pix_colnum = (int)pix_colnumd; */

/*   printf("^^^^ th = %lf, phi = %lf half_pix_size = %lf, colnum = %li rownum = %li\n", 90-180/PI * th, 180/PI * phi, half_pix_size,pix_colnum, pix_rownum); */

  

/*   dbnum = 2*cpeds_log2((double)nside); */
/*   byte_num = (int)dbnum; */
  
/*   /\* reducing the pixel number to the number in a partucular region *\/ */

/*   pix_num = nside*nside*(long int)region; */

/*   /\* setting up the binary representation for finding the pixel x,y position *\/ */
/*   int *tab_num = malloc(byte_num*sizeof(int)); // new int[byte_num]; */

/*   nsteps = byte_num/2; */
/*   for (i = 0; i< nsteps; i++) { */
/*     divrow = floor(pix_rownumd / 2); tab_num[2*(nsteps-i)-1] = (int)(pix_rownumd - 2*divrow); pix_rownumd = divrow; */
/*     divcol = floor(pix_colnumd / 2); tab_num[2*(nsteps-i)-2] = (int)(pix_colnumd - 2*divcol); pix_colnumd = divcol; */
/*     //x    printf("\n   % */
/*   } */

/*   sum = 0;  */
/*   for (i=byte_num-1;i>=0;i--) {    sum = sum+pow((double)2,(double)i)*tab_num[byte_num-1-i];  } */
/*   free(tab_num); */
/*   *pix = pix_num + (long int)sum; */
/*   printf("==========>>>>>>>>>>> pix = %li\n\n",*pix); */

/*   return 0; */
/* } */
/**********************************************************************************/



int main(int ARGC, char* ARGV[]) {
  long int num,i;
  double s,m;
  double *tab;
  char tmpch[100];
  FILE* f;

  if (ARGC == 1) { printf("usage: gaussian_numbers_test arg1: number of numbers to generate arg2: sigma of distribution arg3: mean of distr arg4: outfile.\n"); exit(0); }

    strcpy(tmpch,ARGV[1]);
    num = strtol(&tmpch[0],NULL,10);
    strcpy(tmpch,ARGV[2]);
    s = strtod(&tmpch[0],NULL);
    strcpy(tmpch,ARGV[3]);
    m = strtod(&tmpch[0],NULL);
    strcpy(tmpch,ARGV[4]);

    tab = cpeds_random_gauss_numbers(m,s,num,2);

    f = fopen(&(tmpch[0]),"w");
    for (i=0;i<num;i++) { fprintf(f,"%lf\n",tab[i]); }
    fclose(f);

  return 0;

}
