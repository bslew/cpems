#include <cpgplot.h>
#include "Mscs-map_plot.h"
#include "Mscs-colormap.h"
#include "cpeds-smooth.h"
#include "cpeds-point_set.h"
#include "Mscs-map_proj.h"
//#include <QtGui/QPainter>

//************************************************************************
// plots the flat map in a requested projection
// output : 0 - ineractive do not open and close the pgplot window
// projection: 1 - moll - Mollweide
// output: 1 - XWINDOW
// output: 2 - postscript
// db,dl = offset in galactic lattitude and galactic longitude
// window_size - pgplot window size in inches
long mscsMapPlot::plot(long sizeX, long sizeY, string projection, double l, double b, double Tmin, double Tmax, long color_num, string color_scheme, string title, string bgcolor, double window_size, int output) {
//plot_map(int projection, double dl, double db, int output, long color_num, string color_scheme, strarg title,double window_size) { 
  float *tab;
  float xmin,xmax,ymin,ymax,minTloc,maxTloc;
  float tr[6];
  string logch;
  float x,y; // world coordinates in the PG window
  char ch='\0';
  long i,j;
  long ret;
  bool is_log=false;
  double fmX,fmY; // matrix coordinates (should be integer)
  cpedsDirection reference_dir;
  double projX,projY,fml,fmb,fmT; // from map values 
  double factor, factorx, factory, tmpd, rotation_unit=90;
  cpedsDirection n;
  double zoomLevel=2.0;

  // filenamestr projtype,eastwest;
  // filenamestr command_str;
  projUV p;
  projPJ pj;
  FILE* logf;
  string logs,taken_point_comment,incmd;
  //  double * backup = new double[pix_num];   for (i=0;i<pix_num;i++) { backup[i] = map->T[i]; }
  // const char txtside='B';
  // const char txttext='dupa';

/* EXPERIMENTALL STUFF */
/* 	initscr(); */
/* 	noecho(); */
/* 	nonl(); */
/* 	cbreak(); */
/* EXPERIMENTALL STUFF */

  if (sizeX==0) sizeX=4*nside();
  if (sizeY==0) sizeY=2*nside();
  msgs->say("plotting map of size: "+msgs->toStr(sizeX)+" x "+msgs->toStr(sizeY)+" pix. num: "+msgs->toStr(sizeX*sizeY),High);
  calculate_map_stats(1);
  msgs->say("setting coordinates",Medium);
  set_map_coord();
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // cpedsPointSet3D ps=cpedsProject(get_n(), "moll").projectOnPlane(cpedsDirection(l,b));																																				      //
  // for (long i=0;i<pixNum();i++) { ps[i].setZ(get_T(i)); }																																								      //
  // matrix<long> mask;																																													      //
  // cpedsPointSet3D psint;																																												      //
  // matrix<double> m = ps.exportAsInterpolatedField(sizeX,sizeY,mask,psint);																																						      //
  // cpeds_matrix_save(mask,"ilc5.mask.mat2");																																										      //
  // psint.save("ilc5.coordint");																																											      //
  // // exit(0);																																													      //
  // // inverse project directions to get right temperatures																																								      //
  // string cmdstr="+proj="+projection+" +lat_0="+msgs->toStr(b)+" +lon_0="+msgs->toStr(l)+"w +R="+msgs->toStr(double(sqrt(2.0)/4.0));																															      //
  // pj = pj_init_plus(cmdstr.c_str());																																											      //
  // 																																															      //
  // tab = new float[sizeX*sizeY];																																											      //
  // long k=0;																																														      //
  // for (int j=0;j<sizeY;j++) {																																											      //
  //   for (int i=0;i<sizeX;i++) {																																											      //
  //     																																														      //
  //     if (mask(j,i)>0) { 																																												      //
  // 	tab[k]=m(j,i);																																													      //
  //     }																																														      //
  //     else {																																														      //
  // 	p.u = psint[j*sizeX+i].x();    p.v = psint[j*sizeX+i].y(); // conversion from phi in (0,2pi) --> east and west longitude from (-pi,pi) -- here we inverse the direction of rising l to be consistent with the standards in publications that the l rises form the center of the figure left-wise, jumps to the right edge and rises from pi->2pi towards the picture center	      //
  // 	p.u=-1; p.v=0.5;																																												      //
  // 	p = pj_inv(p,pj); 																																												      //
  // 	if (pj_errno!=0) { tab[k]=m(j,i) = HUGE_VAL;   }																																								      //
  // 	else {																																														      //
  // 	  printf("val ok\n");																																												      //
  // 	  fmb = p.v;      fml = -p.u; if (fml < 0) fml=twoPI+fml;  																																							      //
  // 	  tab[k]=m(j,i) = get_T(fml,fmb);																																										      //
  // 	  mask(j,i)=k;																																													      //
  // 	}																																														      //
  //     }																																														      //
  //     k++;																																														      //
  //   }																																														      //
  // }																																															      //
  // pj_free(pj); 																																													      //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // cpeds_matrix_save(m,"ilc5.mat");
  // cpeds_matrix_save(mask,"ilc5.mask.mat");
  // exit(0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // matrix<double> m = ps.exportAsField(sizeX,sizeY,mask);																																								      //
  // long n;																																														      //
  // double* dtab= cpeds_matrix2array(m, n,true);																																									      //
  // tab=float_double_array(dtab,n);																																											      //
  // delete dtab;																																													      //
  // 																																															      //
  // pj = pj_init_plus(cmdstr.c_str());																																											      //
  // 																																															      //
  // long k=0;																																														      //
  // for (int j=0;j<sizeY;j++) {																																											      //
  //   for (int i=0;i<sizeX;i++) {																																											      //
  // 																																															      //
  //     if (mask(j,i)==-1) {																																												      //
  // 	p.u = psint[j*sizeX+i].x();    p.v = psint[j*sizeX+i].y(); // conversion from phi in (0,2pi) --> east and west longitude from (-pi,pi) -- here we inverse the direction of rising l to be consistent with the standards in publications that the l rises form the center of the figure left-wise, jumps to the right edge and rises from pi->2pi towards the picture center	      //
  // 	p = pj_inv(p,pj); 																																												      //
  // 	if (p.v==HUGE_VAL or p.u==HUGE_VAL) { mask(j,i)=-1; tab[k]=m(j,i) = HUGE_VAL; } 																																				      //
  // 	else {																																														      //
  // 	  fmb = p.v;      fml = -p.u; if (fml < 0) fml=twoPI+fml;  																																							      //
  // 	  tab[k]=m(j,i) = get_T(fml,fmb);																																										      //
  // 	  mask(j,i)=k;																																													      //
  // 	}																																														      //
  // 	k++;																																														      //
  //     }																																														      //
  //   }																																														      //
  // }																																															      //
  // pj_free(pj); 																																													      //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  matrix<long> mask;
  matrix<double> m;
  msgs->say("projecting map on plane",Medium);
  
  mscsMapProj flat(*this,sizeX,sizeY,projection,twoPI-l*PI180,b*PI180);
  m=flat.getFlatMap(mask);

  long N;
  double* dtab= cpeds_matrix2array(m, N,true);																																									      //
  tab=float_double_array(dtab,N);	
  delete [] dtab;

  //
  // map colors
  //
  mscsColormap *mcol;
  if (Tmin==Tmax and Tmax==0) {
    msgs->say("setting Tmin and Tmax to extremal values in the map: "+msgs->toStr(Tmin)+" "+msgs->toStr(Tmax),Medium);
    Tmin=getMinT();
    Tmax=getMaxT();
  }


  //
  // plotting part
  //
  if (projection!="ortho" && projection!="stere")  {
    xmin = 0; xmax = 360; ymin = -90; ymax = 90;
  }
  else {
    xmin = 0; xmax = 180; ymin = -90; ymax = 90;
  }

  tr[0]=xmin;
  tr[1]=(xmax-xmin)/(float)sizeX;
  tr[2]=0;
  tr[3]=ymin;
  tr[4]=0;
  tr[5]=(ymax-ymin)/(float)sizeY;

/*   tr[0]=xmin; */
/*   tr[1]=(xmax-xmin)/(float)xpix_num_flat; */
/*   tr[2]=0; */
/*   tr[3]=ymax; */
/*   tr[4]=0; */
/*   tr[5]=-(ymax-ymin)/(float)ypix_num_flat; */

  if (output == 1) { cpgbeg(0,"/XWINDOW",1,1); cpgask(0);}
  // if (output == 2) { sprintf(file_and_device,"%s/cps",flatmap_file_name);  cpgbeg(0,file_and_device,1,1);}
  // if (output == 3) { sprintf(file_and_device,"%s/GIF",flatmap_file_name);  cpgbeg(0,file_and_device,1,1);}
  // if (output == 4) { sprintf(file_and_device,"%s/jpg",flatmap_file_name);  cpgbeg(0,file_and_device,1,1);}
    //if (output == 10) { }
  do {
    if (output == 1 ) {
      if (projection!="ortho" && projection!="stere") cpgpap(window_size,0.7);
      else cpgpap(window_size,1.0);
      cpgenv(xmin, xmax, ymin, ymax, 1, 2);  //coordinates
      cpglab("180-l [deg]", "b [deg]", title.c_str()); //axes
      cpgsci(10);
/*       cpgsci(0); cpgmtxt(&txtside,7,0.0,0,titlestr3);  cpgsci(10); sprintf(titlestr3,"(l,b)[deg]: %lf %lf, T(l,b)=%lE, ang(ref)[deg]=%lf",fml*180/PI,fmb*180/PI,fmT,180/PI*cpeds_ang_n1n2(PIsnd-reference_dir.b,reference_dir.l,PIsnd-fmb,fml));  cpgmtxt(&txtside,7,0.0,0,titlestr3); */
/*       cpgsci(0); cpgmtxt(&txtside,8,0.0,0,titlestr2); cpgsci(10); sprintf(titlestr2,"ref. dir: (l,b)[deg]: %lf %lf",reference_dir.l*180/PI,reference_dir.b*180/PI);  cpgmtxt(&txtside,8,0.0,0,titlestr2); */

    }

    //
    // generate map colors
    //
    generateMapColors(color_scheme,color_num,Tmin,Tmax);
    
    

    cpgsch(0.5);



    //cpgeras();
    // if (flatmap_change_zrange) {minTloc = flatmap_minT; maxTloc = flatmap_maxT;} 
    // else  {
    //   minTloc = minT; maxTloc = maxT;
/*       cpeds_find_minmax_value(flatmap,pix_num_flat,&minTloc,&maxTloc,&i,&j); */
/*       if (minTloc==-1e10) minTloc = minT; */
  // }

    minTloc=Tmin;
    maxTloc=Tmax;
    cpgimag(tab, sizeX, sizeY, 1, sizeX, 1, sizeY, minTloc, maxTloc, tr);

    // QPainter painter;
    // painter.setPen(Qt::blue);
    // painter.setFont(QFont("Arial", 30));
    // painter.drawText(rect(), Qt::AlignCenter, "Qt");

    // { 
    //   cpgscr(0,0,0,0);
    //   cpgscr(1,1,1,1);
    //   cpgscr(2,0.5,0.5,0.5);
    //   cpgslw(100);
    //   float x,y;
    //   for (long i=0;i<m.RowNo();i++) {
    // 	for (long j=0;j<m.ColNo();j++) {
    // 	  x=tr[0] + tr[1]*float(i) + tr[2]*float(j);
    // 	  y=tr[3] + tr[4]*float(i) + tr[5]*float(j);
    // 	  if (mask(i,j)==-1) { cpgsci(100); cpgpt(1,&x,&y,-1);  }
    // 	  if (mask(i,j)==0) { cpgsci(200); cpgpt(1,&x,&y,-1);  }
    // 	  // cpgsci(0);
    // 	  // cpgpt(1,&y,&x,-1);
    // 	  printf("-------- %f %f\n",x,y);
    // 	}
    //   }
    // }

/*     cpgwedg("BI",1,4,(float)minT,(float)maxT,"K"); */
/*     cpgwedg("BI",1,4,(float)minT,(float)maxT,"(K-<K>)/\\gs"); */
/*     cpgwedg("BI",1,4,(float)minT,(float)maxT,"(S-<S>)/\\gs"); */
    cpgwedg("BI",2,6,(float)Tmin,(float)Tmax,"data ID"); // 2 - displacement, 6 - width
    

    if (output == 1) { 
      ch='A';
      while ( ch=='A') {
    	cpgcurs(&x,&y,&ch);
	
    	// fmX=(double)((x-xmin)*(float)sizeX/(xmax-xmin));
    	// fmY=(double)((y-ymin)*(float)sizeY/(ymax-ymin));
        //       factor = 1.00;//(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed
        // factorx = factor*(maxx_flat-minx_flat)/(double)sizeX; /// how many times the map was squeezed
        // factory = factor*(maxy_flat-miny_flat)/(double)sizeY; /// how many times the map was squeezed
    	// projX=fmX*factorx+minx_flat;
    	// projY=fmY*factory+miny_flat;

    	// strcpy(eastwest,"e");     strcpy(projtype,"moll");    sprintf(command_str,"+proj=%s +lat_0=%lf +lon_0=%lf%s +R=1",projtype,0.0,0.0,eastwest);

    	// pj = pj_init_plus(command_str);

    	// p.u = projX;    p.v = projY; // conversion from phi in (0,2pi) --> east and west longitude from (-pi,pi) -- here we inverse the direction of rising l to be consistent with the standards in publications that the l rises form the center of the figure left-wise, jumps to the right edge and rises from pi->2pi towards the picture center
    	// p = pj_inv(p,pj); fmb = p.v;      fml = -p.u; if (fml < 0) fml=twoPI+fml;  
    	// fmT = get_T(fml,fmb);  pj_free(pj); 
    	// cpgsci(0); cpgmtxt(&txtside,7,0.0,0,titlestr3);  cpgsci(10); sprintf(titlestr3,"(l,b)[deg]: %lf %lf, T(l,b)=%lE, ang(ref)[deg]=%lf",fml*180/PI,fmb*180/PI,fmT,180/PI*cpeds_ang_n1n2(PIsnd-reference_dir.b,reference_dir.l,PIsnd-fmb,fml));  cpgmtxt(&txtside,7,0.0,0,titlestr3);
    	// cpgsci(0); cpgmtxt(&txtside,8,0.0,0,titlestr2); cpgsci(10); sprintf(titlestr2,"ref. dir: (l,b)[deg]: %lf %lf",reference_dir.l*180/PI,reference_dir.b*180/PI);  cpgmtxt(&txtside,8,0.0,0,titlestr2);
      }
    }
    printf("rotation unit: %lf\n",rotation_unit);
    if (ch == 's') { Tmin = Tmin+0.1*Tmin; }     if (ch == 'x') { Tmin = Tmin-0.1*Tmin; }
    if (ch == 'S') { Tmin = Tmin+0.1*get_varianceT(); }     if (ch == 'X') { Tmin = Tmin-0.1*get_varianceT(); }
    if (ch == 'e') { Tmax = Tmax+0.1*Tmax; }     if (ch == 'd') { Tmax = Tmax-0.1*Tmax; }
    if (ch == 'E') { Tmax = Tmax+0.1*get_varianceT(); }     if (ch == 'D') { Tmax = Tmax-0.1*get_varianceT(); }
    if (ch == 'r') { flat.mapSph()=*this; rotation_unit=90; flat.reset();  Tmin=getMinT(); Tmax=getMaxT(); tab=replot(flat,Tmin,Tmax,m,mask); } //for (j=0;j<pix_num;j++) map->T[j] = backup[j];  calculate_map_stats(1); makekill_space_manager("kill","C",2); makekill_space_manager("kill","T",2); calculate_flat_coord(projection,dl-180,db);  create_resized_flat_map(0); }
    if (ch == 'R') { flat.mapSph().rotate_map((PIsnd-reference_dir.b())*180/PI,(reference_dir.l()-PI)*180/PI,0,true,"T"); flat.mapSph().calculate_map_stats(); tab=replot(flat,Tmin,Tmax,m,mask); }
    if (ch == 'W') { printf("give file name: "); cin >> logs; flat.mapSph().savebinT(logs.c_str()); }
    if (ch == '`') { printf("give rotation unit: "); cin >> rotation_unit; }
    if (ch == 'l') { flat.mapSph().rotate_map(0,0,rotation_unit,true,"T"); flat.mapSph().calculate_map_stats(1); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'b') { flat.mapSph().rotate_map(rotation_unit,0,0,true,"T"); flat.mapSph().calculate_map_stats(1); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'L') { flat.setLon0(flat.lon0()+rotation_unit*PI180); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'B') { flat.setLat0(flat.lat0()+rotation_unit*PI180); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'm') { printf("  -- map angular momentum: %lE\n",calculate_map_momentum()); }
    if (ch == 'p') { printf("  -- map asymmetry coefficient: NOT IMPLEMENTED YET\n"); }
    if (ch == 'M') { /* printf("give nside start: "); cin >> i; */ reference_dir=flat.mapSph().maximize_map_momentum(1,nside(),&tmpd,&tmpd,true); flat.mapSph().calculate_map_stats(1);  tab=replot(flat,Tmin,Tmax,m,mask); }//  if (flat_coord_loaded == 0) calculate_flat_coord(projection,dl-180,db); create_resized_flat_map(0); }
    if (ch == 'N') { /* printf("give nside start: "); cin >> i; */ reference_dir=flat.mapSph().maximize_map_momentum(1,nside(),&tmpd,&tmpd,false); flat.mapSph().calculate_map_stats(1); tab=replot(flat,Tmin,Tmax,m,mask); } //   if (flat_coord_loaded == 0) calculate_flat_coord(projection,dl-180,db); create_resized_flat_map(0); }
    if (ch == 'z') { zoomLevel*=1.3; flat.zoomIn((long)getMatrixCoord(x,xmin,xmax,sizeX),(long)getMatrixCoord(y,ymin,ymax,sizeY),zoomLevel); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'Z') { zoomLevel/=1.3; flat.zoomIn((long)getMatrixCoord(x,xmin,xmax,sizeX),(long)getMatrixCoord(y,ymin,ymax,sizeY),1.0/zoomLevel); tab=replot(flat,Tmin,Tmax,m,mask); } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == 'C') { color_scheme=mscsColormap::cycleColorSchemes(color_scheme); generateMapColors(color_scheme,color_num,Tmin,Tmax);  } //makekill_space_manager("kill","C",2); calculate_flat_coord(projection,dl-180,db);   create_resized_flat_map(0); }
    if (ch == '1') { flat.setLonLat1(  flat.getCoordLB( (long)getMatrixCoord(x,xmin,xmax,sizeX),(long)getMatrixCoord(y,ymin,ymax,sizeY) )  ); } 
    if (ch == '2') { flat.setLonLat2(  flat.getCoordLB( (long)getMatrixCoord(x,xmin,xmax,sizeX),(long)getMatrixCoord(y,ymin,ymax,sizeY) )  ); } 

// /*     if (ch == 'Q') {  } */
// /*     if (ch == 'l') {  */
// /*       if (!is_log) for (i=0;i<pix_num;i++) { if (map->T[i] > 0) map->T[i] = log(map->T[i]); if (map->T[i] < 0) map->T[i] = log(-map->T[i]); } */
// /*       else for (i=0;i<pix_num;i++) { if (map->T[i] < 0) map->T[i] = pow(10,map->T[i]); is_log = !is_log; } */

// /*       is_log = !is_log; */
// /*     } */

    if (ch == '?') { // print help
      printf("|%s> * PLOT METHOD HELP MENU:\n",object_name.c_str());
      printf("|%s>   s,x,e,d - in/de/crease the upper/lower z range of the map by 10 %% of the min/max values:\n",object_name.c_str());
      printf("|%s>   S,X,E,D - in/de/crease the upper/lower z range of the map by 10 %% of the variance value:\n",object_name.c_str());
      printf("|%s>   r - reset the plot, and calculate the statistics and replot \n",object_name.c_str()); 
      printf("|%s>   R - rotate the map to the position stored in the reference direction\n",object_name.c_str());
      printf("|%s>   W - dump map to file (will ask for file name)\n",object_name.c_str()); 
      printf("|%s>   ` - change the rotation unit in deg (will ask for value)\n",object_name.c_str()); 
      printf("|%s>   l - rotate map in l by rotation unit\n",object_name.c_str()); 
      printf("|%s>   b - rotate map in b by rotation unit\n",object_name.c_str()); 
      printf("|%s>   L - increase projection angle in galactic longitude by rotation unit\n",object_name.c_str()); 
      printf("|%s>   B - increase projection angle in galactic latitude by rotation unit\n",object_name.c_str()); 
      printf("|%s>   m - calculate the map angular momentum\n",object_name.c_str()); 
      printf("|%s>   p - calculate the map asymmetry correlation coefficient. TO BE IMPLEMENTED INTO THE Mscs-map object... \n",object_name.c_str()); 
      printf("|%s>   M - find the maximal angular momentun axis of the map in real space\n",object_name.c_str()); 
      printf("|%s>   N - find the minimal angular momentun axis of the map in real space\n",object_name.c_str()); 
      printf("|%s>   A,t - store the reference direction from the map over the coursor position, print the direction and map value under the coursor\n",object_name.c_str()); 
      printf("|%s>   T - store the reference direction and value from the map over the coursor position to file (will ask for comment)\n",object_name.c_str()); 
      printf("|%s>   D - store the reference direction from the map over the coursor position to file (will not ask for comment)\n",object_name.c_str()); 
      printf("|%s>   c - read the command from the keyboard\n\tsave_non_zero_value_dirs\n",object_name.c_str()); 
      printf("|%s>   C - cycle color scheme\n",object_name.c_str()); 
      printf("|%s>   z - zoom in region\n",object_name.c_str()); 
      printf("|%s>   Z - zoom out region\n",object_name.c_str()); 
      printf("|%s>   Q - quit the program now\n",object_name.c_str()); 
      printf("|%s>   ? - print this help info\n",object_name.c_str()); 
    }

    if (ch == 'c') {
      printf("|plot> give command name:\n");
      cin>>incmd;

      if (incmd == "save_non_zero_value_dirs") {
	set_map_coord();
	logs = title; logs+=".from_map_NZ_value_dirs";   logf = fopen(logs.c_str(),"w");
	for (i=0;i<pixNum();i++) {
	  if (flat.mapSph().get_T(i)!=0) { 
	    n=flat.mapSph().get_C(i);
	    fprintf(logf,"%lE %lE\n",n.l()*180/PI,n.b()*180/PI); 
	  }
	}
	fclose(logf);
      } 
    }


    if (ch == 't' || ch == 'T' || ch == 'y' || ch == 'A') { // get the coordinates from map and print them out or save to file
      // transform from world coords to matrix coords.
      fmX=getMatrixCoord(x,xmin,xmax,sizeX); //(double)((x-xmin)*(float)sizeX/(xmax-xmin));
      fmY=getMatrixCoord(y,ymin,ymax,sizeY);//      fmY=(double)((y-ymin)*(float)sizeY/(ymax-ymin));
      printf("xw: %f, yw: %f||| fmX:%f fmY:%f xpix_num_flat=%i ypix_num_flat:%i\n",x,y,fmX,fmY,sizeX,sizeY);
      // convert to proj projected coordinates
      n=flat.getCoordLB(long(fmX),long(fmY));
      n.setVal(flat.get_T(long(fmX),long(fmY)));
      n.print_direction("pointer at");
      // // !!! copied from create_resized_flat_map !!! 
      // factor = 1.00;//(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed
      // factorx = factor*(flat.getMaxX()-flat.getMinX())/(double)sizeX; /// how many times the map was squeezed
      // factory = factor*(flat.getMaxY()-flat.getMinY())/(double)sizeY; /// how many times the map was squeezed
      // projX=fmX*factorx+minx_flat;
      // projY=fmY*factory+miny_flat;

// /*       sprintf(command_str,"+proj=%s +lat_0=%lf +lon_0=subroutines.html%lf%s",projtype,0,0,eastwest); */
//       strcpy(eastwest,"e");     strcpy(projtype,"moll");    sprintf(command_str,"+proj=%s +lat_0=%lf +lon_0=%lf%s +R=1",projtype,0.0,0.0,eastwest);

//       pj = pj_init_plus(command_str);

//       p.u = projX;    p.v = projY; // conversion from phi in (0,2pi) --> east and west longitude from (-pi,pi) -- here we inverse the direction of rising l to be consistent with the standards in publications that the l rises form the center of the figure left-wise, jumps to the right edge and rises from pi->2pi towards the picture center
//       p = pj_inv(p,pj); fmb = p.v;      fml = -p.u; if (fml < 0) fml=twoPI+fml;  
//       fmT = get_T(fml,fmb);  pj_free(pj); 
//       printf("fml: %lf, fmb: %f, T[fml,fmb] = %lE\n",fml*180/PI,fmb*180/PI,fmT);
//       printf("relative (l,b)=(%.3lf,%.3lf) angular distance [deg]: %lf\n",reference_dir.l*180/PI,reference_dir.b*180/PI,180/PI*cpeds_ang_n1n2(PIsnd-reference_dir.b,reference_dir.l,PIsnd-fmb,fml));

//       sprintf(titlestr,"(l,b)[deg]: %lf %lf",fml,fmb);  cpgmtxt(&txtside,7,0.0,0,titlestr);

      reference_dir=n;
//       if (ch == 'T') { // save append values to file
// 	cin >> taken_point_comment;
// 	logs = title; logs+=".from_map_taken_vals";   logf = fopen(logs.c_str(),"a"); fprintf(logf,"%lE %lE %lE %s\n",fml*180/PI,fmb*180/PI,fmT,taken_point_comment.c_str()); fclose(logf);      }
//       if (ch == 'y') { // save append values to file
// 	logs = title; logs+=".from_map_taken_dirs";   logf = fopen(logs.c_str(),"a"); fprintf(logf,"%lE %lE\n",fml*180/PI,fmb*180/PI); fclose(logf);      }

//       //printf("--saving to file %s\n",logs.c_str()); 
      long idx=-1;
      get_T(reference_dir,&idx);
      printf("The corresponding pixel index in the healpix map is: %li\n",idx);
      
    }
    

  
  
      
    


// /*     TO BE IMPLEMENTED */
// /*     if (ch == 'e') { calculate_map_stats(1);} // exports the current map to the matrix file */
// /*     if (ch == 'g') { calculate_map_stats(1);} // prints the position of the cursor and l, */
// /*     if (ch == '[') { maxT = maxT+0.1*maxT; }     if (ch == ']') { maxT = maxT-0.1*maxT; } // shifts coordinates in l direction */
// /*     if (ch == 'M') { calculate_map_stats(1);} // sets the circular mask size
// /*     if (ch == 'm') { calculate_map_stats(1);} // puts the mask in the current cursor position of size M
// /*     if (ch == ',') { calculate_map_stats(1);} // saves the mask to mask file
// /*     if (ch == 'h') { calculate_map_stats(1);} // prints help */
// /* ---------------------------------------------------------------------------------------------------- */
//     //for (i=1;i<=multi_mask_reg_num;i++)   {  plot_mm_reg_number(i); }
//     //if (ch == 'd') { maxT = 0.5*sqrt(varianceT); }     if (ch == 'x') { maxT = maxT-0.5*sqrt(varianceT); }
//     if (ch != 't') printf("minT = %lE, maxT = %lE\n",minT,maxT);
//     if (output != 1) ch=' ';
  } while (ch != ' ' && ch != 'Q');

  if (output != 0) {     cpgclos();   }

  if (ch == 'Q') ret = 1; else ret = 0;


  delete [] tab;
  return ret;
}

//************************************************************************
float* mscsMapPlot::replot(mscsMapProj& flat, double& Tmin, double& Tmax, matrix<double>& m, matrix<long>& mask ) {
  // Tmin=flat.mapSph().getMinT(); 
  // Tmax=flat.mapSph().getMaxT(); 
  //  flat.setMap(*this); 
  //flat.resetZoom(); 
  flat.project_map(); 
  long n;
  double* dtab=flat.getFlatArray(mask,n);   
   // dtab= cpeds_matrix2array(m, n,true); 
  float *tab=float_double_array(dtab,n);  
  delete [] dtab;
  return tab;
}

//************************************************************************
double mscsMapPlot::getMatrixCoord(double v,double min,double max, long size) {
  return (v-min)*(double)size/(max-min);
}


//************************************************************************
void mscsMapPlot::generateMapColors(string color_scheme,long color_num, double Tmin, double Tmax) {

    mscsColormap* mcol = new mscsColormap(color_scheme,color_num,Tmin,Tmax);

    cpgscir(0,color_num-1);

    float *lev=mcol->export_thresf_norm();
    float *r=mcol->export_rf();
    float *b=mcol->export_bf();
    float *g=mcol->export_gf();
    float contra=1,bright=0.5;
    cpgctab(lev,r,g,b,color_num,contra,bright);

    // mcol->print_colors();
    delete r;  delete g;   delete b;   delete lev;
    delete mcol;

}

//************************************************************************
// void map_class::plot_mm_reg_number(long reg) {
//   double factor;
//   double factorx, factory;

//   long i,num=0;
//   float x,y;
//   filenamestr tmpch;
//   for (i=0;i<pix_num;i++) { if (map->m[i] == reg) { num=i; i=pix_num; } }
//   i=num;
//   factor = 1.00;//(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed
//   factorx = factor*(maxx_flat-minx_flat)/(double)xpix_num_flat; /// how many times the map was squeezed
//   factory = factor*(maxy_flat-miny_flat)/(double)ypix_num_flat; /// how many times the map was squeezed
// /*   ypix_num_flat = (int)(ceil((maxy_flat-miny_flat)/factor))+1; */
//   pix_num_flat = (long int)xpix_num_flat*(long int)ypix_num_flat;

//   x = (float)ceil((flatcoord[i].x-minx_flat)/factorx);
//   y = (float)ceil((flatcoord[i].y-miny_flat)/factory);

//   sprintf(tmpch,"%li",reg);
// /*   cpgmtxt("B",0,tmpch); */
//   cpgptxt(x,y,0,0,tmpch);

// }

// //************************************************************************
// // this makes the histogram of temperatuers (values in general) in the map->T structure
// // excluding the map->m area (masked areas) and excluding the values in the map that are
// // listed in exclude_values_list array if such exists. -- not implemented yet
// int map_class::plot_temperature_histogram(string filename) {
//   long int i,j=0,num,imin,imax;
//   float * hist, *xbin; 
//   double * tmpmap = new double[pix_num],*subtmpmap; 
//   float min,max,x,y;
//   char ch;
//   int retval='a';
//   FILE *f;

//   if (map->m != NULL) for (i=0;i<pix_num;i++) { if (map->m[i] != 0) { tmpmap[j] = map->T[i]; j++; } }
//   else for (i=0;i<pix_num;i++) {  tmpmap[j] = map->T[i]; j++; }
// /*   printf("*******control pix_num=%li, j=%li\n",pix_num,j);  */
//   subtmpmap = new double[j]; for (i=0;i<j;i++) { subtmpmap[i] = tmpmap[i]; } 
//   delete [] tmpmap; 

//   calculate_map_stats(1);
//   num = 500;
//   xbin = new float[num]; 
//   hist = cpeds_bin_data_float(j,subtmpmap,num,xbin); 
//   cpeds_find_minmax_value(hist,num,&min,&max,&imin,&imax);  
//   if (filename.size() != 0) { f = fopen(filename.c_str(),"w");  
//     if (f==NULL) printf("ERROR: could not open file for writting\n"); 
//     else { 
//       for (i=0;i<num;i++) { fprintf(f,"%lE %lE\n",xbin[i],hist[i]);} 
//     //printf("%E %E\n",xbin[i],hist[i]);
//       fclose(f); 
//     }
//   } // save to file if requested
  
//   printf("|%s>  -- min_freq = %E, max_freq=%E\n",object_name.c_str(),min,max);
//   cpgbeg(0,"/XWINDOW",1,1);
//   cpgask(false);
//   cpgenv(minT, maxT, min, max, 0, 2); //coordinates
//   cpglab("temp", "freq", "map temperature histogram"); //axes
//   cpgsci(9);
//   cpgbin((int)num,xbin,hist,false);
//   do {
//     cpgcurs(&x,&y,&ch);
//   } while ((ch != ' ') && (ch != 'p'));
//   delete [] hist; delete [] xbin;   
//   cpgclos();
//   if (ch == ' ') retval = 0; 
//   if (ch == 'p') retval = -2;

//   return retval;
// }
