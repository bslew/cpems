#include <complex>
#include "Mscs-map_proj.h"

// ****************************************************************************************************
mscsMapProj::mscsMapProj() {
  sizeX=sizeY=pix_num=0;
  map=NULL;
  ctl=NULL;
  setProjection("");
  initiated=false;
  reset();
  initiate_range_variables();
}

// ****************************************************************************************************
mscsMapProj::mscsMapProj(long X, long Y) {
  sizeX=X;
  sizeY=Y;
  pix_num=X*Y;
  map = new matrix<double>(sizeY,sizeX);
  ctl = new matrix<long>(sizeY,sizeX);
  clear();
  setProjection("");
  initiated=false;
  reset();
  initiate_range_variables();
}

// ****************************************************************************************************
mscsMapProj::mscsMapProj(const mscsMap& sphT, long X, long Y, string proj, double lonOff, double latOff) {
  sizeX=X;
  sizeY=Y;
  pix_num=X*Y;
  map = new matrix<double>(sizeY,sizeX);
  ctl = new matrix<long>(sizeY,sizeX);
  clear();
  setProjection(proj);
  initiate_range_variables();
  sphMap=sphT;
  reset();
  set_LonLatOffset(lonOff,latOff);
  project_map();
  initiated=true;
}

// ****************************************************************************************************
mscsMapProj::~mscsMapProj() {
  delete map;
  delete ctl;
}

// ****************************************************************************************************
string mscsMapProj::getProjInitString() { 
  filenamestr tmpch;
  sprintf(tmpch,"+proj=%s  +lon_0=%lfe  +lat_0=%lf +R=%lf",projection().c_str(), lon0()*PI180inv, lat0()*PI180inv,  sqrt(2)/4);
  _projInit=tmpch;
  return _projInit; 
} 
// ****************************************************************************************************
void mscsMapProj::initiate_range_variables() {
  double pixfrac=0.01;
  double avpixSize= 1.0/( cpeds_get_max(sizeX,sizeY) );
  //  double avpixSize= 1.0;

  acc=pixfrac * avpixSize;
  
  if (projection() == "moll" ) { maxX=1.0; maxY=0.5; }
  if (projection() == "merc" ) { maxX=1.11072076405995; maxY=2.10808868761894; }
  if (projection() == "cc" ) { maxX=1.11072076405995; maxY=2.10808868761894; }
  if (projection() == "mill" ) { maxX=1.11072076405995; maxY=0.814379316293293; }
  if (projection() == "cea" ) { maxX=1.11072076405995; maxY=0.3535534; }
  if (projection() == "robin" ) { maxX=0.942668712457677; maxY=0.478110252103702; }
  if (projection() == "ortho" ) { maxX=0.3535534; maxY=0.3535534; }
  if (projection() == "stere" ) { maxX=0.707106799968253; maxY=0.707106799968253; }

  minX=-maxX; 
  minY=-maxY;
  getProjInitString();
}



// ****************************************************************************************************
void mscsMapProj::set_LonLatOffset(double lon, double lat) {
  setLon0(lon);
  setLat0(lat);
}

// ****************************************************************************************************
void mscsMapProj::setZoomWindow(long fromX,long fromY,long toX,long toY) {
  zoomFromX=fromX;
  zoomFromY=fromY;
  zoomToX=toX;
  zoomToY=toY;

}

// ****************************************************************************************************
void mscsMapProj::project_map() {
  projUV p;
  projPJ pj;
  double x,y,dx,dy; // variables in the projected plane domain 
  double startX, startY, endX, endY; // ranges in the projected plane domain defining the current zoom box
  double l,b,th,phi;
  long i,j;

  // sprintf(tmpch,"+proj=%s +lat_0=%lE +lon_0=%lEw +R=%lE",projection().c_str(), lat0(), lon0(), sqrt(2)/4);
  
  pj = pj_init_plus(getProjInitString().c_str());
  msgs->say("initiating proj with: "+getProjInitString(),Medium);
  
  initiate_range_variables();

  // startX=(double)zoomFromX/(double)(sizeX-1) * 2.0*maxX - maxX; // start from this X coordinate in the projected plane 
  // startY=(double)zoomFromY/(double)(sizeY-1) * 2.0*maxY - maxY; // start from this Y coordinate in the projected plane 
  // endX=(double)zoomToX/(double)(sizeX-1) * 2.0*maxX - maxX; // end at this X coordinate in the projected plane 
  // endY=(double)zoomToY/(double)(sizeY-1) * 2.0*maxY - maxY; // end at this Y coordinate in the projected plane 


  // this is the zoom into window in arbitrary location in the projection plane
  // startX=(double)zoomFromX/(double)(sizeX-1) * (getMaxX()-getMinX()) + getMinX(); // start from this X coordinate in the projected plane 
  // startY=(double)zoomFromY/(double)(sizeY-1) * (getMaxY()-getMinY()) + getMinY(); // start from this X coordinate in the projected plane 
  // endX=(double)zoomToX/(double)(sizeX-1) * (getMaxX()-getMinX()) + getMinX(); // start from this X coordinate in the projected plane 
  // endY=(double)zoomToY/(double)(sizeY-1) * (getMaxY()-getMinY()) + getMinY(); // start from this X coordinate in the projected plane 

  // this is for eg. orthographic projections with zoom always in the center of the initial field.
  startX=getMinX()/zoom();
  startY=getMinY()/zoom();
  endX=getMaxX()/zoom();
  endY=getMaxY()/zoom();

  dx=fabs(endX-startX)/(double)(sizeX); 
  dy=fabs(endY-startY)/(double)(sizeY); 

  printf("maxX %lE sizeX %li zoomFromX %li zoomToX %li\n",maxX,sizeX,zoomFromX, zoomToX);
  printf("maxY %lE sizeY %li zoomFromY %li zoomToY %li\n",maxY,sizeY,zoomFromY, zoomToY);
  printf("X: (%lE, %lE), Y: (%lE %lE), dx: %lE, dy: %lE\n",startX,endX, startY,endY, dx,dy);

  x=dx/2.0+startX;    
  for (i=0;i<sizeX;i++) { 
    y=dy/2.0+startY;
    for (j=0;j<sizeY;j++) { 
      x=cpeds_get_min(x,endX);      x=cpeds_get_max(x,startX);
      y=cpeds_get_min(y,endY);      y=cpeds_get_max(y,startY);
      
      p.u=x;
      p.v=y; 
      p = pj_inv(p,pj);
      if (pj_errno==0) {
	l=p.u; phi=p.u;
	b=p.v; th=PIsnd-p.v;
	p = pj_fwd(p,pj);
      
      
      
	if (fabs(p.u-x) > acc || fabs(p.v-y) > acc || pj_errno!=0) { set_ctl(i,j,bad); set_T(i,j,0);  }
	else { 
	  cpeds_check_thphi(&th,&phi);
	  l=-phi; b=PIsnd-th;
	  if (sphMap.get_m(l,b) == 0.0) { set_ctl(i,j,masked); set_T(i,j,0); }
	  else { set_ctl(i,j,ok); set_T(i,j,sphMap.get_T(l,b)); }
	}
      }
      
      y+=dy;
    }
    x+=dx;
  }
  


/*   if (projection == "merc") { */


/*     dx=2.0/(double)(sizeX); // the X size of the map in this coordinate system will span from -1 to 1 */
/*     dy=1.0/(double)(sizeY); // the Y size of the map in this coordinate system will span from -0.5 to 0.5 */
/*     x=dx/2.0-1.0; */
    
/*     for (i=0;i<sizeX;i++) {  */
/*       y=dy/2.0-0.5; */
/*       for (j=0;j<sizeY;j++) {  */
/* 	x=cpeds_get_min(x,1.0);      x=cpeds_get_max(x,-1.0); */
/* 	y=cpeds_get_min(y,0.5);      y=cpeds_get_max(y,-0.5); */
	
/* 	p.u=x; */
/* 	p.v=y;  */
/* 	p = pj_inv(p,pj); */
/* 	l=p.u; phi=p.u; */
/* 	b=p.v; th=PIsnd-p.v; */
/* 	p = pj_fwd(p,pj); */
	
	
/* 	if (fabs(p.u-x) > acc || fabs(p.v-y) > acc) { set_ctl(i,j,bad); set_T(i,j,0);  } */
/* 	else {  */
/* 	  cpeds_check_thphi(&th,&phi); */
/* 	  l=-phi; b=PIsnd-th; */
/* 	  if (sphT->get_m(l,b) == 0.0) { set_ctl(i,j,masked); set_T(i,j,0); } */
/* 	  else { set_ctl(i,j,ok); set_T(i,j,sphT->get_T(l,b)); } */
/* 	} */


/* 	y+=dy; */
/*       } */
/*       x+=dx; */
/*     } */

/*   } */


  pj_free(pj);


}


// ****************************************************************************************************
double mscsMapProj::get_T(long x, long y) {
  x=cpeds_get_min(x,sizeX-1);  x=cpeds_get_max(x,0);
  y=cpeds_get_min(y,sizeY-1);  y=cpeds_get_max(y,0);
  return (*map)(y,x);
}

// ****************************************************************************************************
cpedsDirection mscsMapProj::getCoordLB(long x,long y) {
  projUV p;
  projPJ pj;
  double X,Y,th,phi;

  cpedsDirection n;

  x=cpeds_get_min(x,sizeX-1);  x=cpeds_get_max(x,0);
  y=cpeds_get_min(y,sizeY-1);  y=cpeds_get_max(y,0);
  //printf ("x,y = %li,%li    ",x,y);
  
  if ((*ctl)(y,x) != bad) {
    pj = pj_init_plus(getProjInitString().c_str());
    
    X=2*x*maxX/(double)sizeX-maxX;
    Y=2*y*maxY/(double)sizeY-maxY;
    X=cpeds_get_min(X,maxX);      X=cpeds_get_max(X,-maxX);
    Y=cpeds_get_min(Y,maxY);      Y=cpeds_get_max(Y,-maxY);
    
    //printf("%lE %lE\n",X,Y);
  
    p.u=X;
    p.v=Y; 
    p = pj_inv(p,pj);
    n.lon()=p.u; phi=p.u;
    n.lat()=p.v; th=PIsnd-p.v;
    cpeds_check_thphi(&th,&phi);
    n.lon()=twoPI-phi; n.lat()=PIsnd-th;
    
    

    pj_free(pj);
  }
  else { n.lon()=0; n.lat()=0; }
  return n;

}

// ****************************************************************************************************
void mscsMapProj::set_T(long x, long y, double T) {
  (*map)(y,x)=T;
}

// ****************************************************************************************************
void mscsMapProj::set_ctl(long x, long y, long v) {
  (*ctl)(y,x)=v;
}

// ****************************************************************************************************
bool mscsMapProj::isMasked(long x, long y) {
  if ((*ctl)(y,x) == masked) return true; else return false;
}

// ****************************************************************************************************
bool mscsMapProj::isSet(long x, long y) {
  if ((*ctl)(y,x) > bad) return true; else return false;
}

// ****************************************************************************************************
void mscsMapProj::clear() {
  long i,j;
  for (i=0;i<sizeX;i++) { 
    for (j=0;j<sizeY;j++) { (*map)(j,i)=0.0; (*ctl)(j,i)=0.0; }
  }
}
// ****************************************************************************************************
void mscsMapProj::zoomIn(long x,long y, double zoomLevel) {
  // derive the direction corresponding to x,y
  cpedsDirection n=getCoordLB(x,y);
  msgs->say("setting new zoom box centered at:",High);
  n.print_direction("zoom field center");
  n.lon()=twoPI-n.l();
  setLonLat0(n);
  
  // set new zoom window
  long delta;
  delta=long(double(zoomToX-zoomFromX)/zoomLevel/2.0);
  zoomFromX=x-delta; zoomToX=x+delta;
  zoomFromX=cpeds_get_max(0,zoomFromX);   zoomToX=cpeds_get_min(getSizeX()-1,zoomToX);

  delta=long(double(zoomToY-zoomFromY)/zoomLevel/2.0);
  zoomFromY=y-delta; zoomToY=y+delta;
  zoomFromY=cpeds_get_max(0,zoomFromY);   zoomToY=cpeds_get_min(getSizeY()-1,zoomToY);
  zoom(zoomLevel);
  msgs->say("zoomFromX: "+msgs->toStr(zoomFromX)+" zoomToX: "+msgs->toStr(zoomToX)+" zoomFromY: "+msgs->toStr(zoomFromY)+" zoomToY: "+msgs->toStr(zoomToY),High);
  project_map();

}

