#include "cpeds-project.h"
#include <proj_api.h>

// ****************************************************************************************************
// cpedsProject::cpedsProject(QList<cpedsDirection> &dirs, string projection) {
  
// }
// ****************************************************************************************************
// cpedsProject::cpedsProject(QList<cpedsPoint2D> &points, string projection);
// ****************************************************************************************************
cpedsProject::cpedsProject(const cpedsDirectionSet &dirs, string projection) {
  _ds=dirs;
  _projection=projection;
  _cpeds_project.invertLongitudeOnProject=true;
}
// ****************************************************************************************************
cpedsProject::cpedsProject(const cpedsPointSet3D &points, string projection) {
  _ps=points;
  _projection=projection;
  _cpeds_project.invertLongitudeOnProject=true;
}
// ****************************************************************************************************
cpedsProject::~cpedsProject() {}
// ****************************************************************************************************
cpedsPointSet3D cpedsProject::projectOnPlane(const cpedsDirection& n, string projFormat) {
	projUV p;
	projPJ pj;
	char tmpch[1000];
	string s;
	long i;
	projectionDirection()=n*PI180;
	
	string format="+proj=%s +lat_0=%lE +lon_0=%lEW +R=%lE";
	if (projFormat=="") { // default case typically for Mollweide projections
		if (n.lon()>=0 and n.lon()<180) {
			format="+proj=%s +lat_0=%lE +lon_0=%lEW +R=%lE";
			sprintf(tmpch,format.c_str(),getProjection().c_str(), n.lat(), n.lon(), sqrt(2)/4);
		}
		else {
			format="+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE";	  
			sprintf(tmpch,format.c_str(),getProjection().c_str(), n.lat(), 360.0-n.lon(), sqrt(2)/4);
		}
	}
	else {
		sprintf(tmpch,projFormat.c_str(),getProjection().c_str(), n.lat(), n.lon(), sqrt(2)/4);	  
	}
	
	
	//  printf("proj string: %s\n",tmpch);
	s=tmpch;
	pj = pj_init_plus(s.c_str());
	
	// if (projection == "stere" ) { maxX=0.707106799968253; maxY=0.707106799968253; }
	// x=cpeds_get_min(x,endX);      x=cpeds_get_max(x,startX);
	// y=cpeds_get_min(y,endY);      y=cpeds_get_max(y,startY);
	
	_ps.clear();
	long N=directionsCount();
	for (i=0;i<N;i++) {
		
		if (_cpeds_project.invertLongitudeOnProject)
			p.u=-_ds.at(i).lon(); // conversion to proj lon. convention for maps plotting
		else
			p.u=_ds.at(i).lon(); // some other usages
		p.v=_ds.at(i).lat(); 
		p = pj_fwd(p,pj);
		// if (pj_errno==0) {
		//   l=p.u; phi=p.u;
		//   b=p.v; th=PIsnd-p.v;
		//   p = pj_fwd(p,pj);
		// }
		_ps.append(cpedsPoint3D(p.u, p.v, _ds.at(i).val() ));
		
		//    printf("i=%li pu %lE pv %lE l %lE b %lE\n",i, p.u,p.v,-_ds.at(i).lon(),_ds.at(i).lat());
		
		// if (fabs(p.u-x) > acc || fabs(p.v-y) > acc || pj_errno!=0) { set_ctl(i,j,bad); set_T(i,j,0);  }
		// else { 
		//   cpeds_check_thphi(&th,&phi);
		//   l=-phi; b=PIsnd-th;
		//   if (map_src->get_m(l,b) == 0.0) { set_ctl(i,j,masked); set_T(i,j,0); }
		//   else { set_ctl(i,j,ok); set_T(i,j,map_src->get_T(l,b)); }
		// }
	}
	pj_free(pj);

	return _ps;  
}

// ****************************************************************************************************
cpedsDirectionSet cpedsProject::projectOnSphere(const cpedsDirection& n, string projFormat,bool usePointVals) {
	
	projUV p;
	projPJ pj;
	char tmpch[1000];
	string s;
	long i;
	
	projectionDirection()=n*PI180;
	if (projFormat=="") { // default case typically for Mollweide projections

		if (n.lon()>=0 and n.lon()<180) 
			sprintf(tmpch,"+proj=%s +lat_0=%lE +lon_0=%lEW +R=%lE",getProjection().c_str(), n.lat(), n.lon(), sqrt(2)/4);
		else
			sprintf(tmpch,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE",getProjection().c_str(), n.lat(), 360.0-n.lon(), sqrt(2)/4);
	}
	else {
		sprintf(tmpch,projFormat.c_str(),getProjection().c_str(), n.lat(), n.lon(), sqrt(2)/4);	  		
	}
		
	//  sprintf(tmpch,"+proj=%s +lat_0=%lE +lon_0=%lEw +R=%lE",getProjection().c_str(), projectionDirection().lat(), PI-projectionDirection().lon(), sqrt(2)/4);
	s=tmpch;
	pj = pj_init_plus(s.c_str());
	
	long N=pointsCount();
	_ds.clear();
	for (i=0;i<N;i++) {
		if (_cpeds_project.invertLongitudeOnProject)
			p.u=-_ps.at(i).x();
		else 
			p.u=_ps.at(i).x();
		p.v=_ps.at(i).y(); 
		p = pj_inv(p,pj);
		if (usePointVals)
			_ds.append(cpedsDirection(p.u,p.v,_ps.val(i)));
		else
			_ds.append(cpedsDirection(p.u,p.v,_ps.at(i).z()));
	}
	pj_free(pj);
	
	return _ds;
}
/***************************************************************************************/
double cpedsProject::getLengthScaleFromAng(double ang) {
	double ls;
	cpedsDirectionSet ds;
	ds.append(cpedsDirection(-ang/2.0,0.0));
	ds.append(cpedsDirection( ang/2.0,0.0));
	ds.append(cpedsDirection(0.0,-ang/2.0));
	ds.append(cpedsDirection(0.0, ang/2.0));
	cpedsPointSet3D ps = cpedsProject(ds,getProjection()).projectOnPlane(cpedsDirection(0,0));
	ls=( fabs(ps[1].x()-ps[0].x()) +  fabs(ps[3].y()-ps[2].y()) )/2;  
	return ls;
}
/* ******************************************************************************************** */
