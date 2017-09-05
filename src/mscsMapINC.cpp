/*
 * mscsMapINC.cpp
 *
 *  Created on: Apr 23, 2010
 *      Author: blew
 */

#include "mscsMapINC.h"

//****************************************************************
mscsMapINC::mscsMapINC() {
}

//****************************************************************
mscsMapINC::~mscsMapINC() {
	clearMaps();
}

//****************************************************************
void mscsMapINC::clearMaps() {
	for (int i = 0; i < maps.count(); ++i) {
		delete maps.value(i);
	}
	maps.clear();
}

//****************************************************************
void mscsMapINC::useMap(const mscsMap* map, double sigma0) {
	maps.append(new mscsMap(*map));
	sigmas.append(sigma0);
}
//****************************************************************
const mscsMap* mscsMapINC::makeINC(long DAst, long DAen, string weight, bool precalib, double precalibrate, bool calib, double calibrate) {
//	void mscsMap::mk_INC_DA(long DAst,long DAen,mscsMap* DAs,string weight) {
	  long i,j,l;
	  double w;

	  msgs->say("Generating inverse noise co-added map",High);
	  msgs->say("Weightning: "+weight,Medium);
	  // prepare the INC map space
	  set_nside(maps.at(0)->nside());
	  long pix_num=pixNum();
	  makekill_space_manager("make","T");
	  clear_map();
	  if (DAst==-1) DAst=0;
	  if (DAen==-1) DAen=maps.count()-1;
	  if (DAst>=DAen) return NULL;
	  msgs->say("DA start: "+msgs->toStr(DAst)+" DA end: "+msgs->toStr(DAen),Medium);

	  if (precalib) { for (j=DAst;j<=DAen;j++) { maps.value(j)->import_map_data(*this,"m",1); maps.value(j)->shift_mean_to(precalibrate,true); maps.value(j)->makekill_space_manager("kill","m");    	}	  }


	  if (weight == "inv_noise") {  // inverse noise weighting including number of observations of a given pixel
	    for (i=0;i<pix_num;i++) { // loop for pixel number
	    	for (j=DAst;j<=DAen;j++) { // loop for DA number
	    	  w = 0; for (l=DAst;l<=DAen;l++) { w += maps.value(l)->N(i) / (sigmas.at(l)*sigmas.at(l)); }
	    	  T(i) += maps.value(j)->get_T(i) * (maps.value(j)->get_N(i)/(sigmas.value(j)*sigmas.value(j))) / w;
//	    	  printf("%lE\n", (maps.value(j)->get_N(i)/(sigmas.value(j)*sigmas.value(j))));
	    	}
	    }
	  }

	  if (weight == "simple") { // inverse noise weighting without the Nobs in each pixel information
	    for (i=0;i<pix_num;i++) { // loop for pixel number
	      for (j=DAst;j<=DAen;j++) { // loop for DA number
	    	  w = 0; for (l=DAst;l<=DAen;l++) { w += 1.0 / (sigmas.value(l)*sigmas.value(l)); }
	    	  T(i) += maps.value(j)->get_T(i) / (sigmas.value(j)*sigmas.value(j)) / w;
	      }
	    }
	  }
	  if (calib) { shift_mean_to(calibrate,true); }

	  return this;
}
