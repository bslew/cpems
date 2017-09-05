/*!
 * \file Mscs-WMAPspecifications.cpp
 *
 *  Created on: Apr 24, 2010
 *      Author: blew
 */

#include "Mscs-WMAPspecifications.h"

/* **************************************************************************************************** */
mscsWMAPspecifications::mscsWMAPspecifications() {

  C_2WMAP = 0.12844E-09; // normalization factor;

  // definition of temperature uncertainties in different frequency channels
  // WMAP I yr data.
  sigma0_K1 = 1.439E-3; sigma0_K2 = 1.464E-3; // K1=K1 K2=Ka1 these values are actually from wmap3 TO BE VERIFIED
  sigma0_Q1 = 2.26677E-3; sigma0_Q2 = 2.15567E-3;
  sigma0_V1 = 3.28789E-3; sigma0_V2 = 2.93683E-3;
  sigma0_W1 = 5.85196E-3; sigma0_W2 = 6.53276E-3; sigma0_W3 = 6.88032E-3; sigma0_W4 = 6.72537E-3;

  sigma0_WMAP[0] = sigma0_K1;   sigma0_WMAP[1] = sigma0_K2;
  sigma0_WMAP[2] = sigma0_Q1;   sigma0_WMAP[3] = sigma0_Q2;
  sigma0_WMAP[4] = sigma0_V1;   sigma0_WMAP[5] = sigma0_V2;
  sigma0_WMAP[6] = sigma0_W1;   sigma0_WMAP[7] = sigma0_W2;
  sigma0_WMAP[8] = sigma0_W3;   sigma0_WMAP[9] = sigma0_W4;

// sigma0 in [K] for various DAs of the WMAP
  // WMAP III yr data. from astro-ph/0603451 and astro-ph/0603452
  sigma0_K1 = 1.439E-3; sigma0_K2 = 1.464E-3; // K1=K1 K2=Ka1
  sigma0_Q1 = 2.245E-3; sigma0_Q2 = 2.135E-3;
  sigma0_V1 = 3.304E-3; sigma0_V2 = 2.946E-3;
  sigma0_W1 = 5.883E-3; sigma0_W2 = 6.532E-3; sigma0_W3 = 6.885E-3; sigma0_W4 = 6.744E-3;

  sigma0_WMAP[10] = sigma0_K1;   sigma0_WMAP[11] = sigma0_K2;
  sigma0_WMAP[12] = sigma0_Q1;   sigma0_WMAP[13] = sigma0_Q2;
  sigma0_WMAP[14] = sigma0_V1;   sigma0_WMAP[15] = sigma0_V2;
  sigma0_WMAP[16] = sigma0_W1;   sigma0_WMAP[17] = sigma0_W2;
  sigma0_WMAP[18] = sigma0_W3;   sigma0_WMAP[19] = sigma0_W4;

// sigma0 in [K] for various DAs of the WMAP
  // WMAP 5 yr data. from http://lambda.gsfc.nasa.gov/product/map/dr3/pub_papers/fiveyear/basic_results/wmap5basic.pdf
  sigma0_K1 = 1.436E-3; sigma0_K2 = 1.470E-3; // K1=K1 K2=Ka1
  sigma0_Q1 = 2.254E-3; sigma0_Q2 = 2.141E-3;
  sigma0_V1 = 3.314E-3; sigma0_V2 = 2.953E-3;
  sigma0_W1 = 5.899E-3; sigma0_W2 = 6.565E-3; sigma0_W3 = 6.926E-3; sigma0_W4 = 6.761E-3;

  sigma0_WMAP[20] = sigma0_K1;   sigma0_WMAP[21] = sigma0_K2;
  sigma0_WMAP[22] = sigma0_Q1;   sigma0_WMAP[23] = sigma0_Q2;
  sigma0_WMAP[24] = sigma0_V1;   sigma0_WMAP[25] = sigma0_V2;
  sigma0_WMAP[26] = sigma0_W1;   sigma0_WMAP[27] = sigma0_W2;
  sigma0_WMAP[28] = sigma0_W3;   sigma0_WMAP[29] = sigma0_W4;

  // sigma0 in [K] for various DAs of the WMAP for Stokes I
    // WMAP 7 yr data. from http://lambda.gsfc.nasa.gov/product/map/dr4/pub_papers/sevenyear/basic_results/wmap_7yr_basic_results.pdf p.5
    sigma0_K1 = 1.437E-3; sigma0_K2 = 1.470E-3; // K1=K1 K2=Ka1
    sigma0_Q1 = 2.254E-3; sigma0_Q2 = 2.140E-3;
    sigma0_V1 = 3.319E-3; sigma0_V2 = 2.955E-3;
    sigma0_W1 = 5.906E-3; sigma0_W2 = 6.572E-3; sigma0_W3 = 6.941E-3; sigma0_W4 = 6.778E-3;

    sigma0_WMAP[30] = sigma0_K1;   sigma0_WMAP[31] = sigma0_K2;
    sigma0_WMAP[32] = sigma0_Q1;   sigma0_WMAP[33] = sigma0_Q2;
    sigma0_WMAP[34] = sigma0_V1;   sigma0_WMAP[35] = sigma0_V2;
    sigma0_WMAP[36] = sigma0_W1;   sigma0_WMAP[37] = sigma0_W2;
    sigma0_WMAP[38] = sigma0_W3;   sigma0_WMAP[39] = sigma0_W4;


  // WMAP Beam sizes - the Full Width at Half Maximum in degrees for different DA in degrees
  thFWHM_Q = 0.49;  thFWHM_V = 0.33; thFWHM_W = 0.21;
  //************************************************************************
}

/* **************************************************************************************************** */
double mscsWMAPspecifications::get_FWHM(DAnames da) {
  if (da==Q1 || da==Q2 || da==Qall) { return thFWHM_Q; }
  if (da==V1 || da==V2 || da==Vall) { return thFWHM_V; }
  if (da==W1 || da==W2 || da==W3 || da==W4 || da==Wall) { return thFWHM_W; }
  return -1;
}

/* **************************************************************************************************** */
double mscsWMAPspecifications::get_sigma0(DAnames da, WMAPversion v) {
	switch (v) {
		case WMAP1yr:
		    if (da==K1) { return sigma0_WMAP[0]; }
		    if (da==K2) { return sigma0_WMAP[1]; }
		    if (da==Q1) { return sigma0_WMAP[2]; }
		    if (da==Q2) { return sigma0_WMAP[3]; }
		    if (da==V1) { return sigma0_WMAP[4]; }
		    if (da==V2) { return sigma0_WMAP[5]; }
		    if (da==W1) { return sigma0_WMAP[6]; }
		    if (da==W2) { return sigma0_WMAP[7]; }
		    if (da==W3) { return sigma0_WMAP[8]; }
		    if (da==W4) { return sigma0_WMAP[9]; }
			break;
		case WMAP3yrs:
		    if (da==K1) { return sigma0_WMAP[10]; }
		    if (da==K2) { return sigma0_WMAP[11]; }
		    if (da==Q1) { return sigma0_WMAP[12]; }
		    if (da==Q2) { return sigma0_WMAP[13]; }
		    if (da==V1) { return sigma0_WMAP[14]; }
		    if (da==V2) { return sigma0_WMAP[15]; }
		    if (da==W1) { return sigma0_WMAP[16]; }
		    if (da==W2) { return sigma0_WMAP[17]; }
		    if (da==W3) { return sigma0_WMAP[18]; }
		    if (da==W4) { return sigma0_WMAP[19]; }
		    break;
		case WMAP5yrs:
		    if (da==K1) { return sigma0_WMAP[20]; }
		    if (da==K2) { return sigma0_WMAP[21]; }
		    if (da==Q1) { return sigma0_WMAP[22]; }
		    if (da==Q2) { return sigma0_WMAP[23]; }
		    if (da==V1) { return sigma0_WMAP[24]; }
		    if (da==V2) { return sigma0_WMAP[25]; }
		    if (da==W1) { return sigma0_WMAP[26]; }
		    if (da==W2) { return sigma0_WMAP[27]; }
		    if (da==W3) { return sigma0_WMAP[28]; }
		    if (da==W4) { return sigma0_WMAP[29]; }
		    break;
		case WMAP7yrs:
			if (da==K1) { return sigma0_WMAP[30]; }
			if (da==K2) { return sigma0_WMAP[31]; }
			if (da==Q1) { return sigma0_WMAP[32]; }
			if (da==Q2) { return sigma0_WMAP[33]; }
			if (da==V1) { return sigma0_WMAP[34]; }
			if (da==V2) { return sigma0_WMAP[35]; }
			if (da==W1) { return sigma0_WMAP[36]; }
			if (da==W2) { return sigma0_WMAP[37]; }
			if (da==W3) { return sigma0_WMAP[38]; }
			if (da==W4) { return sigma0_WMAP[39]; }
			break;
		default:
			break;
	}

  return -1;
}
/* **************************************************************************************************** */
mscsMap mscsWMAPspecifications::get_Nobs(DAnames da, WMAPversion v) {
	mscsMap Nobs;
	Nobs.loadfits(MSCS_WMAP_DATA_DIR+get_Nobs_fileName(da,v),"N_OBS","N");
	return Nobs;
}
/* **************************************************************************************************** */
const string mscsWMAPspecifications::get_Nobs_fileName(DAnames da, WMAPversion v) {

	switch (v) {
		case WMAP1yr:
		    if (da==K1) { return MSCS_WMAP_NOBS_FILE_K1_1YR; }
		    if (da==K2) { return MSCS_WMAP_NOBS_FILE_K2_1YR; }
		    if (da==Q1) { return MSCS_WMAP_NOBS_FILE_Q1_1YR; }
		    if (da==Q2) { return MSCS_WMAP_NOBS_FILE_Q2_1YR; }
		    if (da==V1) { return MSCS_WMAP_NOBS_FILE_V1_1YR; }
		    if (da==V2) { return MSCS_WMAP_NOBS_FILE_V2_1YR; }
		    if (da==W1) { return MSCS_WMAP_NOBS_FILE_W1_1YR; }
		    if (da==W2) { return MSCS_WMAP_NOBS_FILE_W2_1YR; }
		    if (da==W3) { return MSCS_WMAP_NOBS_FILE_W3_1YR; }
		    if (da==W4) { return MSCS_WMAP_NOBS_FILE_W4_1YR; }
			break;
		case WMAP3yrs:
		    if (da==K1) { return MSCS_WMAP_NOBS_FILE_K1_3YR; }
		    if (da==K2) { return MSCS_WMAP_NOBS_FILE_K2_3YR; }
		    if (da==Q1) { return MSCS_WMAP_NOBS_FILE_Q1_3YR; }
		    if (da==Q2) { return MSCS_WMAP_NOBS_FILE_Q2_3YR; }
		    if (da==V1) { return MSCS_WMAP_NOBS_FILE_V1_3YR; }
		    if (da==V2) { return MSCS_WMAP_NOBS_FILE_V2_3YR; }
		    if (da==W1) { return MSCS_WMAP_NOBS_FILE_W1_3YR; }
		    if (da==W2) { return MSCS_WMAP_NOBS_FILE_W2_3YR; }
		    if (da==W3) { return MSCS_WMAP_NOBS_FILE_W3_3YR; }
		    if (da==W4) { return MSCS_WMAP_NOBS_FILE_W4_3YR; }
		    break;
		case WMAP5yrs:
		    if (da==K1) { return MSCS_WMAP_NOBS_FILE_K1_5YR; }
		    if (da==K2) { return MSCS_WMAP_NOBS_FILE_K2_5YR; }
		    if (da==Q1) { return MSCS_WMAP_NOBS_FILE_Q1_5YR; }
		    if (da==Q2) { return MSCS_WMAP_NOBS_FILE_Q2_5YR; }
		    if (da==V1) { return MSCS_WMAP_NOBS_FILE_V1_5YR; }
		    if (da==V2) { return MSCS_WMAP_NOBS_FILE_V2_5YR; }
		    if (da==W1) { return MSCS_WMAP_NOBS_FILE_W1_5YR; }
		    if (da==W2) { return MSCS_WMAP_NOBS_FILE_W2_5YR; }
		    if (da==W3) { return MSCS_WMAP_NOBS_FILE_W3_5YR; }
		    if (da==W4) { return MSCS_WMAP_NOBS_FILE_W4_5YR; }
		    break;
		case WMAP7yrs:
		    if (da==K1) { return MSCS_WMAP_NOBS_FILE_K1_7YR; }
		    if (da==K2) { return MSCS_WMAP_NOBS_FILE_K2_7YR; }
		    if (da==Q1) { return MSCS_WMAP_NOBS_FILE_Q1_7YR; }
		    if (da==Q2) { return MSCS_WMAP_NOBS_FILE_Q2_7YR; }
		    if (da==V1) { return MSCS_WMAP_NOBS_FILE_V1_7YR; }
		    if (da==V2) { return MSCS_WMAP_NOBS_FILE_V2_7YR; }
		    if (da==W1) { return MSCS_WMAP_NOBS_FILE_W1_7YR; }
		    if (da==W2) { return MSCS_WMAP_NOBS_FILE_W2_7YR; }
		    if (da==W3) { return MSCS_WMAP_NOBS_FILE_W3_7YR; }
		    if (da==W4) { return MSCS_WMAP_NOBS_FILE_W4_7YR; }
			break;
		default:
			break;
	}

	return "";
}

/* **************************************************************************************************** */
string mscsWMAPspecifications::get_beamtf_fileName(DAnames da, WMAPversion v) {

	switch (v) {
		case WMAP1yr:
		    if (da==K1) { return MSCS_WMAP_BEAMTF_FILE_K1_1YR; }
		    if (da==K2) { return MSCS_WMAP_BEAMTF_FILE_K2_1YR; }
		    if (da==Q1) { return MSCS_WMAP_BEAMTF_FILE_Q1_1YR; }
		    if (da==Q2) { return MSCS_WMAP_BEAMTF_FILE_Q2_1YR; }
		    if (da==V1) { return MSCS_WMAP_BEAMTF_FILE_V1_1YR; }
		    if (da==V2) { return MSCS_WMAP_BEAMTF_FILE_V2_1YR; }
		    if (da==W1) { return MSCS_WMAP_BEAMTF_FILE_W1_1YR; }
		    if (da==W2) { return MSCS_WMAP_BEAMTF_FILE_W2_1YR; }
		    if (da==W3) { return MSCS_WMAP_BEAMTF_FILE_W3_1YR; }
		    if (da==W4) { return MSCS_WMAP_BEAMTF_FILE_W4_1YR; }
			break;
		case WMAP3yrs:
		    if (da==K1) { return MSCS_WMAP_BEAMTF_FILE_K1_3YR; }
		    if (da==K2) { return MSCS_WMAP_BEAMTF_FILE_K2_3YR; }
		    if (da==Q1) { return MSCS_WMAP_BEAMTF_FILE_Q1_3YR; }
		    if (da==Q2) { return MSCS_WMAP_BEAMTF_FILE_Q2_3YR; }
		    if (da==V1) { return MSCS_WMAP_BEAMTF_FILE_V1_3YR; }
		    if (da==V2) { return MSCS_WMAP_BEAMTF_FILE_V2_3YR; }
		    if (da==W1) { return MSCS_WMAP_BEAMTF_FILE_W1_3YR; }
		    if (da==W2) { return MSCS_WMAP_BEAMTF_FILE_W2_3YR; }
		    if (da==W3) { return MSCS_WMAP_BEAMTF_FILE_W3_3YR; }
		    if (da==W4) { return MSCS_WMAP_BEAMTF_FILE_W4_3YR; }
		    break;
		case WMAP5yrs:
		    if (da==K1) { return MSCS_WMAP_BEAMTF_FILE_K1_5YR; }
		    if (da==K2) { return MSCS_WMAP_BEAMTF_FILE_K2_5YR; }
		    if (da==Q1) { return MSCS_WMAP_BEAMTF_FILE_Q1_5YR; }
		    if (da==Q2) { return MSCS_WMAP_BEAMTF_FILE_Q2_5YR; }
		    if (da==V1) { return MSCS_WMAP_BEAMTF_FILE_V1_5YR; }
		    if (da==V2) { return MSCS_WMAP_BEAMTF_FILE_V2_5YR; }
		    if (da==W1) { return MSCS_WMAP_BEAMTF_FILE_W1_5YR; }
		    if (da==W2) { return MSCS_WMAP_BEAMTF_FILE_W2_5YR; }
		    if (da==W3) { return MSCS_WMAP_BEAMTF_FILE_W3_5YR; }
		    if (da==W4) { return MSCS_WMAP_BEAMTF_FILE_W4_5YR; }
		    break;
		case WMAP7yrs:
		    if (da==K1) { return MSCS_WMAP_BEAMTF_FILE_K1_7YR; }
		    if (da==K2) { return MSCS_WMAP_BEAMTF_FILE_K2_7YR; }
		    if (da==Q1) { return MSCS_WMAP_BEAMTF_FILE_Q1_7YR; }
		    if (da==Q2) { return MSCS_WMAP_BEAMTF_FILE_Q2_7YR; }
		    if (da==V1) { return MSCS_WMAP_BEAMTF_FILE_V1_7YR; }
		    if (da==V2) { return MSCS_WMAP_BEAMTF_FILE_V2_7YR; }
		    if (da==W1) { return MSCS_WMAP_BEAMTF_FILE_W1_7YR; }
		    if (da==W2) { return MSCS_WMAP_BEAMTF_FILE_W2_7YR; }
		    if (da==W3) { return MSCS_WMAP_BEAMTF_FILE_W3_7YR; }
		    if (da==W4) { return MSCS_WMAP_BEAMTF_FILE_W4_7YR; }
			break;
		default:
			break;
	}

	return "";
}

/* **************************************************************************************************** */
const mscsWindowFunction mscsWMAPspecifications::get_beamtf(DAnames da, WMAPversion v) {
	mscsWindowFunction b;
	b.load(MSCS_WMAP_DATA_DIR+get_beamtf_fileName(da,v),true);
	return b;
}

/* **************************************************************************************************** */
double mscsWMAPspecifications::get_Cl_calibration(WMAPversion v) {
	switch (v) {
		case WMAP1yr:
			return 23.7e-10; // update this
		case WMAP3yrs:
			return 23.7e-10;
		case WMAP5yrs:
			return 24.1e-10;
		case WMAP7yrs:
			return 24.1e-10; // update this
		default:
			return 0;
	}
	return 0;

}
/***************************************************************************************/
mscsWMAPspecifications::WMAPversion mscsWMAPspecifications::getVersion(string v) {
	if (v=="WMAP1yr") { return WMAP1yr; }
	if (v=="WMAP3yrs") { return WMAP3yrs; }
	if (v=="WMAP5yrs") { return WMAP5yrs; }
	if (v=="WMAP7yrs") { return WMAP7yrs; }
	return WMAP7yrs;
}

/***************************************************************************************/
string mscsWMAPspecifications::getVersion(WMAPversion v) {
	if (v==WMAP1yr) { return "WMAP1yr"; }
	if (v==WMAP3yrs) { return "WMAP3yrs"; }
	if (v==WMAP5yrs) { return "WMAP5yrs"; }
	if (v==WMAP7yrs) { return "WMAP7yrs"; }
	return "WMAP7yrs";
}

/***************************************************************************************/
mscsWMAPspecifications::DAnames mscsWMAPspecifications::getDA(string da) {
	if (da=="K1") return K1;
	if (da=="K2") return K2;
	if (da=="Q1") return Q1;
	if (da=="Q2") return Q2;
	if (da=="V1") return V1;
	if (da=="V2") return V2;
	if (da=="W1") return W1;
	if (da=="W2") return W2;
	if (da=="W3") return W3;
	if (da=="W4") return W4;

	if (da=="Qinc") return Qinc;
	if (da=="Vinc") return Vinc;
	if (da=="Winc") return Winc;

	if (da=="VWinc") return VWinc;
	if (da=="QVWinc") return QVWinc;

	if (da=="Q_V_W_inc") return Q_V_W_inc;
	if (da=="Q_V_W_QVW_inc") return Q_V_W_QVW_inc;
	if (da=="V_W_QVW_inc") return V_W_VW_inc;

	return unknownDA;
}
/***************************************************************************************/
string mscsWMAPspecifications::getDA(DAnames da) {
	string s;
	switch (da) {
		case K1:   { s="K1"; break; }
		case K2:   { s="K2"; break; }
		case Q1:   { s="Q1"; break; }
		case Q2:   { s="Q2"; break; }
		case V1:   { s="V1"; break; }
		case V2:   { s="V2"; break; }
		case W1:   { s="W1"; break; }
		case W2:   { s="W2"; break; }
		case W3:   { s="W3"; break; }
		case W4:   { s="W4"; break; }

		case Kinc:   { s="Kinc"; break; }
		case Qinc:   { s="Qinc"; break; }
		case Vinc:   { s="Vinc"; break; }
		case Winc:   { s="Winc"; break; }
		case QVWinc: { s="QVWinc"; break; }
		case VWinc:  { s="VWinc"; break; }

		case Q_V_W_inc:     { s="Q_V_W_inc"; break; }
		case Q_V_W_QVW_inc: { s="Q_V_W_QVW_inc"; break; }
		case V_W_inc:       { s="V_W_inc"; break; }
		case V_W_VW_inc:    { s="V_W_VW_inc"; break; }
		default: break;
	}

	return s;


}
