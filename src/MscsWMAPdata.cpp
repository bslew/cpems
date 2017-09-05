/*!
 * \file MscsWMAPdata.cpp
 *
 *  Created on: Apr 27, 2010
 *      Author: blew
 */

#include "MscsWMAPdata.h"

string MSCS_WMAP_NOBS_FILE_K1_1YR;
string MSCS_WMAP_NOBS_FILE_K2_1YR;
string MSCS_WMAP_NOBS_FILE_Q1_1YR;
string MSCS_WMAP_NOBS_FILE_Q2_1YR;
string MSCS_WMAP_NOBS_FILE_V1_1YR;
string MSCS_WMAP_NOBS_FILE_V2_1YR;
string MSCS_WMAP_NOBS_FILE_W1_1YR;
string MSCS_WMAP_NOBS_FILE_W2_1YR;
string MSCS_WMAP_NOBS_FILE_W3_1YR;
string MSCS_WMAP_NOBS_FILE_W4_1YR;

string MSCS_WMAP_NOBS_FILE_K1_3YR; //!< NAME OF THE FILE HOLDING N_OBS this should be fits file from the yearly data release maps per DA
string MSCS_WMAP_NOBS_FILE_K2_3YR;
string MSCS_WMAP_NOBS_FILE_Q1_3YR;
string MSCS_WMAP_NOBS_FILE_Q2_3YR;
string MSCS_WMAP_NOBS_FILE_V1_3YR;
string MSCS_WMAP_NOBS_FILE_V2_3YR;
string MSCS_WMAP_NOBS_FILE_W1_3YR;
string MSCS_WMAP_NOBS_FILE_W2_3YR;
string MSCS_WMAP_NOBS_FILE_W3_3YR;
string MSCS_WMAP_NOBS_FILE_W4_3YR;

string MSCS_WMAP_NOBS_FILE_K1_5YR = "wmap_imap_r9_5yr_K1_v3.fits";
string MSCS_WMAP_NOBS_FILE_K2_5YR = "wmap_imap_r9_5yr_Ka1_v3.fits";
string MSCS_WMAP_NOBS_FILE_Q1_5YR = "wmap_imap_r9_5yr_Q1_v3.fits";
string MSCS_WMAP_NOBS_FILE_Q2_5YR = "wmap_imap_r9_5yr_Q2_v3.fits";
string MSCS_WMAP_NOBS_FILE_V1_5YR = "wmap_imap_r9_5yr_V1_v3.fits";
string MSCS_WMAP_NOBS_FILE_V2_5YR = "wmap_imap_r9_5yr_V2_v3.fits";
string MSCS_WMAP_NOBS_FILE_W1_5YR = "wmap_imap_r9_5yr_W1_v3.fits";
string MSCS_WMAP_NOBS_FILE_W2_5YR = "wmap_imap_r9_5yr_W2_v3.fits";
string MSCS_WMAP_NOBS_FILE_W3_5YR = "wmap_imap_r9_5yr_W3_v3.fits";
string MSCS_WMAP_NOBS_FILE_W4_5YR = "wmap_imap_r9_5yr_W4_v3.fits";

string MSCS_WMAP_NOBS_FILE_K1_7YR = "wmap_da_imap_r9_7yr_K1_v4.fits";
string MSCS_WMAP_NOBS_FILE_K2_7YR = "wmap_da_imap_r9_7yr_Ka1_v4.fits";
string MSCS_WMAP_NOBS_FILE_Q1_7YR = "wmap_da_imap_r9_7yr_Q1_v4.fits";
string MSCS_WMAP_NOBS_FILE_Q2_7YR = "wmap_da_imap_r9_7yr_Q2_v4.fits";
string MSCS_WMAP_NOBS_FILE_V1_7YR = "wmap_da_imap_r9_7yr_V1_v4.fits";
string MSCS_WMAP_NOBS_FILE_V2_7YR = "wmap_da_imap_r9_7yr_V2_v4.fits";
string MSCS_WMAP_NOBS_FILE_W1_7YR = "wmap_da_imap_r9_7yr_W1_v4.fits";
string MSCS_WMAP_NOBS_FILE_W2_7YR = "wmap_da_imap_r9_7yr_W2_v4.fits";
string MSCS_WMAP_NOBS_FILE_W3_7YR = "wmap_da_imap_r9_7yr_W3_v4.fits";
string MSCS_WMAP_NOBS_FILE_W4_7YR = "wmap_da_imap_r9_7yr_W4_v4.fits";


/* beam transfer function file name definitions */

string MSCS_WMAP_BEAMTF_FILE_K1_1YR; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
string MSCS_WMAP_BEAMTF_FILE_K2_1YR;
string MSCS_WMAP_BEAMTF_FILE_Q1_1YR;
string MSCS_WMAP_BEAMTF_FILE_Q2_1YR;
string MSCS_WMAP_BEAMTF_FILE_V1_1YR;
string MSCS_WMAP_BEAMTF_FILE_V2_1YR;
string MSCS_WMAP_BEAMTF_FILE_W1_1YR;
string MSCS_WMAP_BEAMTF_FILE_W2_1YR;
string MSCS_WMAP_BEAMTF_FILE_W3_1YR;
string MSCS_WMAP_BEAMTF_FILE_W4_1YR;

string MSCS_WMAP_BEAMTF_FILE_K1_3YR = "wmap_K1_ampl_bl_3yr_v2.txt"; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
string MSCS_WMAP_BEAMTF_FILE_K2_3YR = "wmap_Ka1_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_Q1_3YR = "wmap_Q1_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_Q2_3YR = "wmap_Q2_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_V1_3YR = "wmap_V1_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_V2_3YR = "wmap_V2_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_W1_3YR = "wmap_W1_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_W2_3YR = "wmap_W2_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_W3_3YR = "wmap_W3_ampl_bl_3yr_v2.txt";
string MSCS_WMAP_BEAMTF_FILE_W4_3YR = "wmap_W4_ampl_bl_3yr_v2.txt";

string MSCS_WMAP_BEAMTF_FILE_K1_5YR = "wmap_K1_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_K2_5YR = "wmap_Ka1_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_Q1_5YR = "wmap_Q1_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_Q2_5YR = "wmap_Q2_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_V1_5YR = "wmap_V1_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_V2_5YR = "wmap_V2_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_W1_5YR = "wmap_W1_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_W2_5YR = "wmap_W2_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_W3_5YR = "wmap_W3_ampl_bl_5yr_v3.txt";
string MSCS_WMAP_BEAMTF_FILE_W4_5YR = "wmap_W4_ampl_bl_5yr_v3.txt";

string MSCS_WMAP_BEAMTF_FILE_K1_7YR = "wmap_K1_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_K2_7YR = "wmap_Ka1_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_Q1_7YR = "wmap_Q1_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_Q2_7YR = "wmap_Q2_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_V1_7YR = "wmap_V1_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_V2_7YR = "wmap_V2_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_W1_7YR = "wmap_W1_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_W2_7YR = "wmap_W2_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_W3_7YR = "wmap_W3_ampl_bl_7yr_v4.txt";
string MSCS_WMAP_BEAMTF_FILE_W4_7YR = "wmap_W4_ampl_bl_7yr_v4.txt";


/************************************************/
/* masks file name definitions */
/************************************************/
string MSCS_WMAP_MASK_FILE_KP12_1YR = "map_kp12_mask_yr1_v1.fits";
string MSCS_WMAP_MASK_FILE_KP2_1YR = "map_kp2_mask_yr1_v1.fits";
string MSCS_WMAP_MASK_FILE_KP0_1YR = "map_kp0_mask_yr1_v1.fits";
string MSCS_WMAP_MASK_FILE_KQ75_5YR = "KQ75-wmap_ext_temperature_analysis_mask_r9_5yr_v3.fits";
string MSCS_WMAP_MASK_FILE_KQ85_5YR = "KQ85-wmap_temperature_analysis_mask_r9_5yr_v3.fits";


/************************************************/
/* WMAP fits file storage conventions */
/************************************************/

string MSCS_WMAP_MASK_FITS_FILE_COLUMN_NAME = "N_OBS";


