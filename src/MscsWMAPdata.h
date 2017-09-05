/*!
 * \file MscsWMAPdata.h - defines the file names of the WMAP data
 *
 *  Created on: Apr 27, 2010
 *      Author: blew
 */
#include "Mscs-common.h"
#include "Mscs-global-defs.h"

#ifndef MSCS_WMAP_DATA
#define MSCS_WMAP_DATA

/****************************************************************** */
/* declaration of variables that hold location of various WMAP data */
/****************************************************************** */
/*
 * these variables define the location of the WMAP data files such as:
 * - beam transfer function
 * - sky masks
 * - Nobs
 * - etc...
 *
 * These data are assumed to be held in directory
 * MSCS_WMAP_DATA_DIR
 *
 */


extern string MSCS_WMAP_NOBS_FILE_K1_1YR; //!< NAME OF THE FILE HOLDING N_OBS this should be fits file from the yearly data release maps per DA
extern string MSCS_WMAP_NOBS_FILE_K2_1YR;
extern string MSCS_WMAP_NOBS_FILE_Q1_1YR;
extern string MSCS_WMAP_NOBS_FILE_Q2_1YR;
extern string MSCS_WMAP_NOBS_FILE_V1_1YR;
extern string MSCS_WMAP_NOBS_FILE_V2_1YR;
extern string MSCS_WMAP_NOBS_FILE_W1_1YR;
extern string MSCS_WMAP_NOBS_FILE_W2_1YR;
extern string MSCS_WMAP_NOBS_FILE_W3_1YR;
extern string MSCS_WMAP_NOBS_FILE_W4_1YR;

extern string MSCS_WMAP_NOBS_FILE_K1_3YR; //!< NAME OF THE FILE HOLDING N_OBS this should be fits file from the yearly data release maps per DA
extern string MSCS_WMAP_NOBS_FILE_K2_3YR;
extern string MSCS_WMAP_NOBS_FILE_Q1_3YR;
extern string MSCS_WMAP_NOBS_FILE_Q2_3YR;
extern string MSCS_WMAP_NOBS_FILE_V1_3YR;
extern string MSCS_WMAP_NOBS_FILE_V2_3YR;
extern string MSCS_WMAP_NOBS_FILE_W1_3YR;
extern string MSCS_WMAP_NOBS_FILE_W2_3YR;
extern string MSCS_WMAP_NOBS_FILE_W3_3YR;
extern string MSCS_WMAP_NOBS_FILE_W4_3YR;

extern string MSCS_WMAP_NOBS_FILE_K1_5YR; //!< NAME OF THE FILE HOLDING N_OBS this should be fits file from the yearly data release maps per DA
extern string MSCS_WMAP_NOBS_FILE_K2_5YR;
extern string MSCS_WMAP_NOBS_FILE_Q1_5YR;
extern string MSCS_WMAP_NOBS_FILE_Q2_5YR;
extern string MSCS_WMAP_NOBS_FILE_V1_5YR;
extern string MSCS_WMAP_NOBS_FILE_V2_5YR;
extern string MSCS_WMAP_NOBS_FILE_W1_5YR;
extern string MSCS_WMAP_NOBS_FILE_W2_5YR;
extern string MSCS_WMAP_NOBS_FILE_W3_5YR;
extern string MSCS_WMAP_NOBS_FILE_W4_5YR;

extern string MSCS_WMAP_NOBS_FILE_K1_7YR; //!< NAME OF THE FILE HOLDING N_OBS this should be fits file from the yearly data release maps per DA
extern string MSCS_WMAP_NOBS_FILE_K2_7YR;
extern string MSCS_WMAP_NOBS_FILE_Q1_7YR;
extern string MSCS_WMAP_NOBS_FILE_Q2_7YR;
extern string MSCS_WMAP_NOBS_FILE_V1_7YR;
extern string MSCS_WMAP_NOBS_FILE_V2_7YR;
extern string MSCS_WMAP_NOBS_FILE_W1_7YR;
extern string MSCS_WMAP_NOBS_FILE_W2_7YR;
extern string MSCS_WMAP_NOBS_FILE_W3_7YR;
extern string MSCS_WMAP_NOBS_FILE_W4_7YR;


/************************************************/
/* beam transfer function file name definitions */
/************************************************/
extern string MSCS_WMAP_BEAMTF_FILE_K1_1YR; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
extern string MSCS_WMAP_BEAMTF_FILE_K2_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q1_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q2_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_V1_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_V2_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_W1_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_W2_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_W3_1YR;
extern string MSCS_WMAP_BEAMTF_FILE_W4_1YR;

extern string MSCS_WMAP_BEAMTF_FILE_K1_3YR; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
extern string MSCS_WMAP_BEAMTF_FILE_K2_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q1_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q2_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_V1_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_V2_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_W1_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_W2_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_W3_3YR;
extern string MSCS_WMAP_BEAMTF_FILE_W4_3YR;

extern string MSCS_WMAP_BEAMTF_FILE_K1_5YR; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
extern string MSCS_WMAP_BEAMTF_FILE_K2_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q1_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q2_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_V1_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_V2_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_W1_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_W2_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_W3_5YR;
extern string MSCS_WMAP_BEAMTF_FILE_W4_5YR;

extern string MSCS_WMAP_BEAMTF_FILE_K1_7YR; //!< NAME OF THE FILE HOLDING beam transfer function this should be a txt file from the yearly data release maps per DA
extern string MSCS_WMAP_BEAMTF_FILE_K2_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q1_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_Q2_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_V1_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_V2_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_W1_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_W2_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_W3_7YR;
extern string MSCS_WMAP_BEAMTF_FILE_W4_7YR;

/************************************************/
/* masks file name definitions */
/************************************************/
extern string MSCS_WMAP_MASK_FILE_KP12_3YR;
extern string MSCS_WMAP_MASK_FILE_KP2_3YR;
extern string MSCS_WMAP_MASK_FILE_KP0_3YR;
extern string MSCS_WMAP_MASK_FILE_KQ75_5YR;
extern string MSCS_WMAP_MASK_FILE_KQ85_5YR;

extern string MSCS_WMAP_MASK_FITS_FILE_COLUMN_NAME; //!< this defines the name of the column that holds the actual mask data in the fits file.


#endif /* MSCS_WMAP_DATA */
