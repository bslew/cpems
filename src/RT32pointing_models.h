/*!
  \file pointing_models.h - contains procedures that calculate RT-32 pointing models
*/


#ifndef RT32_POINTING_MODELS_H_
#define RT32_POINTING_MODELS_H_

/*!
	\brief implements Model4c as implemented in the fast_track control system and earlier versions
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	\date Nov 3, 2016, 9:58:35 AM
	\author Bartosz Lew
*/
void cpeds_RT32_Model4(double AZ, double ZD, double *dAZ, double *dZD);

/*!
	\brief implements Model4e (see 2016/2 tech rep. for details)
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	\date Nov 3, 2016, 10:14:09 AM
	\author Bartosz Lew
*/
void cpeds_RT32_Model4e(double AZ, double ZD, double *dAZ, double *dZD);

/*!
	\brief implements Model5 as implemented in the COCONUT version of the control system
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	(see 2016/2 tech rep. for details)
	\date Nov 3, 2016, 10:14:09 AM
	\author Bartosz Lew
*/
void cpeds_RT32_Model5(double AZ, double ZD, double *dAZ, double *dZD);

/* ******************************************************************************************** */
/*!
	\brief implements Model4g (see 2017/1 tech rep. for details)
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	Temperature independent model.
	Although this model is mathematically the same as Model4e for AZ corrections (for the same data)
	and similar to Model4e for ZD corrections, the processing of the pointing data to 
	derive this model is different in processing the weather dependent refraction. 
	Now the full history of the refraction angles is used rather than optical refraction.
	
	The model is calibrated against DR2017sum pointing data. See 2017/1 tech rep. for details.

	\date Mar 6, 2018, 7:16:07 PM
	\author Bartosz Lew
*/
void Model4g_DR2017sum(double AZ, double ZD, double* dAZ, double *dZD);

/*!
	\brief implements Model4gr 
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	The model adds rail track corrections atop model 4g: 
	hence the name has "r" - according to the new naming convention
	
	Calibrated against DR2017sum pointing data. See 2017/1 tech rep. for details.

	\date Mar 6, 2018, 7:24:13 PM
	\author Bartosz Lew
*/
void Model4gr_DR2017sum(double AZ, double ZD, double *dAZ, double *dZD);

/***************************************************************************************/
/*!
	\brief implements Model6a (see 2017/1 tech rep. for details)
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction
	@param Tmeteo [degC] - temperature from weather station IRDAM WST-7000

	Temperature dependent model for ZD corrections.
	This model is mathematically the same as Model4g for AZ corrections.
	
	The model is calibrated against DR2017 pointing data. See 2017/1 tech rep. for details.

	\date Mar 6, 2018, 7:16:07 PM
	\author Bartosz Lew
*/
void Model6a_DR2017(double AZ, double ZD, double *dAZ, double *dZD, double Tmeteo);

/* ******************************************************************************************** */
/*!
	\brief implements Model4gr - rail track for model 4g (hence the name has "r"- new naming convention)
	\details 
	@param AZ [deg] - from South
	@param ZD [deg]
	@param dAZ [deg] - model position correction
	@param dZD [deg] - model position correction

	Calibrated against DR2017sum pointing data. See 2017/1 tech rep. for details.

	\date Mar 6, 2018, 7:24:13 PM
	\author Bartosz Lew
*/
void Model6ar_DR2017(double AZ, double ZD, double *dAZ, double *dZD, double Tmeteo);



#endif /* RT32_POINTING_MODELS_H_ */ 



