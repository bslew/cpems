/*!
 \file Cthwisdom.h - 
 */

#ifndef SRC_CTHWISDOM_H_
#define SRC_CTHWISDOM_H_

#include <vector>
#include <string>
#include <map>
#include "cpeds-list.h"
#include <unordered_map>


/*!
 \class Cthwisdom
 \brief Encapsulates structures to keep information on how to calculate Cth fast
 \details 
 	 
 
 
 \date created: Apr 29, 2019, 2:53:02 PM 
 \author Bartosz Lew
 */

class Cthwisdom {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef std::unordered_map<long, std::vector<long> > pairs_t;
		typedef struct {
//			long i=0;
//			long mask=0;
//			std::vector<long> j;
			pairs_t pairs;
		} element_t;
		typedef struct {
			std::vector< element_t > idx;
			std::vector<double> ang;
			std::vector<double> hits;
			double min,max,res,bin,halfbin;
			long Npix;
		} Cthwisdom_t;
		
		double use_wisdom_fraction;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		/*!
			\brief default constructor
			\details 
			@param theta_min - 
			@param theta_max - 
			@param step - 
			@param bin - full bin width

			The wisdom for the correlation function will be defined for the following angles:
			theta_i = theta_min+i*step, where
			i=0,1,2,...nmax, where
			theta_min+nmax*step is the maximal angle that is smaller or equal than theta_max.
			
			Hence the size of the wisdom is given by
			
			size = (theta_max-theta_min)//step+1 (where // is the integer division operator)
			
			Angles between pairs of directions are attributed to a bin according to:
			
			[theta_i-bin/2, theta_i+bin/2).
			
			The bin size should be small enough to assure that neighboring bins do not overlap.
		
			\date May 9, 2019, 1:51:16 PM
		*/
		Cthwisdom(double theta_min=1, double theta_max=180, double step=1, double bin=1);
		
		Cthwisdom(cpedsList<double> th, cpedsList<double> bin);
//		Cthwisdom(cpedsList<double> th=cpedsList<double>(), cpedsList<double> bin=cpedsList<double>());
		~Cthwisdom();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		void initialize(cpedsList<double> th, cpedsList<double> bin);
		
		/* ---------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------- */
		
		/*!
			\brief save wisdom to file
			\details 
			@param fname - wisdom file name. If not given, default name will be created for this wisdom 
			object.
			
			\date May 7, 2019, 4:43:36 PM
		*/
		void save(std::string fname="");
		void saveHDF(std::string fname="");
		
		void setWisdomFraction(double frac) { use_wisdom_fraction=frac; }
		
		/*!
			\brief load wisdom from file
			\details 
			@param fname - wisdom file name. If not given, default name will be used for this wisdom 
			object.
			
			@return returns 0 if loaded correctly, and something else otherwise
			\date May 7, 2019, 4:43:36 PM
		*/
		int load(std::string fname="");
		long size() { return wisdom.ang.size(); }
		long get_size(double theta_min, double theta_max, double step);
		void add_wisdom(double angle, long i, long j);
		
		/*!
			\brief perform wisdom calibration and optimization actions 
			\details 

			This should be called before saving to file.
		
			\date May 14, 2019, 11:59:22 AM
		*/
		void compile();
		
		std::string get_wisdom_file_name();
		
		Cthwisdom::element_t& operator[](long i) { return wisdom.idx[i]; }
		double& operator()(long i) { return wisdom.ang[i]; }
		
	protected:
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		long get_ang_idx(double angle);

		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */
		
		Cthwisdom_t wisdom;
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */

		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		Cthwisdom_t _Cthwisdom_data;
		
};
#endif /* SRC_CTHWISDOM_H_ */ 

