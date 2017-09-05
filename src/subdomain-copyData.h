/*!
 * \file subDomain.h
 *
 *  Created on: Jan 26, 2012
 *      Author: blew
 */

#ifndef SUBDOMAIN_H_
#define SUBDOMAIN_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include <vector>
#include <string>
//#include "Mscs-function3dregc.h"
#include "cpeds-point3d.h"
#include "cpeds-list.h"
#include "mscsVector.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

typedef struct {
	double xmin,xmax,ymin,ymax,zmin,zmax;
	long subx,suby,subz; // number of subdivisions on each of the dimensions per one subdivision level
} subDomain_region_t;

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates 
 \details 
 
 \date Jan 26, 2012, 4:53:03 PM 
 \author Bartosz Lew
 */

class subDomain {
	public:
//		static const int subregionsCount1d=2;
		typedef struct {
			// subDomain definition
			subDomain_region_t region;
			// structuring
			unsigned long minObjectsCount; //!< number of objects (particles) below which no more subdivisions are made
			subDomain **sub; //!< points to array of sub-domains
			int subDomainsCount; //!< number of sub-domains in this domain 1 level down
			vector<long>* idx; //!< pointer to an array with indexes of objects located within this subDomain
			bool root; //! true if the domain is the root domain; false otherwise
			subDomain* rootDomain; //!< pointer to the root subdomain
		} subDomainDef_t;
		
	public:
		//! this is the constructor to be used for making only the top domain.
		subDomain(const vector<cpedsPoint3D>* points=NULL, const subDomain_region_t *domain=NULL, long minPoints=1);
		//! sub-domain constructor (do not use it for the initialization of this class. It is used internally recursively)
		subDomain(const subDomain* parent, const subDomain_region_t &newSubDomain);
		virtual ~subDomain();
		void initiate(const vector<cpedsPoint3D>* points);

		subDomain(const subDomain& rhs);
		subDomain(const subDomain& rhs, subDomain* rootP);

//		void tree(const subDomain *domain=NULL);
		void tree();
		
		const subDomainDef_t& getSubdomainInfo() const { return _def; }
		
		void print_domain_range() { print_domain_range(_def.region); }
		void print_domain_range(subDomain_region_t r);
		/*!
			\brief saves the domain to a file for plotting
			\details 
			@param fileName - file name

			In each row a single domain is saved defined by two points coordinates (6 columns) in 7th column number of particles is saved.
			In the 8'th column 0 indicates that is isn't the root sub-domain. 1 - indicates the root sub-domain.
		
			\date Jan 27, 2012, 1:31:26 AM
			\author Bartosz Lew
		*/
		void saveDomains(string fileName);
		bool inDomain(const cpedsPoint3D& p) const;
		
		/*!
			\brief checks if the parallel-pipe shaped region centered at p and given side half-sizes is fully inside of this domain
			\details 
			@param p - center of the cubic region
			@param rx - half size of the side of the parallel-piped region along X
			@param ry - half size of the side of the parallel-piped region along Y
			@param rz - half size of the side of the parallel-piped region along Z
			@return true if the region is fully inside of the domain, false - otherwise
		
			\date Jan 27, 2012, 11:13:07 AM
			\author Bartosz Lew
		*/
		bool inDomain(const cpedsPoint3D& p, double rx,double ry,double rz) const;
		bool inDomain(const subDomain_region_t& r, const subDomain* D=NULL) const;
		/*!
			\brief checks if the point p is inside of the region r
			\details 
			@param p - point to check
			@param r - definition of the region
			@return returns true in p is inside of r and false otherwise
		
			\date Jan 27, 2012, 10:19:45 PM
			\author Bartosz Lew
		*/
		bool inRegion(const subDomain_region_t& r,const cpedsPoint3D& p);
		bool inRange(double r,const cpedsPoint3D& p0, const cpedsPoint3D& p);
		
		/*!
			\brief get a list of Nmin to Nmax particles from domain D that lie within sphere centered at p0
			\details 
			@param p0 - center of the spherical region to search particles around
			@param Nmin - minimal number of particles inside of the spherical region
			@param Nmax - maximal number of particles inside of the spherical region
			@param hsml - pointer to the smoothing length - the radius containing the found particles
			@param acc - accuracy parameter defining how many iterations will be possible in trying to set the smoothing length
			by means of divisinon by two search, in order to find the one that would grant number of particles within the sphere around 
			p0 from within the range [Nmin,Nmax]. The number of iterations is defined by rMax/acc where rMax the diagonal of the domain size.
			@return list of Nmin to Nmax particles.
			
			If the particles cound not have been found within the domain or its subdomains then 
		
			\date Jan 30, 2012, 9:37:37 AM
			\author Bartosz Lew
		*/
		mscsVector<long> getNeighbors(cpedsPoint3D p0, long Nmin, double* hsml=NULL, double acc=1e-5);

		/*!
			\brief returns the indexes of points located within the requested parallel-pipe region
			\details 
			@param r - region definition
			@return list of indexes of the points within that region
		
			\date Jan 27, 2012, 5:40:32 PM
			\author Bartosz Lew
		*/
//		cpedsList<long> getContainingDomainPointIdx(const subDomain_region_t& r);
//		vector<long>* getContainingDomainPointIdx(const subDomain_region_t& r, int level=0);
//		vector<long>* getContainingDomainPointIdx(const subDomain_region_t& r, subDomain** sd=NULL, bool ignoreRegionSize=false);
		vector<long>* getContainingDomainPointIdx(const subDomain_region_t& r, subDomain** sd=NULL);
		cpedsList<long> getPointsIdx(const subDomain_region_t& r);
		/*!
			\brief return the list of indexes of points located within the requested distance from the point p0 in this or in the pointed subDomain
			\details 
			@param p0 - center of the spherical region to search points in
			@param r - radii of the region
			@param D - subdomain to seach within; if null then this domain is searched; if not null then the D subDomain is used first and if 
			the region is too big then this subDomain is searched instead 
			@return list of indexes of the points within specified region
			
			If D=NULL then the ignoreRegionSize parameter is not used. Then first the smallest containing domain is found and of those points
			in the domain there are selected only these points that match the criteria of being inside of the sphere S(p0,r).
			
			If D!=NULL then the ignoreRegionSize parameter is used. If it is false (default) then first the smallest containing domain is found starting
			with the domain D and of those points in the domain there are selected only these points that match the criteria of being inside of the 
			sphere S(p0,r). if the ignoreRegionSize==true then the D domain is used and its points are tested against the criteria of
			being inside of the sphere S(p0,r).
		
			\date Jan 30, 2012, 9:23:46 AM
			\author Bartosz Lew
		*/
		mscsVector<long> getPointsIdx(const cpedsPoint3D& p0, double r, subDomain* D=NULL, bool ignoreRegionSize=false);
		
	protected: 
		
		void saveDomain(FILE** f);
		/*!
			\brief finds sub-set of points belonging to this domain from the set of points contained in the parent domain
			\details 
			@param parent - pointer to the parent subDomain 

			The function generates an index of points that are located in this domain selecting 
			from the list of points defined (listed) in the parent subDomain.
		
			\date Jan 26, 2012, 11:22:36 PM
			\author Bartosz Lew
		*/
		void findPointsInDomain(const subDomain* parent);
		
		
		subDomainDef_t _def; //!< definition of this subDomain
//		const vector<cpedsPoint3D>* _points; //!< scattered data points
		vector<cpedsPoint3D>* _points; //!< a copy of the scattered data points
		
//		const subDomain_region_t* _topDomainRegion; //! definition of the domain to be split into subDomains
};


#endif /* SUBDOMAIN_H_ */
