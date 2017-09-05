//! This class defines the containers and procedures for Monte Carlo Marcov Chain based multidimentional likelihood maximization in a pre-defined
//! parameter space. This class allows for starting multiple chains simultaneously - implemented via MPI interface.
//! It  also allows for prior definition and likelihood or posterior integration using previously interpolated likelihood surface based on the 
//! information probed by the  MC Markov Chains. This class also allows for marginalization over requested parameters after the interpolation was done.

/*! The interpolation is done in multidimentional space using spline algorithm. The class stores the whole information about the trace of the 
MCMC also the stetes in which the likelihood is lower than the current state. This increases the amount of information needed for posterior integration.
The class includes the lists  from the QT library as these lists allow for management of lists consisting of objects of arbitrary type. 

The default implementation is the parallel implementation where the parallelization is done such that each MCMC dwells on separate core. Hence if
eg. 8 chains are requested then the program will automatically run 8 processes in parallel.. This probably makes sense on multi-processor machines.
If the physical number of cores is less then the number of requested MCMCs then of course there will be some task sharring on the same core. 

The parameter space is predefined but also the user defined configuration can be read in. 
The information about the user defined parameter space can be introduced to the program by including the Mscs-MCMC-parameter_space.h file
in which the apropriate definitions are given. These definitions are then read into objects of this class at the run-time level.
The definitions must include the requested number of MCMCs 
*/
 

#include <qptrlist.h> 
#include "Mscs-common.h"


class Mscs-MCMC {

 public:
  
  // definition of the parameter space


 protedted:


 private:


};

