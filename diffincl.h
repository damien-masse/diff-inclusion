///////////////////////////////////////////////////////////////////////////
//  diffincl.h : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////

#ifndef __DIFFINCL_H__
#define __DIFFINCL_H__

#include <iostream>
#include <vector>
#include <codac.h>

using namespace codac;

namespace diffincl {

class DiffInclusion 
{
   public :
      DiffInclusion(unsigned int dim, const TFunction* tfun,
		    Interval& time, unsigned int nbsteps);

      DiffInclusion(unsigned int dim, const Function* fun,
		    Interval& time, unsigned int nbsteps);

      // reading from input 
      DiffInclusion(std::istream& input);

      int get_dim() const;
     
      // evaluate function on box + time (time dependant or not)
      const IntervalVector teval_vector(const Interval &tim,
                 const IntervalVector& cdom);
      const IntervalVector teval_vector(double tim,
                 const IntervalVector& cdom);
      const IntervalVector teval_vector(double tim,
                 const Vector& cdom);

      // compute jacobian ; in case of time dependency, tvec is the first 
      // column of the jacobian
      IntervalMatrix jacobian(const IntervalVector& codom,
                 const Interval& tdom,
                 IntervalVector& tvec);

      ///// basic enclosing of a time frame
      IntervalVector extend_box_basic
		(const IntervalVector& frame, const IntervalVector& Xbox,
                 const Interval& tim);

	// compute next step, from linear approximation and uncertainty
        // centzone : center of approxbox
	// fcent = f(center of the zone, half time)
	// tVect = mid df/dt
	// jac : mid jacobian
	// prodExpM : cumulated product of exponential (without time)
	// uncert : U value (uncertity on f)
	// double tsteps : time step
	// f(x,t) \in fcent + (t-tmid)*tvect + jac*(x-x0) + uncert
      void compute_next_step(const Vector& centzone,
			 const Vector& fcent,
                         const Vector& tVect,
                         const Matrix& jac,
                         const IntervalVector& uncert,
                         double tsteps,
                         IntervalVector& Cent,
                         IntervalMatrix& prodExpM,
                         IntervalVector& tauCent,
                         IntervalMatrix& tauProdExpM);


      ///// differential inclusion        X' = f(X)
      // frame : global frame (no computation outside)
      // X0 : initial position
      TubeVector fwd_inclusion(const IntervalVector& frame, 
	      	               const IntervalVector& X0);

    private:

      unsigned int dim; // number of variables (without time)
      const TFunction* tfun=NULL;
      const Function* fun=NULL;
      bool time_dependent;
          // if time_dependant is true, then the first variable is time
      Interval time;
      unsigned int nbsteps;

      IntervalMatrix* proj=NULL;
};

}


namespace diffincl {
      inline int DiffInclusion::get_dim() const { return this->dim; }
}

#endif
