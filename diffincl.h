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
      /**
       * constructor for non-autonomous function
       */
      DiffInclusion(unsigned int dim, const TFunction* tfun,
		    Interval& time, unsigned int nbsteps);

      /**
       * constructor for autonomous function
       */
      DiffInclusion(unsigned int dim, const Function* fun,
		    Interval& time, unsigned int nbsteps);

      /** 
       * construction from an input 
       * syntax :
       * dim:=<dimension>
       * timedependent:=<non-autonomous?>
       * time:=[<start>,<end>]
       * function:= (<v1>,<v2>...) : <function>
       * nbsteps:= <minsteps>
       */
      DiffInclusion(std::istream& input);

      /**
       * number of variables (without time)
       */
      int get_dim() const;
     
      /** evaluate function on box + time interval
       */
      const IntervalVector teval_vector(const Interval &tim,
                 const IntervalVector& cdom);
      /** evaluate function on box + time
       */
      const IntervalVector teval_vector(double tim,
                 const IntervalVector& cdom);
      /** evaluate function on point + time
       */
      const IntervalVector teval_vector(double tim,
                 const Vector& cdom);

      /** compute the jacobian ; 
       * in case of time dependency, tvec is the first 
       * column of the jacobian	(dfun / dt), the return is the rest
       */
      IntervalMatrix jacobian(const IntervalVector& codom,
                 const Interval& tdom,
                 IntervalVector& tvec);

      /** basic estimation of the set of states in the following time interval
       * two steps : 
       *  1) compute f1=fun(Xbox), f2=fun(Xbox+0.5*dt*f1), 
       *             f3=fun(Xbox+dt*(f1|f2)) and B = Xbox+[0,dt]*(f1|f2|f3)
       *  2) ``inflate'' B, intersect with frame, recompute
       *                 B = Xbox+[0,dt]*fun(B) and intersect with frame
       * @param frame : constraint on the set
       * @param Xbox : initial set of tates
       * @param tim : time interval
       * TODO : how to customize the policy?
       */
      IntervalVector extend_box_basic
		(const IntervalVector& frame, const IntervalVector& Xbox,
                 const Interval& tim, 
		 double inflation_factor = default_inflation_factor,
		 bool last_contraction=true);

	// compute next step, from linear approximation and uncertainty
	// fcent = f(center of the zone, half time)
	// tVect = mid df/dt
	// jac : mid jacobian
	// prodExpM : cumulated product of exponential (without time)
	// uncert : U value (uncertity on f)
	// double tsteps : time step
	// f(x,t) \in fcent + (t-tmid)*tvect + jac*(x-x0) + uncert
      void compute_next_step(const Vector& fcent,
                         const Vector& tVect,
                         const Matrix& jac,
                         const IntervalVector& uncert,
                         double tsteps,
                         IntervalVector& evolCenter,
                         IntervalMatrix& ExpM,
                         IntervalVector& tauCent,
                         IntervalMatrix& tauProdExpM);


      void fwd_step(const IntervalVector& constraint,
                const IntervalVector& Xbox,
                const IntervalVector& approxbox,
                const Vector& center_Approxbox,
                const Interval& timeslice,
                double tsteps,
                IntervalMatrix& ExpM,
                IntervalVector& evolCenter,
                IntervalVector& tauXbox);


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

      /** customizable fields */
      int debug_level = 3;

      /** inflation factor in extend_box_basic */
      constexpr static double default_inflation_factor = 1.2; 

      IntervalMatrix* proj=NULL;
};

}


namespace diffincl {
      inline int DiffInclusion::get_dim() const { return this->dim; }
}

#endif
