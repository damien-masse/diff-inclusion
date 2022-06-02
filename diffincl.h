///////////////////////////////////////////////////////////////////////////
//  diffincl.h : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////

#ifndef __DIFFINCL_H__
#define __DIFFINCL_H__

#include <iostream>
#include <vector>
#include <codac.h>
#include <codac-capd.h>
#include "IVparals.h"

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

      /**
       * predefined time interval
       */
      const Interval& get_time() const;
     
      /** evaluate function on box + time interval
       */
      const IntervalVector teval_vector(const Interval &tim,
                 const IntervalVector& cdom);
      const IntervalVector teval_vector(const Interval &tim,
                 const IVparals& cdom);
      /** evaluate function on box + time
       */
      const IntervalVector teval_vector(double tim,
                 const IntervalVector& cdom);
      const IntervalVector teval_vector(double tim,
                 const IVparals& cdom);
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
      IntervalMatrix jacobian(const IVparals& codom,
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
       * @inflation_factor : inflation of the box
       * @nb_tries : number of tries : if >0, extends even more...
       * TODO : how to customize the policy?
       */
      IVparals extend_box_basic
		(const IntervalVector& frame, const IVparals& startIV,
                 const Interval& tim, 
		 double inflation_factor = default_inflation_factor,
		 int nb_tries = 0);

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
                         const IntervalMatrix* coord,
                         const IntervalMatrix* invcoord,
                         double tsteps,
                         IntervalVector& evolCenter,
                         IntervalMatrix& ExpM,
                         IntervalMatrix& invExpM,
                         IntervalVector& tauCent,
                         IntervalMatrix& tauProdExpM);


      /** try one step of computation, may fail
       *  constraint : external frame (cannot grow bigger)
       *  actState : actual point 
       *  tauState : actual approxBox (may be modified)
       *  finState : final states (to be computed)
       *  tsteps : elapsed time
       *  outDoors : to compute exits from the box
       *  returns false : tauState fails, need a new one. finState not modified
       *  returns true :  success, tauState minimized, finState modified
       */ 
      bool compute_step(const IntervalVector& frame,
			const IVparals& actState,
			IVparals& tauState,
			IVparals& finState,
			const Interval& timeslice,
			vector<IVparals> *outDoors);


      void fwd_step(const IntervalVector& constraint,
                const IntervalVector& approxbox,
                const Vector& center_Approxbox,
                const Interval& timeslice,
                const IntervalMatrix* coord,
                const IntervalMatrix* invcoord,
                double tsteps,
                IntervalMatrix& ExpM,
                IntervalMatrix& invExpM,
                IntervalMatrix& tauExpM,
                IntervalVector& evolCenter,
                IntervalVector& tauCent);


      ///// differential inclusion        X' = f(X)
      // frame : global frame (no computation outside)
      // X0 : initial position
      TubeVector fwd_inclusion(const IntervalVector& frame, 
	      	               const IntervalVector& X0,
			       VIBesFig *fig=NULL,
                               const Matrix *pr=NULL,
			       int nbsl=0);

      ///// differential inclusion        X' = f(X)
      // frame : global frame (no computation outside)
      // already defined tube
      // and precise time computation
      IVparals fwd_inclusion(const IntervalVector& frame, 
                               const IVparals& X0,
			       TubeVector &Result,
			       const Interval &ptime,
			       unsigned int pnbsteps,
                               VIBesFig *fig=NULL,
                               const Matrix *pr=NULL,
                               int nbsl=0);

      ///// CAPD differential inclusion        X' = f(X)
      // frame : global frame (no computation outside)
      // X0 : initial position
      TubeVector capd_fwd_inclusion(const IntervalVector& frame, 
	      	               const IntervalVector& X0);

      ///// ``room'' traversal (generating doors from initial input door,
      // time interval, time steps...
      // X0 : initial position (input door)
      // intTime : time interval to compute the exits (e.g. not starting
      //      from 0 => the initial door will not appear as an output door)
      // timestep : time steps
      // outDoors : output (2*dim doors)
      // timeshift (optional) : for time-dependant function, the value of t 
      // corresponding to the initial position
      vector<IVparals> room_fwd_computation(const IntervalVector& box,
			const IVparals& X0,
			const Interval& intTime,
			const double timesteps,
			double timeShift = 0.);


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
      constexpr static double default_inflation_factor = 1.01; 
      constexpr static int narrow_factor = 3;

      IntervalMatrix* proj=NULL;
};

}


namespace diffincl {
      inline int DiffInclusion::get_dim() const { return this->dim; }
      inline const Interval& DiffInclusion::get_time() const { return this->time; }
}

#endif
