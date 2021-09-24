///////////////////////////////////////////////////////////////////////////
//  IVdouble.h : representation of interval vectors as A V_0+B (V_0 centered)
///////////////////////////////////////////////////////////////////////////

#ifndef __DOUBLEREP_H__
#define __DOUBLEREP_H__

#include <iostream>
#include <vector>
#include <codac.h>

using namespace codac;

namespace diffincl {

/** basic class, with one representation
 */
class IVdouble
{
   public :
      /**
       * constructor from initial intervalVector
       */
      IVdouble(const IntervalVector& iv);

      /**
       * copy
       */
      IVdouble(const IVdouble& iv);

      /**
       * number of variables
       */
      int get_dim() const;
     
      /** bounding box
       */
      IntervalVector bounding_box() const;

      /** multiply by a new matrix, and add 
       *  D' = MD + V
       *  (mA',vB') = (M*mA , M*vB + V)
       */
      void mult_and_add(const IntervalMatrix& M, const IntervalVector& V);

      /** centered multiply and add 
       *  D' = c + M(D-c) + V
       *  (mA', vB') = (M*mA, M*(vB-c) + c + V)
       */
      void cmult_and_add(const Vector& center,
		const IntervalMatrix& M, const IntervalVector& V);

      /** ``simplification''
       *  for each ``quasi-straight'' column vector of mA
       *    associated to a significant coefficient in vB, include
       * the coefficient in x0
       * threshold :     max(OtherC.mag) / (Pivot.mig * vB.c) < threshold
       */
      void simplify(double threshold1=1e-4, double threshold2=1e-3);

    private:

      unsigned int dim; // number of variables
      IntervalMatrix mA;
      IntervalVector x0;
      IntervalVector vB;

      /** customizable fields */
      int debug_level = 3;

};

}


namespace diffincl {
      inline int IVdouble::get_dim() const { return this->dim; }
}

#endif
