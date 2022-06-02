///////////////////////////////////////////////////////////////////////
//  IVparals.h : representation of polyhedron as C + inter M_i X_i
///////////////////////////////////////////////////////////////////////

#ifndef __IVPARAL_H__
#define __IVPARAL_H__

#include <iostream>
#include <vector>
#include <codac.h>

using namespace codac;

namespace diffincl {

/** basic class, with one representation
 */
class IVparals
{
   public :
      /**
       * constructor from dimension (empty domain, ID matrix)
       */
      IVparals(int dim);

      /**
       * constructor from an initial box (id matrix)
       */
      IVparals(const IntervalVector& iv);

      /**
       * constructor from an initial box and matrices from an existing IVparal
       */
      IVparals(const IVparals& iv, const IntervalVector& box);

      /**
       * copy
       */
      IVparals(const IVparals& iv);

      /**
       * from a slanted box  ( M.V )
       */
      IVparals(const IntervalMatrix& M, const IntervalMatrix &rM, 
			const IntervalVector &V);

      /**
       * number of variables
       */
      int get_dim() const;

      /**
       * number of matrices
       */
      int get_nbmat() const;

      /**
       * ``basic'' simplification
       */
      void simplify(bool bwd);

      /**
       * ``strong'' simplification (if possible) 
       */
      void strong_simplify();


      /**
       * recenter 
       */
      void recenter(const Vector& nCent);


      /**
       * ``middle'' (e.g. center of the bounding box ?
       */
      Vector mid() const;

      /**
       * union (keep the actual matrices)
       */
      void join_with(const IVparals& iv);
     
      /**
       * intersection (keeping the actual matrices)
       */
      void intersect_with(const IntervalVector& iv);
      void intersect_with(/* const */IVparals& iv);

       /**
        * intersection with quasi-linear constraints [M][X]=[Y] 
        */
      void intersect_with(const IntervalMatrix& M,
		const IntervalVector& Y);

       /**
        * intersection with quasi-linear constraints [M][X-c]=[Y] 
        */
      void intersect_with(const IntervalMatrix& M,
		const Vector& c,
		const IntervalVector& Y);

      /**
       * union with an intersection of an IVparals and the border of
       * a box (keep the matrices)
       * returns true if modified
       */
      bool join_intersect_with
		(const IVparals& iv, 
		 const IntervalVector& box, int d, double val);

      /**
       * relative distance, fast (with same base) 
       * max of the distance of each box
       */
       double rel_distance_fast(const IVparals& iv) const;
 
      /**
       * union with 
			box[d->val] 
				intersected with 
			iv + ]0,1]*(M (iv-center)+V) */
      bool join_intersect_with_tau
		(const IVparals& iv, const Vector&center,
                const IntervalMatrix& M,
                const IntervalVector& V,
		 const IntervalVector& box, int d, double val);
      /**
       * intersection with the bound of an interval vector
       */
      bool intersect_with(/* const */ IntervalVector& iv, int d, double val);

      /** is_empty
       */
      bool is_empty() const;

      /** bounding box
       */
      IntervalVector bounding_box() const;

      /** contains point
       * (approximate...)
       */
      bool contains(const Vector& iv) const;

      /** is subset, **with the same base** */
      bool is_subset_fast(const IVparals& iv) const;

      /** expansion form a bigger set, **with the same base** */
      /* this |= this+(iv-this)*fact */
      void inflate_from_base_fast(const IVparals& iv, double fact);

      /** multiply by a new matrix, and add 
       *  M : the new matrix, IM its inverse
       */
      void mult_and_add(const IntervalMatrix& M, const IntervalMatrix &invM,
 			const IntervalVector& V);

      /** centered multiply and add 
       */
      void cmult_and_add(const Vector&center,
		const IntervalMatrix& M, const IntervalMatrix &invM,
		const IntervalVector& V);

      /** centered multiply and add 2
       */
      void cmult_and_add2(const Vector&center,
		const IntervalMatrix& M, const IntervalMatrix &invM,
		const IntervalVector& V);

      /** cone-addition, with retur (keep matrices and center) */
      friend IVparals tau_add(const IVparals& iv, const IntervalVector& V);

      /** tau-centered multiply and add (not modifying the matrices)
       * IV = IV + [0,1]*(M (IV-center) + V)
       */
      void ctau_mult_and_add(const Vector&center,
		const IntervalMatrix& M,
		const IntervalVector& V);

      /** centered multiply and add3
       *  D' = c + M(D-c + mA*V) + V2 (V centered on 0)
       *  (mA', x0', vB') = 
                (M*mA,x0 + V, M*(vB-c) + c + V2)
       */
      void cmult_and_add3(const Vector& center,
                const IntervalMatrix& M,
                const IntervalMatrix& invM,
                const IntervalVector& V,
                const IntervalVector& V2);

      /** centered multiply and add4
       *  D' = c + M(D-c) + V2 + (M mA V)inter(Vb) (with V centered on 0)
       *  (mA', x0', vB') = 
                (M*mA,x0 + Vm, M*(vB-c) + c + V2)
       */
      void cmult_and_add4(const Vector& center,
                const IntervalMatrix& M,
                const IntervalMatrix& invM,
                const IntervalVector& V,
                const IntervalVector& V2,
                const IntervalVector& Vb);

      /** generate a ConvexPolygon which is an overapproximation of the
       * projection of the polyhedron (basic form)
       */
      ConvexPolygon over_polygon(const Matrix& M) const;

      /** operations */
      friend IVparals operator+(const IVparals& ivp, 
				const IntervalVector& iv);
      friend IVparals operator-(const IVparals& ivp, 
				const IntervalVector& iv);
      friend IVparals operator+(const IVparals& iv, const Vector& v);
      friend IVparals operator-(const IVparals& iv, const Vector& v);
      friend IntervalVector operator*(const IntervalMatrix& M, 
				const IVparals& iv);

      /** display
       */
      friend std::ostream& operator<<(std::ostream& str, const IVparals& iv);


    private:

      unsigned int dim; // number of variables
      unsigned int nbmat; // number of intersections, without BBox
      unsigned int nbNcst; // number of ``constant'' matrices
      bool isempty;
      Vector center; // ``center''
      std::vector<IntervalMatrix> mats; // nbmat matrix (first the variables)
      std::vector<IntervalMatrix> Imats; // nbmat matrix
      std::vector<IntervalVector> rhs; // nbmat+1 rhs, the last one is BBox

      /** customizable fields */
      int debug_level = 3;

};


}


namespace diffincl {
      inline int IVparals::get_dim() const { return this->dim; }
      inline int IVparals::get_nbmat() const { return this->nbmat; }
      inline bool IVparals::is_empty() const { return this->isempty; }
}

#endif
