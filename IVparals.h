///////////////////////////////////////////////////////////////////
//  IVparals.h : representation of polyhedron as inter M_i X_i
///////////////////////////////////////////////////////////////////

#ifndef __IVPARAL_H__
#define __IVPARAL_H__

#include <iostream>
#include <vector>
#include <codac.h>

using namespace codac;

namespace diffincl {

class IVparals;

  IVparals sum_tau(const IVparals& iv, const IntervalVector& V,
				 bool keep=false);

class IVparals
{
   public :

      /******************** CONSTRUCTORS *******************/
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
       * from a slanted box  ( M.V, with reverse known )
       */
      IVparals(const IntervalMatrix& M, const IntervalMatrix &rM, 
			const IntervalVector &V);
      /**
       * from a slanted box  ( M.V, reverse computed )
       */
      IVparals(const IntervalMatrix& M, const IntervalVector &V);

      /***************** ACCESS **************************/
      /**
       * number of variables
       */
      int get_dim() const;
      /**
       * number of matrices
       */
      int get_nbmat() const;
      /** is_empty
       */
      bool is_empty() const;
      /** bounding box
       */
      IntervalVector bounding_box() const;
      /**
       * ``middle'' (e.g. center of the bounding box ?
       */
      Vector mid() const;
      /** contains point (approximate...)
       */
      bool contains(const Vector& iv) const;
      /** get i-th matrice (i<nbmat) */
      const IntervalMatrix& getMat(int i) const;
      /** get i-th vector, not modified */
      const IntervalVector& getVec(int i) const;
      /** get (M*Mi) Vi */
      IntervalVector getPar(const IntervalMatrix& M, int i) const;
      /** get inter_i (M*Mi) Vi */
      IntervalVector getPar(const IntervalMatrix& M) const;
      /**
       * relative distance, fast (with same base) 
       * max of the distance of each box
       */
       double rel_distanceFast(const IVparals& iv) const;
 

      /************* Modification **********************/
      /** empty */
      void set_empty();
      /** return to 0, not modifying the matrices */
      void clear();
      /** inflation by a cube
       *  this <- this + [-r,r]^d
       *  keep the matrices
       *  return this
       */
      IVparals& inflate(double rad);
      /** inflation by a ball
       *  this <- this + ball(rad)
       *  keep the matrices
       *  return this
       */
      IVparals& inflateBall(double rad);
      /** centered homothety
       *  x <- [c] + delta*([x]-[c]) ( (1-delta)[c] + delta*[x] )
       *  keep the matrices
       *  return this
       */
      IVparals& homothety(IntervalVector c, double delta);
      /** set this to x
       *  keep the matrices
       */
      IVparals& operator=(const IntervalVector& x);
      /** "project" x on this, get an overapproximation
       *  keep the matrices
       */
      IVparals& assign(const IVparals& iv);

      /** expansion form a bigger set, **with the same base** */
      /* this |= this+(iv-this)*fact */
      void inflate_from_baseFast(const IVparals& iv, double fact);


      /** apply a (non-singular) endomorphism ;
          change the matrices.
          this <- M*this ;    IM : inverse of M */
      IVparals& linMult(const IntervalMatrix& M, const IntervalMatrix& IM);

      /***** Intersections *****/

      /** intersection with a box, keep the matrices
       */
      IVparals& operator&=(const IntervalVector& x);
      /** intersection with a box, new value */
      friend IVparals operator&(const IVparals& iv, const IntervalVector& x);
      /** intersection with IVparals, equals matrices (fast) */
      IVparals& meetFast(const IVparals& iv);
      /** intersection with IVparals, keep matrices for now */
      IVparals& meet(const IVparals& iv);
      /** intersection with a interval linear constraint,
       *      M x \in b 
       *      keep = false => modify a constraint in last matrice */
      bool meetLN(const IntervalVector& V, const Interval& b, bool keep);
      /** intersection with a centered interval linear constraint,
       *      M (x-c) \in b 
       *      keep = false => modify a contraint in last matrice */
      bool meetLN(const IntervalVector& V, const IntervalVector& C,
				const Interval& b, bool keep);
      /** with a set of linear constraints */
      bool meetLM(const IntervalMatrix& S, const IntervalVector& b, bool keep);
      /** with a centered set of linear constraints */
      bool meetLM(const IntervalMatrix& S, const IntervalVector& C,
				const IntervalVector& b, bool keep);
      /** with a centered set of linear constraints */
      bool meetLM(const IntervalMatrix& S, const Vector& C,
				const IntervalVector& b, bool keep);
      
      /** union with a box 
       */
      IVparals& operator|=(const IntervalVector& x);
      friend IVparals operator|(const IVparals& iv, const IntervalVector& x);

      /***** operations  *****/

      /** sum, difference */
      IVparals& operator+=(const IntervalVector& V);
      IVparals& operator-=(const IntervalVector& V);
      friend IVparals operator+(const IVparals& iv, const IntervalVector& V);
      friend IVparals operator-(const IVparals& iv, const IntervalVector& V);
      friend IVparals sum_tau(const IVparals& iv, const IntervalVector& V,
						bool keep);


      /** product : compared with linMult :
           M is (generally) small and contains singular matrices 
           we just want an IntervalVector */
      friend IntervalVector operator*(const IntervalMatrix& M, 
				const IVparals& iv);

      IVparals& sumFast(const IVparals& iv);
      IVparals& diffFast(const IVparals& iv);

      /*** specific function from diffincl */
 
      /** centered multiply and add 
          this <- center + M *(this-center) + V
       */
      void cmult_and_add(const Vector&center,
		const IntervalMatrix& M, const IntervalMatrix &invM,
		const IntervalVector& V);
      /** tau-centered multiply and add (not modifying the matrices)
       * IV = IV + [0,1]*(M (IV-center) + V)
       * quick algorithm, maybe not precise (generally M is thick)
       */
      void ctau_mult_and_add(const Vector&center,
		const IntervalMatrix& M,
		const IntervalVector& V);
      /**
       * union with     box[d->val] 
				intersected with 
			iv + ]0,1]*(M (iv-center)+V) */
      bool join_intersect_with_tau
		(const IVparals& iv, const Vector&center,
                const IntervalMatrix& M,
                const IntervalVector& V,
		 const IntervalVector& box, int d, double val);

      /** comparison **/
      bool is_subset(const IntervalVector& V) const;
      bool is_subsetFast(const IVparals& iv) const;


      /** replace the matrices with almost-points matrices */
      IVparals& toPointMatrices();
      /** "orthogonalise" a vector of a matrix : change the
        * other M-generators s.t. they become orthogonal to the
        * vector (numM : number of the matrix, ncol : column of the vector */
      IVparals& orthogonalise (int numM, int ncol);
      /** replace a M-generator (column) in a matrix by a new vector
       * numM= -1 : last matrix  ;   ncol=-1 : find best column 
       * modify numM and ncol to get the matrix and column changed 
       * ortho : "orthogonalise" afterwards */
      IVparals& replaceVectorMat
                (const IntervalVector& nl, int& numM, int& ncol,
                 bool ortho);
      /** replace a constraint (line of inverse) in a matrix by another one 
       * numM= -1 : last matrix  ;   nlig=-1 : find best line 
       * modify numM and ncol to get the matrix and line changed */
      IVparals& replaceVectorImat
                (const IntervalVector& nl, int& numM, int& nlig);


      /**
       * simplification
       */
      void simplify(double ratio = 1e-4, int nbit = 2);

      /** generate a ConvexPolygon which is an overapproximation of the
       * projection of the polyhedron (basic form)
       */
      ConvexPolygon over_polygon(const Matrix& M) const;

      /** display
       */
      friend std::ostream& operator<<(std::ostream& str, const IVparals& iv);

#if 0


      /**
       * union (keep the actual matrices)
       */
      void join_with(const IVparals& iv);
     

      /**
       * union with an intersection of an IVparals and the border of
       * a box (keep the matrices)
       * returns true if modified
       */
      bool join_intersect_with
		(const IVparals& iv, 
		 const IntervalVector& box, int d, double val);

      /**
       * intersection with the bound of an interval vector
       */
      bool intersect_with(/* const */ IntervalVector& iv, int d, double val);

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


      /** centered multiply and add 2
       */
      void cmult_and_add2(const Vector&center,
		const IntervalMatrix& M, const IntervalMatrix &invM,
		const IntervalVector& V);

      /** cone-addition, with retur (keep matrices and center) */
      friend IVparals tau_add(const IVparals& iv, const IntervalVector& V);


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

      /** display
       */
      friend std::ostream& operator<<(std::ostream& str, const IVparals& iv);

#endif

    private:
      static uint_fast64_t matIdCnt; // used to identify identical matrices
 

      unsigned int dim; // number of variables
      unsigned int nbmat; // number of intersections, without BBox
      unsigned int nbNcst; // number of ``constant'' matrices
      uint_fast64_t matId; // id of matrices, to identify identical matrices
      bool empty;
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
      inline bool IVparals::is_empty() const { return this->empty; }
}

#endif
