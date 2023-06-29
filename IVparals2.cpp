////////////////////////////////////////////////////////////////////////////
//  IVparals.cpp : representation of polyhedrons as intersections of M_i X_i
////////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
// #include "expIMat.h"
#include "IVparals.h"

using namespace codac;
using namespace codac2;

namespace diffincl {
  uint_fast64_t IVparals::matIdCnt=0;

  IVparals::IVparals(int dim) :
      dim(dim), empty(true), nbmat(1), nbNcst(1), mats(1,IntervalMatrix(dim,dim)), Imats(1,IntervalMatrix(dim,dim)), matId(0), rhs(2)
   {
       IntervalMatrix mA(dim,dim);
       mA = Matrix::eye(dim);
       this->mats[0] = mA;
       this->Imats[0] = mA;
       this->rhs[0] = this->rhs[1] = IntervalVector::empty(dim);
   }
  IVparals::IVparals(const IntervalVector& iv) :
      dim(iv.size()), empty(iv.is_empty()), nbmat(1), nbNcst(1), mats(1,IntervalMatrix(iv.size(),iv.size())), Imats(1,IntervalMatrix(iv.size(),iv.size())), matId(0), rhs(2)
   {
       IntervalMatrix mA(dim,dim);
       mA = Matrix::eye(dim);
       this->mats[0] = mA;
       this->Imats[0] = mA;
       this->rhs[0] = this->rhs[1] = iv;
   }
   IVparals::IVparals(const IVparals& iv, const IntervalVector& box) :
       dim(iv.dim), empty(box.is_empty()),
       nbmat(iv.nbmat), nbNcst(iv.nbNcst), mats(iv.mats), 
       Imats(iv.Imats), matId(iv.matId), rhs(iv.rhs) { 
       assert(iv.dim == box.size());
       if (box.is_empty()) {
           for (int i=0;i<=nbmat;i++) this->rhs[i].set_empty();
           return;
       }
       for (int i =0;i<this->nbmat;i++) {
          this->rhs[i] = this->Imats[i] * box;
       }
       this->rhs[this->nbmat] = box;
   }
   IVparals::IVparals(const IVparals& iv) :
       dim(iv.dim), empty(iv.empty),
       nbmat(iv.nbmat), nbNcst(iv.nbNcst), matId(iv.matId), mats(iv.mats), 
       Imats(iv.Imats), rhs(iv.rhs) { }

   IVparals::IVparals(const IntervalMatrix& M, const IntervalMatrix& rM,
			const IntervalVector& V) :
       dim(V.size()), empty(V.is_empty()),
       nbmat(1), nbNcst(1), matId(++matIdCnt), mats(1,M), 
       Imats(1,rM), rhs(2) { 
       if (V.is_empty()) {
           for (int i=0;i<=nbmat;i++) this->rhs[i].set_empty();
           return;
       }
       rhs[0] = V;
       rhs[1] = M*V;
   }
   IVparals::IVparals(const IntervalMatrix& M, const IntervalVector& V) :
       dim(V.size()), empty(V.is_empty()),
       nbmat(1), nbNcst(1), matId(++matIdCnt), mats(1,M), 
       Imats(1,inv_IntervalMatrix(M)), rhs(2) { 
       assert(!Imats[0].is_empty()); /* FIXME : approach for singular M */
       if (V.is_empty()) {
           for (int i=0;i<=nbmat;i++) this->rhs[i].set_empty();
           return;
       }
       rhs[0] = V;
       rhs[1] = M*V;
   }


   /****** Access ******/
   IntervalVector IVparals::bounding_box() const {
       return this->rhs[this->nbmat];
   }
   Vector IVparals::mid() const {
       return this->bounding_box().mid(); /* FIXME : better ? */
   }

   bool IVparals::contains(const Vector& iv) const {
       if (this->empty) return false;
       if (!this->rhs[this->nbmat].contains(iv)) return false;
       for (int i=0;i<this->nbmat;i++) {
          if (!this->rhs[i].intersects(this->Imats[i]*iv)) return false;
       }
       return true;
   }
   const IntervalMatrix& IVparals::getMat(int i) const {
       assert (i>=0 && i<nbmat);
       return this->mats[i];
   }
   const IntervalVector& IVparals::getVec(int i) const {
       assert (i>=0 && i<=nbmat);
       return this->rhs[i];
   }
   IntervalVector IVparals::getPar(const IntervalMatrix& M, int i) const {
       if (i==nbmat) return M*this->rhs[nbmat];
       return (M*this->mats[i])*this->rhs[i];
   }
   IntervalVector IVparals::getPar(const IntervalMatrix& M) const {
       IntervalVector Res = M*this->rhs[nbmat];
       for (int i=0;i<this->nbmat;i++) {
           Res &= (M*this->mats[i])*this->rhs[i];
       }
       return Res;
   }
   
   double IVparals::rel_distanceFast(const IVparals& iv) const {
        assert(this->matId==iv.matId);
        double a=0;
	for (int i=0;i<this->nbmat;i++) {
           double b= this->rhs[i].rel_distance(iv.rhs[i]);
           if (a<b) a=b;
        }
	return a;
   }

   /*** modification ***/

   void IVparals::set_empty() {
        this->empty=true;
        for (int i=0;i<=this->nbmat;i++) (this->rhs[i]).set_empty();
   }
   void IVparals::clear() {
        this->empty=false;
        for (int i=0;i<=this->nbmat;i++) (this->rhs[i]).clear();
   }
   IVparals& IVparals::inflate(double rad) {
        if (this->empty) return *this;
        IntervalVector V(dim,rad);
        return (*this += V);
   }
   IVparals& IVparals::inflateBall(double rad) {
        if (this->empty) return *this;
        for (int i=0;i<=this->nbmat;i++) {
           for (int j=0;j<dim;j++) {
             double a = (i<this->nbmat ? (rad/this->Imats[i][j].norm2()).ub()
				: rad);
	     this->rhs[i][j].inflate(a);
           }
        }
        return *this;
   }
   void IVparals::inflate_from_baseFast(const IVparals& iv, double fact) {
        assert(this->matId==iv.matId);
        for (int i=0;i<=this->nbmat;i++) {
            IntervalVector& a = this->rhs[i];
            const IntervalVector& b = iv.rhs[i];
            for (int j=0;j<dim;j++) {
                double mn = fact*(b[j].lb() - a[j].lb());
                if (mn>0) mn=0.0;
                double mx = fact*(b[j].ub() - a[j].ub());
                if (mx<0) mx=0.0;
                Interval ev(mn,mx);
                a[j] += ev;
            }
        }
    }

   IVparals& IVparals::homothety(IntervalVector c, double delta) {
        if (this->empty) return *this;
        IntervalVector nc = (1-delta)*c;
        for (int i=0;i<=this->nbmat;i++) {
	     this->rhs[i] *= delta;
        }
        return (*this += nc);
   }
   IVparals& IVparals::operator=(const IntervalVector& x) {
        if (x.is_empty()) { this->set_empty(); return *this; }
        this->empty=false;
        this->rhs[nbmat] = x;
        for (int i=0;i<this->nbmat;i++) {
	     this->rhs[i] = this->Imats[i]*x;
        }
        return *this;
   }
   IVparals& IVparals::assign(const IVparals& iv) {
        if (iv.empty) { this->set_empty(); return *this; }
        this->empty=false;
        this->rhs[nbmat] = iv.rhs[iv.nbmat];
        for (int i=0;i<this->nbmat;i++) {
	     this->rhs[i] = this->Imats[i]*iv.rhs[iv.nbmat];
        }
        for (int i=0;i<this->nbmat;i++) {
            for (int j=0;i<this->nbmat;i++) {
	       this->rhs[i] &= (this->Imats[i]*this->mats[j])*iv.rhs[j];
            }
        }
        this->simplify();
        return *this;
   }

   IVparals& IVparals::linMult(const IntervalMatrix& M,
			const IntervalMatrix& IM) {
        if (this->is_empty()) return *this;
        /* not an optimal algorithm, we may need the 
           piecewise-affine interval functions to get a better
           (but possibly non-optimal application */
       for (int i=this->nbNcst;i<this->nbmat;i++) { /* for each unchanging
                                                        matrices */
           IntervalMatrix ImpM = this->Imats[i]*M;
           IntervalVector newRhs((ImpM*this->mats[i])*this->rhs[i]);
                        /* basic modification. */
           for (int j=0;j<this->nbNcst;j++) {
              /* computing  Imats[i]*M*mats[j] / Imats[i]*mat[j]
                                (componentwise division) */
              IntervalMatrix R = ImpM*this->mats[j];
              IntervalMatrix denom = this->Imats[i]*this->mats[j];
              for (int k=0;k<dim;k++) 
              for (int l=0;l<dim;l++) {
                if (denom[k][l].contains(0.0)) { R[k][l]=1.0; continue; }
                R[k][l]/=denom[k][l];
              } 
              R = R.mid();
              /* foreach column V of R */
              for (int k=0;k<dim;k++) {
                /* compute X = V*rhs[i] */
                IntervalVector X(dim);
                for (int l=0;l<dim;l++) X[l] = R[l][k]*this->rhs[i][l];
                /* compute K = ImpM- V*Imats[i] (line by line) */
                IntervalMatrix K(ImpM);
                for (int l=0;l<dim;l++) K[l]-=R[l][k]*this->Imats[i][l];
                X += (K*this->mats[j])*this->rhs[j];
                newRhs &= X;
              }
           }
           this->rhs[i]=newRhs;
       }
       /* case i=nbmat */
       { 
           IntervalVector newRhs(M*this->rhs[nbmat]);
           for (int j=0;j<this->nbNcst;j++) {
               IntervalMatrix R = M*this->mats[j];
               IntervalMatrix& denom = this->mats[j]; 
               for (int k=0;k<dim;k++) 
               for (int l=0;l<dim;l++) {
                 if (this->mats[j][k][l].contains(0.0)) 
                        { R[k][l]=1.0; continue; }
                 R[k][l]/=denom[k][l];
               } 
               R = R.mid();
               /* foreach column V of R */
               for (int k=0;k<dim;k++) {
                 /* compute X = V*rhs[nbmat] */
                 IntervalVector X(dim);
                 for (int l=0;l<dim;l++) X[l] = R[l][k]*this->rhs[nbmat][l];
                 /* compute K = M-V*Id (line by line) */
                 IntervalMatrix K(M);
                 for (int l=0;l<dim;l++) K[l][l]-=R[l][k];
                 X += (K*this->mats[j])*this->rhs[j];
                 newRhs &= X;
               }
            }
            this->rhs[nbmat]=newRhs;
       }
       this->matId=(++matIdCnt);
       for (int i=0;i<this->nbNcst;i++) {
           IntervalMatrix MProd(this->Imats[i]);
           this->Imats[i] = this->Imats[i] * IM;
           this->mats[i] = M * this->mats[i];
       }
       this->simplify(); // absolutely needed
       return *this;
   }


   /**** intersections ******/
   
   IVparals& IVparals::operator&=(const IntervalVector& x) {
        if (this->empty) return *this;
        this->rhs[this->nbmat] &= x;
        if (this->rhs[this->nbmat].is_empty()) { 
           this->set_empty(); return *this; 
        }
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] &= this->Imats[i]*x;
            if (this->rhs[i].is_empty()) { 
              this->set_empty(); return *this; 
            }
        }
        this->simplify();
        return *this;
   }
   IVparals operator&(const IVparals& iv, const IntervalVector& x) {
        IVparals res(iv);
        res &= x;
        return res;
   }
   IVparals& IVparals::meetFast(const IVparals& iv) {
        assert(this->matId==iv.matId);
        if (this->empty) return *this;
        if (iv.empty) { this->set_empty(); return *this; }
        this->rhs[this->nbmat] &= iv.rhs[iv.nbmat];
        if (this->rhs[this->nbmat].is_empty()) { 
           this->set_empty(); return *this; 
        }
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] &= iv.rhs[i];
            if (this->rhs[i].is_empty()) { 
              this->set_empty(); return *this; 
            }
        }
        this->simplify();
        return (*this);
   }
   IVparals& IVparals::meet(const IVparals& iv) {
        if (this->empty) return *this;
        if (iv.empty) { this->set_empty(); return *this; }
        this->rhs[this->nbmat] &= iv.rhs[iv.nbmat];
        if (this->rhs[this->nbmat].is_empty()) { 
           this->set_empty(); return *this; 
        }
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] &= this->Imats[i]*iv.rhs[iv.nbmat];
            for (int j=0;j<this->nbmat;j++) {
                this->rhs[i] &= (this->Imats[i]*iv.mats[j])*iv.rhs[iv.nbmat];
            }
            if (this->rhs[i].is_empty()) { 
              this->set_empty(); return *this; 
            }
        }
        this->simplify();
        return *this;
   }

   bool IVparals::meetLN
      (const IntervalVector& V, const Interval& b, bool keep) {
        /* WARNING : use of bwd_mul */
        if (this->empty) return false;
        Interval bcopy(b);
        IntervalVector VCopy(V);
        if (!bwd_mul(bcopy,VCopy,this->rhs[this->nbmat])) {
            this->set_empty(); return false;
        }
    	for (int i=0;i<this->nbmat;i++) {
            IntervalVector VM = VCopy*this->mats[i];
	    if (!bwd_mul(bcopy,VM,this->rhs[i])) {
		this->set_empty(); return false;
	    }
        }
        if (keep) {
	   this->simplify();
           return !(this->empty);
        }
        IntervalVector Vmid = VCopy.mid();
        int numM=-1; int nlig=-1;
        this->replaceVectorImat (Vmid, numM, nlig);
	/* VCopy.x \in b 
	   => Vmid.x + Veps.x \in b
           => Vmid.x \in b - Veps.x
           => Vmid.x \in (b - Veps rhs[nbmat]) */
        VCopy -= Vmid;
        this->rhs[numM][nlig] &= bcopy - (VCopy*this->rhs[nbmat]);
        if (this->rhs[numM][nlig].is_empty()) {
	   this->set_empty(); return false;
        }
        this->simplify();
        return !(this->empty);
   }

	
   bool IVparals::meetLN(const IntervalVector& V, const IntervalVector& C,
			const Interval& b, bool keep) {
        return this->meetLN(V,(b-V*C),keep);
   }

   /* FIXME : maybe we can do better... though I'm not sure */
   bool IVparals::meetLM(const IntervalMatrix& S, const IntervalVector& b, bool keep) {
        assert(S.nbrows()==b.size());
        /* WARNING : use of bwd_mul */
        if (this->empty) return false;
        IntervalVector bcopy(b);
        IntervalMatrix SCopy(S);
        if (!bwd_mul(bcopy,SCopy,this->rhs[this->nbmat],1e-3)) {
            this->set_empty(); return false;
        }
    	for (int i=0;i<this->nbmat;i++) {
            IntervalMatrix SM = SCopy*this->mats[i];
	    if (!bwd_mul(bcopy,SM,this->rhs[i],1e-3)) {
		this->set_empty(); return false;
	    }
        }
        if (keep) {
	   this->simplify();
           return !(this->empty);
        }
        /* not keep : we apply replaceVectorImat sequentially
	   I'm not sure that's a good idea nevertheless... */
        IntervalMatrix Smid = SCopy.mid();
        SCopy -= Smid;
        for (int i=0;i<S.nb_rows();i=i+1) {
           int numM=-1; int nlig=-1;
           this->replaceVectorImat (Smid[0], numM, nlig);
	   /* VCopy.x \in b 
	      => Vmid.x + Veps.x \in b
              => Vmid.x \in b - Veps.x
              => Vmid.x \in (b - Veps Mi rhs[i]) */
           this->rhs[numM][nlig] &= bcopy[i] - 
			((SCopy[i]*this->mats[numM])*this->rhs[numM]);
           if (this->rhs[numM][nlig].is_empty()) {
	      this->set_empty(); return false;
           }
	   this->simplify();
           if (this->empty) return false;
        }
        return true;
   }
   /* FIXME : maybe we can do better... though I'm not sure */
   bool IVparals::meetLM(const IntervalMatrix& S, const IntervalVector& C,
			const IntervalVector& b, bool keep) {
	  return this->meetLM(S,b+S*C,keep);
   }
   bool IVparals::meetLM(const IntervalMatrix& S, const Vector& C,
			const IntervalVector& b, bool keep) {
	  return this->meetLM(S,b+S*C,keep);
   }

   
   /** union with a box 
    */
   IVparals& IVparals::operator|= (const IntervalVector& x) {
       if (x.is_empty()) return *this;
       if (this->empty) {
           return (*this = x);
       }
       for (int i=1;i<this->nbmat;i++) {
	  this->rhs[i] |= this->Imats[i]*x;
       }
       this->rhs[this->nbmat] |= x;
       return *this;
   }
   IVparals operator|(const IVparals& iv, const IntervalVector& x) {
       IVparals Res(iv);
       Res |= x;
       return Res;
   }


   IVparals& IVparals::operator+=(const IntervalVector& V) {
        if (this->empty) return *this;
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] += this->Imats[i]*V;
        }
        this->rhs[this->nbmat] += V;
        return *this; 
	/* simplification not needed : max C(x+v) = max Cx + max Cv */
   }
   IVparals& IVparals::operator-=(const IntervalVector& V) {
        if (this->empty) return *this;
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] -= this->Imats[i]*V;
        }
        this->rhs[this->nbmat] -= V;
        return *this; 
	/* simplification not needed : max C(x+v) = max Cx + max Cv */
   }
   IVparals operator+(const IVparals& iv, const IntervalVector& V) {
        IVparals res(iv);
        res += V;
        return res;
   }
   IVparals operator-(const IVparals& iv, const IntervalVector& V) {
        IVparals res(iv);
        res -= V;
        return res;
   }

      /** product : compared with linMult :
           M is (generally) small and contains singular matrices 
           we just want an IntervalVector */
   IntervalVector operator*(const IntervalMatrix& M, const IVparals& iv) {
       IntervalVector res = M*iv.rhs[iv.nbmat];
       for (int i=0;i<iv.nbmat;i++) {
           IntervalMatrix MP = M*iv.mats[i];
           res &= MP*iv.rhs[i];
       }
       return res;
   }

   IVparals sum_tau(const IVparals& iv, const IntervalVector& V, bool keep) {
       IVparals res(iv);
       Interval Tau(0.0,1.0);
       for (int i=0;i<res.nbmat;i++) {
	  res.rhs[i] += Tau*(res.Imats[i]*V);
			/* better than res.Imats[i]*(Tau*V) */
       }  
       res.rhs[res.nbmat] += Tau*V;
       res.simplify(true);
       return res;
   }

   void IVparals::cmult_and_add (const Vector& center,
			const IntervalMatrix& M, 
			const IntervalMatrix& IM,
			const IntervalVector& V)
   {
       if (this->empty) return;
       this->linMult(M,IM);
       (*this) += (-M*center + center +V);
   }

   /* quick algorithm :
      1) we want to keep the matrices
      2) more precise would require (much) more work
      FIXME : can we adapt a bit nevertheless, from the algorithm of 
              linMult ? */
   void IVparals::ctau_mult_and_add
		(const Vector& center,
			const IntervalMatrix& M, 
			const IntervalVector& V)
   {
       if (this->empty) return;
       Interval Tau(0.0,1.0);
       for (int i=0;i<this->nbmat;i++) {
         this->rhs[i] += Tau*((this->Imats[i]*M*this->mats[i])*this->rhs[i]
			-(this->Imats[i]*M)*center
			+this->Imats[i]*V);
       }
       this->rhs[this->nbmat] += Tau*(M*(this->rhs[this->nbmat]-center)+V);
       this->simplify(true);
   }
   
   /* FIXME : check and simplify this algorithm */
   bool IVparals::join_intersect_with_tau(const IVparals& iv,
			const Vector& center, const IntervalMatrix& M,
			const IntervalVector& V,
			const IntervalVector& box, int d, double val) {
        if (iv.empty) return false;
        Interval Tau(0.0,1.0);
        bool ret = false;
        IntervalVector cbox(box);
        cbox[d]=val;
        IntervalVector bbox = iv.bounding_box();
        Interval v1(val-bbox[d]);
        Interval v2((M*(bbox-center)+V)[d]);
        if (!bwd_mul(v1,Tau,v2)) return false;
        if (Tau.ub()<=0.0) return false;
        v1 &= Tau*v2;
        bbox[d] = val-v1;
        
        IVparals ivInter(*this, bbox);
        IntervalMatrix mId(dim,dim);
        mId = (Matrix::eye(dim) + Tau*M);
        ivInter.rhs[this->nbmat] &= mId*(iv.rhs[this->nbmat])
			+ Tau * (V-M*center);
        if (ivInter.rhs[this->nbmat].is_empty()) return false;
        for (int j=0;j<iv.nbmat;j++) {
            IntervalVector tmpVect= iv.mats[j]*iv.rhs[j];
            ivInter.rhs[this->nbmat] &= mId*tmpVect
			+ Tau * (V-M*center);
            if (ivInter.rhs[this->nbmat].is_empty()) return false;
            ivInter.rhs[j] &= (ivInter.Imats[j]*mId*ivInter.mats[j])*iv.rhs[j]
			+ ivInter.Imats[j]*(Tau*(V-M*center));
            if (ivInter.rhs[j].is_empty()) return false;
        }
/*
        for (int i=0;i<this->nbmat;i++) {
            IntervalVector Vi(ivInter.Imats[i]*V);
            for (int j=0;j<iv.nbmat;j++) {
                IntervalMatrix tmpM(ivInter.Imats[i] * M);
                ivInter.rhs[i] &= (ivInter.Imats[i] * iv.mats[j]) * iv.rhs[j]
		   + Tau * ((tmpM * iv.mats[j]) * iv.rhs[j] + 
			tmpM*(iv.center - center)+Vi)
			+ ivInter.Imats[i] * (iv.center - ivInter.center);
                if (ivInter.rhs[i].is_empty()) return false;
            }
        }
*/
        ivInter.simplify(true);
        if (ivInter.is_empty()) return false;
        if (this->is_empty()) {
            for (int i=0;i<=this->nbmat;i++) this->rhs[i] = iv.rhs[i];
	    return true;
        }
        for (int i=0;i<=this->nbmat;i++) {
            if (!ivInter.rhs[i].is_subset(this->rhs[i])) {
		this->rhs[i] |= ivInter.rhs[i];
                ret = true;
            }
        }
        if (ret) { this->simplify(true); return true; }
        return false;
   }
 

   IVparals& IVparals::sumFast(const IVparals& iv) {
        if (this->empty) return *this;
        if (iv.empty) { this->set_empty(); return *this; }
        assert(this->matId==iv.matId);
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] += iv.rhs[i];
        }
        this->rhs[this->nbmat] += iv.rhs[this->nbmat];
        return *this; 
   }
   IVparals& IVparals::diffFast(const IVparals& iv) {
        if (this->empty) return *this;
        if (iv.empty) { this->set_empty(); return *this; }
        assert(this->matId==iv.matId);
        for (int i=0;i<this->nbmat;i++) {
            this->rhs[i] -= iv.rhs[i];
        }
        this->rhs[this->nbmat] -= iv.rhs[this->nbmat];
        return *this; 
   }

   
   bool IVparals::is_subset(const IntervalVector& V) const {
       if (this->empty) return true;
       return this->rhs[this->nbmat].is_subset(V);
   }
   bool IVparals::is_subsetFast(const IVparals& iv) const {
       if (this->empty) return true;
       if (iv.empty) return false;
       assert(this->matId==iv.matId);
       for (int i=0;i<this->nbmat;i++) {
           if (!this->rhs[i].is_subset(iv.rhs[i])) return false;
       }
       return this->rhs[this->nbmat].is_subset(iv.rhs[iv.nbmat]);
   }
  

   IVparals& IVparals::toPointMatrices() {
       this->matId = (++matIdCnt);
       for (int i=0;i<this->nbmat;i++) {
           IntervalMatrix M1 = this->mats[i].mid();
           this->Imats[i] = inv_IntervalMatrix(M1); 
           this->rhs[i] = (this->Imats[i]*this->mats[i])*this->rhs[i];
		/* FIXME : is there a different approach ? */
           this->mats[i] = M1;
       }
       this->simplify();
       return *this;
   }

   /* "orthogonalise" a vector of a matrix:
    *  modify the other vectors to make them orthogonal */
   IVparals& IVparals::orthogonalise
		(int numM, int ncol) {
       this->matId = (++matIdCnt);
       IntervalVector nl = this->mats[numM].col(ncol);
       Interval sqnormNL = nl.sqnorm2();
       /* constructing new inverse matrix :
	    M''^-1 = T2 M'^-1
            with T2 = (Id except nlig which is nl.ci/nl^2 (or 1) )
            and rhs :
               V' = T2 V
       */
       IntervalVector prod = nl*this->mats[numM]; /* dot products */
       for (int i=0;i<dim;i=i+1) {
          if (i==ncol) continue;
          this->Imats[numM][ncol] += (prod[i]/sqnormNL)*this->Imats[numM][i];
          this->rhs[numM][ncol] += (prod[i]/sqnormNL)*this->rhs[numM][i];
       }
       /* constructing new matrix (just change nth column) */
       for (int i=0;i<dim;i=i+1) {
          if (i==ncol) continue;
          Interval md = prod[i]/sqnormNL;
          for (int j=0;i<dim;i=i+1) {
             this->mats[numM][j][i] -= md*nl[ncol];
          }
       }
       return *this;
   }
   /* replace a generator (column) in a matrix by another one */
   /* numM= -1 : last matrix  ;   ncol=-1 : find best column */
   /* modify numM and ncol to get the matrix and column changed */
   IVparals& IVparals::replaceVectorMat
		(const IntervalVector& nl, int& numM, int& ncol,
		 bool ortho) {
     assert(numM>=-1 && numM<this->nbmat);
     if (numM==-1) numM=this->nbmat-1;
     assert(ncol>=-1 && ncol<dim);
     bool tryPoint=false;
     IntervalVector u = this->Imats[numM] * nl;
     while (ncol==-1) {
        /* look for a replaceable column : we consider
	   M^(-1) * nl = (u1,u2,...) and select u_i s.t. 0 \notin u_i
           et mag(rhs[i]/ui) minimal
	   -> pb if all u_i has 0 */
	double valbest=0.0;
        for (int i=0;i<dim;i=i+1) {
            if (!u[i].contains(0.0)) {
               double val = (this->rhs[numM][i]/u[i]).mag();
	       if (ncol==-1 || val<valbest) {
	  	  ncol=i; valbest=val;
	       }
            }
        }
        if (ncol==-1) {
             assert(!tryPoint); /* FIXME : other approach? */
             this->toPointMatrices();
             u = this->Imats[numM] * nl;
	     tryPoint=true;
	}
     }
     this->matId = (++matIdCnt);
     /* constructing new inverse matrix :
	M'^-1 = T M^-1
        with T = (Id except ncol which is -uk/ui (or 1/ui) )
        and rhs :
            V' = T V
      */
     this->Imats[numM][ncol] *= (1.0/u[ncol]);
     this->rhs[numM][ncol] /= u[ncol];
     for (int i=0;i<dim;i=i+1) {
        if (i==ncol) continue;
        this->Imats[numM][i] -= u[i]*this->Imats[numM][ncol];
        this->rhs[numM][i] -= u[i]*this->rhs[numM][ncol];
     }
     /* constructing new matrix (just change nth column) */
     for (int i=0;i<dim;i=i+1) {
        this->mats[numM][i][ncol] = nl[i];
     }
     if (ortho) { /* warning : quasi-point matrices recommended... */
         this->orthogonalise(numM,ncol);
     }
     return *this;
   }
   /* replace a constraint (line of inverse) in a matrix by another one */
   /* numM= -1 : last matrix  ;   nlig=-1 : find best line */
   /* modify numM and ncol to get the matrix and line changed */
   IVparals& IVparals::replaceVectorImat
		(const IntervalVector& nl, int& numM, int& nlig) {
     assert(numM>=-1 && numM<this->nbmat);
     if (numM==-1) numM=this->nbmat-1;
     assert(nlig>=-1 && nlig<dim);
     bool tryPoint=false;
     IntervalVector u = nl * this->mats[numM];
     while (nlig==-1) {
        /* look for a replaceable column : we consider
	   nl * M= (u1,u2,...) and select u_i s.t. diam(rhs[i])*mig(ui) maximal
		(diam(rhs[i]) : we don't want to lose "strong" constraint)
	   -> pb if all u_i has 0 */
	double valbest=-1.0;
        for (int i=0;i<dim;i=i+1) {
            if (!u[i].contains(0.0)) {
               double val = this->rhs[numM][i].diam() * u[i].mig();
	       if (val>valbest) {
	  	  nlig=i; valbest=val;
	       }
            }
        }
        if (nlig==-1) {
             assert(!tryPoint); /* FIXME : other approach? */
             this->toPointMatrices();
             u = nl * this->mats[numM];
	     tryPoint=true;
	}
     }
     this->matId = (++matIdCnt);
     /* constructing new inverse matrix :
	M'^-1 = replacement of nl in ith line
        M' = M T with T = (Id except nlig which is -uk/ui (or 1/ui)
        and rhs :
            V' = T^-1 V ( T^-1 = (Id except nlig which is u) )
      */
     this->Imats[numM][nlig] = nl;
     for (int i=0;i<dim;i=i+1)
          this->mats[numM][i][nlig] /= u[i]; 
     for (int i=0;i<dim;i=i+1)
     for (int j=0;j<dim;j=j+1) {
        if (j==nlig) continue;
        this->mats[numM][i][j] -= u[j]*this->mats[numM][i][nlig];
     }
     this->rhs[numM][nlig] = u*this->rhs[numM];
     return *this;
   }

   void IVparals::simplify(double ratio, int nbit) {
     if (this->empty) return;
     /* FIXME : use ratio, rel_distance and a stack ??? */
     for (int i=0;i<=nbit;i++) {
       this->rhs[0] &= this->Imats[0] * this->rhs[nbmat];
       this->rhs[nbmat] &= this->mats[0] * this->rhs[0];
       if (!bwd_mul(this->rhs[0],this->Imats[0],
                    this->rhs[nbmat],ratio)) {
            this->set_empty(); return; }
       if (!bwd_mul(this->rhs[nbmat],this->mats[0],
                    this->rhs[0],ratio)) {
            this->set_empty(); return; }
       for (int j=0;j<nbmat-1;j++) {
         IntervalMatrix Jp1J = this->Imats[j+1] * this->mats[j];
         IntervalMatrix JJp1 = this->Imats[j] * this->mats[j+1];
         this->rhs[j+1] &= Jp1J * this->rhs[j];
         this->rhs[j] &= JJp1 * this->rhs[j+1];
         if (this->rhs[j].is_empty() || this->rhs[j+1].is_empty())
		  { this->set_empty(); return; }
         if (!bwd_mul(this->rhs[j+1],Jp1J,this->rhs[j],ratio)) {
              this->set_empty(); return; }
         if (!bwd_mul(this->rhs[j],JJp1,this->rhs[j+1],ratio)) {
              this->set_empty(); return; }
       }
       this->rhs[nbmat] &= this->mats[nbmat-1] * this->rhs[nbmat-1];
       this->rhs[nbmat-1] &= this->Imats[nbmat-1] * this->rhs[nbmat];
       if (!bwd_mul(this->rhs[nbmat],this->mats[nbmat-1],
                    this->rhs[nbmat-1],ratio)) {
            this->set_empty(); return; }
       if (!bwd_mul(this->rhs[nbmat-1],this->Imats[nbmat-1],
                    this->rhs[nbmat],ratio)) {
            this->set_empty(); return; }
     }
   }


   /** generate a list of (2D) points, the convex hull of which is an
    * (over)-approximation of the projection of the polyhedron
    */
   ConvexPolygon IVparals::over_polygon(const Matrix& M) const {
        /* first we generate a projection of the parallelotope */
        if (this->empty) return ConvexPolygon();
        /* just the first polygon (not manage intersection) */
        Vector V1(this->dim);
        ConvexPolygon res;
        /* compute the projection for large dimension is a bit complex
           (but interesting), will do it dirty */
        for (int k=0;k<=this->nbmat;k++) {
          bool val[this->dim];
          vector<ThickPoint> lpoints;
          for (int i=0;i<this->dim;i++) {
             val[i]=false;
             V1[i] = this->rhs[k][i].lb(); 
          }
          while (true) {
             if (k<this->nbmat) {
	            lpoints.push_back(ThickPoint(M*(this->mats[k]*V1)));
             } else {
	            lpoints.push_back(ThickPoint(M*V1));
             }
             int j=dim-1;
             while (j>=0 && val[j]==true) {
                  V1[j]=this->rhs[k][j].lb();
                  val[j]=false; 
                  j--;
             }
             if (j<0) break;
             val[j]=true;
             V1[j] = this->rhs[k][j].ub(); 
          }
          ConvexPolygon a(lpoints);
          if (k==0) res=a; else res= res & a;
       }
       return res; 
   }


   std::ostream& operator<<(ostream& str, const IVparals& iv) {
       if (iv.empty) { str << "IVparals : empty\n" << flush; return str; }
       str << "IVparals : box " << iv.rhs[iv.nbmat] << "\n";
       for (int i=0;i<iv.nbmat;i++) {
            str << " /\\ " << iv.mats[i] << "\n     X " << iv.rhs[i];
       }
       str << "\n" << flush;
       return str;
   }

#if 0
     /**** old code ****/


   void IVparals::join_with(const IVparals& iv) {
        if (iv.is_empty()) return;
        if (this->isempty) {
	    this->center = iv.center;
            this->rhs[this->nbmat] = iv.rhs[iv.nbmat];
            for (int i=0;i<this->nbmat;i++) {
               this->rhs[i] = this->Imats[i] * iv.rhs[iv.nbmat];
               for (int j=0;j<this->nbmat;i++) {
                  this->rhs[i] &= (this->Imats[i] * iv.mats[j]) * iv.rhs[j];
	       }
            }
        } else {
           this->rhs[this->nbmat] |= (iv.rhs[iv.nbmat] + iv.center - this->center);
           Vector nCent = this->rhs[this->nbmat].mid() + this->center;
           this->recenter(nCent);
           Vector bCent = iv.center - this->center;
           for (int i=0;i<this->nbmat;i++) {
 	      IntervalVector nM = 
		this->Imats[i] * (iv.rhs[iv.nbmat] + bCent);
              for (int j=0;j<this->nbmat;i++) {
                 nM &= (this->Imats[i] * iv.mats[j]) * iv.rhs[j] + 
			(this->Imats[i] * bCent);
	      }
              this->rhs[i] |= nM;
           }
           this->simplify(true);
        }
   }



/*
   void IVparals::intersect_with(const IntervalVector& iv, int d, double val) {
        if (this->isempty) return;
        Vector newCenter = iv.mid();
        newCenter[d]=val;
        this->recenter(newCenter);
	this->rhs[this->nbmat] &= iv - newCenter[d];
        this->rhs[this->nbmat][d]= Interval::zero();
        this->simplify(true);
        /* FIXME : make a cuboid */ /*
        if (!this->isempty) {
           this->rhs[
        }
   }
*/
   void IVparals::intersect_with(/* const */ IVparals& ivp) {
        // TODO : intersection with the different elements
	this->intersect_with(ivp.bounding_box());
   }

  /* intersection with quasi-linear constraints,
     keeping the actual "constraints boxes"
   */
   void IVparals::intersect_with(const IntervalMatrix& M,
			const IntervalVector& Y) {
        if (this->isempty) return;
        IntervalMatrix M1 = M;
        IntervalVector cRhs = this->rhs[this->nbmat] + this->center;
        bwd_mul(Y,M1,cRhs,0.001);
        if (cRhs.is_empty()) { this->isempty=true; 
			this->rhs[this->nbmat].set_empty(); return; }
       this->rhs[nbmat] &= cRhs - this->center;

       for (int j=0;j<this->nbmat;j++) {
           IntervalMatrix M1 = M * this->mats[j];
           IntervalVector cRhs = this->rhs[j] + this->Imats[j] * this->center;
           bwd_mul(Y,M1,cRhs,0.001);
           if (cRhs.is_empty()) { 
		this->isempty=true; 
		this->rhs[this->nbmat].set_empty(); return; }
           this->rhs[j] &= cRhs - this->Imats[j] * this->center;
        }
        this->simplify(true);
   }

  /* intersection with quasi-linear constraints,
     keeping the actual "constraints boxes"
   */
   void IVparals::intersect_with(const IntervalMatrix& M, const Vector& c,
			const IntervalVector& Y) {
        if (this->isempty) return;
        IntervalMatrix M1 = M;
        IntervalVector cRhs = this->rhs[this->nbmat] + (this->center - c);
        bwd_mul(Y,M1,cRhs,0.001);
        if (cRhs.is_empty()) { this->isempty=true; 
			this->rhs[this->nbmat].set_empty(); return; }
       this->rhs[nbmat] &= cRhs - (this->center - c);

       for (int j=0;j<this->nbmat;j++) {
           IntervalMatrix M1 = M * this->mats[j];
           IntervalVector cRhs = this->rhs[j] + this->Imats[j] * (this->center-c);
           bwd_mul(Y,M1,cRhs,0.001);
           if (cRhs.is_empty()) { 
		this->isempty=true; 
		this->rhs[this->nbmat].set_empty(); return; }
           this->rhs[j] &= (cRhs - this->Imats[j] * (this->center-c));
        }
        this->simplify(true);
   }

   void IVparals::inflate_from_base_fast(const IVparals& iv, double fact) {
        assert(this->nbmat==iv.nbmat);
        for (int i=0;i<=this->nbmat;i++) {
            IntervalVector& a = this->rhs[i];
            const IntervalVector& b = iv.rhs[i];
            for (int j=0;j<dim;j++) {
	        double mn = fact*(b[j].lb() - a[j].lb());
                if (mn>0) mn=0.0;  
	        double mx = fact*(b[j].ub() - a[j].ub());
                if (mx<0) mx=0.0;
	        Interval ev(mn,mx);
	        a[j] += ev;
            }
        }
    }


   bool IVparals::join_intersect_with(const IVparals& iv,
				const IntervalVector& box, int d, double val) {
        if (iv.is_empty()) return false;
        bool ret = false;
        IntervalVector bbox(iv.bounding_box());
        bbox &= box;
        if (!bbox[d].contains(val)) return false;
        bbox[d] = val;
        IVparals ivInter(*this,bbox);
        for (int j=0;j<iv.nbmat;j++) {
            ivInter.rhs[this->nbmat] &= iv.mats[j] * iv.rhs[j] + 
					iv.center - ivInter.center;
            if (ivInter.rhs[this->nbmat].is_empty()) return false;
        }
        for (int i=0;i<this->nbmat;i++) {
            for (int j=0;j<iv.nbmat;j++) {
                ivInter.rhs[i] &= (ivInter.Imats[i] * iv.mats[j]) * iv.rhs[j]
			+ ivInter.Imats[i] * (iv.center - ivInter.center);
                if (ivInter.rhs[i].is_empty()) return false;
            }
        }
        ivInter.simplify(true);
        if (ivInter.is_empty()) return false;
        if (this->is_empty()) {
            this->center = iv.center;
            for (int i=0;i<=this->nbmat;i++) this->rhs[i] = iv.rhs[i];
	    return true;
        }
        for (int i=0;i<=this->nbmat;i++) {
            if (!ivInter.rhs[i].is_subset(this->rhs[i])) {
		this->rhs[i] |= ivInter.rhs[i];
                ret = true;
            }
        }
        if (ret) { this->simplify(true); return true; }
        return false;
   }



   bool IVparals::contains(const Vector& iv) const {
   }

   bool IVparals::is_subset_fast(const IVparals& iv) const {
       assert(this->nbmat==iv.nbmat);
       for (int i=0;i<=this->nbmat;i++) {
          if (!this->rhs[i].is_subset(iv.rhs[i])) return false;
       }
       return true;
   }




   // mult and add 
   // principle for kept matrices:
   // we consider rhs[i] = Md*rhs[i] + inter_j (MNd * IMi Mj) rhs[j]
   //             Md = diagonal of IMi * M * Mi
   //                           
   void IVparals::mult_and_add(const IntervalMatrix& M,
 		const IntervalMatrix &IM, const IntervalVector& V)
   {
       if (this->isempty) return;
       // the center 
       IntervalVector mCenter = V+M*this->center;
       this->center = mCenter.mid();
       IntervalVector nV = mCenter - this->center;
       // the multiplication
       vector<IntervalVector> vRes(this->nbmat+1-this->nbNcst);
       IntervalVector Z(this->dim,Interval::zero());
       vector<IntervalVector> MD(this->nbmat+1-this->nbNcst,Z);
       for (int i=this->nbNcst;i<this->nbmat;i++) {
           IntervalMatrix MProd(this->Imats[i]*M*this->mats[i]);
           for (int j=0;j<MProd.nb_rows();j++) {
                MD[i-this->nbNcst][j]=MProd[j][j];
                MProd[j][j] = Interval::zero();
           }
	   vRes[i-this->nbNcst] =
			 (MProd*this->Imats[i]) * this->rhs[this->nbmat];
             // [rhs] + intersect (IMi (M-Id) Mj [rhs(j)])
             for (int j=0;j<this->nbmat;j++) {
	        vRes[i-this->nbNcst] &= 
			(MProd *this->Imats[i]*this->mats[j]) * this->rhs[j];
	    }
           
       }
       // case of the Id matrix, always kept
       {
           IntervalMatrix MProd(M);
           IntervalMatrix MD(this->dim,this->dim,Interval::zero());
           for (int j=0;j<MProd.nb_rows();j++) {
                MD[this->nbmat-this->nbNcst][j]=MProd[j][j];
                MProd[j][j] = Interval::zero();
           }
           vRes[this->nbmat-this->nbNcst] = MProd * this->rhs[this->nbmat];
           // [rhs] + intersect (IMi (M-Id) Mj [rhs(j)])
           for (int j=0;j<this->nbmat;j++) {
                vRes[this->nbmat-this->nbNcst] &= 
			(MProd * this->mats[j]) * this->rhs[j];
           }
       }
       vRes[this->nbmat-this->nbNcst] += nV;
       this->rhs[this->nbmat] = 
	hadamard_product(MD[this->nbmat-this->nbNcst],this->rhs[this->nbmat])
		+vRes[this->nbmat-this->nbNcst];
       for (int i=0;i<this->nbNcst;i++) {
           this->Imats[i] = this->Imats[i] * IM;
           this->mats[i] = M * this->mats[i];
           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }
       for (int i=this->nbNcst;i<this->nbmat;i++) {
            vRes[i-this->nbNcst] += this->Imats[i] * nV;
            this->rhs[i] = 
	      hadamard_product(MD[i-this->nbNcst],this->rhs[i])
		+vRes[i-this->nbNcst];
       }
       // this->simplify(); // needed ?
   }




   
   IVparals operator-(const IVparals& iv, const Vector& v) {
      IVparals res(iv);
      res.center -=v; /* FIXME : not exact? */
      return res;
   }

      

   std::ostream& operator<<(ostream& str, const IVparals& iv) {
       if (iv.isempty) { str << "IVparals : empty\n" << flush; return str; }
       str << "IVparals : (c) " << iv.center << " box " << (iv.center + iv.rhs[iv.nbmat]) << "\n";
       for (int i=0;i<iv.nbmat;i++) {
            str << " /\\ " << iv.mats[i] << "\n     X " << iv.rhs[i];
       }
       str << "\n" << flush;
       return str;
   }

#endif

}

