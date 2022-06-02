////////////////////////////////////////////////////////////////////////////
//  IVparals.cpp : representation of polyhedrons as intersections of M_i X_i
////////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include "expIMat.h"
//#include "IVdouble.h"
#include "IVparals.h"

using namespace codac;

namespace diffincl {

  IVparals::IVparals(int dim) :
      dim(dim), isempty(true), nbmat(1), nbNcst(1), mats(1,IntervalMatrix(dim,dim)), Imats(1,IntervalMatrix(dim,dim)), rhs(2),
      center(dim)
   {
       IntervalMatrix mA(dim,dim);
       mA = Matrix::eye(dim);
       Vector v(dim);
       this->mats[0] = mA;
       this->Imats[0] = mA;
       this->center = v;
       this->rhs[0] = this->rhs[1] = IntervalVector::empty(dim);
   }


  IVparals::IVparals(const IntervalVector& iv) :
      dim(iv.size()), isempty(iv.is_empty()), nbmat(1), nbNcst(1), mats(1,IntervalMatrix(iv.size(),iv.size())), Imats(1,IntervalMatrix(iv.size(),iv.size())), rhs(2),
      center(iv.size())
   {
       IntervalMatrix mA(dim,dim);
       mA = Matrix::eye(dim);
       Vector v = iv.mid();
       this->mats[0] = mA;
       this->Imats[0] = mA;
       this->center = v;
       this->rhs[0] = this->rhs[1] = iv-v;
   }

   IVparals::IVparals(const IVparals& iv) :
       dim(iv.dim), isempty(iv.isempty),
       nbmat(iv.nbmat), nbNcst(iv.nbNcst), mats(iv.mats), 
       Imats(iv.Imats), rhs(iv.rhs), center(iv.center) { }

   IVparals::IVparals(const IVparals& iv, const IntervalVector& box) :
       dim(iv.dim), isempty(box.is_empty()),
       nbmat(iv.nbmat), nbNcst(iv.nbNcst), mats(iv.mats), 
       Imats(iv.Imats), rhs(iv.rhs), center(iv.dim) { 
       assert(iv.dim == box.size());
       if (box.is_empty()) return;
       this->center = box.mid();
       IntervalVector nBox = box - box.mid();
       for (int i =0;i<this->nbmat;i++) {
          this->rhs[i] = this->Imats[i] * nBox;
       }
       this->rhs[this->nbmat] = nBox;
   }

   IntervalVector IVparals::bounding_box() const {
       if (this->isempty) return this->rhs[this->nbmat];
       return (this->center+this->rhs[this->nbmat]);
   }

   void IVparals::recenter(const Vector& nCent) {
       Vector vect = this->center - nCent; /* FIXME : interval ? */
       for (int i=0;i<this->nbmat;i++) {
          this->rhs[i] += this->Imats[i] * vect;
       }
       this->rhs[this->nbmat] += vect;
   }

   Vector IVparals::mid() const {
       /* FIXME : can we get something better? */
       return (this->bounding_box()).mid();
   }

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

   void IVparals::intersect_with(const IntervalVector& iv) {
        if (this->isempty) return;
	this->rhs[this->nbmat] &= (iv-this->center);
        this->simplify(true);
   }

   double IVparals::rel_distance_fast(const IVparals& iv) const {
        double a=0;
	for (int i=0;i<this->nbmat;i++) {
           double b= this->rhs[i].rel_distance(iv.rhs[i]);
           if (a<b) a=b;
        }
	return a;
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

   bool IVparals::join_intersect_with_tau(const IVparals& iv,
			const Vector& center, const IntervalMatrix& M,
			const IntervalVector& V,
			const IntervalVector& box, int d, double val) {
        if (iv.is_empty()) return false;
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
        ivInter.rhs[this->nbmat] &= mId*(iv.rhs[this->nbmat] + iv.center)
			+ Tau * (V-M*center)
			- ivInter.center;
        if (ivInter.rhs[this->nbmat].is_empty()) return false;
        for (int j=0;j<iv.nbmat;j++) {
            IntervalVector tmpVect= iv.mats[j]*iv.rhs[j] + iv.center;
            ivInter.rhs[this->nbmat] &= mId*tmpVect
			+ Tau * (V-M*center)
			- ivInter.center;
            if (ivInter.rhs[this->nbmat].is_empty()) return false;
            ivInter.rhs[j] &= (ivInter.Imats[j]*mId*ivInter.mats[j])*iv.rhs[j]
			+ (ivInter.Imats[j]*mId)*iv.center
			+ ivInter.Imats[j]*(Tau*(V-M*center) - ivInter.center);
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
       if (this->isempty) return false;
       Vector ivc = iv - this->center;
       if (!this->rhs[this->nbmat].contains(ivc)) return false;
       for (int i=0;i<this->nbmat;i++) {
          if (!this->rhs[i].intersects(this->Imats[i]*ivc)) return false;
       }
       return true;
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

   void IVparals::cmult_and_add
		(const Vector& center,
			const IntervalMatrix& M, 
			const IntervalMatrix& IM,
			const IntervalVector& V)
   {
       if (this->isempty) return;
       // the center 
       IntervalVector mCenter = center+V+M*(this->center-center);
       this->center = mCenter.mid();
       IntervalVector nV = mCenter - this->center;
#if 0
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
       this->rhs[this->nbmat] = 
	hadamard_product(MD[this->nbmat-this->nbNcst],this->rhs[this->nbmat])
		+vRes[this->nbmat-this->nbNcst];
//        this->rhs[this->nbmat] = M*this->rhs[this->nbmat];
    //   this->rhs[this->nbmat] &= Z;
       for (int i=0;i<this->nbNcst;i++) {
           IntervalMatrix MProd(this->Imats[i]);
           this->Imats[i] = this->Imats[i] * IM;
           IntervalVector Mbla(dim);
           for (int j=0;j<MProd.nb_rows();j++) {
                Mbla[j]=MProd[j][j];
                MProd[j][j] = Interval::zero();
           }
           this->rhs[i] &= hadamard_product(Mbla,this->mats[i]*this->rhs[i]) + MProd*this->rhs[nbmat];
           this->mats[i] = M * this->mats[i];
//           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }
       for (int i=this->nbNcst;i<this->nbmat;i++) {
            this->rhs[i] = 
	      hadamard_product(MD[i-this->nbNcst],this->rhs[i])
		+vRes[i-this->nbNcst];
       }
#else    /** new algorithm **/
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
       /* modified matrices */
       for (int i=0;i<this->nbNcst;i++) {
           IntervalMatrix MProd(this->Imats[i]);
           this->Imats[i] = this->Imats[i] * IM;
/*
           IntervalVector Mbla(dim);
           for (int j=0;j<MProd.nb_rows();j++) {
                Mbla[j]=MProd[j][j];
                MProd[j][j] = Interval::zero();
           }
           this->rhs[i] &= hadamard_product(Mbla,this->mats[i]*this->rhs[i]) + MProd*this->rhs[nbmat];
*/
           this->mats[i] = M * this->mats[i];
//           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }
#endif
       this->simplify(true); // absolutely needed
//       std::cout << nV << "\n" << vRes[this->nbmat-this->nbNcst] << "\n" << this->rhs[this->nbmat] << "\n";

       this->rhs[this->nbmat] +=  nV;
       for (int i=0;i<this->nbmat;i++) {
           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }

//       this->simplify(true); // needed
//       std::cout << this->rhs[1] << "\n" << this->mats[0] << "\n" << this->Imats[0] << "\n" << this->rhs[0] << "\n";
   }


   IVparals tau_add(const IVparals& iv, const IntervalVector& V) {
       IVparals res(iv);
       Interval Tau(0.0,1.0);
       for (int i=0;i<res.nbmat;i++) {
	  res.rhs[i] += Tau*(res.Imats[i]*V);
       }  
       res.rhs[res.nbmat] += Tau*V;
       res.simplify(true);
       return res;
   }

   void IVparals::ctau_mult_and_add
		(const Vector& center,
			const IntervalMatrix& M, 
			const IntervalVector& V)
   {
       if (this->isempty) return;
       Interval Tau(0.0,1.0);
       for (int i=0;i<this->nbmat;i++) {
         this->rhs[i] += Tau*((this->Imats[i]*M*this->mats[i])*this->rhs[i]
			+(this->Imats[i]*M)*(this->center-center)
			+this->Imats[i]*V);
       }
       this->rhs[this->nbmat] += Tau*(M*(this->rhs[this->nbmat]+this->center-center)+V);
       this->simplify(true);
   }
 
   void IVparals::simplify(bool bwd) {
     if (this->isempty) return;
     /* FIXME : make something for more than one dimension */
     for (int i=0;i<=3;i++) {
       this->rhs[1] &= this->mats[0] * this->rhs[0];
       this->rhs[0] &= this->Imats[0] * this->rhs[1];
       if (this->rhs[0].is_empty() ||
           this->rhs[1].is_empty())  { this->rhs[1].set_empty(); 
				     this->isempty=true; }

     if (bwd) {
       bwd_mul(this->rhs[0],this->Imats[0],this->rhs[1],1e-4);
       bwd_mul(this->rhs[1],this->mats[0],this->rhs[0],1e-4);
     }
     }

   }

   /** generate a list of (2D) points, the convex hull of which is an
    * (over)-approximation of the projection of the polyhedron
    */
   ConvexPolygon IVparals::over_polygon(const Matrix& M) const {
        /* first we generate a projection of the parallelotope */
        if (this->isempty) return ConvexPolygon();
        /* just the first polygon (not manage intersection) */
        Vector V1(this->dim);
        ConvexPolygon res;
	Vector cent(M*this->center);
        /* compute the projection for large dimension is a bit complex
           (but interesting), will do it dirty */
        for (int k=0;k<=this->nbmat;k++) {
          bool val[this->dim];
          vector<Point> lpoints;
          for (int i=0;i<this->dim;i++) {
             val[i]=false;
             V1[i] = this->rhs[k][i].lb(); 
          }
          while (true) {
             if (k<this->nbmat) {
	            lpoints.push_back(Point(cent+M*(this->mats[k]*V1)));
             } else {
	            lpoints.push_back(Point(cent+M*V1));
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
   
   IVparals operator-(const IVparals& iv, const Vector& v) {
      IVparals res(iv);
      res.center -=v; /* FIXME : not exact? */
      return res;
   }

   IntervalVector operator*(const IntervalMatrix& M, const IVparals& iv) {
       IntervalVector res = M*(iv.center + iv.rhs[iv.nbmat]);
       for (int i=0;i<iv.nbmat;i++) {
           IntervalMatrix MP = M*iv.mats[i];
           res &= MP*(iv.rhs[i]+iv.Imats[i]*iv.center);
       }
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


}

