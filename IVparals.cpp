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
        IntervalVector bb = this->bounding_box();
	bb &= (iv-this->center);
	this->rhs[this->nbmat] = bb;
        this->simplify(true);
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
   void IVparals::intersect_with(const IVparals& ivp) {
        // TODO : intersection with the different elements
	this->intersect_with(ivp.bounding_box());
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
        for (int j=0;j<iv.nbmat;j++) {
            IntervalVector tmpVect= iv.mats[j]*iv.rhs[j] + iv.center;
            ivInter.rhs[this->nbmat] &= tmpVect
			+ Tau * (M * (tmpVect - center)+V)
			- ivInter.center;
            if (ivInter.rhs[this->nbmat].is_empty()) return false;
        }
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
    //   this->rhs[this->nbmat] &= Z;
       for (int i=0;i<this->nbNcst;i++) {
           this->Imats[i] = this->Imats[i] * IM;
           this->mats[i] = M * this->mats[i];
//           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }
       for (int i=this->nbNcst;i<this->nbmat;i++) {
            this->rhs[i] = 
	      hadamard_product(MD[i-this->nbNcst],this->rhs[i])
		+vRes[i-this->nbNcst];
       }
       this->simplify(false); // needed
//       std::cout << nV << "\n" << vRes[this->nbmat-this->nbNcst] << "\n" << this->rhs[this->nbmat] << "\n";

       this->rhs[this->nbmat] +=  nV;
       for (int i=0;i<this->nbNcst;i++) {
           this->rhs[i] = this->rhs[i]+(this->Imats[i]*nV);
       }

       this->simplify(true); // needed
//       std::cout << this->rhs[1] << "\n" << this->mats[0] << "\n" << this->Imats[0] << "\n" << this->rhs[0] << "\n";
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
			+(this->Imats[i]*M)*(this->center-center+V));
       }
       this->rhs[this->nbmat] += Tau*(M*(this->rhs[this->nbmat]+this->center-center+V));
       this->simplify(true);
   }
 
   void IVparals::simplify(bool bwd) {
     if (this->isempty) return;
     /* FIXME : make something for more than one dimension */
     this->rhs[1] &= this->mats[0] * this->rhs[0];
     this->rhs[0] &= this->Imats[0] * this->rhs[1];
     if (bwd) {
       bwd_mul(this->rhs[1],this->mats[0],this->rhs[0],0.001);
       bwd_mul(this->rhs[0],this->Imats[0],this->rhs[1],0.001);
     }
     if (this->rhs[0].is_empty() ||
         this->rhs[1].is_empty())  this->isempty=true;
   }
}

