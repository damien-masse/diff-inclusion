///////////////////////////////////////////////////////////////////////////
//  IVdouble.cpp : Interval Vector as double interval
///////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include "expIMat.h"
#include "IVdouble.h"

using namespace codac;

namespace diffincl {

  IVdouble::IVdouble(const IntervalVector& iv) :
      dim(iv.size()), mA(iv.size(),iv.size(), Interval::zero()),
      mAinv(iv.size(),iv.size(),Interval::zero()),
      x0(iv.size()), vB(iv.size())
   {
       for (int i=0;i<this->dim;i++) this->mA[i][i] = Interval::ONE;
       for (int i=0;i<this->dim;i++) this->mAinv[i][i] = Interval::ONE;
       this->vB = iv.mid();
       this->x0 = iv - this->vB;
   }

   IVdouble::IVdouble(const IVdouble& iv) :
       dim(iv.dim), mA(iv.mA), mAinv(iv.mAinv), x0(iv.x0), vB(iv.vB) { }

   IntervalVector IVdouble::bounding_box() const {
       return (this->mA)*(this->x0) + (this->vB);
   }

   void IVdouble::mult_and_add(const IntervalMatrix& M, 
		const IntervalMatrix& Minv, const IntervalVector& V)
   {
       this->vB = M*(this->vB) + V;
       this->mA = M*this->mA;
       this->mAinv = this->mA*Minv;
   }

   void IVdouble::cmult_and_add
		(const Vector& center,
		 const IntervalMatrix& M, 
		 const IntervalMatrix& Minv,
		 const IntervalVector& V)
   {
       this->vB = center + M*(this->vB - center) + V;
       this->mA = M*this->mA;
       this->mAinv = this->mAinv*Minv;
       this->simplify(0.1,0);
   }

   void IVdouble::cmult_and_add2
		(const Vector& center,
		 const IntervalMatrix& M, 
		 const IntervalMatrix& Minv,
		 const IntervalVector& V)
   {
       Vector Vc = V.mid();
       IntervalVector NV = V - Vc;
       this->vB = center + M*(this->vB - center + this->mA * Vc);
       this->x0 = this->x0 + NV;
       this->mA = M*this->mA;
       this->mAinv = this->mAinv*Minv;
   }

   void IVdouble::cmult_and_add3
		(const Vector& center,
		 const IntervalMatrix& M, 
		 const IntervalMatrix& Minv,
		 const IntervalVector& V,
		 const IntervalVector& V2)
   {
       this->vB = center + M*(this->vB - center) + V2;
       this->x0 = this->x0 + V;
       this->mA = M*this->mA;
       this->mAinv = this->mAinv*Minv;
   }

   void IVdouble::cmult_and_add4
		(const Vector& center,
		 const IntervalMatrix& M, 
		 const IntervalMatrix& Minv,
		 const IntervalVector& V, // unused
		 const IntervalVector& V2,
		 const IntervalVector& Vb) 
   {
       this->vB = center + M*(this->vB - center) + V2;
       this->x0 = this->x0 + V;
       this->mA = M*this->mA;
       this->mAinv = this->mAinv*Minv;
   }

   void IVdouble::simplify(double threshold1, double threshold2) {
       int i,j;
       for (i=0;i<this->dim;i++) {
         double maxmig = this->mA[0][i].mig();
         double maxmagO = 0.0;
	 int line_maxmig =0;
         for (j=1;j<this->dim;j++) {
           double mig = this->mA[j][i].mig();
           double mag;
	   if (mig>maxmig) {
              mag = this->mA[line_maxmig][i].mag();
              if (mag>maxmagO) {
		 maxmagO = mag;
              }
              maxmig = mig; 
              line_maxmig = j;
           } else {
              mag = this->mA[j][i].mag();
              if (mag>maxmagO) {
		 maxmagO = mag;
              }
           }
         }
         double prod = (maxmig*this->vB[line_maxmig].rad());
         if (prod==0.0) continue;
//         std::cout << this->mA << "\n";
//         std::cout << "mm " << maxmagO << " mg " << maxmig << 
//                       "  mm/prod " << maxmagO/prod << "\n";
         if (maxmagO/maxmig>=threshold1) continue;
         if (this->vB[line_maxmig].rad()<=threshold2) continue;
         if (this->debug_level>2) {
	   std::cout << "IVdouble simplification (column : " << i 
                     << " line : " << line_maxmig << ")\n";
           std::cout << "bBox before " << this->bounding_box() << "\n";
           std::cout << "XO before " << this->x0 << "\n";
           std::cout << "vB before " << this->vB << "\n";
         }
         double center= this->vB[line_maxmig].mid();
         Interval Idiam = this->vB[line_maxmig] - center;
         Interval pivInter = Idiam / this->mA[line_maxmig][i];
          // pivInter est normalement centré sur 0
         if (this->debug_level>2) 
             std::cout << "pivInter " << pivInter << "\n";
	 this->x0[i] += pivInter;
         for (j=0;j<this->dim;j++) {
	     if (j==line_maxmig) continue;
             this->vB[j] -= this->mA[j][i] * pivInter;
         }
         this->vB[line_maxmig] = center;
         if (this->debug_level>2)
             std::cout << "bBox after " << this->bounding_box() << "\n";
           std::cout << "XO after " << this->x0 << "\n";
           std::cout << "vB after " << this->vB << "\n";
       }
   }
}
