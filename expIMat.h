///////////////////////////////////////////////////////////////////////////
//  expIMat.h : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////

#ifndef __EXPIMAT_H__
#define __EXPIMAT_H__

#include <ibex_IntervalMatrix.h>
#include <ibex_IntervalVector.h>
#include <ibex_Matrix.h>
#include <iostream>

using namespace codac;

/*** A) Utility functions **/

/**  0) basic functions **/

/* no function here */

/** 1) quadratic expressions (optimal computation of M(M+A)  ) **/

/** 
*  returns [x]([x]+a) 
*/
Interval quad_I(const Interval& x, double a);
/**
*   returns [x]([x]+[a]) (optimal if [x] and [a] independant)
*/
Interval quad_II(const Interval& x, const Interval& a);

/**
*    returns [M]([M]+a I)
*/
IntervalMatrix quad_M(const IntervalMatrix& M, double a);

/**
*  returns [M]([M]+[A]) (optimal if [M] and [A] independant)
*/
IntervalMatrix quad_MM(const IntervalMatrix& M, const IntervalMatrix& A);

/** 2) power functions **/

/**
*  returns [M]^2 (special case of quad_M(M,0))
*/
IntervalMatrix square_IntervalMatrix (const IntervalMatrix& M);
/**
*  power of M : compute [M]^n (with n>=0), using product and squares 
*  ( M^(2n) = (M^n)^2 and M^(2n+1)= (M^n)^2 *M )
*/
IntervalMatrix power_IntervalMatrix (const IntervalMatrix& M, unsigned int n);
/**
* squaring (Id+M) : compute Phi(M,k) = [Id+M]^(2^l)-Id , k>=0
*  Phi(M,l+1) = quad_M(Phi(M,l),2)    
*/
void squaring_IntervalMatrix(IntervalMatrix& M, unsigned int l);

/** "primitive" squaring Phi(S,M,0) = S
*   Phi(S,M,k+1) = Phi(S,M,k) + 2^(k-1) M Phi(S,M,k)^2
*/
void primSquaring_IntervalMatrix(IntervalMatrix& S, 
     const IntervalMatrix& M, unsigned int l);


/** B) Matrix inversion **/

/** 1) Inversion **/

/**
*    approximative matrix inversion. (Res := A^(-1) Res )
*    modify the (approximate) inverse matrix, returns true if 
*    successful false if failed (singular?)
*/
bool inv_Matrix(const Matrix& A, Matrix& Res);

/**
*    matrix inversion, crude Gaussian elimination using Rows
*    returns empty if [M] may contains a singular matrix (finding a subset
*    of M would be useful...)
*/
IntervalMatrix inv_IntervalMatrix(const IntervalMatrix& M);
IntervalMatrix inv_IntervalMatrix2(const IntervalMatrix& M);

/** 2) Res :=  A^(-1)Res   ( empty if A contains a singular matrix) 
*   (still crude Gaussian elimination using Rows)
*/
void inv_IntervalMatrix(const IntervalMatrix& A, 
        IntervalMatrix& B);

/** C) Matrix exponential (minus Id) **/

/** 1) Horner scheme **/

/**
*  compute the exponential using Horner scheme, including or not error term,
*  minus Id. Condition : k>=2, use (k/2)+1 matrix products
*/
IntervalMatrix exp_Horner(const IntervalMatrix& M,
		 unsigned int k, bool error);
/**
* adapted Horner scheme for degre 4 (minus Id)
*/
IntervalMatrix exp_Horner4(const IntervalMatrix& M, bool error);

/**
* exp(tau A) - ID with tau \in [0,1]
*/
IntervalMatrix exp_HornerTau(const IntervalMatrix& M,
                 unsigned int k, bool error);

/** 2) Padé approximant **/

/**
* Padé 3,2 approximation (simplified numerator), no error term (FIXME)
* using Id + Den^(-1)*(M*Num)
*/
IntervalMatrix exp_pade32(const IntervalMatrix& M);

/** 3) Scaling and squaring **/

/**
* using generic Horner scheme with error term
*/
IntervalMatrix exp_SSHorner(const IntervalMatrix& M,
           unsigned int k, unsigned int l);
IntervalMatrix exp_SSHornerTau(const IntervalMatrix& M,
           unsigned int k, unsigned int l);
/**
* using degre 4 Horner scheme, and bound for infinity norm instead of l
*/
IntervalMatrix exp_SSHorner4(const IntervalMatrix& M, double norm);

/** D) Integral of matrix exponential **/

/** Compute S(M,t) = Int_{0}^{t} e^Ms ds = Int_{0}{t} e^M(t-s) ds
* with t>0 
*/

/** S(M,t) = t*S0(M,t)
*  S0(M,t) = Id + (tM)/(2!) + (tM)²/(3!) + ... + (tM)^k/((k+1)!) + err
*           with err = ||M||^(k+1) t^(k+1) / ((k+1)!(k+2-||M|| |t|)) E
*        ( TODO : prove it )
*/

// furthermore, we have :
// S((2^k)M,t) = 1/(2^k) (S(M,t))^(2^k)   (idem with S0)
// and if Phi(M,t,k) = S0((2^k)M,t) - Id = 1/(2^k) (S0(M,t)^(2^k) - Id)
//    then Phi(M,t,k+1) = 1/(2^(k+1)) (S0(M,t)^(2^(k+1)) - Id)
//                      = 1/(2^(k+1)) ((2^k Phi(M,t,k) + Id)^2 - Id)
//          = 1/(2^(k+1)) (2^(2k) Phi(M,t,k)^2 + 2^(k+1) Phi(M,t,k) + Id - Id)
//          = 2^(k-1) Phi(M,t,k)^2 + Phi(M,t,k)
//          = 2^(k-1) [ Phi(M,t,k)* [Phi(M,t,k) + 1/(2^(k-1)) Id]]


// compute S0(M,t,k) using Horner scheme, including or not error term,
// minus Id. Condition : k>=2, use (k/2)+1 matrix products
IntervalMatrix expI_Horner(const IntervalMatrix& M,
		 unsigned int k, bool error);


/**
* using generic Horner scheme with error term
*/
IntervalMatrix expI_SSHorner(const IntervalMatrix& M,
           const Interval& t, unsigned int k, unsigned int l);
IntervalMatrix expI_SSHorner4(const IntervalMatrix& M,
            const Interval &t, double norm);

/** E) General computation **/
/** basic version without squaring or equivalent, M (or Mt) is very small,
*  k is number of iterations and must be odd, >=5
*  We need : 
*  a) ExpA = exp(A)-Id =~ sum(1,k+1) A^n/(n!)
*          = 1/2A(A+2Id)+1/6A²(A+ Base1) + ET0
*     error term : nET0 = NA^(k+2) / ((k+2)!(1-NA/(k+3)))
*     Base1 = A²/4(Id+1/5(A+A²/6(Id+1/7(A+...+A²/(k+1)))))
*  b) tauExpA = { 1/tau * (exp(tau A)-Id) | tau \in U } avec u = [0,1]
*             = 1/2A(UA+2Id)+UA²/6(A+UBase1) + ET0
*     UBase1 = UA²/4(Id+U/5(A+UA²/6(Id+U/7(A+...+UA²/(k+1))))) 
*  c) IexpA =~ sum(1,k) A^n/((n+1)!)
*           = 1/6A(A+3Id + Base1) + ET1
*     error term : nET1 = NA^(k+1) / ((k+1)!(k+2-NA))
*  d) tauIexpA =~ sum(1,k) tau^(n-1)A^n/((n+1)!)
*           = 1/6A(UA+3Id + UBase1) + ET1
*  e) VexpA = sum(1,k) -nA^n/(2(n+2)!)
*       = -1/24A(A+2Id+1/2 A²/5(3Id+1/6(4A+A²/7(5Id + 
*                           1/8(6A+A²/9(...+A²/(k+2))))))) + ET2
*     error term : nET2 = (k+5) NA^(k+1)/(2 (k+1)!(k+2-NA)(k+3))
*  f) tauVexpA = { V_tau(A,2) | tau \in U } 
*              = sum(0,k) (-n-2(1-tau))(tau A)^n/(2(n+1)!)
*       = -Id/2 - U/2 (-Id+1/6 A([1,3]Id+ U/4([2,4]A+U/2 A²/5([3,5]Id +
*                 U/6([4,6]A + U A^2/7([5,7]Id + ... U A²/(k+2))))))) + ET2
*  tauIexpA is needed only with slices, VexpA with non-autonomous eqdiff,
*  tauVexpA with both
*/
void global_exp_base(const IntervalMatrix& A,
          int k, bool with_slices, bool with_time,
          IntervalMatrix& ExpA,
          IntervalMatrix& tauExpA,
          IntervalMatrix& IexpA,
          IntervalMatrix& tauIexpA,
          IntervalMatrix& VexpA,
          IntervalMatrix& tauVexpA);

/**
*  Input : matrix M, time t, 
*  Compute simultaneously 
*     a) exp(Mt)  (expM)
*     a') exp(-Mt) (invExpM)
*     b) { 1/a (exp(a t M) - Id), a \in [0,1] }  ( tauExpM )
*     c) int_0^t exp(tau M) dtau  (IexpM)
*     d) { 1/a int_0^(at) exp(tau M) dtau, a \in [0,1] }
*                  (tauIexpM) (with_slices)
*     e) int_0^t int_0^tau exp(theta M) dtheta (VexpM)  (with_time)
*     d) { VexpM(at)/a , a \in [0,1] } (tauVexpM) (with_time, with_slices)
*     e) ub of int_0^t | exp(tau m) | dtau  (intAbsEM)
*/
void global_exp(const IntervalMatrix& M,
          double tim, bool with_slices, bool with_time,
          IntervalMatrix& expM,
	  IntervalMatrix& invExpM,
          IntervalMatrix& tauExpM,
          IntervalMatrix& IexpM,
          IntervalMatrix& tauIexpM,
          IntervalMatrix& VexpM,
          IntervalMatrix& tauVexpM,
          Matrix& intAbsEM);

 
#endif
