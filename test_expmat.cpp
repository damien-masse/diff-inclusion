///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
// #include "expIMat.h"


// tests

using namespace codac;
using namespace codac2;

int main(int argc, char *argv[]) {
    IntervalMatrix X(4,4,Interval::zero());
    X[0][0]=2; 
    X[0][1]=0; 
    X[0][2]=-1; 
    X[0][3]=0; 
    X[1][0]=1; 
    X[1][1]=-1; 
    X[1][2]=0;
    X[1][3]=2;
    X[2][0]=-1; 
    X[2][1]=0;
    X[2][2]=1;
    X[2][3]=0;
    X[3][0]=0; 
    X[3][1]=-1;
    X[3][2]=1;
    X[3][3]=-1;
    IntervalMatrix A(4,4,Interval::zero());
    A[0][0]=0.3; 
    A[0][1]=0.3; 
    A[0][2]=0; 
    A[0][3]=0; 
    A[1][0]=0; 
    A[1][1]=0.3; 
    A[1][2]=0;
    A[1][3]=0;
    A[2][0]=0; 
    A[2][1]=0;
    A[2][2]=0;
    A[2][3]=0.3;
    A[3][0]=0; 
    A[3][1]=0;
    A[3][2]=-0.3;
    A[3][3]=0;
/*
    mnorm=5;
    nbitbase=15;
*/
#if 0
    std::cout << "CoefMat,nbitbase=14,mnorm=5.0,radExpM,radIexpM\n";

    nbitbase=14;
	    global_exp(X, 1, true, false, A1,A2,A3,A4,A5,A6,A7,A8);
	    std::cout << infinite_norm(A1.diam()) << "\n" << infinite_norm(A4.diam()) << "\n";
            std::cout << A1 << "\n" << A4 << "\n";
	    std::cout << infinite_norm(A3.diam()) << "\n" << infinite_norm(A5.diam()) << "\n";
            std::cout << A3 << "\n" << A5 << "\n";
    std::cout << "CoefMat,nbitbase=2,mnorm=0.02,radExpM,radIexpM\n";
#endif
#if 0
    nbitbase=2;
    mnorm = 0.02;
	    global_exp(X, 1, true, false, A1,A2,A3,A4,A5,A6,A7,A8);
	    std::cout << infinite_norm(A1.diam()) << "\n" << infinite_norm(A4.diam()) << "\n";
            std::cout << A1 << "\n" << A4 << "\n";
//	    std::cout << infinite_norm(A3.diam()) << "\n" << infinite_norm(A5.diam()) << "\n";
//            std::cout << A3 << "\n" << A5 << "\n";
   IntervalMatrix SaveA1(A1);
    X[0][1]=1.0;
    X[1][0]=-1.0;
    X[1][1]=-2.0;
	    global_exp(X, 1, true, false, A1,A2,A3,A4,A5,A6,A7,A8);
	    std::cout << infinite_norm(A1.diam()) << "\n" << infinite_norm(A4.diam()) << "\n";
            std::cout << A1 << "\n" << A4 << "\n";
//	    std::cout << infinite_norm(A3.diam()) << "\n" << infinite_norm(A5.diam()) << "\n";
//            std::cout << A3 << "\n" << A5 << "\n";
   IntervalMatrix A1a(A1);
    X[0][1]=0.0;
    X[1][0]=0.0;
    X[1][1]=Interval(-0.01,0.01);
	    global_exp(X, 1, true, false, A1,A2,A3,A4,A5,A6,A7,A8);
	    std::cout << infinite_norm(A1.diam()) << "\n" << infinite_norm(A4.diam()) << "\n";
            std::cout << A1 << "\n" << A4 << "\n";
//	    std::cout << infinite_norm(A3.diam()) << "\n" << infinite_norm(A5.diam()) << "\n";
//            std::cout << A3 << "\n" << A5 << "\n";
   IntervalMatrix A1b(A1);
   std::cout << A1a*A1b << "\n";
#endif
   std::cout << "X : " << X << "\n";
   IntervalMatrix InvX = inv_IntervalMatrix(X);
   std::cout << "InvX " << InvX << "\n";
    IntervalMatrix AP(4,4,Interval::zero());
    AP[0][0]=1.2; 
    AP[0][1]=-0.3; 
    AP[0][2]=1.5; 
    AP[0][3]=-0.3; 
    AP[1][0]=-0.3; 
    AP[1][1]=0; 
    AP[1][2]=-0.9;
    AP[1][3]=0;
    AP[2][0]=-0.6; 
    AP[2][1]=0.2;
    AP[2][2]=-0.7;
    AP[2][3]=0.1;
    AP[3][0]=0; 
    AP[3][1]=0.2;
    AP[3][2]=0.2;
    AP[3][3]=0.1;
//   IntervalMatrix AP = X * A *InvX;
   std::cout << "AP - XAX-1 " << AP - X * A * InvX  << "\n";
   IntervalMatrix A1(4,4),A2(4,4),A3(4,4),A4(4,4),A5(4,4),A6(4,4),A7(4,4);
   Matrix A8(4,4);
   global_exp(A, 1, true, true, A1,A2,A3,A4,A5,A6,A7,A8);
   std::cout << "A1 : " << A1 << "\n"
             << "A2 : " << A2 << "\n"
             << "A3 : " << A3 << "\n"
             << "A4 : " << A4 << "\n"
             << "A5 : " << A5 << "\n"
             << "A6 : " << A6 << "\n"
             << "A7 : " << A7 << "\n"
             << "A8 : " << A8 << "\n";
   global_exp(AP, 1, true, true, A1,A2,A3,A4,A5,A6,A7,A8);
   std::cout << "AP1 : " << A1 << "\n"
             << "AP2 : " << A2 << "\n"
             << "AP3 : " << A3 << "\n"
             << "AP4 : " << A4 << "\n"
             << "AP5 : " << A5 << "\n"
             << "AP6 : " << A6 << "\n"
             << "AP7 : " << A7 << "\n"
             << "AP8 : " << A8 << "\n";
   std::cout << "PAP1 : " << InvX * A1 * X << "\n"
             << "PAP2 : " << InvX * A2 * X << "\n"
             << "PAP3 : " << InvX * A3 * X << "\n"
             << "PAP4 : " << InvX * A4 * X << "\n"
             << "PAP5 : " << InvX * A5 * X << "\n"
             << "PAP6 : " << InvX * A6 * X << "\n"
             << "PAP7 : " << InvX * A7 * X << "\n"
             << "PAP8 : " << InvX * A8 * X << "\n";

/*
    while (X[1][0].ub()>1e-6) {
      for (int i=1;i<10;i++) {
	nbitbase=i;
        mnorm=1e-1;
	for (int j=1;j<20;j++) {
            std::cout << X[1][0].ub() << ",";
            std::cout << nbitbase << "," << mnorm << ",";
          
            IntervalMatrix A1(2,2),A2(2,2),A3(2,2),A4(2,2),A5(2,2),A6(2,2),
			A7(2,2);
            Matrix A8(2,2);
	    global_exp(X, 1, false, false, A1,A2,A3,A4,A5,A6,A7,A8);
            std::cout << infinite_norm(A1.diam()) << "," <<
		infinite_norm(A4.diam()) << "\n";
            mnorm = mnorm/2.0;
	}
      }
      X *= 0.25;
    }
*/
    return 0;
}

