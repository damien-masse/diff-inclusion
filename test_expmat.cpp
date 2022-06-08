///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "expIMat.h"


// tests

using namespace codac;

int main(int argc, char *argv[]) {
    IntervalMatrix X(2,2,Interval::zero());
    X[1][0]=X[0][1]=0.5;
    std::cout << "CoefMat,nbitbase,mnorm,radExpM,radIexpM\n";
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
    return 0;
}

