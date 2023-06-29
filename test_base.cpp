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
    vibes::beginDrawing();

    TFunction f("x","y","z","u","v","w","(-y*z+u,x*z+v,-y*z+w)");
    TFunction f2("(sin(1.732*t)*sin(t)+[-0.001,0.001],cos(1.732*t),sin(1.732*t)*cos(t))");
    Interval a(0.0,6.5);
    codac2::CtcDiffInclusion *ctc = new codac2::CtcDiffInclusion(f);
    std::shared_ptr<codac2::TDomain> tdom = codac2::create_tdomain(a,0.001,true);
    codac2::Tube<codac2::IParals>* ufun = new codac2::Tube<codac2::IParals>(tdom,f2);
    codac2::Tube<codac2::IParals>* Res= new codac2::Tube<codac2::IParals>(tdom,3);

     Interval b(-10.0,10.0);
     IntervalVector frame(3,b);
     Res->set(frame);
     Interval c(0.0,0.0);
     IntervalVector X0(3,c);
     Res->first_slice().next_slice_ptr()->set(X0);
     ctc->contract(*Res,ufun,TimePropag::FORWARD);

    codac::TubeVector Result1 = to_codac1(*Res);
    codac::TubeVector Result2 = to_codac1(*ufun);

    std::cout << Result2  << " ok\n";
    VIBesFigTubeVector VB("essai");
    VB.add_tube(&Result1,"results");
    VB.add_tube(&Result2,"donnees");
    VB.show(true);

//    delete f2;
//    delete ufun;
    vibes::endDrawing();
    return 0;
}

