///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include "vibes.h"
#include "expIMat.h"
#include "diffincl.h"
#include "IVparals.h"


// tests

using namespace diffincl;
using namespace codac;

// jacobian of phi = atan(y,x)-z-alpha (alpha is constant)
//    dphi/dx = -y/(x^2+y^2)   dphi/dy = x/(x^2+y^2)
IntervalVector jac_atan2(const IntervalVector& box) {
    IntervalVector res(3,Interval::empty_set());
    res[2] = Interval::one();
    /* base : extrema of the box */
    Interval xmin(box[0].lb());
    Interval xmax(box[0].ub());
    Interval ymin(box[1].lb());
    Interval ymax(box[1].ub());
    res[0] |= -ymin/(sqr(box[0])+sqr(ymin));
    res[0] |= -ymax/(sqr(box[0])+sqr(ymax));
    res[1] |= xmin/(sqr(box[1])+sqr(xmin));
    res[1] |= xmax/(sqr(box[1])+sqr(xmax));
    /* now cases with |x|=|y| */
    Interval p = Interval::pos_reals() & box[0] & box[1];
    if (!p.is_empty()) { Interval a = 1.0/(2.0*p);
			 res[0] |= a; res[1] |= -a; }
    p = Interval::pos_reals() & (-box[0]) & box[1];
    if (!p.is_empty()) { Interval a = 1.0/(2.0*p);
			 res[0] |= a; res[1] |= a; }
    p = Interval::pos_reals() & box[0] & (-box[1]);
    if (!p.is_empty()) { Interval a = (-1.0)/(2.0*p);
			 res[0] |= a; res[1] |= a; }
    p = Interval::neg_reals() & box[0] & box[1];
    if (!p.is_empty()) { Interval a = 1.0/(2.0*p);
			 res[0] |= a; res[1] |= -a; }
    return res;
}

/* goniometric localisation */
void goniometric(IVparals& actval, 
		const Vector& pos,
		const vector<Vector>& beacons, 
		double uncert) {

    IntervalMatrix M(beacons.size(),3);
    IntervalVector av = actval.bounding_box();

    IntervalVector V(beacons.size());
    for (int i=0;i<beacons.size();i++) {
	const Vector& B = beacons[i];
	
        IntervalVector box = B-av;
        M[i] = jac_atan2(box);
        Vector center = B-av.mid();
        double alpha = atan2(B[1]-pos[1],B[0]-pos[0])-pos[2];
        Interval alphaI = Interval(alpha-uncert,alpha+uncert);
        V[i] = atan2(center[1],center[0]) - (alphaI - center[2]);
						// center[2] = av.mid()[2]
        while (V[i].mid()>M_PI) V[i] -= Interval::two_pi();
        while (V[i].mid()<-M_PI) V[i] += Interval::two_pi();
    }
    std::cout << "M : " << M << "\n";
    std::cout << "V : " << V << "\n";
    actval.meetLM(M,av.mid(),V,true);

}
		


int main(int argc, char *argv[]) {
    vibes::beginDrawing();
    std:ifstream inpf;
    char name[50];
    sprintf(name,"%s-fun",argv[1]);
    inpf.open(name);
    DiffInclusion diffincl(inpf); 
    inpf.close();

    IntervalVector frame(diffincl.get_dim());
    IntervalVector X0(diffincl.get_dim());
    
    /* init : [-0.2,0.2] [-0.2,0.2] [0,0] */
    {
      for (int i=0;i<=2;i++) X0[i] = Interval::zero();
    }
    /* frame : [-1000,1000] ... */
    {
       Interval Z1(-1000.0,1000.0);
       for (int i=0;i<=2;i++) frame[i] = Z1;
    }
    TubeVector Result(diffincl.get_time(),frame);

    IVparals actval(X0);
    std::cout << "time:0 " << actval;
    VIBesFig fig("Polygons");
    Matrix projec(2,3);
    projec[0][0]=1.0; projec[0][1]=0.0; projec[0][2]=0.0;
    projec[1][0]=0.0; projec[1][1]=1.0; projec[1][2]=0.0;
    fig.draw_polygon(actval.over_polygon(projec), "red");

    /* position of beacons : (0,10) (10,0) */
    Vector u(3);
    vector<Vector> B;
    u[0] = 0.0; u[1] = 8.0; u[2]=0.0; B.push_back(u);
    u[0] = 8.0; u[1] = -2.0; u[2]=0.0; B.push_back(u);
    u[0] = -8.0; u[1] = -2.0; u[2]=0.0; B.push_back(u);
    Vector pos(3);

    double t=0.0; /* de 1 en 1 */
    double deltat=0.5;
    while (t<30.0) {
      Interval t1(t,t+deltat);
      actval = diffincl.fwd_inclusion(frame,actval,Result,t1,16);
      t=t+deltat;
//      std::cout << "time:" << t << " " << actval;

      fig.draw_polygon(actval.over_polygon(projec), t<=15 ? "red" : "blue");
/*
    // possible position t=10 : (-4,3 ; 1,7 ; 4,5 ) 
    pos[0] = -4.3; pos[1] = 1.9; pos[2]=4.52;
    goniometric(actval,pos,B,0.01);
    std::cout << "time(R):10 " << actval;
*/

    // position time 20 : 0.65, 3.3, 9.27
    
      if (t==15.0) {
        std::cout << "time:15 " << actval;
        pos[0] = -0.5; pos[1] = -0.7; pos[2]=7;
        if (actval.contains(pos)) 
            std::cout << "OK : in\n";
        else
            std::cout << "pas OK : out\n";
        for (int i=0;i<10;i++) goniometric(actval,pos,B,0.001);
        std::cout << "time(R):15 " << actval;
        fig.draw_polygon(actval.over_polygon(projec), "blue");
      }
    }
/*
    // position time 20 : -1.5, -0.4, 12.9.
    
    pos[0] = -1.5; pos[1] = -0.4; pos[2]=12.9;
    goniometric(actval,pos,B,0.01);
    std::cout << "time(R):30 " << actval;
*/

    std::cout << "expo : " << (Result(Result.nb_slices()-1)) << "\n";
    VIBesFigTubeVector VB(name);
    VB.add_tube(&Result, "resultat");
    VB.show(true);

    vibes::endDrawing();
    return 0;
}

