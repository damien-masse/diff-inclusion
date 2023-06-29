///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include "vibes.h"
#include "diffincl.h"


// tests

using namespace codac;
using namespace diffincl;


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
    Matrix projec(2,3);
    projec[0][0]=1.0; projec[0][1]=0.0; projec[0][2]=0.0;
    projec[1][0]=0.0; projec[1][1]=1.0; projec[1][2]=0.0;
    
    /* init : [-0.2,0.2] [-0.2,0.2] [0,0] */
    {
      for (int i=0;i<=2;i++) X0[i] = Interval::zero();
      X0[0].inflate(0.3);
      X0[1].inflate(0.3);
    }
    /* frame : [-1000,1000] ... */
    {
       Interval Z1(-1000.0,1000.0);
       for (int i=0;i<=2;i++) frame[i] = Z1;
    }


    codac2::TubeVector *Res2 = diffincl.fwd_inclusion(frame,NULL,X0,NULL,50);
    std::cout << "fwd termine " << *Res2 << "\n";

    /* dessin des polygones */
    VIBesFig fig("Polygons");
    {
      double a=0.0;
      codac2::SliceVector *p = &(Res2->first_slice());
      while (a<30.0) {
          while (p!=NULL && (!p->is_gate() || p->t0_tf().lb()<a))
		 p = p->next_slice();
          if (p==NULL) break;
          fig.draw_polygon(p->codomainI().over_polygon(projec), "black");
          a=a+3.0;
      }
    }          

    /* fonction polar : atan2(x(1)-u) */
    Function gonio("x[3]","u[2]","atan2(u[1]-x[1],u[0]-x[0])-x[2]");
    codac2::CtcBwdFun bwdFun(gonio);


/*
    IVparals actval(X0);
    std::cout << "time:0 " << actval;
    VIBesFig fig("Polygons");
    fig.draw_polygon(actval.over_polygon(projec), "red");
*/

    /* position of beacons : (0,10) (10,0) */
    std::vector<std::pair<const IntervalVector*, const IntervalVector>> r;
    
    IntervalVector u(2);
    
    u[0] = 0.0; u[1] = 8.0; 
    IntervalVector rs(1,-7.05+atan2(8.0+0.65,0.0+0.5));
    rs.inflate(0.002);
    std::pair<const IntervalVector*, const IntervalVector> 
		pr1(new IntervalVector(u),rs);
    r.push_back(pr1);

/*
    u[0] = 8.0; u[1] = -2.0;
    rs[0]=-7.05+atan2(-2.0+0.65,8.0+0.5);
    rs.inflate(0.002);
    std::pair<const IntervalVector*, const IntervalVector> 
		pr2(new IntervalVector(u),rs);
    r.push_back(pr2);
*/

/*

    u[0] = -8.0; u[1] = -2.0;
    rs[0]=-7.05+atan2(-2.0+0.65,-8.0+0.5);
    rs.inflate(0.002);
    std::pair<const IntervalVector*, const IntervalVector> 
		pr3(new IntervalVector(u),rs);
    r.push_back(pr3);

*/ 
    // look for time 15
    codac2::SliceVector *p = &(Res2->first_slice());
    while (!p->is_gate() || p->t0_tf().lb()<15.0) p = p->next_slice();
    std::cout << "avant contraction : " << (*p) << " " << p->codomainI() << "\n";
    bwdFun.contract(*p,r);
    std::cout << "apres contraction : " << (*p) << " " << p->codomainI() << "\n";
   
    diffincl.fwdbwd_contract(*Res2,*p,NULL);
    {
      double a=0.0;
      codac2::SliceVector *p = &(Res2->first_slice());
      while (a<30.0) {
          while (p!=NULL && (!p->is_gate() || p->t0_tf().lb()<a))
		 p = p->next_slice();
          if (p==NULL) break;
          fig.draw_polygon(p->codomainI().over_polygon(projec), "blue");
          a=a+3.0;
      }
    }          
    fig.draw_polygon(p->codomainI().over_polygon(projec), "red");
    
    TubeVector Result=Res2->to_codac1();
    VIBesFigTubeVector VB(name);
    VB.add_tube(&Result, "resultat");
    VB.show(true);

    vibes::endDrawing();

    // position time 15 : -0.5, -0.7, 7.0
    
/*
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
*/

#if 0
    double t=0.0; /* de 1 en 1 */
    double deltat=0.5;
    while (t<30.0) {
      Interval t1(t,t+deltat);
      actval = diffincl.fwd_inclusion(frame,actval,Result,t1,16);
      t=t+deltat;
//      std::cout << "time:" << t << " " << actval;

      fig.draw_polygon(actval.over_polygon(projec), t<=15 ? "red" : "blue");
#endif
/*
    // possible position t=10 : (-4,3 ; 1,7 ; 4,5 ) 
    pos[0] = -4.3; pos[1] = 1.9; pos[2]=4.52;
    goniometric(actval,pos,B,0.01);
    std::cout << "time(R):10 " << actval;
*/

/*
    // position time 20 : -1.5, -0.4, 12.9.
    
    pos[0] = -1.5; pos[1] = -0.4; pos[2]=12.9;
    goniometric(actval,pos,B,0.01);
    std::cout << "time(R):30 " << actval;
*/
/*
    std::cout << "expo : " << (Result(Result.nb_slices()-1)) << "\n";
*/
    return 0;
}

