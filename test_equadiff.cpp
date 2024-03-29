///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "vibes.h"
// #include "expIMat.h"
#include "diffincl.h"


// tests

using namespace diffincl;
using namespace codac;

void parse_command(const char *name,DiffInclusion& di, std::ifstream& input) {
    char str[200], *token;
    IntervalVector frame(di.get_dim());
    IntervalVector X0(di.get_dim());
    IntervalVector uval(di.get_udim());
    TFunction *ufun=NULL;
    
    while (!input.eof()) {
      str[0] = '\0';
      input >> str;
       // FIXME : check errors
      if (strncasecmp(str,"frame",5)==0) {
           for (int i=0;i<di.get_dim();i++) {
              input >> str;
              double v1=0.0,v2=0.0;
              sscanf(str,"[%lg,%lg]",&v1,&v2);
              Interval a(v1,v2);
              frame[i] = a;
           }
      } else if (strncasecmp(str,"init",4)==0) {
           for (int i=0;i<di.get_dim();i++) {
              input >> str;
              double v1=0.0,v2=0.0;
              sscanf(str,"[%lg,%lg]",&v1,&v2);
              Interval a(v1,v2);
              X0[i] = a;
           }
      } else if (strncasecmp(str,"uval",4)==0) {
           for (int i=0;i<di.get_udim();i++) {
              input >> str;
              double v1=0.0,v2=0.0;
              sscanf(str,"[%lg,%lg]",&v1,&v2);
              Interval a(v1,v2);
              uval[i] = a;
           }
      } else if (strncasecmp(str,"ufun",4)==0) {
             assert(di.get_udim()>0);
             input >> str;
             std::cout << "ufun : " << str << "\n";
             ufun = new TFunction(str);
      } else if (strcasecmp(str,"compute_fwd")==0) {
            std::cout << "compute_fwd\n";

            codac2::Tube<codac2::IParals> * Res2;
            if (ufun!=NULL) {
               Res2 = di.fwd_inclusion(frame,*ufun,X0,NULL,50);
            } else { 
               Res2 = di.fwd_inclusion(frame,(di.get_udim()==0?NULL:&uval),X0,NULL 
	    /* ,&projec */
   	    ,50);
            }
            std::cout << "termine " << *Res2 << "\n";
            codac::TubeVector Result = to_codac1(*Res2);
    // FIXME : bounding boxes of tubes !!
#if 0
    for (int j=0;j<=1;j++) {
    for (int i=0;i<Result[0].nb_slices();i++) {
      ((Result[j]).slice(i))->set_envelope
		(((Result[j]).slice(i))->input_gate() | ((Result[j]).slice(i))->output_gate());
    }
    }
#endif
           {
              std::shared_ptr<codac2::Slice<codac2::IParals>> p = Res2->last_slice_ptr();
              while (p!=NULL && !p->is_gate())
                       p = p->prev_slice_ptr();
              std::cout << "expo : " << *p << "\n";
           }
           std::cout << "expo : " << (Result(di.get_time().ub())) << "\n";
           VIBesFigTubeVector VB(name);
           VB.add_tube(&Result, "resultat");

           VIBesFig fig("Polygons");

           Matrix projec(2,di.get_dim());
           projec[0][0]=1.0; projec[0][1]=0.0;
           projec[1][0]=0.0; projec[1][1]=1.0;
           for (int i=2;i<di.get_dim();i++) projec[0][i]=projec[1][i]=0.0;
           {
               double a=di.get_time().lb();
               std::shared_ptr<codac2::Slice<codac2::IParals>> p = Res2->first_slice_ptr();
               while (a<=di.get_time().ub()) {
                     while (p!=NULL && (!p->is_gate() || p->t0_tf().lb()<a))
                       p = p->next_slice_ptr();
                     if (p==NULL) break;
                     fig.draw_polygon(p->codomain().over_polygon(projec), "blue");
                     a=a+di.get_time().diam()/75.0;
              }
           }
    	   VB.show(true);
      } else if (strcasecmp(str,"capd_compute_fwd")==0) {
#if 0
           TubeVector Result = di.capd_fwd_inclusion(frame,X0);
    // FIXME : bounding boxes of tubes !!
#if 0
    for (int j=0;j<=1;j++) {
    for (int i=0;i<Result[0].nb_slices();i++) {
      ((Result[j]).slice(i))->set_envelope
		(((Result[j]).slice(i))->input_gate() | ((Result[j]).slice(i))->output_gate());
    }
    }
#endif
           std::cout << "CAPD : " << (Result(Result.nb_slices()-1)) << "\n";
           char namep[50];
           strcpy(namep,name);
           strcat(namep,"-capd");
           VIBesFigTubeVector VB(namep);
           VB.add_tube(&Result, "resultat_capd");
    	   VB.show(true);
#endif
      }
    }
//    if (ufun!=NULL) delete(ufun);
}

int main(int argc, char *argv[]) {
    vibes::beginDrawing();
    std:ifstream inpf;
    char name[50];
    sprintf(name,"%s-fun",argv[1]);
    inpf.open(name);
    DiffInclusion diffincl(inpf); 
    inpf.close();
    sprintf(name,"%s-rec",argv[1]);
    inpf.open(name);
    parse_command(argv[1],diffincl,inpf);
 
/*    TFunction *fun = new TFunction("x","y","(y+[-0.001,0.001],-x+[-0.001,0.001])");
    Interval time(0.0,10.0); 
    DiffInclusion diffincl(2,fun,true,time,50);  */

#if 0
    std::cout.precision(7);

    Interval base(-10.0,10.0);
    IntervalVector frame(2,base);
    Interval d1(-2.0,-2.0);
    Interval d2(-2.0,-2.0);
    IntervalVector  X0(2);
    X0[0] = d1; X0[1]=d2;

    TubeVector Result = diffincl.fwd_inclusion(frame,X0);

    std::cout << "sortie\n" << *((Result[0]).slice(0)) << "\n";
    for (int j=0;j<=1;j++) {
    for (int i=0;i<Result[0].nb_slices();i++) {
      ((Result[j]).slice(i))->set_envelope
		(((Result[j]).slice(i))->input_gate() | ((Result[j]).slice(i))->output_gate());
    }
   }
    VIBesFigTubeVector VB("test1");
    std::cout << "ok1 " << Result.size() << " " << VB.size() << " " << VB.subfigs_number() << "\n";
    VB.add_tube(&Result, "resultat");
    std::cout << "ok2\n";
    VB.show();
#endif
    vibes::endDrawing();
    return 0;
}

