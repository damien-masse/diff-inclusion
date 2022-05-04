///////////////////////////////////////////////////////////////////////////
//  test_equadiff.cpp : equadiff application of exponential computation
///////////////////////////////////////////////////////////////////////////

#include <codac.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "vibes.h"
#include "expIMat.h"
#include "diffincl.h"


// tests

using namespace diffincl;
using namespace codac;

void parse_command(const char *name,DiffInclusion& di, std::ifstream& input) {
    char str[50], *token;
    IntervalVector frame(di.get_dim());
    IntervalVector X0(di.get_dim());
    
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
      } else if (strcasecmp(str,"compute_fwd")==0) {
           TubeVector Result = di.fwd_inclusion(frame,X0);
    // FIXME : bounding boxes of tubes !!
#if 0
    for (int j=0;j<=1;j++) {
    for (int i=0;i<Result[0].nb_slices();i++) {
      ((Result[j]).slice(i))->set_envelope
		(((Result[j]).slice(i))->input_gate() | ((Result[j]).slice(i))->output_gate());
    }
    }
#endif
           std::cout << "expo : " << (Result(Result.nb_slices()-1)) << "\n";
           VIBesFigTubeVector VB(name);
           VB.add_tube(&Result, "resultat");
    	   VB.show(true);
      } else if (strcasecmp(str,"capd_compute_fwd")==0) {
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
      }
    }
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

