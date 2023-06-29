///////////////////////////////////////////////////////////////////////////
//  diffincl.cpp : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <array>
#include <cstdio>
#include <cstring>
// #include "expIMat.h"
#include "vibes.h"
#include "diffincl.h"

using namespace codac;

namespace diffincl {

  DiffInclusion::DiffInclusion(unsigned int dim, unsigned int udim,
                      const TFunction* tfun,
                      Interval& time, unsigned int nbsteps) :
       dim(dim), udim(udim), ctc(NULL), 
       time(time), nbsteps(nbsteps),
       with_delay(false), with_init(false)
  {
       ctc = new codac2::CtcDiffInclusion(*tfun);
  }
  
  DiffInclusion::DiffInclusion(std::istream& input) {
       char str[256], *token;
       bool fini = false;
       this->ctc = NULL;
       this->dim=-1;
       this->udim=0;
       this->nbsteps=0;
       this->with_delay=false;
       this->with_init=false;
       while (!input.eof()) {
          input.getline(str,255);
          std::cout << str << "\n";
          token=strtok(str," :=");
          if (token==NULL) continue;
          if (strncasecmp(token,"dim",3)==0) {
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"%d",&(this->dim));
             if (this->dim<=0) {
               std::cerr << "Invalid dimension.\n";
             }
          } else if (strncasecmp(token,"udim",4)==0) {
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"%d",&(this->udim));
             if (this->udim<=0) {
               std::cerr << "Invalid dimension.\n";
             }
          } else if (strncasecmp(token,"time",4)==0) {
             double v1=0.0,v2=0.0;
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"[%lg,%lg]",&v1,&v2);
             if (v1>=v2) {
               std::cerr << "Invalid time : " << v1 << " " << v2 << "\n";
             }
             ibex::Interval iv(v1,v2);
               std::cerr << "Interval : " << v1 << " " << v2 << "\n";
             this->time = iv;
          } else if (strncasecmp(token,"delay",5)==0) {
             double v1=0.0,v2=0.0;
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"[%lg,%lg]",&v1,&v2);
             if (v1>=v2 || v1<=0.0) {
               std::cerr << "Invalid delay : " << v1 << " " << v2 << "\n";
             }
             ibex::Interval iv(-v2,-v1);
               std::cerr << "Delay : " << v1 << " " << v2 << "\n";
             this->delay = iv;
             this->with_delay=true;
          } else if (strncasecmp(token,"withinit",8)==0) {
             double v1=0.0,v2=0.0;
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             this->with_init=strchr("tTyYoO",token[0]);
             std::cerr << "With_init : " << this->with_init << "\n";
          } else if (strncasecmp(token,"function",8)==0) {
             assert(this->dim>0);
             const char *xdyn[this->dim+this->udim];
             for (int i=0;i<this->dim+this->udim;i=i+1) {
                 token=strtok(NULL,"(,) :=");
                 xdyn[i] = token;
             } 
             token=strtok(NULL," :");
             for (int i=0;i<this->dim+this->udim;i=i+1) {
                 std::cout << xdyn[i] << "\n";
             } 
             std::cout << token << "\n";
             TFunction *tfun = new TFunction(dim+udim,xdyn,token);
             this->ctc = new codac2::CtcDiffInclusion(*tfun);
          } else if (strncasecmp(token,"nbsteps",7)==0) {
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"%d",&(this->nbsteps));
             assert(this->nbsteps>0);
          } else {
              std::cerr << "Invalid line!";
          }
       }
       assert(this->dim>0);
       assert(this->fun!=NULL || this->tfun!=NULL);
       assert(this->nbsteps>0);
       assert(!(this->time.is_degenerated()));
  }

  ///// forward differential inclusion X' = f(X,t)
  codac2::Tube<codac2::IParals> *DiffInclusion::fwd_inclusion
			(const IntervalVector& frame, 
		  	  const IntervalVector* u_val,
                               const IntervalVector& X0,
                               VIBesFig *fig,
                               int nbsl)
  {
       std::shared_ptr<codac2::TDomain> tdom = codac2::create_tdomain(this->time,this->time.diam()/nbsteps,true);
     codac2::Tube<codac2::IParals>* Res = new codac2::Tube<codac2::IParals>(tdom,dim);
//       codac2::TDomain *tdom = new codac2::TDomain(this->time,this->time.diam()/nbsteps,true);
//       codac2::Tube<codac2::IParals>* Res = new codac2::Tube<codac2::IParals>(*tdom,dim);
       Res->set(frame);
       Res->first_slice().next_slice_ptr()->set(X0);
       std::cout << "first slice " << *Res << "\n";
       assert(X0.dim()==dim);

       assert(ctc!=NULL);
       codac2::Tube<codac2::IParals> *UDom=NULL;
       if (u_val && this->udim>0) {
          UDom = new codac2::Tube<codac2::IParals>(tdom,udim);
          UDom->set(*u_val);
       }
       ctc->contract(*Res,UDom,TimePropag::FORWARD);
       if (UDom!=NULL) delete UDom;
       return Res;
  }

  ///// forward differential inclusion X' = f(X,t)
  codac2::Tube<codac2::IParals> *DiffInclusion::fwd_inclusion
			(const IntervalVector& frame, 
		  	  const TFunction& ufun,
                               const IntervalVector& X0,
                               VIBesFig *fig,
                               int nbsl)
  {
       std::shared_ptr<codac2::TDomain> tdom = codac2::create_tdomain(this->time,this->time.diam()/nbsteps,true);
     codac2::Tube<codac2::IParals>* Res = new codac2::Tube<codac2::IParals>(tdom,dim);
//       codac2::TDomain *tdom = new codac2::TDomain(this->time,this->time.diam()/nbsteps,true);
//       codac2::Tube<codac2::IParals>* Res = new codac2::Tube<codac2::IParals>(*tdom,dim);
       Res->set(frame);
       Res->first_slice().next_slice_ptr()->set(X0);
       std::cout << "first slice " << *Res << "\n";
       assert(X0.dim()==dim);

       assert(ctc!=NULL);
       assert(this->udim>0);
//       codac2::Tube<codac2::IParals> *UDom = new codac2::Tube<codac2::IParals>(tdom,ufun);
       codac2::Tube<codac2::IParals> UDom(tdom,ufun);
       ctc->contract(*Res,&UDom,TimePropag::FORWARD);
       //delete(UDom);
       return Res;
  }

  ///// forward and backward from a slice inclusion X' = f(X,t)
  void DiffInclusion::fwdbwd_contract(codac2::Tube<codac2::IParals> &current, 
                       std::shared_ptr<codac2::Slice<codac2::IParals>> gate,
				const IntervalVector* u_val) {
       codac2::Tube<codac2::IParals> *UDom=NULL;
       if (u_val && this->udim>0) {
          UDom = new codac2::Tube<codac2::IParals>(current.tdomain(),udim);
          UDom->set(*u_val);
       }
       ctc->contract_from_slice(current,UDom,gate,TimePropag::FORWARD | 
					TimePropag::BACKWARD);
       if (UDom!=NULL) delete UDom;
  }

		 
#if 0

  ///// differential inclusion        X' = f(X)
  // dim : number of variables
  // fun : function f of dim variables
  // frame : global frame (no computation outside)
  // X0 : initial position
  // time : time interval [0,time]
  // steps : number of steps
  // result : list of each steps
  // proj : projection for display ( matrix 2*dim )
  void diffInclusion(unsigned int dim, const ibex::Function& fun,
          const ibex::IntervalVector& frame, const ibex::IntervalVector& X0,
          double time, unsigned int nbsteps,
          std::vector<ibex::IntervalVector>& res,
          const ibex::IntervalMatrix& proj) {
       // time step
       const double tsteps = time/nbsteps;
       ibex::IntervalVector fp = proj * frame;
       vibes::drawBox(fp[0].lb(),fp[0].ub(),fp[1].lb(),fp[1].ub(),"blue[cyan]");
  
       ibex::IntervalVector Cent(X0.mid());
       const ibex::IntervalVector X0c = X0-Cent; // centered "initial" place
       const ibex::Matrix Id = ibex::Matrix::eye(dim);
       ibex::IntervalMatrix prodExpM(Id); // products of exponential
       ibex::IntervalVector Xbox(X0); // recomputed box
  
       res.push_back(X0);
       for (int currentstep=1;currentstep<=nbsteps;currentstep++) {
           double currenttime=currentstep*tsteps;
           ibex::IntervalVector approxbox = 
  		this->extend_box_basic(fun,frame,Xbox,tsteps);
           ibex::Vector center_Approxbox = approxbox.mid();
           ibex::IntervalMatrix jac = fun.jacobian(approxbox);
           ibex::Matrix center_jac = jac.mid();
           jac = jac - center_jac;
           ibex::IntervalVector uncert = jac*(approxbox-center_Approxbox);
            
           ibex::IntervalVector fun_eval = fun.eval_vector(center_Approxbox)
                           - center_jac*center_Approxbox; // FIXME : check
           ibex::Vector fun_evalc = fun_eval.mid();
           std::cout << "uncert : " << uncert << "\n";
           uncert += (fun_eval-fun_evalc);
  
           ibex::IntervalVector u2 = fun.eval_vector(approxbox)-fun.eval_vector(center_Approxbox) - (approxbox-center_Approxbox)*center_jac;
           std::cout << "uncert2 : " << u2 << "\n";
           uncert &= u2;
  
           std::cout << "uncert : " << uncert << "\n";
           compute_next_step(dim,fun_evalc,Cent,prodExpM,
  			center_jac,uncert,tsteps);
  //         std::cout << center_jac << "\n" << prodExpM << "\n\n";
           Xbox = Cent +prodExpM * X0c ;
           if (!Xbox.is_subset(approxbox)) {
              std::cerr << "Erreur Xbox//approxbox !!!\n" <<
                  Xbox << "  " << approxbox << "\n\n";
              exit(-1);
           }
           display_state(Cent,Xbox,proj,"black[magenta]","red[magenta]");
           res.push_back(Xbox);
       }
  }
#endif

}
