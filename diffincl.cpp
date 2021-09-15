///////////////////////////////////////////////////////////////////////////
//  diffincl.cpp : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include "expIMat.h"
#include "vibes.h"
#include "diffincl.h"

using namespace codac;

namespace diffincl {

  DiffInclusion::DiffInclusion(unsigned int dim,
                      const TFunction* tfun,
                      Interval& time, unsigned int nbsteps) :
       dim(dim), tfun(tfun),
       time(time), nbsteps(nbsteps) {
       this->time_dependent=true;
  }
  
  DiffInclusion::DiffInclusion(unsigned int dim,
                      const Function* fun, 
                      Interval& time, unsigned int nbsteps) :
       dim(dim), tfun(tfun),
       time(time), nbsteps(nbsteps) {
       this->time_dependent=false;
  }

  DiffInclusion::DiffInclusion(std::istream& input) {
       char str[256], *token;
       bool fini = false;
       this->dim=-1;
       this->nbsteps=0;
       this->time_dependent = false;
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
          } else if (strncasecmp(token,"timedependent",12)==0) {
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             if ((token[0]=='1') || (token[0]=='y') || (token[0]=='Y')) {
                 this->time_dependent = true;
             }
          } else if (strncasecmp(token,"time",4)==0) {
             double v1=0.0,v2=0.0;
             token=strtok(NULL," :=");
             if (token==NULL) continue;
             sscanf(token,"[%lg,%lg]",&v1,&v2);
             if (v1>=v2) {
               std::cerr << "Invalid time.\n";
             }
             ibex::Interval iv(v1,v2);
             this->time = iv;
          } else if (strncasecmp(token,"function",8)==0) {
             assert(this->dim>0);
             const char *xdyn[this->dim];
             for (int i=0;i<this->dim;i=i+1) {
                 token=strtok(NULL,"(,) :=");
                 xdyn[i] = token;
             } 
             token=strtok(NULL," :");
             for (int i=0;i<this->dim;i=i+1) {
                 std::cout << xdyn[i] << "\n";
             } 
             std::cout << token << "\n";
             if (this->time_dependent) 
                this->tfun = new TFunction(dim,xdyn,token);
             else
                this->fun = new Function(dim,xdyn,token);
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

  // evaluate function on box + time (time dependant or not)
  const IntervalVector DiffInclusion::teval_vector(const Interval &tim,
                 const IntervalVector& cdom) {
        if (fun) {
           return fun->eval_vector(cdom);
        } 
        assert(tfun);
	const Function& func = tfun->getFunction();
        IntervalVector box(dim+1);
        box[0]=tim;
        box.put(1,cdom);
        return func.eval_vector(box);
  }
  const IntervalVector DiffInclusion::teval_vector(double tim,
                 const IntervalVector& cdom) {
        const Interval timI(tim);
        return this->teval_vector(timI,cdom);
  }
  const IntervalVector DiffInclusion::teval_vector(double tim,
                 const Vector& cdom) {
        const Interval timI(tim);
        const IntervalVector cdomI(cdom);
        return this->teval_vector(timI,cdomI);
  }
  
  // compute jacobian ; in case of time dependency, tvec is the first 
  // column of the jacobian
  IntervalMatrix DiffInclusion::jacobian(const IntervalVector& codom,
                 const Interval& tdom,
                 IntervalVector& tvec) {
        if (this->fun) {
           tvec.clear();
           return fun->jacobian(codom);
        }
        assert(this->tfun);

	const Function& func = tfun->getFunction();
        IntervalVector box(dim+1);
        box[0]=tdom;
        box.put(1,codom);
        IntervalMatrix bjac = func.jacobian(box);
//        std::cout << "jacob : " << bjac << "\n";
        for (int i=0;i<=dim-1;i++) { tvec[i]=bjac[i][0]; }
        return bjac.cols(1,dim);
  }
  
  // very basic enclosing for evolution
  IntervalVector DiffInclusion::extend_box_basic(const IntervalVector& frame,
                 	const IntervalVector& Xbox,
                 	const Interval& tim) {
      double tstep = tim.diam();
      Interval btime(0.0,tstep);
      /* estimation des pentes */
      IntervalVector k1 = this->teval_vector(tim.lb(),Xbox);
      k1 |= this->teval_vector(tim.mid(), Xbox+tstep/2.0*k1);
      k1 |= this->teval_vector(tim, Xbox+btime*k1);
      /* compute approximation of the result */
      IntervalVector Res = Xbox + btime*k1;
      /* inflation */
      Res.inflate(1.2,0.0);
      Res &= frame;
      Res = Xbox + btime * this->teval_vector(tim,Res);
      Res &= frame;
      return Res;
  }
  
        
        

  // compute next step, from linear approximation and uncertainty
  // centzone : center of the zone
  // fcent = f(center of the zone, half time)
  // tVect = mid of temporal differential
  // jac : mid of jacobian
  // prodExpM : cumulated product of exponential (without time)
  // uncert : U value (uncertity on f)
  // double tsteps : time step
  // f(x,t) \in fcent + (t-tmid)*tvect + jac*(x-x0) + uncert
  void DiffInclusion::compute_next_step(const Vector& centzone,
			 const Vector& fcent,
                         const Vector& tVect,
                         const Matrix& jac,
                         const IntervalVector& uncert,
                         double tsteps,
  	                 IntervalVector& Cent, 
                         IntervalMatrix& prodExpM,
			 IntervalVector& tauCent,
			 IntervalMatrix& tauProdExpM) {
     // global computation
     IntervalMatrix ExpM(dim,dim);
     IntervalMatrix IExpM(dim,dim);
     IntervalMatrix VExpM(dim,dim);
     IntervalMatrix tauExpM(dim,dim);
     IntervalMatrix tauIExpM(dim,dim);
     IntervalMatrix tauVExpM(dim,dim);
     Matrix IntAbs(dim,dim);
     const Matrix Id = Matrix::eye(dim);
     // computation of exponentials and other values
     global_exp(jac,tsteps,true,time_dependent,ExpM,tauExpM,IExpM,tauIExpM,VExpM,tauVExpM,IntAbs);
//     std::cout << "\n\nExpM" << ExpM << "\ntauExpM" << tauExpM << "\nIExpM" << IExpM << "\ntauIExpM" << tauIExpM << "\nVExpM" << VExpM << "\ntauVExpM" << tauVExpM << "\nIntAbs" << IntAbs << "\n";
     // intermediates values
     Vector vuncert = IntAbs * uncert.ub();
     IntervalVector evolCenter = tauIExpM * fcent;
     if (time_dependent) {
           evolCenter += tauVExpM * tVect;
     }
//     std::cout << "tauIExpM * fcent : " << evolCenter << "\n";
     tauCent = evolCenter;
     for (int i=0;i<dim;i=i+1) (tauCent[i]).inflate(vuncert[i]); // FIXME !!!
     tauProdExpM = tauExpM;

     // final values
     evolCenter = IExpM * fcent;
     if (time_dependent) {
           evolCenter += VExpM * tVect;
     }
     Cent = centzone + ExpM*(Cent-centzone) + evolCenter;
     for (int i=0;i<dim;i=i+1) (Cent[i]).inflate(vuncert[i]);
     prodExpM = ExpM * prodExpM;

  }
  
  ///// forward differential inclusion X' = f(X,t)
  TubeVector DiffInclusion::fwd_inclusion(const IntervalVector& frame, 
                               const IntervalVector& X0) {
       assert(fun != NULL);
       assert(X0.dim()==dim);
       const double tsteps = time.diam()/nbsteps;
  
       IntervalVector Cent(X0.mid());
       const Interval Unt(0.0,1.0);
       const IntervalVector X0c = X0-Cent; // centered "initial" place
       const Matrix Id = Matrix::eye(dim);
       IntervalMatrix prodExpM(Id); // products of exponential
       IntervalVector Xbox(X0); // recomputed box

       // result
       TubeVector Result(time,frame);
  
       Result.sample(time.lb(),X0);
       double oldtime = time.lb();
       for (int currentstep=1;currentstep<=nbsteps;currentstep++) {
           double currenttime=time.lb()+currentstep*tsteps;
           if (currenttime>time.ub()) currenttime=time.ub();
           Interval timeslice(oldtime,currenttime);
           IntervalVector approxbox = 
  		this->extend_box_basic(frame,Xbox,timeslice);
           Vector center_Approxbox = approxbox.mid();
           IntervalVector tdiff(dim);
           IntervalMatrix jac = this->jacobian(approxbox,timeslice,tdiff);
           Matrix center_jac = jac.mid();
           jac = jac - center_jac;
           Vector center_tdiff = tdiff.mid();
           tdiff = tdiff - center_tdiff;

           // compute uncert on f : 
           // f(x,t) in f(xc,tc) + [df/dt(x,t)] (t-tc) + [df/dx(x,t)](x-c)
           // in mid(f(xc,tc)) + err(f(xc,tc))
           //             + mid(df/dt) (t-tc) + mid(df/dx)(x-c)
           //             + |err(df/dt)| |rad(time)| + |err(df/dx)| |rad(x)|
//           std::cout << "tdiff : " << center_tdiff << "\ntdiffuncert : " << tdiff << "\n";
           IntervalVector uncert = jac*(approxbox-center_Approxbox);
           uncert += timeslice.rad() * tdiff;
           IntervalVector fun_eval = 
			this->teval_vector(timeslice.mid(),center_Approxbox);
           ibex::Vector fun_evalc = fun_eval.mid();
           uncert += (fun_eval-fun_evalc);
  
           // FIXME : most probably useless
#if 0
           ibex::IntervalVector u2 = fun->eval_vector(approxbox)-
		fun->eval_vector(center_Approxbox) -
		(approxbox-center_Approxbox)*center_jac;
           uncert &= u2;
#endif
  
           IntervalVector tauCent(dim);
           IntervalVector tauXbox(dim);
           IntervalMatrix tauprodExpM(dim,dim);
           this->compute_next_step(center_Approxbox,
		fun_evalc,center_tdiff,center_jac,
		uncert, tsteps, Cent,prodExpM, tauCent, tauprodExpM);
           tauXbox = Xbox + 
		Unt*(tauprodExpM*(Xbox-center_Approxbox) + tauCent); // FIXME
           Xbox = Cent + prodExpM * X0c ;
           if (!tauXbox.is_subset(approxbox)) {
              std::cerr << "Erreur Xbox//tauXbox//approxbox !!!\n" <<
                  Xbox << "  " << tauXbox << "  " << approxbox << "\n\n";
           } 
           std::cerr << "time : " << currenttime << " box : " << Xbox << "\n";
           Result.set(tauXbox,timeslice);
           Result.set(Xbox,currenttime);
           oldtime=currenttime;
       }
       return Result;
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
