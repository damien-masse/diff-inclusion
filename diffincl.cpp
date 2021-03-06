///////////////////////////////////////////////////////////////////////////
//  diffincl.cpp : exponential of interval matrices 
///////////////////////////////////////////////////////////////////////////


#include <codac.h>
#include <iostream>
#include <vector>
#include <array>
#include <cstdio>
#include <cstring>
#include "expIMat.h"
#include "vibes.h"
#include "diffincl.h"
#include "IVdouble.h"
#include "IVparals.h"

#define USE_IVPARALS	1

using namespace codac;

namespace diffincl {

  DiffInclusion::DiffInclusion(unsigned int dim,
                      const TFunction* tfun,
                      Interval& time, unsigned int nbsteps) :
       dim(dim), tfun(tfun), fun(NULL),
       time(time), nbsteps(nbsteps) {
       this->time_dependent=true;
  }
  
  DiffInclusion::DiffInclusion(unsigned int dim,
                      const Function* fun, 
                      Interval& time, unsigned int nbsteps) :
       dim(dim), fun(fun), tfun(NULL),
       time(time), nbsteps(nbsteps) {
       this->time_dependent=false;
  }

  DiffInclusion::DiffInclusion(std::istream& input) {
       char str[256], *token;
       bool fini = false;
       this->fun = NULL;
       this->tfun=NULL;
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
  const IntervalVector DiffInclusion::teval_vector(const Interval &tim,
                 const IVparals& cdom) {
     return this->teval_vector(tim, cdom.bounding_box());
  }
  const IntervalVector DiffInclusion::teval_vector(double tim,
                 const IntervalVector& cdom) {
        const Interval timI(tim);
        return this->teval_vector(timI,cdom);
  }
  const IntervalVector DiffInclusion::teval_vector(double tim,
                 const IVparals& cdom) {
        const Interval timI(tim);
        return this->teval_vector(timI,cdom.bounding_box());
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
  IntervalMatrix DiffInclusion::jacobian(const IVparals& codom,
                 const Interval& tdom,
                 IntervalVector& tvec) {
        return this->jacobian(codom.bounding_box(), tdom, tvec);
  }
  
  // basic enclosing for evolution
  IVparals DiffInclusion::extend_box_basic(const IntervalVector& frame,
                 	const IVparals& startIV,
                 	const Interval& tim, 
			double inflation_factor,
			int nb_tries) {
      double tstep = tim.diam();
      Interval btime(0.0,tstep);
      /* estimation des pentes */
      IntervalVector k1 = this->teval_vector(tim.lb(),startIV);
      IVparals XB1 = sum_tau(startIV,(tstep/2.0)*k1);
      IntervalVector k2 = this->teval_vector(tim.mid(), XB1);
      XB1 = sum_tau(XB1,(tstep/2.0)* k2);
      IntervalVector k3 = this->teval_vector(tim.ub(), XB1);
      /* compute approximation of the result */
      IVparals Res = sum_tau(startIV,
			(tstep*0.505)* ((k1|k2) + (k2|k3)));

//      cout << Res << "\n" << tim << "\n";
      Res&=frame;
//    cout << Res << "   ";
      if (nb_tries==0) return Res;
      /* inflation */
//      Res.inflate(inflation_factor);
//      Res &= frame;
//      if (nb_tries<3) {
        double ifact = nb_tries==0 ? 1.0 : inflation_factor;
        if (nb_tries==2) ifact *= inflation_factor;
        Res = sum_tau(startIV, (ifact * btime) * (this->teval_vector(tim,Res)));
        Res&=frame;
 //     } 
      return Res;
  }
  
        
  bool DiffInclusion::compute_step(const IntervalVector& frame,
				const IVparals& actState,
				IVparals& tauState,
				IVparals& finState,
				const Interval& timeslice,
				vector<IVparals> *outDoors) {

//       cout << timeslice << " " << tauState << "\n";
       IntervalVector tdiff(dim);
       IntervalMatrix jacM = this->jacobian(tauState,timeslice,tdiff);
       const Matrix jac = jacM.mid();
       jacM -= jac;
//       cout << "jac " << jac << "\njacM " << jacM << "\n";

       IntervalMatrix ExpM(dim,dim);
       IntervalMatrix invExpM(dim,dim);
       IntervalMatrix tauExpM(dim,dim);
       IntervalMatrix IExpM(dim,dim);
       IntervalMatrix tauIExpM(dim,dim);
       IntervalMatrix VExpM(dim,dim);
       IntervalMatrix tauVExpM(dim,dim);
       Matrix IntAbs(dim,dim);

       double tsteps = timeslice.diam();

       global_exp(jac,tsteps,true,time_dependent,
	   	ExpM,invExpM,tauExpM,IExpM,tauIExpM,
		VExpM,tauVExpM,IntAbs);
        
       /* other variables which needs to be kept */
       Vector cent_tdiff(dim);


       bool ok=false;
       bool reducing=false;
       int nb_red=0;
       int nb_enl=0;
       while(!ok) {
           /* computing uncert and other values */
           Vector cent_tauState = tauState.mid();
           IVparals ctauState = tauState - cent_tauState;
           IntervalVector uncert = jacM * ctauState;
           if (time_dependent) {
              cent_tdiff = tdiff.mid();
              IntervalVector ctdiff = tdiff-cent_tdiff;
              uncert += timeslice.rad() * ctdiff;
           }
           IntervalVector fun_eval = 
	  	  this->teval_vector(timeslice.mid(),cent_tauState);
           Vector fun_evalc = fun_eval.mid();
           uncert += (fun_eval-fun_evalc);
           
//           cout << "uncert : " << uncert << "\n";
           /* now computing the new "tau" states */
           Vector vuncert = IntAbs * uncert.ub();
           // ``derivation'' of the center ( f(cent) * int exp(M\tau) )
           IntervalVector tauCent = tauIExpM * fun_evalc;
           // non-autonomous factor
           if (time_dependent) {
                 tauCent += tauVExpM * cent_tdiff;
           }
           /* evolution of the center */
           for (int i=0;i<dim;i=i+1) (tauCent[i]).inflate(vuncert[i]); 
          
           /* next tauState, do NOT change the underlying matrices */
           IVparals ntauState(actState);
//           cout << "avant ctau" << ntauState << "\ncent " << cent_tauState <<
//			"\ntauExpM " << tauExpM << "\ntauCent " << tauCent << "\n";
           ntauState.ctau_mult_and_add(cent_tauState,tauExpM,tauCent);

           if (ntauState.is_subsetFast(tauState)) {
	      reducing=true;
//              cout << "plus petit " << nb_red << " " << ntauState << "\t" << tauState << "\n";
              if (tauState.rel_distanceFast(ntauState)<1e-5) ok=true;
				/* FIXME : select a value for rel_distance */
				/* maybe dependent on the computation */
	      tauState=ntauState;
              nb_red++; 
           } else {
              if (reducing) {
                  if (outDoors!=NULL) {
                    for (int i=0;i<dim;i++) {
	  	        (*outDoors)[2*i].join_intersect_with_tau(actState,
		      	      cent_tauState, tauExpM,tauCent,
			      frame,i,frame[i].lb());
		        (*outDoors)[2*i+1].join_intersect_with_tau(actState,
			      cent_tauState, tauExpM,tauCent,
			      frame,i,frame[i].ub());
                    }  
		  }
                  return true; /* return last computation */
              }
              /* problem of computation... should we compute from a larger
                 state? */
//             cout << "plus gros " << nb_enl << " " << ntauState << "\t" << tauState << "\n";
              nb_enl++; 
              if (nb_enl==narrow_factor) return false; 	
			/* FIXME : narrow_factor is arbitrary */
              tauState.inflate_from_baseFast(ntauState,default_inflation_factor);
           }
           if (reducing) {
	      /* compute doors */
              if (ok && outDoors!=NULL) {
                for (int i=0;i<dim;i++) {
	  	    (*outDoors)[2*i].join_intersect_with_tau(actState,
		    	    cent_tauState, tauExpM,tauCent,
			    frame,i,frame[i].lb());
		    (*outDoors)[2*i+1].join_intersect_with_tau(actState,
			    cent_tauState, tauExpM,tauCent,
			    frame,i,frame[i].ub());
                }  
              }
              /* computing ``arrival'' states */
              IntervalVector evolCenter = IExpM * fun_evalc;
              // non-autonomous factor
              if (time_dependent) {
                  evolCenter += VExpM * cent_tdiff;
              }
              for (int i=0;i<dim;i=i+1) (evolCenter[i]).inflate(vuncert[i]);
	      finState = actState;
	      finState.cmult_and_add(cent_tauState,ExpM,invExpM,evolCenter);
           }
           if (!ok) {
	 	/* recompute jacM and tdiff */
 		jacM = this->jacobian(tauState,timeslice,tdiff);
                jacM = jacM - jac; /* maybe not centered */
           } 
	}
        return true;
  }
        

  // compute next step, from linear approximation and uncertainty
  // fcent = f(center of the zone, half time)
  // tVect = mid of temporal differential (dfun/dt (mid_t))
  // jac : mid of jacobian (dfun/dx (zone))
  // coord : actual abstract domain representation
  // invcoord : its inverse
  // ExpM (out) : exponential (not cumulated)
  // evolCenter (out) : evolution of the ``center'' 
  // tauProdExpM (out) : intermediate exponential/[0,1] (!not! cumulated)
  // tauCent (out) : evolution of the center / [0,1] (with uncertainty)
  // uncert : U value (uncertity on f)
  // double tsteps : time step
  // f(x,t) \in fcent + (t-tmid)*tvect + jac*(x-x0) + uncert
  void DiffInclusion::compute_next_step(const Vector& fcent,
                         const Vector& tVect,
                         const Matrix& jac,
                         const IntervalVector& uncert,
                         const IntervalMatrix* coord,
			 const IntervalMatrix* invcoord,
                         double tsteps,
  	                 IntervalVector& evolCenter, 
                         IntervalMatrix& ExpM,
                         IntervalMatrix& invExpM,
			 IntervalVector& tauCent,
			 IntervalMatrix& tauProdExpM) {
     // global computation of the exponentials and their integral
     IntervalMatrix IExpM(dim,dim);
     IntervalMatrix VExpM(dim,dim);
     IntervalMatrix tauExpM(dim,dim);
     IntervalMatrix tauIExpM(dim,dim);
     IntervalMatrix tauVExpM(dim,dim);
     Matrix IntAbs(dim,dim);
     const Matrix Id = Matrix::eye(dim);

       // with coord and invcoord
         IntervalMatrix M(dim,dim);
         IntervalMatrix Minv(dim,dim);
     // computation of exponentials and other values
     if (coord!=NULL) {
         IntervalMatrix dummy1(dim,dim);
         IntervalMatrix dummy2(dim,dim);
         IntervalMatrix dummy3(dim,dim);
         // first computation, with the matrix and its inverse
         global_exp(jac,tsteps,false,false,
		ExpM,invExpM,tauExpM,IExpM,tauIExpM,
		VExpM,tauVExpM,IntAbs);
//	 std::cout << "ExpM " << ExpM << "\ninvExpM " << invExpM << "\nproduct " << (ExpM*invExpM) << "\n"; 
         // new computation with jac = invcoord * jac * coord
         M = ExpM*(*coord);
         Minv = (*invcoord)*invExpM;
//	 std::cout << "M " << M << "\nMinv " << Minv << "\nproduct " << (M*Minv) << "\n"; 
         IntervalMatrix modifiedJac = 
		Minv * jac * M;
//	 std::cout << "jac " << jac << "\nmodified jac " << modifiedJac << "\n"; 
         global_exp(modifiedJac,tsteps,true,time_dependent,
		dummy1,dummy2,dummy3,IExpM,tauIExpM,
		VExpM,tauVExpM,IntAbs);
     } else {
         global_exp(jac,tsteps,true,time_dependent,
	   	ExpM,invExpM,tauExpM,IExpM,tauIExpM,
		VExpM,tauVExpM,IntAbs);
     }
     if (debug_level>5)
	std::cout << "\n\nExpM" << ExpM << "\ntauExpM" 
		  << tauExpM << "\nIExpM" << IExpM << "\ntauIExpM" 
		  << tauIExpM << "\nVExpM" << VExpM << "\ntauVExpM" 
		  << tauVExpM << "\nIntAbs" << IntAbs << "\n";

     ////// ``intermediates'' values  (between 0 and tsteps) 

     if (coord==NULL) {
       // uncertainty
       Vector vuncert = IntAbs * uncert.ub();
       // ``derivation'' of the center ( f(cent) * int exp(M\tau) )
       tauCent = tauIExpM * fcent;
       // non-autonomous factor
       if (time_dependent) {
             tauCent += tauVExpM * tVect;
       }
       if (debug_level>5)
          std::cout << "derivation without uncertainty : " << tauCent << "\n";
       // adding the uncertainty to evolCenter
       for (int i=0;i<dim;i=i+1) (tauCent[i]).inflate(vuncert[i]); 
       // new product
       tauProdExpM = tauExpM;
  
       ///// ``final'' values (at tsteps)
       // ``derivation'' of the center ( f(cent) * int exp(M\tau) )
       evolCenter = IExpM * fcent;
       // non-autonomous factor
       if (time_dependent) {
             evolCenter += VExpM * tVect;
       }
       for (int i=0;i<dim;i=i+1) (evolCenter[i]).inflate(vuncert[i]);
    } else {

       // uncertainty
       Vector vuncert = IntAbs * (Minv * uncert).ub();
       // ``derivation'' of the center ( f(cent) * int exp(M\tau) )
       tauCent = tauIExpM * (Minv * fcent);
       // non-autonomous factor
       if (time_dependent) {
             tauCent += tauVExpM * (Minv * tVect);
       }
       if (debug_level>5)
          std::cout << "derivation without uncertainty : " << tauCent << "\n";
       // adding the uncertainty to evolCenter
       for (int i=0;i<dim;i=i+1) (tauCent[i]).inflate(vuncert[i]); 
       // new product
       tauProdExpM = tauExpM;
  
       ///// ``final'' values (at tsteps)
       // ``derivation'' of the center ( f(cent) * int exp(M\tau) )
       evolCenter = IExpM * (Minv * fcent);
       // non-autonomous factor
       if (time_dependent) {
             evolCenter += VExpM * (Minv * tVect);
       }
       for (int i=0;i<dim;i=i+1) (evolCenter[i]).inflate(vuncert[i]);
    }
}

  ///// possible step 
  void DiffInclusion::fwd_step(const IntervalVector& constraint,
		const IntervalVector& approxbox, 
                const Vector& center_Approxbox,
		const Interval& timeslice,
                const IntervalMatrix* coord,
		const IntervalMatrix* invcoord,
                double tsteps,
		IntervalMatrix& ExpM,
		IntervalMatrix& invExpM,
		IntervalMatrix& tauExpM,
		IntervalVector& evolCenter,
		IntervalVector& tauCent) {
        /* computation of the Jacobian , and temporal differentiation */
        const Interval Unt(0.0,1.0);
        IntervalVector tdiff(dim);
        IntervalMatrix jac = this->jacobian(approxbox,timeslice,tdiff);
        Matrix center_jac = jac.mid();
        jac = jac - center_jac;
        Vector center_tdiff = tdiff.mid();
        tdiff = tdiff - center_tdiff;

        // compute cumulated uncertainty on f (parallel linearization) :
        // f(x,t) in f(xc,tc) + [df/dt(x,t)] (t-tc) + [df/dx(x,t)](x-c)
        // in mid(f(xc,tc)) + err(f(xc,tc))
        //             + mid(df/dt) (t-tc) + mid(df/dx)(x-c)
        //             + |err(df/dt)| |rad(time)| + |err(df/dx)| |rad(x)|
        IntervalVector uncert = jac*(approxbox-center_Approxbox);
        uncert += timeslice.rad() * tdiff;
        IntervalVector fun_eval = 
		this->teval_vector(timeslice.mid(),center_Approxbox);
        ibex::Vector fun_evalc = fun_eval.mid();
        uncert += (fun_eval-fun_evalc);
  
        // most probably useless
#if 0
        ibex::IntervalVector u2 = fun->eval_vector(approxbox)-
  	   fun->eval_vector(center_Approxbox) -
	       (approxbox-center_Approxbox)*center_jac;
        uncert &= u2;
#endif
  
        this->compute_next_step(fun_evalc,center_tdiff,center_jac,
	    uncert, coord, invcoord,
	    tsteps, evolCenter, ExpM, invExpM, tauCent, tauExpM);

#if 0
        if (coord==NULL) {
	// imprecise computation of intermediate box from
        // Xbox (instead of Cent + prodExpM X0c)
	// probably not useful to do better (?)
           tauXbox = Xbox + 
		Unt*(tauprodExpM*(Xbox-center_Approxbox) + tauCent);
           tauXbox &= constraint;
        } else {
           IntervalMatrix M = ExpM * (*coord);
           IntervalMatrix Minv = (*invcoord) * invExpM;
        // very imprecise computation ?
//           tauXbox = Xbox + Unt*
//		(tauprodExpM*(Xbox-center_Approxbox) + M*tauCent);
//          tauXbox &= constraint;
 	   tauXbox = Xbox;
        }
#endif
  }

  IVparals DiffInclusion::fwd_inclusion(const IntervalVector& frame, 
                                const IVparals& X0,
				TubeVector& Result,
				const Interval& ptime,
				unsigned int pnbsteps,
                               VIBesFig *fig,
                               const Matrix *pr,
                               int nbsl) {

       assert(fun != NULL || tfun != NULL);
       const double cst_tsteps = ptime.diam()/pnbsteps;
       double fig_step=0;
       double next_fig=ptime.lb()-1.0;
       if (nbsl>0) {
           fig_step = ptime.diam()/nbsl;
           next_fig = ptime.lb()+fig_step;
       }
       double tsteps = cst_tsteps;
       double inflation_f = default_inflation_factor;
       int nb_tries = 0;

/*
#ifndef USE_IVPARALS
       IVdouble actStates(X0); // actuel set of states
#else
*/
       IVparals actStates(X0); // actuel set of states
       IVparals nextStates(X0); // actuel set of states
/*
#endif
*/
//       IntervalVector Cent(X0.mid());
//       const IntervalVector X0c = X0-Cent; // centered "initial" place
       const Matrix Id = Matrix::eye(dim);
//       IntervalMatrix prodExpM(Id); // products of exponential
       IntervalVector Xbox(actStates.bounding_box());
       Result.sample(ptime.lb(),Xbox);
       double oldtime = ptime.lb();
       if (fig) {
//             std::cout << "drawing " << *pr << actStates << "\n";
             fig->draw_polygon(actStates.over_polygon(*pr),"blue");
       }
       bool okreduced=false; 
       IVparals approxIV(X0);
       while (oldtime<ptime.ub()) {
           double currenttime=oldtime+tsteps;
	/// Reinitialisation (FIXME : customize)
/*
#ifndef USE_IVPARALS
	   if ((int) (oldtime*7/time.ub()) < (int) (currenttime*7/time.ub()))
                   actStates = IVdouble(Xbox);
#endif
*/
           if (currenttime>ptime.ub()) currenttime=ptime.ub();
           Interval timeslice(oldtime,currenttime);
           if (!okreduced) {
              approxIV = 
  		this->extend_box_basic(frame,actStates,timeslice,inflation_f,
				nb_tries);
           }
	   bool safe = compute_step(frame,
		actStates, approxIV, nextStates, timeslice,NULL);
           if (!safe) {
              if (debug_level>=5)
                std::cerr << "Erreur tauXbox//approxbox !!!\n";
//                    << bbtauXbox 
//		    << "  " << approxbox << " tsteps " <<
//                      tsteps <<  " inf_f " << inflation_f << "\n\n";
	      // trying a smaller time slice
              if (nb_tries>2) {
	        // tsteps/=1.4;
                inflation_f *= 1.5;
              }
              nb_tries++;
              continue;
           } 
           if (safe && nb_tries<5) { /* FIXME : narrowing strategy */
               okreduced=true;
               nb_tries++;
               continue;
           }

	   actStates=nextStates;
           oldtime=currenttime;
	   // Computation of Xbox (for the tube)
           Xbox = actStates.bounding_box();
	   Xbox &= frame;
           if (debug_level>=5) 
              std::cerr << "time : " << currenttime 
		        << " box : " << Xbox << "\n";
           Result.set(approxIV.bounding_box(),timeslice);
           Result.set(Xbox,currenttime);
           // resetting the time step 
	   tsteps=cst_tsteps;
	   
 	   // display
	   if (currenttime>= next_fig) {
               next_fig+=fig_step;
               if (fig) {
//                  std::cout << "drawing " << *pr << actStates << "\n";
                  fig->draw_polygon(actStates.over_polygon(*pr),"blue");
	       }
               std::cout << "Time : " << currenttime << "  " << actStates << "\n";
               actStates.toPointMatrices();
	   }
	   inflation_f=default_inflation_factor;
           okreduced=false;
           nb_tries=0;
       }
       return actStates;
  }
  
  ///// forward differential inclusion X' = f(X,t)
  TubeVector DiffInclusion::fwd_inclusion(const IntervalVector& frame, 
                               const IntervalVector& X0,
                               VIBesFig *fig,
                               const Matrix *pr,
                               int nbsl)
  {
       assert(X0.dim()==dim);
  
#ifndef USE_IVPARALS
       IVdouble actStates(X0); // actuel set of states
#else
       IVparals actStates(X0); // actuel set of states
#endif

       // result
       TubeVector Result(time,frame);
  
       this->fwd_inclusion(frame, actStates,
				Result,
				time,
				nbsteps,fig,pr,nbsl);
       return Result;
  }



  TubeVector DiffInclusion::capd_fwd_inclusion(const IntervalVector& frame, 
                               const IntervalVector& X0) {
      if (this->time_dependent) {
       return CAPD_integrateODE(this->time,*(this->tfun),X0,
		this->time.diam()/this->nbsteps,2,
		this->time.diam()/this->nbsteps);
      } else {
       return CAPD_integrateODE(this->time,*(this->fun),X0,this->time.diam()/this->nbsteps,2,
		this->time.diam()/this->nbsteps);
      }
   }

		 
   vector<IVparals>
		 DiffInclusion::room_fwd_computation(const IntervalVector& box,
                        const IVparals& X0,
                        const Interval& intTime,
                        double timesteps,
                        double timeShift) {
       assert ((fun != NULL) || (tfun != NULL));
       assert(X0.dim()==dim);
       IVparals initParals(dim);
       vector<IVparals> outDoors(2*dim,initParals);
       double tsteps = timesteps;
       double inflation_f = default_inflation_factor;
       int nb_tries;
  
       IVparals actStates(X0); // actuel set of states
       IVparals nextStates(X0); // actuel set of states
//       IntervalVector Cent(X0.mid());
//       const IntervalVector X0c = X0-Cent; // centered "initial" place
       const Matrix Id = Matrix::eye(dim);
//       IntervalMatrix prodExpM(Id); // products of exponential
       IntervalVector Xbox(actStates.bounding_box()); // recomputed box (= actStates.bounding_box())

       double oldtime = 0.0;
       bool okreduced=false; 
       IVparals approxIV(X0);
       while (oldtime<time.ub()) {
           double currenttime=oldtime+tsteps;
	/// Reinitialisation (FIXME : customize)
/*
#ifndef USE_IVPARALS
	   if ((int) (oldtime*7/time.ub()) < (int) (currenttime*7/time.ub()))
                   actStates = IVdouble(Xbox);
#endif
*/
           if (currenttime>time.ub()) currenttime=time.ub();
           Interval timeslice(oldtime+timeShift,currenttime+timeShift);
           if (!okreduced) {
              approxIV = 
  	  	  this->extend_box_basic(box,actStates,timeslice,inflation_f,
				nb_tries);
           }

	   bool safe = compute_step(box,
		actStates, approxIV, nextStates, timeslice, &outDoors);
           if (!safe) {
              if (debug_level>=5)
                std::cerr << "Erreur tauXbox//approxbox !!!\n";
//                    << bbtauXbox 
//		    << "  " << approxbox << " tsteps " <<
//                      tsteps <<  " inf_f " << inflation_f << "\n\n";
	      // trying a smaller time slice
              if (nb_tries>2) {
//	        tsteps/=1.4;
                inflation_f *= 1.5;
              }
              nb_tries++;
              continue;
           } 
/*
           if (safe && nb_tries==0) {
               okreduced=true;
               nb_tries++;
               continue;
           }
*/

	   actStates=nextStates;
           oldtime=currenttime;
 	   actStates &= box;
	   // Computation of Xbox (for the tube)
           Xbox = actStates.bounding_box();
	   Xbox &= box;
           if (debug_level>=5) 
              std::cerr << "time : " << currenttime 
		        << " box : " << Xbox << "\n";
           // resetting the time step 
	   tsteps=timesteps;
	   inflation_f=default_inflation_factor;
           nb_tries=0;
           okreduced=false;
       }
       return outDoors;
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
