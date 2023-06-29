///////////////////////////////////////////////////////////////////
//  PW_affine.cpp : piecewise affine function (and interval?)
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>
#include "PW_affine.h"
//#include <pair>
//#include <codac.h>

using namespace codac;


class PW_affine
{
      PW_affine::PW_affine():
	nb_points(1), points(1,0.0), values(1,0.0), slopes(2,0.0) { }
      private PW_affine::PW_affine(int nbp):
	nb_points(0), points(), values(), slopes() { 
          points.reserve(nbp); values.reserve(nbp); slopes.reserve(nbp+1);
      }
      PW_affine::PW_affine(double a, double b):
	nb_points(1), points(1,0.0), values(1,b), slopes(2,a) { }
      PW_affine::PW_affine(const PW_affiche& pw):
        nb_points(pw.nb_points), points(pw.points), values(pw.values),
	slopes(pw.slopes) { }

      /** utility **/
      void PW_affine::simplify() {
         int i=0;
         while (i<nb_points+1) {
           if (slopes[i]==slopes[i+1] && nb_points>1) {
              points.erase[i];
              values.erase[i];
              slopes.erase[i];
              nb_points--:
           } else {
              i=i+1;
           }
         }
      }
  
      /***** access *****/
      int PW_affine::nb_points() const { return nb_points; }

      double PW_affine::value(double a) const {
         int i=0;
         while (i<nb_points && a<=points[i]) i++;
         if (i==0) {
            return slopes[0]*(a-points[0])+values[0]; /* FIXME */
         } else {
            return slopes[i]*(a-points[i-1])+values[i-1]; /* FIXME */
         }
      }
      double PW_affine::min() const {
         if (slopes[0]>0.0 || slopes[nb_points]<0.0) return NEG_INFINITY;
         double a = values[0];
         for (int i=1;i<nb_points;i++) 
             if (values[i]<a) a=values[i];
         }
         return a;
      }
      double PW_affine::max() const {
         if (slopes[0]<0.0 || slopes[nb_points]>0.0) return POS_INFINITY;
         double a = values[0];
         for (int i=1;i<nb_points;i++) 
             if (values[i]>a) a=values[i];
         }
         return a;
      }

      /***** operators ****/
      PW_affine operator-(const PW_affine &pw) {
         PW_affine res(*this);
         for (int i=0;i<res.nb_points;i++) {
             res.values[i]=-res.values[i];
         }
         for (int i=0;i<=res.nb_points;i++) {
             res.slopes[i]=-res.slopes[i];
         }
         return res;
      }
      PW_affine operator+(const PW_affine &pw, double x) {
         PW_affine res(*this);
         for (int i=0;i<res.nb_points;i++) {
             res.values[i]+=x;  /* FIXME */
         }
         return res;
      }
      PW_affine operator-(const PW_affine &pw, double x) {
         PW_affine res(*this);
         for (int i=0;i<res.nb_points;i++) {
             res.values[i]-=x;  /* FIXME */
         }
         return res;
      }
      PW_affine operator*(const PW_affine &pw, double x) {
         if (x==0.0) {
            PW_affine res(0.0);
            return res;
         }
         PW_affine res(*this);
         for (int i=0;i<res.nb_points;i++) {
             res.values[i]*=x;  /* FIXME */
         }
         for (int i=0;i<=res.nb_points;i++) {
             res.slopes[i]*=x;  /* FIXME */
         }
         return res;
      }
      PW_affine operator/(const PW_affine &pw, double x) {
         assert(x!=0.0);
         PW_affine res(*this);
         for (int i=0;i<res.nb_points;i++) {
             res.values[i]/=x;  /* FIXME */
         }
         for (int i=0;i<=res.nb_points;i++) {
             res.slopes[i]/=x;  /* FIXME */
         }
         return res;
      }
      PW_affine operator+(const PW_affine &pw1, const PW_affine &pw2) {
         PW_affine res(pw1.nbpoints+pw2.nbpoints);
         int i=0; int j=0; int k=0;
         res.slopes.push_back(pw1.slopes[0]+pw2.slopes[0]);
         while (i<pw1.nb_points || j<pw2.nb_points) {
            if (i==pw1.nb_points) {
               res.points.push_back(pw2.points[j]);
               res.values.push_back(pw2.values[j]+
		     (pw1.slopes[i]*(pw2.points[j]-pw1.points[i-1])
			   +pw1.values[i-1]));
	       j++;
            } else if (j==pw2.nb_points) {
               res.points.push_back(pw2.points[i]);
               res.values.push_back(pw1.values[i]+
		     (pw2.slopes[j]*(pw1.points[i]-pw2.points[j-1])
			   +pw2.values[j-1]));
	       i++;
            } else {
               if (pw1.points[i]<=pw2.points[j]) {
                   res.points.push_back(pw1.points[i]);
                   res.values.push_back(pw1.values[i]+
		(pw2.slopes[j]*(pw1.points[i]-pw2.points[j])+pw2.values[j]));
                   if (pw1.points[i]==pw2.points[j]) j++;
                   i++;
               } else {
                   res.points.push_back(pw2.points[j]);
                   res.values.push_back(pw2.values[j]+
		(pw1.slopes[i]*(pw2.points[j]-pw1.points[i])+pw1.values[i]));
		   j++;
               }
            }
            res.nbpoints++;
            res.slopes.push_back(pw1.slopes[i]+pw2.slopes[j]);
         }
         return res;
      }
      PW_affine operator-(const PW_affine &pw1, const PW_affine &pw2) {
         PW_affine res(pw1.nbpoints+pw2.nbpoints);
         int i=0; int j=0;
         res.slopes.push_back(pw1.slopes[0]-pw2.slopes[0]);
         while (i<pw1.nb_points || j<pw2.nb_points) {
            if (i==pw1.nb_points) {
               res.points.push_back(pw2.points[j]);
               res.values.push_back(-pw2.values[j]+
		     (pw1.slopes[i]*(pw2.points[j]-pw1.points[i-1])
			   +pw1.values[i-1]));
	       j++;
            } else if (j==pw2.nb_points) {
               res.points.push_back(pw2.points[i]);
               res.values.push_back(pw1.values[i]-
		     (pw2.slopes[j]*(pw1.points[i]-pw2.points[j-1])
			   +pw2.values[j-1]));
	       i++;
            } else {
               if (pw1.points[i]<=pw2.points[j]) {
                   res.points.push_back(pw1.points[i]);
                   res.values.push_back(pw1.values[i]-
		(pw2.slopes[j]*(pw1.points[i]-pw2.points[j])+pw2.values[j]));
                   if (pw1.points[i]==pw2.points[j]) j++;
                   i++;
               } else {
                   res.points.push_back(pw2.points[j]);
                   res.values.push_back(-pw2.values[j]+
		(pw1.slopes[i]*(pw2.points[j]-pw1.points[i])+pw1.values[i]));
		   j++;
               }
            }
            res.nbpoints++;
            res.slopes.push_back(pw1.slopes[i]-pw2.slopes[j]);
         }
         return res;
      }

      PW_affine& PW_affine::operator=(const PW_affine &pw) {
	  this->nb_points=pw.nb_points;
          this->values = pw.values;
          this->points = pw.points;
          this->slopes = pw.slopes.
          return *this;
      }
      PW_affine& PW_affine::operator=(double a) {
          PW_affine res(a);
          this->nb_points=1;
          this->values.swap(a.values);
          this->points.swap(a.points);
          this->slopes.swap(a.slopes);
          return *this;
      }
      PW_affine& PW_affine::operator+=(double x) {
         for (int i=0;i<this->nb_points;i++) {
             this->values[i]+=x;
         }
         return *this;
      }
      PW_affine& PW_affine::operator-=(double x) {
         for (int i=0;i<this->nb_points;i++) {
             this->values[i]-=x;
         }
         return *this;
      }

      PW_affine& PW_affine::operator*=(double x) {
         if (x==0.0) {
            *this = 0.0;
            return *this;
         } else {
            for (int i=0;i<this->nb_points;i++) {
                this->values[i]*=x;
            }
            for (int i=0;i<=this->nb_points;i++) {
                this->slopes[i]*=x;
            }
         }
         return *this;
      }
      PW_affine& PW_affine::operator/=(double x) {
         assert(x!=0.0);
         for (int i=0;i<this->nb_points;i++) {
             this->values[i]/=x;
         }
         for (int i=0;i<=this->nb_points;i++) {
             this->slopes[i]/=x;
         }
      }

      PW_affine& PW_affine::operator+=(const PW_affine &pw) {
          std::vector<double>::iterator itP=this->points.begin();
          std::vector<double>::iterator itV=this->values.begin();
          std::vector<double>::iterator itS=this->slopes.begin();
          int j=0;
          double lastV = *itV;
          double lastP = *itP;
          while (itP!=this->points.end() && j<pw2.nb_points) {
            if (itP==this->points.end()) {
               itP = this->points.insert(itP,pw2.points[j]);
               itV = this->values.insert(itV,pw2.values[j]+
		     ((*itS)*(pw2.points[j]-lastP)+lastV));
               itP++; itV++;
               itS = this->slopes.insert(itS,(*itS)+slopes[j]);
               itS++;
	       j++;
               this->nb_points++;
            } else if (j==pw2.nb_points) {
               *(itS++) += slopes[j];
               lastV = *itV;
               lastP = *itP;
               *itV += (pw2.slopes[j]*((*itP)-pw2.points[j-1])
                           +pw2.values[j-1]);
               itV++; itP++;
            } else {
               if (*itP<=pw2.points[j]) {
                   *(itS++) += slopes[j];
		   lastP=*itP; lastV=*itV;
                   *itV += (pw2.slopes[j]*((*itP)-pw2.points[j])
                           +pw2.values[j]);
                   itV++; itP++;
                   if (pw1.points[i]==pw2.points[j]) j++;
               } else {
                   itP = this->points.insert(itP,pw2.points[j]);
                   itV = this->values.insert(itV,pw2.values[j]+
		        (itS*(pw2.points[j]-lastP)+lastV));
                   itP++; itV++;
                   itS = this->slopes.insert(itS,(*itS)+slopes[j]);
                   itS++;
	           j++;
               }
            }
         }
         *itS += slopes[j];
         return res;
      }
      PW_affine& PW_affine::operator-=(const PW_affine &pw) {
          std::vector<double>::iterator itP=this->points.begin();
          std::vector<double>::iterator itV=this->values.begin();
          std::vector<double>::iterator itS=this->slopes.begin();
          int j=0;
          double lastV = *itV;
          double lastP = *itP;
          while (itP!=this->points.end() && j<pw2.nb_points) {
            if (itP==this->points.end()) {
               itP = this->points.insert(itP,pw2.points[j]);
               itV = this->values.insert(itV,-pw2.values[j]+
		     ((*itS)*(pw2.points[j]-lastP)+lastV));
               itP++; itV++;
               itS = this->slopes.insert(itS,(*itS)-slopes[j]);
               itS++;
	       j++;
               this->nb_points++;
            } else if (j==pw2.nb_points) {
               *(itS++) -= slopes[j];
               lastV = *itV;
               lastP = *itP;
               *itV -= (pw2.slopes[j]*((*itP)-pw2.points[j-1])
                           +pw2.values[j-1]);
               itV++; itP++;
            } else {
               if (*itP<=pw2.points[j]) {
                   *(itS++) -= slopes[j];
		   lastP=*itP; lastV=*itV;
                   *itV -= (pw2.slopes[j]*((*itP)-pw2.points[j])
                           +pw2.values[j]);
                   itV++; itP++;
                   if (pw1.points[i]==pw2.points[j]) j++;
               } else {
                   itP = this->points.insert(itP,pw2.points[j]);
                   itV = this->values.insert(itV,pw2.values[j]-
		        (itS*(pw2.points[j]-lastP)+lastV));
                   itP++; itV++;
                   itS = this->slopes.insert(itS,(*itS)-slopes[j]);
                   itS++;
	           j++;
               }
            }
         }
         *itS -= slopes[j];
         return res;

      }

      PW_affine PW_affine::abs() const {
         PW_affine res(this->nb_points);
         res.slopes.push_back(-fabs(this->slopes[0]));
         if (this->slopes[0]*this->values[0]>0) {
            double pos = this->points[0]-this->values[0]/this->slopes[0];
            res.points.push_back(pos);
            res.values.push_back(0.0);
            res.nb_points++;
            res.slopes.push_back(fabs(this->slopes[0]));
         }
         res.points.push_back(this->points[0]);
         res.values.push_back(fabs(this->points[0]));
         res.nb_points++;
         for (int i=1;i<this->nb_points;i++) {
             /* entre point[i-1] (ajouté) et point[i], avec slopes[i] */
             if (this->values[i]*this->values[i-1]<0) {
                double pos = this->points[i]-this->values[i]/this->slopes[i];
                res.points.push_back(pos);
                res.values.push_back(0.0);
                res.nb_points++;
                res.slopes.push_back(-fabs(this->slopes[i]));
                res.slopes.push_back(fabs(this->slopes[i]));
             } else {
                if (this->values[i]<0) 
                    res.slopes.push_back(-this->slopes[i]);
                else if (this->values[i]==0.0) 
                        res.slopes.push_back(-fabs(this->slopes[i]));
                     else res.slopes.push_back(this->slopes[i]);
             }
             res.points.push_back(this->points[i]);
             res.values.push_back(fabs(this->points[i]));
             res.nb_points++;
         }
         int i = this->nb_points;
         if (this->slopes[i]*this->values[i-1]<0) {
            double pos = this->points[i-1]-this->values[i-1]/this->slopes[i];
            res.points.push_back(pos);
            res.values.push_back(0.0);
            res.nb_points++;
            res.slopes.push_back(-fabs(this->slopes[i]));
         }
         res.slopes.push_back(fabs(this->slopes[i]));
         return res;
      }

      PW_affine PW_affine::max(double x) const {
         PW_affine res(this->nb_points);
         double act_pos;
         if (this->slopes[0]>0) {
             res.slopes.push_back(0.0);
             if (this->values[0]>x) {
		double pos = this->points[0]-
				(this->values[0]-x)/this->slopes[0];
                res.points.push_back(pos);
                res.values.push_back(x);
                res.nb_points++;
                res.slopes.push_back(this->slopes[0]);
                act_pos=this->values[0];
             } else act_pos=x;
         } else {
             res.slopes.push_back(this->slopes[0]);
             if (this->values[0]<x && this->slopes[0]<0) {
		double pos = this->points[0]-
				(this->values[0]-x)/this->slopes[0];
                res.points.push_back(pos);
                res.values.push_back(x);
                res.nb_points++;
                res.slopes.push_back(0.0);
                act_pos=x;
             } else act_pos=this->values[0];
         }
         for (int i=1;i<this->nb_points;i++) {
             /* entre point[i-1] (non ajouté) et point[i], avec slopes[i] */
             /* if this->values[x-1]<x, actual slope is 0 */
             if (this->values[i-1]<x) {
		 if (this->values[i]<=x) continue;
                 double pos = this->points[i-1]+
				(x-this->values[i-1])/this->slopes[i];
                 res.points.push_back(pos); res.values.push_back(x); 
		 res.nb_points++;
		 res.slopes.push_back(this->slopes[i]);
             } else { /* si this->values[i-1]==x, on place un point
			 même si il pourrait ne pas être nécessaire */
                 res.points.push_back(this->points[i-1]);
		 res.values.push_back(this->values[i-1]); 
		 res.nb_points++;
                 if (this->values[i]>=x) 
			res.slopes.push_back(this->slopes[i]);
                 else {
		    double pos = this->points[i-1]+
				   (x-this->values[i-1])/this->slopes[i];
		    if (pos>this->points[i-1]) {
                         res.points.push_back(pos);
		         res.values.push_back(x); 
		         res.nb_points++;
		         res.slopes.push_back(this->slopes[i]); 
                    }
		    res.slopes.push_back(0.0); 
	         }
	     }
         }
         int i = this->nb_points;
	 if (this->values[i-1]<x) { /* actual slope is 0 */
	     if (this->slopes[i]<=0.0) {
		/* just check if there is no point */
                if (res.nb_points==0) {
                   res.points.push_back(0.0);
		   res.values.push_back(x);
		   res.nb_points++;
		   res.slopes.push_back(0.0);
		} /* else nothing, just keep slope=0 */
	     } else {
		/* identify growth */
		double pos = this->points[i-1]+
			(x-this->values[i-1])/this->slopes[i];
                res.points.push_back(pos);
		res.values.push_back(x);
		res.nb_points++;
		res.slopes.push_back(this->slopes[i]);
            }
	} else {
            res.points.push_back(this->points[i-1]);
	    res.values.push_back(this->values[i-1]); 
	    res.nb_points++;
            if (this->slopes[i]>=0) res.slopes.push_back(this->slopes[i]); 
	    else {
	       double pos = this->points[i-1]+
		  	 (x-this->values[i-1])/this->slopes[i];
	       if (pos>this->points[i-1]) {
                      res.points.push_back(pos);
	              res.values.push_back(x); 
		      res.nb_points++;
		      res.slopes.push_back(this->slopes[i]); 
               }
	       res.slopes.push_back(0.0); 
            }
         }
         return res;
      }
      
      PW_affine PW_affine::min(double x) const {
         PW_affine res(this->nb_points);
         double act_pos;
         if (this->slopes[0]<0) {
             res.slopes.push_back(0.0);
             if (this->values[0]<x) {
		double pos = this->points[0]-
				(this->values[0]-x)/this->slopes[0];
                res.points.push_back(pos);
                res.values.push_back(x);
                res.nb_points++;
                res.slopes.push_back(this->slopes[0]);
                act_pos=this->values[0];
             } else act_pos=x;
         } else {
             res.slopes.push_back(this->slopes[0]);
             if (this->values[0]>x && this->slopes[0]>0) {
		double pos = this->points[0]-
				(this->values[0]-x)/this->slopes[0];
                res.points.push_back(pos);
                res.values.push_back(x);
                res.nb_points++;
                res.slopes.push_back(0.0);
                act_pos=x;
             } else act_pos=this->values[0];
         }
         for (int i=1;i<this->nb_points;i++) {
             /* entre point[i-1] (non ajouté) et point[i], avec slopes[i] */
             /* if this->values[x-1]>x, actual slope is 0 */
             if (this->values[i-1]>x) {
		 if (this->values[i]>=x) continue;
                 double pos = this->points[i-1]+
				(x-this->values[i-1])/this->slopes[i];
                 res.points.push_back(pos); res.values.push_back(x); 
		 res.nb_points++;
		 res.slopes.push_back(this->slopes[i]);
             } else { /* si this->values[i-1]==x, on place un point
			 même si il pourrait ne pas être nécessaire */
                 res.points.push_back(this->points[i-1]);
		 res.values.push_back(this->values[i-1]); 
		 res.nb_points++;
                 if (this->values[i]<=x) 
			res.slopes.push_back(this->slopes[i]);
                 else {
		    double pos = this->points[i-1]+
				   (x-this->values[i-1])/this->slopes[i];
		    if (pos>this->points[i-1]) {
                         res.points.push_back(pos);
		         res.values.push_back(x); 
		         res.nb_points++;
		         res.slopes.push_back(this->slopes[i]); 
                    }
		    res.slopes.push_back(0.0); 
	         }
	     }
         }
         int i = this->nb_points;
	 if (this->values[i-1]>x) { /* actual slope is 0 */
	     if (this->slopes[i]>=0.0) {
		/* just check if there is no point */
                if (res.nb_points==0) {
                   res.points.push_back(0.0);
		   res.values.push_back(x);
		   res.nb_points++;
		   res.slopes.push_back(0.0);
		} /* else nothing, just keep slope=0 */
	     } else {
		/* identify growth */
		double pos = this->points[i-1]+
			(x-this->values[i-1])/this->slopes[i];
                res.points.push_back(pos);
		res.values.push_back(x);
		res.nb_points++;
		res.slopes.push_back(this->slopes[i]);
            }
	} else {
            res.points.push_back(this->points[i-1]);
	    res.values.push_back(this->values[i-1]); 
	    res.nb_points++;
            if (this->slopes[i]<=0) res.slopes.push_back(this->slopes[i]); 
	    else {
	       double pos = this->points[i-1]+
		  	 (x-this->values[i-1])/this->slopes[i];
	       if (pos>this->points[i-1]) {
                      res.points.push_back(pos);
	              res.values.push_back(x); 
		      res.nb_points++;
		      res.slopes.push_back(this->slopes[i]); 
               }
	       res.slopes.push_back(0.0); 
            }
         }
         return res;
      }

      PW_affine max(const PW_affine &pw1, const PW_affine &pw2) {
         PW_affine res(pw1.nbpoints+pw2.nbpoints);
         int i=0; int j=0; int k=0;
         int act_pos;
         double act_slope;
         if (pw1.slopes[0]<pw2.slopes[0]) {
            act_slope=pw1.slopes[0];
            act_pos=1;
         } else {
            act_slope=pw2.slopes[0];
            act_pos=2;
         }
         res.slopes.push_back(act_slope);

         while (i<pw1.nb_points || j<pw2.nb_points) {
            double v1;
	    double v2;
            double pt;
            bool on1, on2;
            if (j==pw2.nb_points || 
			(i<pw1.nb_points && pw1.points[i]<=pw2.points[j])) {
                v1=pw1.values[i];
		pt=pw1.points[i];
		on1=true;
            } else {
                double a=(i<pw1.nb_point ? i : i-1);
		v1=pw1.slopes[i]*(pw2.points[j]-pw1.points[a])
                                        +pw1.values[a];
		pt=pw2.points[j];
		on1=false;
	    }
            if (i==pw1.nb_points || 
			(j<pw2.nb_points && pw2.points[j]<=pw1.points[i])) {
                v2=pw2.values[j];
		on2=true;
            } else {
                double a=(j<pw2.nb_point ? j : j-1);
		v2=pw2.slopes[j]*(pw1.points[i]-pw2.points[a])
                                        +pw2.values[a];
		on2=false;
	    }
            if (v1>v2 && act_pos==2) {
		
            
               if (pw1.points[i]<=pw2.points[j]) {
                   double v1 = pw1.values[i];
		   double v2 = pw2.slopes[j]*(pw1.points[i]-pw2.points[j])
					+pw2.values[j];
		   if (v1==v2) {
                       res.points.push_back(pw1.points[i]);
	               res.values.push_back(v1); 
		       res.nb_points++;
		       i++; if (pw1.points[i]==pw2.points[j]) j++;
		       if (pw1.slopes[i]>pw2.slopes[j]) {
            		  act_slope=pw1.slopes[i];
            		  act_pos=1;	
		       } else {
            		  act_slope=pw2.slopes[j];
            		  act_pos=2;	
		       }
                       res.slopes.push_back(act_slope);
		   } else 
		   if (v1>v2) {
                      if (act_pos==2) { /* on croise */
			 double pos = 
			     (pw1.values[i]-pw2.values[j]
				-pw1.slopes[i]*pw1.points[i]
				+pw2.slopes[j]*pw2.points[j])/
					(pw2.slopes[j]-pw1.slopes[i]);
 			 if (pos<pw1.points[i]) {
                            double vpos =
				(pw1.values[i]*pw2.slopes[j]
				- pw2.values[j]*pw1.slopes[i]
				+pw1.slopes[i]*pw2.slopes[j]
					*(pw2.points[j]-pw1.points[i]=))/
				(pw2.slopes[j]-pw1.slopes[i]);
                            res.points.push_back(pos);
	                    res.values.push_back(vpos); 
			    res.nb_points++;
			    res.slopes.push_back(pw1.slopes[i]);
                         }
			 act_pos=1;
                      }
		      res.points.push_back(pw1.points[i]);
		      res.values.push_back(pw1.values[i]);
	              res.nb_points++;
		      res.slopes.push_back(pw1.slopes[i+1]);
		    } else if (v1==v2) {
			if (act_pos==1 && pw1.slopes[i+1]>pw1.slop
		    }
			
		      
			    

                   res.points.push_back(pw1.points[i]);
                   res.values.push_back(pw1.values[i]+
		(pw2.slopes[j]*(pw1.points[i]-pw2.points[j])+pw2.values[j]));
                   if (pw1.points[i]==pw2.points[j]) j++;
                   i++;
               } else {
                   res.points.push_back(pw2.points[j]);
                   res.values.push_back(pw2.values[j]+
		(pw1.slopes[i]*(pw2.points[j]-pw1.points[i])+pw1.values[i]));
		   j++;
               }
            }
            res.nbpoints++;
            res.slopes.push_back(pw1.slopes[i]+pw2.slopes[j]);
         }
         return res;
      }
      PW_affine max(const PW_affine& pw) const;
      PW_affine min(const PW_affine& pw) const;

      /** display
       */
      friend std::ostream& operator<<(std::ostream& str, const PW_affine& iv);


    private:
      unsigned int nbPoints;
      std::vector<double> points;
      std::vector<double> values;
      std::vector<double> slopes; /* always 1 more than the points */
      /* naive representation : not really good... is using map or set an
         overkill? */
};


}

#endif
