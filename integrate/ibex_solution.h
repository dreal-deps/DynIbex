/* ============================================================================
 * D Y N I B E X - Definition of the Solution of one simulation step
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_H
#define IBEX_SOL_H

#include <iomanip>
#include <stdlib.h>

const double hmin = 1e-8;
const double hmax = 0.9;

const int ordre_max = 7;

namespace ibex{
  
class solution_j
{
  public:
	IntervalVector* box_j0; //encadrement grossier sur t=[tn,tn+h]
	IntervalVector* box_j1; //encadrement fin sur t=[tn,tn+h]
	Interval time_j; //t=[tn,tn+h]

	IntervalVector* box_jn; //encadrement à tn
	IntervalVector* box_jnh; //encadrement à tn+h
	
	Affine2Vector* box_err_aff;
	
	Affine2Vector* box_jn_aff;
	Affine2Vector* box_jnh_aff;
	
	//accuracy expected
	double atol;
	
	///control of stepsize
	double factor; //factor for the next stepsize : if > 1 the current step was easy with a small error
		       // if < 1 the current stepsize used a small stepsize or the error was too large
	
	//few different picard operator
	//euler
	IntervalVector picard_euler(IntervalVector y0, ivp_ode* _ode)
	{
	  return *box_jn + Interval(0,time_j.diam())*_ode->compute_derivatives_aff(1, Affine2Vector(y0,true)).itv();  
	};
	
	//taylor
	IntervalVector picard_tayl(IntervalVector y0, ivp_ode* _ode, int ordre)
	{
	  double h = time_j.diam();
	  int n=ordre;
	  int fac_i=1;
	  
	  Affine2Vector jn = Affine2Vector(*box_jn,true); //restart of the affine form
	  Affine2Vector int_tayl(jn);
	  for (int i=1;i<n;i++)
	  {
	      fac_i = fac_i*i;
	      Affine2Vector df = _ode->compute_derivatives_aff(i, jn);
	      df*=( 1.0 /fac_i);
	      df*=Interval(0,std::pow(h,i));
	      int_tayl=int_tayl + df; 

	  }
	    
	  IntervalVector err = _ode->compute_derivatives_aff(n, Affine2Vector(y0)).itv(); 
	    
	  fac_i = fac_i*(n);
	  err*=(1.0/fac_i);
	  err*=Interval(0,std::pow(h,n));  
	  
	  return int_tayl.itv()+err;
	};
	
	//rk4
	IntervalVector picard_rk4(IntervalVector y0, ivp_ode* _ode)
	{
	  double h = time_j.diam();
	  Interval time0h = Interval(0,h);  
	  
	  Affine2Vector boxj1 = Affine2Vector(*box_jn,true);
	    
	  Affine2Vector k1 = _ode->compute_derivatives_aff(1, boxj1);
	  
	  Affine2Vector boxj2(k1);
	  boxj2*=(0.5*time0h);
	  
	  Affine2Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj2);
	  
	  Affine2Vector boxj3(k2);
	  boxj3*=(0.5*time0h);
	  
	  Affine2Vector k3 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj3);
	  
	  Affine2Vector boxj4(k3);
	  boxj4*=(time0h);
	  
	  Affine2Vector k4 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj4);
	  
	  k2*=(2.0);
	  k3*=(2.0);
	  
	    
	  Affine2Vector int_rk4 = k1+k2+k3+k4;
	  int_rk4*=(time0h/6.0); 
	  int_rk4+=*box_jn_aff;
	  
	  Affine2Vector err_aff = _ode->computeRK4derivative(Affine2Vector(y0,true));    
	  err_aff*=(pow(time0h,5) / 120);
	    
	  IntervalVector err = err_aff.itv();
	  
	  return int_rk4.itv()+err;
	  
	}
	
	
	//******virtual methods to define to build a new solution scheme********///
	
	//the picard operator
	virtual IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
	{  
	  return picard_tayl(y0,_ode,ordre);
	}
	
	//the LTE
	virtual Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode,double h)
	{
	    Affine2Vector err_aff = _ode->computeRK4derivative(Affine2Vector(y0,true));    
	    err_aff*=(std::pow(h,5) / 120);
	    
	    return err_aff;	  
	}
	
	//the factor for the next stepsize computation
	virtual double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.2)));	  
	}
	
	//compute the sharpest jn+1
	virtual int calcul_jnh(ivp_ode* _ode){    
	   
	  //with RK4 and affine form
	  *box_jnh_aff = remainder_rk4(_ode);   
		  
	  return 1;
	};
	
	
	
	///***algorithm for scheme, general for all scheme***///
	//compute j0 roughly
	int calcul_j0(IntervalVector _box_j, ivp_ode* _ode){   
 
	  int nb=0;
	  //approx of jn+h with rk4
	  IntervalVector yn_p1 = approx_rk4(*box_jn, time_j.diam(), _ode);  	  
	  IntervalVector yn_p0(_box_j);	  	  
	  yn_p0|=yn_p1;
	  //first inflation
	  double infl = yn_p0.diam().max()*0.1;
	  yn_p0.inflate(infl);
	  
	  //picard with euler rectangle rule, with affine evaluation
	  yn_p1 = picard(yn_p0,_ode,3);
	    
	  int iter = 0;
	  while ((!yn_p1.is_subset(yn_p0))&&(iter < 8*_ode->nbvar+1))
	  {
	    iter++;
	    yn_p0=yn_p1;
	    
	    double infl = yn_p0.diam().max()*0.1;
	    //yn_p0.inflate(infl);
	    
	    yn_p1 = picard(yn_p0,_ode,3);
	  }
	  
	  
	  
	  //prune of the step
	  double h;	  
	  if ((yn_p1).is_subset(yn_p0))
	  {
	    //success of picard operator
	    *box_j0 = yn_p1;
	    
	    //contractor based on picard
	    calcul_j1(_ode);
	    yn_p1=*box_j1;
	    
	    Affine2Vector err_aff = LTE(yn_p1,_ode, time_j.diam());
	    IntervalVector err = err_aff.itv();
       
	    double norm_err = infinite_norm(err);
	    
	    //test if truncature error lower than tolerance (or if we have no choice due to the min step reaching)
	    //from Hairer
	    double rtol = atol;
	    double err_test = atol + infinite_norm(yn_p1)*rtol;
	    
	    double test = norm_err/err_test;
	    
	    if ((test <= 1.0)||(time_j.diam()<=hmin))
	    {    	      
	      LOGGER->inc_accepted_picard();
	      LOGGER->update_lte_max(norm_err);
	      LOGGER->update_step(time_j.diam());
	      //computation of the factor for the next step
	      factor = factor_stepsize(test);
	      *box_err_aff = err_aff;
	      return 1;
	    }
	    else
	    {
	      LOGGER->inc_rejected_picard();
	      //truncature error not accepted => step division
	      h = std::max(hmin,time_j.diam()/2.0);
	      time_j = Interval(time_j.lb(),time_j.lb()+h);
	      return 0;
	    }
	  }
	  else
	  {
	    LOGGER->inc_rejected_picard();
	    //picard rejected => step division
	    h = std::max(hmin,time_j.diam()/2.0);
	    time_j = Interval(time_j.lb(),time_j.lb()+h);
	    return 0;
	  }
	};
	
	
	//compute a sharp j1
	int calcul_j1(ivp_ode* _ode){    	  

	  double h = time_j.diam();	    
	  //here : picard accepted on box_j0 => integration rule is now contracting
	  IntervalVector yn_p1(*box_j0);
	  IntervalVector yn_p0(*box_j0);

	  do
	  {
	    yn_p0=yn_p1;
	    
	    //evaluation with affine
	    yn_p1 &= picard(yn_p0,_ode,3);
	  } while(yn_p1.rel_distance(yn_p0) > 1e-18); //till a fix point
	  
	  *box_j1 = yn_p1;
	  
	    
	  return 1;
	};
	
		

	double infinite_norm(const IntervalVector _vec)
	{
	  double res = std::max(std::abs(_vec[0].lb()),std::abs(_vec[0].ub()));
	  for (int i=1;i<_vec.size();i++)
	  {
	   res =  std::max(std::abs(_vec[i].lb()),res);    
	   res =  std::max(std::abs(_vec[i].ub()),res);    
	  }  
	  return res;  
	};
	
	
	void set_atol(double a)
	{
	  atol = a;
	}
	
	//flush affine form to keep memory
	void flush()
	{
	  *box_err_aff = box_err_aff->itv();
	  *box_jn_aff = *box_jn;
	  *box_jnh_aff = *box_jnh;
	};
	
	
	//empty constructor
	solution_j(){};
	
	solution_j(const Affine2Vector _box_jn, double tn, double h, double a)
	{
	    box_jn_aff = new Affine2Vector(_box_jn); 
	    box_jn = new IntervalVector(_box_jn.itv());
	    box_j0 = new IntervalVector(*box_jn);
	    box_j1 = new IntervalVector(*box_jn); 
	    box_err_aff = new Affine2Vector(box_jn->size());
	    box_jnh  = new IntervalVector(*box_jn);
	    
	    //important:
	    time_j = Interval(tn-h,tn);
	    box_jnh_aff = new Affine2Vector(_box_jn);
	    atol = a;
	    factor = 1.0;
	}
	
	//constructor
	solution_j(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac)
	{
	  factor = fac;
  	  atol = a;
	  box_jn_aff = new Affine2Vector(_box_jn); 
	  
	  box_jn = new IntervalVector(_box_jn.itv()); 	  
	  
	  time_j = Interval(tn,tn+std::min(hmax,std::max(hmin,h)));   
	  
	  //compute j0 roughly
	 box_j0 = new IntervalVector(box_jn->size());
	 box_j1 = new IntervalVector(box_jn->size()); 
	 box_err_aff = new Affine2Vector(box_jn->size()); 
	  
	  //compute jnh*/
	  box_jnh = new IntervalVector(_box_jn.size());
	  box_jnh_aff = new Affine2Vector(_box_jn.size());
	  
	}
	
	int compute_oneStep(const Affine2Vector _box_jn, ivp_ode* _ode)
	{
	  
	  //compute j0 roughly
	 int ok=calcul_j0(*box_jn,_ode);
	    
	  int nb=0;
	  while ((ok != 1)&&(nb < 2*_ode->nbvar+1))
	  {
	       IntervalVector _box(*box_j0);
	       ok=calcul_j0(_box, _ode);
	       nb++;
	  }
	 if (ok != 1)
	 {
		cout << "Step:"<<nb<<" refused !" << endl;
		LOGGER->Log("Solution at t=%f : ",time_j.lb());
		LOGGER->Log_sol(_box_jn.itv());
		LOGGER->Log_end();
		//////////**************/////////////
		time_j = Interval(time_j.lb(),1e8);   
		box_j1->init(Interval::ALL_REALS);
		box_jnh->init(Interval::ALL_REALS);
		*box_jnh_aff = *box_jnh;
		return 0;
		//exit(EXIT_FAILURE);
		///////////////**************//////////////
		
	 }
	 else
	 {
		if (box_j0->is_unbounded())
		{
		  cout << "Step:"<<nb<<" refused !" << endl;
		  LOGGER->Log("Solution at t=%f : ",time_j.lb());
		  LOGGER->Log_sol(_box_jn.itv());
		  LOGGER->Log_end();
		 //////////**************/////////////
		time_j = Interval(time_j.lb(),1e8);   
		box_j1->init(Interval::ALL_REALS);
		box_jnh->init(Interval::ALL_REALS);
		*box_jnh_aff = *box_jnh;
		return 0;
		//exit(EXIT_FAILURE);
		///////////////**************//////////////
		}
	 }
	  
	  //compute jnh
	  calcul_jnh(_ode);
	  box_jnh_aff->compact();   ///compactage
	  
	  *box_jnh = box_jnh_aff->itv();  
	  
	  if ((box_jnh->is_unbounded()) || (box_jnh->is_empty()))
	  {
	    cout << "Unacceptable solution!" << endl;
	    LOGGER->Log("Solution at t=%f : ",time_j.lb());
	    LOGGER->Log_sol(_box_jn.itv());
	    LOGGER->Log_end();
	    //////////**************/////////////
		time_j = Interval(time_j.lb(),1e8);   
		box_j1->init(Interval::ALL_REALS);
		box_jnh->init(Interval::ALL_REALS);
		*box_jnh_aff = *box_jnh;
		return 0;
		//exit(EXIT_FAILURE);
		///////////////**************//////////////
	  }  
	  
	}
	
	
	//printer
	void print_soljn(){	  
	    cout << setprecision(20)<< "Solution at t=" << time_j.ub() << " : " << *box_jnh << endl;
	    cout << "affine form : " << *box_jnh_aff << endl;
	}
	
	
	//destructor
	~solution_j(){
	}
	
	
	
    IntervalVector approx_rk4(const IntervalVector yj, double h, ivp_ode* _ode)
	{
	    IntervalVector int_rk4(_ode->nbvar);
	    IntervalVector boxj1(yj);
	    
	    IntervalVector k1 = _ode->compute_derivatives(1, boxj1);
	    
	    IntervalVector boxj2(k1);
	    boxj2*=(0.5*h);
	    
	    IntervalVector k2 = _ode->compute_derivatives(1, yj+boxj2);
	    
	    IntervalVector boxj3(k2);
	    boxj3*=(0.5*h);
	    
	    IntervalVector k3 = _ode->compute_derivatives(1, yj+boxj3);
	    
	    IntervalVector boxj4(k3);
	    boxj4*=(h);
	    
	    IntervalVector k4 = _ode->compute_derivatives(1, yj+boxj4);
	    
	    k2*=(2.0);
	    k3*=(2.0);
	        
	  
	    int_rk4 = k1+k2+k3+k4;
	    int_rk4*=(h/6.0);
	    int_rk4 +=(yj);
	    
	  return int_rk4;

	};

	  
	  
	  //rk4 with remainder
	  Affine2Vector remainder_rk4(ivp_ode* _ode)
	  {
	    double h=time_j.diam();
	    
	    Affine2Vector boxj1(*box_jn_aff);
	    
	    Affine2Vector k1 = _ode->compute_derivatives_aff(1, boxj1);
	    
	    Affine2Vector boxj2(k1);
	    boxj2*=(0.5*h);
	    
	    Affine2Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj2);
	    
	    Affine2Vector boxj3(k2);
	    boxj3*=(0.5*h);
	    
	    Affine2Vector k3 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj3);
	    
	    Affine2Vector boxj4(k3);
	    boxj4*=(h);
	    
	    Affine2Vector k4 = _ode->compute_derivatives_aff(1, *box_jn_aff+boxj4);
	    
	    k2*=(2.0);
	    k3*=(2.0);
	    
	      
	    Affine2Vector int_rk4 = k1+k2+k3+k4;
	    int_rk4*=(h/6.0); 
	    int_rk4+=*box_jn_aff;
	    
	    return int_rk4+*box_err_aff;//df
	    
	  };
	  
	  

};
}

#endif
