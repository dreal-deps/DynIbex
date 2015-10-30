/* ============================================================================
 * D Y N I B E X - Definition of the Lobatto III A Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_LA3_H
#define IBEX_SOL_LA3_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{
  
class solution_j_la3 : public solution_j
{
  public:
	//method to define
    
	IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
	{  
	  return picard_tayl(y0,_ode,ordre);
	}
	
	//the LTE
	Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode, double h)
	{	
	    //computation of the LA3 error
	    Affine2Vector err_aff = _ode->computeLA3derivative_aff(Affine2Vector(y0,true));    
	    err_aff*=(std::pow(h,5) / 120.0);

	    return err_aff;	  
	}
	
	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.2)));
	}
	     
	
	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){    
	  //with LA3 and affine form
	  *box_jnh_aff = remainder_la3(_ode); 
	  return 1;
	};
	
	
	//constructor
	solution_j_la3(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{
	  
	}
	
	
	//destructor
	~solution_j_la3(){

	}
	
	
	
private:
	  
	  //la3 with remainder
	  Affine2Vector remainder_la3(ivp_ode* _ode)
	  {         	    
	    double h=time_j.diam();
	    double tol = 1e-20;
	    
	    IntervalVector k1 = _ode->compute_derivatives(1,*box_j1);   
	    IntervalVector k2(k1);
	    IntervalVector k3(k1);
	    
	    IntervalVector k1_old(k1);
	    IntervalVector k2_old(k2);
	    IntervalVector k3_old(k3);
	    k1 = _ode->compute_derivatives(1,*box_jn);
	    k1_old=k1;

	    do
	    {

	      k2_old=k2;
	      k3_old=k3;
	      
	      
	      k2 &= _ode->compute_derivatives(1,*box_jn + h*((5.0/24.0)*k1 + (1.0/3.0)*k2 + (-1.0/24.0)*k3));
	      
	      k3 &= _ode->compute_derivatives(1,*box_jn + h*((1.0/6.0)*k1 + (2.0/3.0)*k2 + (1.0/6.0)*k3));
	      
	    }while((k2.rel_distance(k2_old) > tol)||(k3.rel_distance(k3_old) > tol));
	    
	    Affine2Vector k1_aff = Affine2Vector(k1,true);
	    Affine2Vector k2_aff = Affine2Vector(k2,true);
	    Affine2Vector k3_aff = Affine2Vector(k3,true);
	    
	    k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff);
	    k2_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((5.0/24.0)*k1_aff + (1.0/3.0)*k2_aff + (-1.0/24.0)*k3_aff));	      
	    k3_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + h*((1.0/6.0)*k1_aff + (2.0/3.0)*k2_aff + (1.0/6.0)*k3_aff));
	    
	      
	    Affine2Vector la3=*box_jn_aff + h*( (1.0/6.0)*k1_aff + (2.0/3.0)*k2_aff + (1.0/6.0)*k3_aff);

	    
	    return la3+*box_err_aff;
	    
	  };
	  
	  

};
}

#endif