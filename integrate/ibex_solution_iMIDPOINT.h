/* ============================================================================
 * D Y N I B E X - Definition of the Implicit Midpoint Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_IMIDPOINT_H
#define IBEX_SOL_IMIDPOINT_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{
  
class solution_j_imidpoint : public solution_j
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
	    //computation of the ieuler error
	    Affine2Vector err_aff = _ode->computeIMIDPOINTderivative(Affine2Vector(y0,true));    
	    err_aff*=(std::pow(h,3) / 6);

	    return err_aff;	  
	}
	
	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.33)));
	}     
	
	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){    
	  //with imidpoint and affine form
	  *box_jnh_aff = remainder_imidpoint(_ode); 
	  return 1;
	};
	
	
	//constructor
	solution_j_imidpoint(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{
	  	  
	}
	
	
	//destructor
	~solution_j_imidpoint(){

	}
	
	
	
private:
	  
	  //imidpoint with remainder
	  Affine2Vector remainder_imidpoint(ivp_ode* _ode)
	  {         	
	    double h=time_j.diam();
	    double tol = atol*0.001;
	    IntervalVector k1 = _ode->compute_derivatives(1,*box_j1); 
	    IntervalVector k1_old(k1);
	    
	    do
	    {

	      k1_old=k1;
	      k1 &= _ode->compute_derivatives(1,*box_jn + 0.5*h*k1);
	      
	      
	    } while (k1.rel_distance(k1_old) > tol);
	      
	    Affine2Vector k1_aff = Affine2Vector(k1,true);
	    
	    k1_aff = _ode->compute_derivatives_aff(1,*box_jn_aff + 0.5*h*k1_aff);
	    
	    Affine2Vector ieuler = *box_jn_aff + h*k1_aff;
	      
	    return ieuler+*box_err_aff;
	    
	  };
	  
	  

};
}

#endif