/* ============================================================================
 * D Y N I B E X - Definition of the Heun Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_HEUN_H
#define IBEX_SOL_HEUN_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{
  
class solution_j_heun : public solution_j
{
  public:
	
	//method to define
    
	IntervalVector picard(IntervalVector y0, ivp_ode* _ode, int ordre)
	{  
      return picard_euler(y0,_ode);//picard_tayl(y0,_ode,2);
	}
	
	//the LTE
	Affine2Vector LTE(IntervalVector y0,ivp_ode* _ode, double h)
	{
	    Affine2Vector err_aff = _ode->computeHEUNderivative(Affine2Vector(y0,true));    
	    err_aff*=(std::pow(time_j.diam(),3) / 6.0);
	    return err_aff;	  
	}
	
	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.33)));	  
	}
	
	
	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){    
	  //with heun and affine form
	  *box_jnh_aff = remainder_heun(_ode); 
	  return 1;
	};
	
	
	//constructor
	solution_j_heun(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{
	 
	}
	
	
	//destructor
	~solution_j_heun(){

	}
	
	
	
private:
	  
	  //heun with remainder
	  Affine2Vector remainder_heun(ivp_ode* _ode)
	  {         	    
	    double h=time_j.diam();

	    Affine2Vector k1 = _ode->compute_derivatives_aff(1, *box_jn_aff);	    
	    Affine2Vector k2 = _ode->compute_derivatives_aff(1, *box_jn_aff+h*k1);
      
	    Affine2Vector heun = k1+k2;
	    heun*=(h/2.0);
	    heun+=*box_jn_aff;
	    
	    return heun+*box_err_aff;
	    
	  };
	  
	  

};
}

#endif
