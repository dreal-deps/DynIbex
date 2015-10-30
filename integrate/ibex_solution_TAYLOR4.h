/* ============================================================================
 * D Y N I B E X - Definition of the Taylor order 4 Method
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SOL_TAYL4_H
#define IBEX_SOL_TAYL4_H

#include <iomanip>
#include <stdlib.h>
#include "ibex_solution.h"


namespace ibex{
  
class solution_j_tayl4 : public solution_j
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
	    Affine2Vector err_aff = _ode->compute_derivatives_aff(5,Affine2Vector(y0,true));    
	    err_aff*=(std::pow(h,5) / 120);

	    return err_aff;	  
	}
	
	//the factor for the next stepsize computation
	double factor_stepsize(double test)
	{
	  return std::min(1.8,std::max(0.4,0.9*std::pow(1.0/test,0.2)));
	}

	//compute the sharpest jn+1
	int calcul_jnh(ivp_ode* _ode){    
	  //with TAYL4 and affine form
	  *box_jnh_aff = remainder_taylor4(_ode); 
	  return 1;
	};
	
	
	//constructor
	solution_j_tayl4(const Affine2Vector _box_jn, double tn, double h, ivp_ode* _ode,double a, double fac) : solution_j(_box_jn, tn, h, _ode, a, fac)
	{
	  
	  
	}
	
	
	//destructor
	~solution_j_tayl4(){

	}
	
	
	
private:
	  
	  //tayl4 with remainder
	  Affine2Vector remainder_taylor4(ivp_ode* _ode)
	  {         	
	    double h=time_j.diam();
	    double tol = atol*0.001;
	   
	    Affine2Vector tayl4(*box_jn_aff);
	    int n=5;
	    int fac_i=1;
	    for (int i=1;i<n;i++)
	    {
		fac_i = fac_i*i;
		Affine2Vector df = _ode->compute_derivatives_aff(i, *box_jn_aff);
		df*=(1.0 /fac_i);
		df*=pow(h,i);

		tayl4+=df; 

	    }

	    
	    return tayl4+*box_err_aff;
	    
	  };
	  
	  

};
}

#endif