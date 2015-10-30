/* ============================================================================
 * D Y N I B E X - Definition of the ODE structure
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_ODE_H
#define IBEX_ODE_H

const int MAX_FCT=60;

namespace ibex{

class ivp_ode
{
public:
  ivp_ode(const Function _ydot, double _t0, const IntervalVector _yinit){
      nbvar = _yinit.size();
   
      yinit = new IntervalVector(_yinit);
      yinit_aff = new Affine2Vector(_yinit,true);
      
      ydot = new Function(_ydot);
      t0 = _t0;
      
  };
  
  ivp_ode(const Function _ydot, double _t0, const Affine2Vector _yinit_aff){
      nbvar = _yinit_aff.size();
   
      yinit = new IntervalVector(_yinit_aff.itv());
      yinit_aff = new Affine2Vector(_yinit_aff);
      
      ydot = new Function(_ydot);
      t0 = _t0;
  };
  
  void frechet_precomputing(int order)
  {
    edtfr = new edtree_frechet(order,ydot,nbvar);    
  }

  
  IntervalVector compute_derivatives(int ordre, IntervalVector yi){
    assert(ordre < MAX_FCT);
    
    if (ordre == 1)
    {
        IntervalMatrix J(yi.size(), yi.size());
	ydot->hansen_matrix(yi,J);
	return ydot->eval_vector(yi.mid()) + J*(yi-yi.mid());
    }
    IntervalVector der(nbvar);

    for (int j=0;j<nbvar;j++)
      der[j]=edtfr->get_derivatives(ordre,Affine2Vector(yi,true),j).itv();
      
    return der;
  };
  
  
    
  Affine2Vector computeRK4derivative(Affine2Vector y)
  {
    Affine2Vector rk4_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      rk4_deriv[i]=edtfr->lteExplicitRK4(i, y);
      
    }
    return rk4_deriv;
    
  }
  
  Affine2Vector computeRADAU3derivative(Affine2Vector y)
  {
    Affine2Vector radau3_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      radau3_deriv[i]=edtfr->lteImplicitRadau3(i, y);
    }
    return radau3_deriv;
    
  }
  
  Affine2Vector computeLC3derivative(Affine2Vector y)
  {
    Affine2Vector lc3_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      lc3_deriv[i]=edtfr->lteImplicitLobbato3c4(i, y);
    }
    return lc3_deriv;
    
  }

  
  Affine2Vector computeHEUNderivative(Affine2Vector y)
  {
    Affine2Vector heun_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      heun_deriv[i]=edtfr->lteExplicitHeun(i, y);
    }
    return heun_deriv;
    
  }
  
  Affine2Vector computeIEULERderivative(Affine2Vector y)
  {
    Affine2Vector ieuler_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      ieuler_deriv[i]=edtfr->lteImplicitEuler(i, y);
      
    }
    return ieuler_deriv;
    
  }
  
  
  Affine2Vector computeIMIDPOINTderivative(Affine2Vector y)
  {
    Affine2Vector imidpoint_deriv(y.size(),0);
    
    for (int i=0;i<y.size();i++)
    {
      imidpoint_deriv[i]=edtfr->lteImplicitMidpoint(i, y);
      
    }
    return imidpoint_deriv;
    
  }
  
   Affine2Vector computeLA3derivative_aff(Affine2Vector y)
  {
    Affine2Vector la3_deriv(y.size());
    
    for (int i=0;i<y.size();i++)
    {
      la3_deriv[i]=edtfr->lteImplicitLobbato3a4(i, y);
      
    }
    return la3_deriv;
    
  }
  
    //with complete affine form
   Affine2Vector compute_derivatives_aff(int ordre, Affine2Vector yi){
    assert(ordre < MAX_FCT);
  
    Affine2Vector der(nbvar);
    
    edtfr->edfr->inc_global_step();

    for (int j=0;j<nbvar;j++)
      der[j]=edtfr->get_derivatives(ordre,yi,j);
    
    return der;
  };
  
  
  
  void compute_derivatives(int ordre){
    assert(ordre < MAX_FCT);

    
    for (int i=0;i<=ordre;i++)
    {	
	cout << "eval "<<i<<"^th diff" << endl;
	for (int j=0;j<nbvar;j++)
	    cout <<  edtfr->get_derivatives(i,Affine2Vector(*yinit,true),j).itv() << endl;
    }
    
    
  };
  
  
    void compute_derivatives2(int ordre){
    assert(ordre < MAX_FCT);
    
    for (int i=0;i<=ordre;i++)
    {	
	cout << "eval "<<i<<"^th diff" << endl;
	for (int j=0;j<nbvar;j++)
	    cout << edtfr->get_derivatives(i,Affine2Vector(*yinit,true),j) << endl;
    }
    
    
  };
  
  
  ~ivp_ode(){
    delete yinit;
    delete yinit_aff;
    delete ydot;
    delete edtfr;
    
  };
  

public:
  int nbvar;
  Function* ydot;
  Affine2Vector* yinit_aff;
  double t0;

private:
  IntervalVector* yinit;
  edtree_frechet* edtfr;
  

};
  

}


#endif