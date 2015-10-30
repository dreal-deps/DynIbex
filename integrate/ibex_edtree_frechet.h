/* ============================================================================
 * D Y N I B E X - Definition of the Trees based on Frechet derivatives
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_ED_TREE_FRECHET_H
#define IBEX_ED_TREE_FRECHET_H

#include "ibex_edfrechet.h"

namespace ibex{

class edtree_frechet
{
public:
  
  int nbvar;
  int order;
  Function* func;
  
  edfrechet* edfr;
  
  edtree_frechet(int _order,Function* _func,int _nbvar){
    //assert(order <= MAX_ORDER);

      nbvar=_nbvar;
      order=_order;
      func = _func;
      edfr = new edfrechet(order,func,nbvar);

  };

  ~edtree_frechet(){
    delete edfr;
  }
  
  
  Affine2 get_derivatives(int order, Affine2Vector y, int j); 
  Affine2 lteImplicitRadau3 (int j, Affine2Vector y);
  Affine2 lteExplicitRK4 (int j, Affine2Vector y);
  Affine2 lteImplicitLobbato3a4 (int j, Affine2Vector y);
  Affine2 lteImplicitEuler (int j, Affine2Vector y);
  Affine2 lteImplicitMidpoint (int j, Affine2Vector y);
  Affine2 lteImplicitLobbato3c4 (int j, Affine2Vector y);
  Affine2 lteExplicitHeun (int j, Affine2Vector y);
};


}

#endif