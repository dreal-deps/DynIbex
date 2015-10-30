/* ============================================================================
 * D Y N I B E X - Example for a first validated simulation : oil reservoir
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#include "ibex.h"

using namespace ibex;

int main(){

const int n= 3;
Variable y(n);
IntervalVector yinit(n);
yinit[0] = Interval(0);
yinit[1] = Interval(0);
yinit[2] = Interval(0);
Interval par0(2.78,2.79);

Function ydot = Function(y,Return( Interval(1.),
				    y[2],
				    (y[1]*y[1])*y[1]/6.0-y[1]+2.0*sin(par0*y[0])));

ivp_ode problem = ivp_ode(ydot,0.0,yinit);

simulation simu = simulation(&problem,8.0,LC3,1e-10);

simu.run_simulation();

//For an export in order to plot
simu.export3d_yn("export_3d.txt",0,1,2);


return 0;

} 
