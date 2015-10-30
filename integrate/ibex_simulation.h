/* ============================================================================
 * D Y N I B E X - Definition of the Simulation
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_SIMU_H
#define IBEX_SIMU_H

#include <fstream>
#include "ibex_integrate.h"
#include <stdlib.h>

namespace ibex{
//namespace std{
  
  enum Method {IEULER,IMIDPOINT,RADAU3,HEUN,TAYLOR4,LA3,LC3,RK4};

class simulation
{
public:
	double time_T;
	list<solution_j> list_solution_j;
	int process;
	ivp_ode* embedded_ode;
	double atol;
	int test;
	Method meth;
	
	
	//compute solutions after solutions till t >= time_T 
	int run_simulation(){
	  double tn ;
	  
	  do
	  {
		
		//reuse of the first (and verified step) to build the followers
		solution_j sol_temp = list_solution_j.back();
		//sol_temp.print_soljn();
		
		//initialize next solution
		tn = sol_temp.time_j.ub();
		double fac = sol_temp.factor;
		double newh = std::abs(sol_temp.time_j.diam())*fac; //we apply the factor estimated the last step
		
		if (tn < time_T)
		{
		  
		  //cout << "\r Current time : " << tn;
		  double delta_t = time_T-tn;
		  
		  if (delta_t < newh)
		    newh = delta_t;
		    
		    switch(meth)
		    {
			case IEULER:
			{
			  solution_j_ieuler u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case IMIDPOINT:
			{
			  solution_j_imidpoint u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case RADAU3:
			{
			  solution_j_radau3 u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case HEUN:
			{
			  solution_j_heun u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case TAYLOR4:
			{
			  solution_j_tayl4 u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case LA3:
			{
			  solution_j_la3 u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case LC3:
			{
			  solution_j_lc3 u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			case RK4:
			{
			  solution_j_rk4 u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);
			  break;
			}
			
			default : 
			{
			  /*solution_j u_j(*sol_temp.box_jnh_aff,tn,newh, embedded_ode,atol,fac);	
			  u_j.compute_oneStep(*sol_temp.box_jnh_aff,embedded_ode);
			  sol_temp.flush();
			  list_solution_j.push_back(u_j);*/
			  std::cout << "This method does not exist !" << std::endl;
			  tn = time_T;
			}
		    
		    }

		  
		}
	    		tn = sol_temp.time_j.ub();
	  }while (tn < time_T);
	  ///////////////**************//////////////
	  if (tn > 1e7)
	    tn = time_T;
	  ///////////////**************//////////////
	  else
	  {
	    //print of the last solution
	    solution_j sol_temp = list_solution_j.back();
	    //sol_temp.print_soljn();
	    cout << "\nSolution at t=" << tn << " : " << sol_temp.box_jnh_aff->itv() << endl;
	    //cout << "Diam : " << sol_temp.box_jnh->diam() << endl;

	    LOGGER->Log("Solution at t=%f : ",tn);
	    LOGGER->Log_sol(sol_temp.box_jnh_aff->itv());
	    LOGGER->Log_end();
	    
	  }
	  
	  return 1;
	};
		
	//find the box of a stack which contains a sub-box (stack size nbvar/2, y size nbvar)
	IntervalVector find_box(list<IntervalVector>* stack, IntervalVector y)
	{
	  
	  list<IntervalVector>::iterator iterator_list;
  
	  IntervalVector suby = y.subvector(0,y.size()/2 - 1);
	  for(iterator_list=stack->begin();iterator_list!=stack->end();iterator_list++)
	  {
	    IntervalVector yv =*iterator_list;
	    
	    
	    if (!(suby & yv).is_empty())
	      return yv;
	  }	
	  
	  return IntervalVector(y.size()/2);
	}
	
	
	int order_of_method(Method _meth)
	{
	  switch(meth)
		    {
			case IEULER:
			{
			  return 2;
			  break;
			}
			case IMIDPOINT:
			{
			  return 3;
			  break;
			}
			case RADAU3:
			{
			  return 4;
			  break;
			}
			case HEUN:
			{
			  return 3;
			  break;
			}
			case TAYLOR4:
			{
			  return 5;
			  break;
			}
			case LA3:
			{
			  return 5;
			  break;
			}
			case LC3:
			{
			  return 5;
			  break;
			}
			case RK4:
			{
			  return 5;
			  break;
			}	
			
			default : 
			{
			  return 5;
			}
		    }
	 return 5;
	}
	
	//constructor
	simulation(ivp_ode* _ode, double T,Method _meth, double a){
	  process = 0;
	  embedded_ode = _ode;
	  time_T = T;
	  atol = a;
	  meth = _meth;
	  
	  
	  double h = 0.01;
	  
	  embedded_ode->frechet_precomputing(order_of_method(meth));
	  
	  solution_j u_j0(*embedded_ode->yinit_aff,embedded_ode->t0,h,atol);
	  list_solution_j.push_back(u_j0);
	}
	
	//constructor
	simulation(ivp_ode* _ode, double T){
	  process = 0;
	  embedded_ode = _ode;
	  time_T = T;
	  atol = 1e-6;	
	  meth = RK4;
	  
	  double h = 0.01;
	  embedded_ode->frechet_precomputing(order_of_method(meth));
	  solution_j u_j0(*embedded_ode->yinit_aff,embedded_ode->t0,h,atol);
	  list_solution_j.push_back(u_j0);
	}
	
	//constructor
	simulation(ivp_ode* _ode, double T,Method _meth){
	  process = 0;
	  embedded_ode = _ode;
	  time_T = T;
	  atol = 1e-6;
	  meth = _meth;
	  
	  double h = 0.01;
	  embedded_ode->frechet_precomputing(order_of_method(meth));
	  solution_j u_j0(*embedded_ode->yinit_aff,embedded_ode->t0,h,atol);
	  list_solution_j.push_back(u_j0);
	}
	
	//constructor
	simulation(ivp_ode* _ode, double T,double a){
	  process = 0;
	  embedded_ode = _ode;
	  time_T = T;
	  atol = a;
	  meth = RK4;
	  
	  double h = 0.01;
	  embedded_ode->frechet_precomputing(order_of_method(meth));
	  solution_j u_j0(*embedded_ode->yinit_aff,embedded_ode->t0,h,atol);
	  list_solution_j.push_back(u_j0);
	}

	

	//destructor
	~simulation(){
	  list_solution_j.clear();  
	}
	

	//functions to test the solution
	IntervalVector get_last()
	{
	    solution_j sol_temp = list_solution_j.back();
	    return sol_temp.box_jnh_aff->itv();
	}


	//final solution is include in a box
	bool finished_in(IntervalVector y_final)
	{
	  if (y_final.size()!=embedded_ode->nbvar)
	    return false;
	  
	  if(!list_solution_j.empty())
	  {
	    solution_j sol_temp = list_solution_j.back();
	    return (sol_temp.box_jnh_aff->itv()).is_subset(y_final);	    
	  }
	  else 
	    return false;	  
	}
	
	//final solution is include in at least one box of a list
	bool finished_in(list<IntervalVector> *stack)
	{
	  list<IntervalVector>::iterator iterator_l;
	  for(iterator_l=stack->begin();iterator_l!=stack->end();iterator_l++)
	  {
	      if (finished_in((IntervalVector)*iterator_l))
		return true;
	  }
	  return false;
	}
	
	//tube crosses a box
	bool has_crossed(IntervalVector y)
	{
	  if (y.size()!=embedded_ode->nbvar)
	    return false;
	    
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		 IntervalVector y_temp = *(iterator_list->box_j1);
		 if (!((y_temp&y).is_empty()))
		   return true;
	      }
	  }
	  else 
	    return false;	  
	}

	
	//tube stays in a box
	bool stayed_in(IntervalVector y_hull)
	{
	  if (y_hull.size()!=embedded_ode->nbvar)
	    return false;
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		 IntervalVector y_temp = *(iterator_list->box_j1);
		 if (!y_temp.is_subset(y_hull))
		   return false;
	      }
	  }
	  return true;	  
	}
	
	//at least one y(t) outside a box
	bool go_out(IntervalVector y_hull)
	{
	  if (y_hull.size()!=embedded_ode->nbvar)
	    return false;
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		 IntervalVector y_temp = iterator_list->box_jnh_aff->itv();
		 if ((y_temp&y_hull).is_empty())
		   return true;
	      }
	  }
	  return false;
	}
	
	
	//tube stays in a box till t
	bool stayed_in_till(IntervalVector y_hull, double t)
	{
	  if (y_hull.size()!=embedded_ode->nbvar)
	    return false;
	  
	  if (t <0)
	    return true;
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		 IntervalVector y_temp = *(iterator_list->box_j1);
		 if (iterator_list->time_j.lb() > t)
		   return true;
		 if (!y_temp.is_subset(y_hull))
		   return false;
	      }
	  }
	  return true;	  
	}
	
	
	//final solution crosses a box
	bool has_reached(IntervalVector y_final)
	{
	  if (y_final.size()!=embedded_ode->nbvar)
	    return false;
	    
	  if(!list_solution_j.empty())
	  {
	    solution_j sol_temp = list_solution_j.back();
	    return (!((sol_temp.box_jnh_aff->itv())&y_final).is_empty());    
	  }
	  else 
	    return false;	
	  
	}
	
	//final solution crosses at least one box of a list
	bool has_reached(list<IntervalVector> *stack)
	{
	  list<IntervalVector>::iterator iterator_l;
	  for(iterator_l=stack->begin();iterator_l!=stack->end();iterator_l++)
	  {
	      if (has_reached(*iterator_l))
		return true;
	  }
	  return false;
	}
	
	
	//one solution is inside at least one box of a list
	double one_in(list<IntervalVector> *stack)
	{
	  list<IntervalVector>::iterator iterator_l;
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		
		IntervalVector y_temp = iterator_list->box_jnh_aff->itv();
		double t= iterator_list->time_j.ub();
		 		 
		for(iterator_l=stack->begin();iterator_l!=stack->end();iterator_l++)
		{
		  if (y_temp.is_subset(*iterator_l))
		    return t;
		}
	      }
	  }
	  return (-1.0);
	}
	
	
	//return a box containing y(t)
	IntervalVector get(double t)
	{
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		if ((iterator_list->time_j).contains(t))
		  return *iterator_list->box_j1;
	      }
	  }
	  return IntervalVector(embedded_ode->nbvar);  
	  
	}
	
	
	
	//return a box contained in the next box (possibly an attractor)
	IntervalVector get_attractor()
	{
	  
	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      iterator_list=list_solution_j.begin();
	      IntervalVector attractor = iterator_list->box_jnh_aff->itv();
	      
	      iterator_list++; //prob 2 val idem
	      
	      for(;iterator_list!=list_solution_j.end();iterator_list++)
	      {
		if (attractor.is_subset(iterator_list->box_jnh_aff->itv()))
		  return attractor;
		else
		  attractor = iterator_list->box_jnh_aff->itv();
	      }
	  }
	  
	  return IntervalVector(embedded_ode->nbvar);  
	  
	}
	
	
	//return the domain covered by the simulation (the hull of y)
	IntervalVector get_domain()
	{	  	  
	  if(!list_solution_j.empty())
	  {
	      list<solution_j>::iterator iterator_list;
	      iterator_list=list_solution_j.begin();
	      
	      IntervalVector hull = iterator_list->box_jnh_aff->itv();
	      
	      for(;iterator_list!=list_solution_j.end();iterator_list++)
	      {
		hull |= iterator_list->box_jnh_aff->itv();
	      }
	      return hull;
	  }
	  
	  return IntervalVector(embedded_ode->nbvar);  
	  
	}
	
	
	
	//export in a file for ploting
	void export2d_yn(const char* filename, int x, int y)
	{
	  assert(embedded_ode->nbvar > 1);
	  cout << "export in progress..." << endl;
	  if(!list_solution_j.empty())
	  {
	      ofstream file(filename, ios::out | ios::trunc);
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		  file << iterator_list->box_jnh->operator[](x) << " ; " << iterator_list->box_jnh->operator[](y) << endl;
	      }
	  
	      file.close();
	  }
	  
	};
	
	//export in a file for ploting
	void export3d_yn(const char* filename, int x, int y, int z)
	{
	  assert(embedded_ode->nbvar > 2);
	  cout << "export in progress..." << endl;
	  if(!list_solution_j.empty())
	  {
	      ofstream file(filename, ios::out | ios::trunc);
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		  file << iterator_list->box_jn->operator[](x) << 
		  " ; " << iterator_list->box_jn->operator[](y) << 
		  " ; " << iterator_list->box_jn->operator[](z) << endl;
	      }
	      
	  
	      file.close();
	  }
	  
	};
	
	//export in a file for ploting
	void export1d_yn(const char* filename, int x)
	{
	  assert(embedded_ode->nbvar > 1);
	  cout << "export in progress..." << endl;
	  if(!list_solution_j.empty())
	  {
	      ofstream file(filename, ios::out | ios::trunc);
	      list<solution_j>::iterator iterator_list;
	      for(iterator_list=list_solution_j.begin();iterator_list!=list_solution_j.end();iterator_list++)
	      {
		  file << iterator_list->box_j1->operator[](x) << 
		  " ; " << iterator_list->time_j << endl;
	      }
	  
	      file.close();
	  }
	  
	};
	
  
};
//}
}

#endif
