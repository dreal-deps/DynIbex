/* ============================================================================
 * D Y N I B E X - Definition of the Frechet derivatives
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#ifndef IBEX_ED_FRECHET_H
#define IBEX_ED_FRECHET_H

#include <map>
#include <iostream>

const int ORDER_MAX=20;
  
namespace ibex{

  typedef std::map<std::vector<int>,Function*> function_map;
  
  typedef struct memo {int step; Affine2 value;} memo; 
  typedef std::map<std::vector<int>,memo> value_map;
  
  
  template<typename T>
    ostream& operator<< (ostream& out, const vector<T> v) {
    int last = v.size() - 1;
    out << "[";
    for(int i = 0; i < last; i++)
        out << v[i] << ", ";
    out << v[last] << "]";
    return out;
}

  
class edfrechet
{
  public:   
    
    int _global_step;
    
    void inc_global_step(){
      _global_step++;
    };
    
    edfrechet(int _order, Function* _func, int _nbvar){
     // assert(_order < ORDER_MAX);
      m_order = _order;
      std::vector<int> vec;
      
      _global_step=1;
      
      memo zero = {0,Affine2(0)};
      
      //sans recurrence
      for (int j=0;j<_nbvar;j++)
      {
	  std::vector<int> temp(1,j);
	  
	  std::pair<std::vector<int>,Function*> elnt(temp,&_func->operator[](j));
	  tab_derivatives[0].insert(elnt);
	  
	  
	  std::pair<std::vector<int>,memo> elnt2(temp,zero);
	  tab_values[0].insert(elnt2);
	  
      }
      for (int i=1;i <= m_order;i++)
      {
	compute_derivatives(i,_nbvar);
      }
      
      //std::cout << "computation of elementary derivatives is done !" << std::endl; 
	
	//print();
    }; 
    
    
    ~edfrechet(){
      
    }
    
    
    Interval eval_frechet(int order, std::vector<int> key, IntervalVector y){  
      //assert(order < ORDER_MAX);
      
      
      return tab_derivatives[order][key]->eval(y);     
    };
    
    
   Affine2 eval_frechet(int order, std::vector<int> key, Affine2Vector y){  
      //assert(order < ORDER_MAX);
      
      memo temp = tab_values[order][key];
      Affine2 temp2;
      
      if (temp.step < _global_step)
      {
	
	  temp2 = tab_derivatives[order][key]->eval_affine2(y);
	  tab_values[order][key].value = temp2;
	  tab_values[order][key].step = _global_step;
      }
      else
      {
	temp2 = temp.value;
      }
      
      return temp2;
    };
    
    void print(){
      for (int j=0;j<m_order;j++)
      {
	function_map mymap = tab_derivatives[j];
	
	cout << "order : " << j << endl;
	
	  for (function_map::iterator it=mymap.begin(); it!=mymap.end(); ++it)
	  {
	      std::cout << " " << it->first << ":" << *(it->second) << " ; ";
	  }
	  
	  std::cout << std::endl;
	}
	std::cout << std::endl;
      };
    

  private:
    
    int m_order;
    function_map tab_derivatives[ORDER_MAX];
    value_map tab_values[ORDER_MAX];
    
    

    void compute_derivatives(int order,const Function* function,int nbvar,std::vector<int> key)
    {
      for (int j=0;j<nbvar;j++)
      {
	  std::vector<int> temp = key;
	  temp.push_back(j);
	  function_map temp_map = tab_derivatives[order];
	  std::pair<std::vector<int>,Function*> elnt(temp,&function->operator[](j));
	  temp_map.insert(elnt);
	  tab_derivatives[order] = temp_map;
	    
	  if (order >= m_order)
	  {
	      return;
	  }
	  else
	  {
	    compute_derivatives(order+1,&(function->operator[](j).diff()),nbvar,temp);	    
	  }   
      }
           
    };
    
    
    void compute_derivatives(int order,int nbvar)
    {
      assert(order > 0);
      function_map temp_map = tab_derivatives[order-1];
      
      function_map::iterator it = temp_map.begin();
      
      for (;it!=temp_map.end();++it)
      {
	std::vector<int> temp_key = it->first;
	Function* temp_func = it->second;
	
	const Function& temp_func2 = temp_func->diff();
	for (int j=0;j<nbvar;j++)
	{
	    std::vector<int> temp_key2 = temp_key;
	    temp_key2.push_back(j);
	    
	    
	    std::pair<std::vector<int>,Function*> elnt(temp_key2,&(temp_func2.operator[](j)));
	    tab_derivatives[order].insert(elnt);
	    
	    memo zero = {0,Affine2(0)};
	    std::pair<std::vector<int>,memo> elnt2(temp_key2,zero);
	    tab_values[order].insert(elnt2);
	  
	}
	
      }
           
    };

};


}

#endif