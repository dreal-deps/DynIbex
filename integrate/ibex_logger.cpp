/* ============================================================================
 * D Y N I B E X - Definition of the Logger
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#include "ibex_logger.h"

namespace ibex{

const std::string CLogger::m_sFileName = "./Log.txt";
CLogger* CLogger:: m_pThis = NULL;
std::ofstream CLogger::m_Logfile;

int CLogger::nb_rejected_picard = 0;
int CLogger::nb_accepted_picard = 0;
double CLogger::norm_lte_max = 0.0;
double CLogger::step_min = 100.0;
double CLogger::step_max = 0.0;
    
CLogger::CLogger()
{
 
}
CLogger* CLogger::getLogger(){
    if(m_pThis == NULL){
        m_pThis = new CLogger();
        m_Logfile.open(m_sFileName.c_str(), std::ios::out | std::ios::trunc );
	
    }
    return m_pThis;
}

void CLogger::Log( const char * format, ... )
{
    char sMessage[256];
    va_list args;
    va_start (args, format);
    vsprintf (sMessage,format, args);
    m_Logfile << sMessage;
    va_end (args);
    m_Logfile << "\n";
}

void CLogger::Log( const std::string& sMessage )
{
    m_Logfile << sMessage << "\n";
}


void CLogger::Log_sol(IntervalVector sol)
{
  m_Logfile << sol << "\n";
  m_Logfile << "Diameter : " << sol.diam() << "\n";
  
}
  
  
void CLogger::Log_end()
{
    m_Logfile << "Rejected picard :" << nb_rejected_picard << "\n";
    m_Logfile << "Accepted picard :" << nb_accepted_picard << "\n";
    m_Logfile << "Step min :" << step_min << "\n";
    m_Logfile << "Step max :" << step_max << "\n";
    m_Logfile << "Truncature error max :" << norm_lte_max << "\n";
}

CLogger& CLogger::operator<<(const std::string& sMessage )
{
    m_Logfile << sMessage << "\n";
    return *this;
}


//methods for logging in simulation
void CLogger::inc_rejected_picard(){
    nb_rejected_picard++;
}
void CLogger::inc_accepted_picard(){
    nb_accepted_picard++;
}

void CLogger::update_lte_max(double lte)
{
    if (norm_lte_max < lte)
      norm_lte_max = lte;
}

void CLogger::update_step(double step){
    step_min = std::min(step_min,step);
    step_max = std::max(step_max,step);
}

}

