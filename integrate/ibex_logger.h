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
#ifndef CUSTOM_CLogger_H
#define CUSTOM_CLogger_H

#include <fstream>
#include <iostream>
#include <cstdarg>
#include <string>
#include "ibex_IntervalVector.h"

namespace ibex{

#define LOGGER CLogger::getLogger()
/**
 *   Singleton Logger Class.
 */
class CLogger
{
public:
  
    /**
     *   Logs a message
     *   @param sMessage message to be logged.
     */
    void Log(const std::string& sMessage);
    
    /**
     *   Variable Length Logger function
     *   @param format string for the message to be logged.
     */
    void Log( const char * format, ... );
    
    //flush of logging members
    void Log_end();
    
    //methods for logging in simulation
    void inc_rejected_picard();
    void inc_accepted_picard();
    void update_lte_max(double lte);
    void Log_sol(IntervalVector sol);
    void update_step(double step);
    
        /**
     *   << overloaded function to Logs a message
     *   @param sMessage message to be logged.
     */
    CLogger& operator<<(const std::string& sMessage );
    /**
     *   Funtion to create the instance of logger class.
     *   @return singleton object of Clogger class..
     */
    static CLogger* getLogger();
private:
    /**
     *    Default constructor for the Logger class.
     */
    CLogger();
    /**
     *   copy constructor for the Logger class.
     */
    CLogger( const CLogger&){};             // copy constructor is private
    /**
     *   assignment operator for the Logger class.
     */
    CLogger& operator=(const CLogger& ){ return *this;};  // assignment operator is private
    /**
     *   Log file name.
     **/
    static const std::string m_sFileName;
    /**
     *   Singleton logger class object pointer.
     **/
    static CLogger* m_pThis;
    /**
     *   Log file stream object.
     **/
    static std::ofstream m_Logfile;
    
    
    //members for logging in simulation
    static int nb_rejected_picard;
    static int nb_accepted_picard;
    static double norm_lte_max;
    static double step_min;
    static double step_max;
    
};

}

#endif