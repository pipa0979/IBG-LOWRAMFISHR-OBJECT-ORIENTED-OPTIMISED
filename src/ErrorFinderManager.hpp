#ifndef ERRORFINDERMANAGER_HPP
#define ERRORFINDERMANAGER_HPP

#include "Consolidator.hpp"
#include "ErrorCalculator.hpp"


#include <iostream>
#include <string>


class ErrorFinderManager	//flr GERMLINE.h
{
 public:
       ErrorFinderManager();
       void performConsolidation(int argc, char *argv[]);
       void displayError(std::string argv);
       void initiateErrorFinder();
       void initiateErrorFinder_PRE();
       void initiateErrorFinder_CURRENT(int);
       void initiateErrorFinder_POST();
       bool T_ERROR;
       bool T_TRIM;
       std::string CERRMESSAGE;
       int SITUATION_NO;



 private:
       int WINDOW, MIN_SNP,GAP, MA_SNP_END, TRUESNP;
       float MIN_CM,MA_ERR_THRESHOLD_START, MA_ERR_THRESHOLD_END,
             PCT_ERR_THRESHOLD,HO_THRESHOLD, TRUECM, PIELENGTH, MA_THRESHOLD,EMPIRICAL_MA_RESULT,EMPIRICAL_PIE_RESULT;
       bool ISMOL, COUNTGAPERR;
       std::string BMIDFILE,BMATCHFILE,BSIDFILE,PEDFILE,
                     HPEDFILE,HMAPFILE,OPTION, HO_MISSING, LOGFILE, TRUEIBDFILE;
       ErrorCalculator eCalculator;
       Consolidator consolidator;
       std::string wrongParam;
};

#endif
