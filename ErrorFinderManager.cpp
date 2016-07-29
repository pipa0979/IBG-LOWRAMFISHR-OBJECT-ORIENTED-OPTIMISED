#include <cstdlib>
#include <iostream>
#include <cstring>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>

#include "ErrorFinderManager.hpp"
#include "ErrorCalculator.hpp"
#include "Consolidator.hpp"




#define GetCurrentDir getcwd

/*This function was copied from consolidator.cpp. This is a standalone function*/
template <typename T>
  std::string NumberToString ( T Number )
  {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }

ErrorFinderManager::ErrorFinderManager():WINDOW(50),	//lrm ErrorFinderManager.hpp
                    MIN_SNP(30),GAP(1), MA_SNP_END(50), TRUESNP( 600 ),
                    MIN_CM(0.4), MA_ERR_THRESHOLD_START(0.08),
                    MA_ERR_THRESHOLD_END(0.08),PCT_ERR_THRESHOLD( 0.90 ),
                    HO_THRESHOLD( 0.98 ), TRUECM( 6 ),PIELENGTH( 3 ),
                    ISMOL( false), COUNTGAPERR( false ),MA_THRESHOLD(0.8),EMPIRICAL_MA_RESULT(-1.0),EMPIRICAL_PIE_RESULT(-1.0),T_ERROR(false),T_TRIM(false),CERRMESSAGE(""),SITUATION_NO(0)

{
}

void ErrorFinderManager::performConsolidation(int argc,char *argv[])
{
       bool goodParam=true;
       bool thresholdError = false; //if the user supplies both -empricial-ma-threshold and -ma-threshold, that is an error. This will be used to detect such an error.
       bool pieThresholdError = false;
       char currentPath[FILENAME_MAX];
       if(!GetCurrentDir(currentPath,sizeof(currentPath)))
       {
std::cerr << "Error reading current directory" << std::endl;
		return;
       }
       for(int i=1;i<argc;i++)
        {

                if(strcmp(argv[i],"-bmatch")==0&&i<argc-1)
                {
                     BMATCHFILE=std::string(argv[++i]);
                }
                else if(strcmp(argv[i],"-bmid")==0&&i<argc-1)
                {
                     BMIDFILE=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-bsid")==0&&i<argc-1)
                {
                     BSIDFILE=std::string(argv[++i]);
                }
                 else if(strcmp(argv[i],"-reduced")==0&&i<argc-2)
                {
                     MIN_SNP=atoi(argv[++i]);
                     MIN_CM=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ped-file")==0&&i<argc-1)
                {
                     PEDFILE=std::string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-ped")==0&&i<argc-1)
                {
                     HPEDFILE=std::string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-map")==0&&i<argc-1)
                {
                     HMAPFILE=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-window")==0&&i<argc-1)
                {
                     WINDOW=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-err-threshold-start")==0&&i<argc-1)
                {
                     MA_ERR_THRESHOLD_START=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-threshold")==0&&i<argc-1)
                {
                     HO_THRESHOLD=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-trueCM")==0&&i<argc-1)
                {
                     TRUECM=atof(argv[++i]);
                }
                else  if( strcmp( argv[i],"-trueSNP" )==0 && i < argc-1 )
                {
                     TRUESNP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-missing")==0&&i<argc-1)
                {
                     HO_MISSING= std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-err-threshold-end")==0&&i<argc-1)
                {
                     MA_ERR_THRESHOLD_END=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-gap")==0&&i<argc-1)
                {
                     GAP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-snp")==0&&i<argc-1)
                {
                     MA_SNP_END=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-pct-err-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
                     std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
                      exit(1);
                     }
                     PCT_ERR_THRESHOLD=atof(argv[++i]);
                     pieThresholdError = true;
                }
                else if(strcmp(argv[i],"-emp-pie-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
                    	 std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
                      exit(1);
                     }
                     EMPIRICAL_PIE_RESULT=atof(argv[++i]); 
                     pieThresholdError = true;               
                }
		            else  if(strcmp(argv[i],"-output.type")==0&&i<argc-1)
                {
                     OPTION=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-log.file")==0&&i<argc-1)
                {
                     LOGFILE=std::string(argv[++i]);
                } 
                else if(strcmp(argv[i], "-trueibd")==0&&i<argc-1)
                {
                     TRUEIBDFILE = std::string(argv[++i]);
                }
            		else if(strcmp(argv[i],"-ma-threshold")==0&&i<argc-1)//adding new -ma-threshold argument
            		{
            		     if(thresholdError == true){ //user has already supplied an empirical-ma-threshold, so exit the program with an error message
            		    	 std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
                       exit(1);
            		     }
            		     MA_THRESHOLD=atof(argv[++i]);
            		     thresholdError = true;
            		}
            		else if(strcmp(argv[i],"-empirical-ma-threshold")==0 && i<argc-1){
            		     if(thresholdError == true){
            		    	 std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
            			exit(1);
            		     }	
            		     EMPIRICAL_MA_RESULT = atof(argv[++i]); //use the user supplied empirical ma threshold, instead of calculating it via true ibd segments	
            		     thresholdError = true;
            		}
                else  if(strcmp(argv[i],"-PIE.dist.length")==0&&i<argc-1)
                {
                     std::string MOL=std::string(argv[++i]);
                     if( MOL.compare( "MOL" ) ==0 )
                     {
                        ISMOL = true;
                     }
                     else
                     {
                        PIELENGTH = atof( MOL.c_str() );
                     }
                }
                else  if(strcmp(argv[i],"-count.gap.errors")==0&&i<argc-1)
                {
                     std::string option=std::string(argv[++i]);
                     if( option.compare( "TRUE" ) ==0 )
                     {
                        COUNTGAPERR = true;
                     }
                }

                else
                {
                        wrongParam += " " + std::string(argv[i]);

                        goodParam=false;
                }
        }
        if((!goodParam)||BMATCHFILE.compare("")==0||BSIDFILE.compare("")==0||BMIDFILE.compare("")==0||PEDFILE.compare("")==0)
        {

                displayError( argv[0] );
                return;
             
        }
        if( OPTION.compare( "" ) == 0 )
        {
        	std::cerr<< " please provide a valid output.type option " <<std::endl;
            exit( -1 );
        }
        if( LOGFILE.compare( "" ) == 0 )
        {
        	std::cerr<< " default log file name is FISH " <<std::endl;
            LOGFILE = "FISH";
        }
        eCalculator.createLogFile( LOGFILE  );
        eCalculator.countGapErrors( COUNTGAPERR ); 
        time_t startTime;

        time (&startTime);
        std::string str1 = " The program started at: " + std::string( ctime ( &startTime ) );
        std::string str = " Program working directory was: " + std::string(currentPath) +
		     " \nProgram version was: " + std::string(argv[0]) +
		     " \nProgram options:\n-bmatch file: " + BMATCHFILE +
                     " \n-bmid file: " + BMIDFILE +
                     " \n-bsid file: " + BSIDFILE +
                     " \n-ped file: " + PEDFILE +
                     " \n-holdout ped file: " + HPEDFILE +
                     " \n-holdoutmap file: " + HMAPFILE +
                     " \n-output type: " + OPTION +
                     " \n- missing SNP representation in pedfile: " + HO_MISSING +
                     " \n-log file: " + LOGFILE;
        str = str + " \nmin snp length : " + NumberToString( MIN_SNP )  +
                    " \nmin cm length : " + NumberToString( MIN_CM  ) +
                    " \ngap to consolidate : " + NumberToString( GAP ) +
                    " \nmoving averages window size : " + NumberToString(  WINDOW ) +
                    " \ndiscard ends to calculate pct err : " + NumberToString( MA_SNP_END ) +
                    "  \npercentage error threshold: " + NumberToString( PCT_ERR_THRESHOLD );
                  eCalculator.log( str );  

      //initiateErrorFinder();
//int pers_count=eCalculator.getNoOfPersons();
initiateErrorFinder_PRE();	//Load bsid bmid ped filesa and set the person count
int pers_count=eCalculator.getNoOfPersons();

SITUATION_NO = 0;	//falls in none of the 4 categories (explained in the header file), if 0 throw error and exit
/*
 * 0 indicates HPED HMAP: TRUE,ISMOL: TRUE
 * 1 indicates HPED HMAP: TRUE,ISMOL: FALSE
 * 2 indicates HPED HMAP: FALSE,ISMOL: TRUE
 * 3 indicates HPED HMAP: FALSE,ISMOL: FALSE
*/
initiateErrorFinder_CURRENT(pers_count);

if (	(SITUATION_NO == 1)	||	(SITUATION_NO == 3)		)
{
//	/consolidator.findTruePctErrors( eCalculator, MA_SNP_END, T_ERROR, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count, TRUESNP, TRUECM);
	ErrorCalculator &e_obj = eCalculator;
	int ma_snp_ends = MA_SNP_END;
	//bool holdOut = T_ERROR;
	int window = WINDOW;
	float ma_threshold = MA_THRESHOLD;
	float empirical_ma_threshold = EMPIRICAL_MA_RESULT;
	std::string path = BMATCHFILE;
	int person_count = pers_count;
	int trueSNP = TRUESNP;
	float trueCM = TRUECM;


	//ErrorCalculator& e_obj,
	//int window = WINDOW;
	//int ma_snp_ends = MA_SNP_END;
	//float ma_threshold = MA_THRESHOLD;
	int min_snp = MIN_SNP;
	float min_cm = MIN_CM;
	float per_err_threshold = PCT_ERR_THRESHOLD;
	std::string option = OPTION;
	float hThreshold = HO_THRESHOLD;
	//bool holdOut,
	float empirical_threshold = EMPIRICAL_MA_RESULT;
	//std::string path = BMATCHFILE;
	float empirical_pie_threshold = EMPIRICAL_PIE_RESULT;


	int removed1 = 0, removed2 = 0, removed3 = 0, removed4 = 0;
	int not_removed = 0;
	int total_count = consolidator.global_initial; //for keeping track of total #SH
	//bool wrongOption = false;
	float per_err_threshold1;
	if(empirical_pie_threshold >= 0.0)
		{
			per_err_threshold1 = empirical_pie_threshold;
		}
	else
		{
			per_err_threshold1 = consolidator.getPctErrThreshold( per_err_threshold );
		}
	str = " \ncorresponding empirical error threshold for " + NumberToString( per_err_threshold ) + " is: " + NumberToString( per_err_threshold1 )  + " \n";
	float hThreshold1 = 0;
	if( T_TRIM )
		{
			hThreshold1 = consolidator.getHoldOutThreshold( hThreshold );
			str = str + " \n corresponding holdout emperical error threshold for " + NumberToString( hThreshold ) + " is: " + NumberToString( hThreshold1 )  + " \n";
		}
	e_obj.log( str );
	per_err_threshold = per_err_threshold1;
	hThreshold = hThreshold1;
	//begin streaming in bmatch data
	unsigned int pid[2];
	unsigned int sid[2];
	unsigned int dif,hom[2];
	std::ifstream file_bmatch(path.c_str(),std::ios::binary);
	int pers_count = person_count;
	SNP_lrf snp;
	if( !file_bmatch )
		{
			std::cerr<<"unable to open the bmatch file, exiting the program" << std::endl;
			exit( -1 );
		}



	while ( !file_bmatch.eof())
		{
			//read in bmatch data

			pid[0] = -1;
			file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
			if ( pid[0] == -1 ) continue;
			file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &dif , sizeof( int ) );
			file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
			file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
			if(pid[0]>=pers_count||pid[1]>=pers_count)
				{
					std::cerr<<"problem with bsid file in findTruePctErrors, check it please"<<std::endl;
					return;
				}

			snp.start=sid[0];
			snp.end=sid[1];
			if(pid[0] > pid[1])
				{
					int temp = pid[0];
					pid[0] = pid[1];
					pid[1] = temp;
				}

			int i = pid[0];
			int j = pid[1];

			if(snp.end==-1)
				{
					continue;
				}

			if( ( e_obj.getCMDistance( sid[ 1 ] ) -  e_obj.getCMDistance( sid[ 0 ] ) ) >= trueCM && ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
				{

					//handle moving averages calculation
					//
					int t1 = snp.start + ( snp.end - snp.start ) * 0.25;
					int t2 = snp.end -  ( snp.end -  snp.start ) * 0.25;
					//now we have the positions of the first and last 25% of the truly ibd SH
					//all that's left to do is to pass them into the moving averages function, and obtain the max ma
					//then store that in a vector, sort them, and find the xth percentile of that vector. That will be
					//the ma that we use later
					//for that "finalErrors" parameters, need to get the number of errors along the truly IBD SH first...
					std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
					std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
					//handles MA calculations
					std::vector<float> av;
					float current_max;
					if(empirical_ma_threshold < 0.0)
						{
							av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
							current_max = av[0];
							for(int q = 1; q < av.size(); q++)
								{
									if(av[q] > current_max)
										{
											current_max = av[q];
										}
								}
							e_obj.addMaxAverage(current_max);
						}
					//
					int temp1 = snp.start +	 ( snp.end -  snp.start ) * 0.15; //Should probably stop doing this
					int temp2 = snp.end -  ( snp.end -	snp.start ) * 0.15;
					int start =0, end =0, fend = ( temp2 -temp1 )  ;

					//since we are using MOL at this point, this will pick out a random SH from the set of non-truly IBD SH
					//and use that length to define the region over which we find PIE. Unless you are changing something with MOL,
					//don't ever read this next block
					/*
					int randPers1, randPers2, pos;
					randPers1 = std::rand() % person_count;
					randPers2 = std::rand() % person_count;
					if( randPers1 > randPers2 )
					{
					randPers1 = randPers1 + randPers2;
					randPers2 = randPers1 - randPers2;
					randPers1 = randPers1 - randPers2;
					}
					while( m_matches[ randPers1 ][ randPers2 ].size() <= 0 )
					{
					randPers1 = std::rand() % person_count;
					randPers2 = std::rand() % person_count;
					if( randPers1 > randPers2 )
					{
					 randPers1 = randPers1 + randPers2;
					 randPers2 = randPers1 - randPers2;
					 randPers1 = randPers1 - randPers2;
					}

					}
					pos = std::rand() % m_matches[ randPers1 ][ randPers2 ].size();
					int len = m_matches[ randPers1 ][ randPers2 ][ pos ].end
						   - m_matches[ randPers1 ][ randPers2 ][ pos ].start;
					if( len >= fend || len <= 0)
					{
					 continue;
					}

					temp1 = temp1;
					temp2 = temp1 + len;
					*/
					//end crazy MOL stuff

					std::vector<std::vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);
					std::vector<int>finalErrors=e_obj.getFinalErrors( errors );
					//                  float per_err = e_obj.getThreshold(finalErrors,temp1,temp2,ma_snp_ends );
					float per_err = e_obj.getThreshold(finalErrors,temp1,temp2);//overload
					consolidator.m_errors.push_back( per_err );
					if( T_ERROR  )
						{
							float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
							consolidator.m_holdOutErrors.push_back( oppHom );
						}

				}
			else
				{
					//continue;
				}
		    int temp1=snp.start;
		    int temp2=snp.end;
		    int pers1 = pid[0], pers2 = pid[1];
		    if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
				{
				  pers1 = std::rand() % e_obj.getNoOfPersons();
				  pers2 = std::rand() % e_obj.getNoOfPersons();
				  if( pers1 > pers2 )
					  {
						pers1 = pers1 + pers2;
						pers2 = pers1 - pers2;
						pers1 = pers1 - pers2;
					  }
				}

		    std::vector<std::vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
		    std::vector<int>finalErrors=e_obj.getFinalErrors(errors);
		    /*Inject implied error at start/end of SH here*/
		    std::vector<int>::iterator it;
		    it = finalErrors.begin(); //go to the start of the vector
		    if(finalErrors[0] != 1)
				{
				  finalErrors.insert(it,1); //inject an error at position 1, if not already there
				}
		    /*End inject implied error section*/
		    std::vector<int>trimPositions;
		    std::vector<float>movingAverages;
		    float threshold;
		    //length 32 should survive this, since the reduced arg is 32 (although the reduced cM may be the culprit)

		    if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) )
				{ //initial drop. Don't calculate MA
				  trimPositions.push_back(temp1);
				  trimPositions.push_back(temp2);
				  trimPositions.push_back(1);
				}
		    else
				{
					movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window);
					if(empirical_threshold < 0.0)
						{
						  threshold = e_obj.getCutoff(); //this is the empirical average threshold for moving averages
						}
					else
						{
						  threshold = empirical_threshold;
						}
					trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm);
				}
		    //-----------------
		    int beforeTrimStart = temp1;
		    int beforeTrimEnd = temp2;
		    snp.end = temp2 = temp1+trimPositions[1];
		    snp.start = temp1 = temp1+trimPositions[0];
		    int del0 = trimPositions[0];
		    int del1 = trimPositions[1];
		    //in order to update PIE calculation, may have to alter the finalErrors vector before passing it into getThreshold(). But the alteration
		    //would need to occur at the start and end of the trimmed SH, not just at the start and end points of the initial SH.
		    float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);
		    //add new weighted option
		    /*
		     For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
		    */
		    if( (option.compare("weightedOutput") == 0) )
				{
					int snp1 = 0, snp2 = 0, hlength = 0;
					float noOfOppHom = 0;
					/*if( holdOut )
					{
					snp1 = e_obj.getNewSnp( temp1 );
					snp2 = e_obj.getNewSnp( temp2 );
					hlength = snp2 - snp1;
					if( hlength <= 0 )
					{
					  hlength = 1;
					}
					noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
					}*/
					if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) )
						{
							removed4++;
							continue;
						}
					if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
						{ //removed2 a tpos.size of 3 indicates trimming due ot cM
							removed2++;
							continue;
						}
					if( per_err > per_err_threshold)
						{
							removed1++;
							continue;
						}
					if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
						{
							removed3++;
							continue;
						} //removed3

					not_removed++;
					consolidator.m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
					continue;
				}//end weghtedOutput
				if( (option.compare("FullPlusDropped") == 0) )
					{
						int snp1 = 0, snp2 = 0, hlength = 0;
						float noOfOppHom = 0;
						if( T_TRIM )
							{
								snp1 = e_obj.getNewSnp( temp1 );
								snp2 = e_obj.getNewSnp( temp2 );
								hlength = snp2 - snp1;
								if( hlength <= 0 )
									{
									  hlength = 1;
									}
								noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
							}


				if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) )
					{
						std::vector<float>movingAverages;
						temp1 = beforeTrimStart;
						temp2 = beforeTrimEnd;
						e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);
						snp.start= snp.end=-1;
						removed4++;
						continue;
					}

				if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
					{
						e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2);
						removed2++;
						snp.start= snp.end=-1;
						continue;
					}
				if( per_err > per_err_threshold)
					{
						e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
						removed1++;
						continue;
					}
				if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
					{
						e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
						removed3++;
						snp.start= snp.end=-1;
						continue;
					}
				e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
				continue;
					} //end FullPlusDropped

		    //add finalOutput here
				if(option.compare("finalOutput") == 0)
					{
						//only print if it passes all of the tests.
						int snp1 = 0, snp2 = 0, hlength = 0;
						float noOfOppHom = 0;
						if( T_TRIM )
							{
								snp1 = e_obj.getNewSnp( temp1 );
								snp2 = e_obj.getNewSnp( temp2 );
								hlength = snp2 - snp1;
								if( hlength <= 0 )
									{
									  hlength = 1;
									}
								noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
							}
						if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) )
							{
								std::vector<float>movingAverages;
								temp1 = beforeTrimStart;
								temp2 = beforeTrimEnd;
								snp.start= snp.end=-1;
								removed4++;
								continue;
							}
						if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
							{
								removed2++;
								snp.start= snp.end=-1;
								continue;
							}
						if( per_err > per_err_threshold)
							{
								removed1++;
								continue;
							}
						if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
							{
								removed3++;
								snp.start= snp.end=-1;
								continue;
							}
					  e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
					  continue;
					}//end finaloutput

		    //Calculate Error1
				if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) )
					{

						if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) )
							{ //dropped before trimming
								//don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
								std::vector<float>movingAverages;//null
								//trying something special in this case. This can be removed once idrops aren't being trimmed
								//test code
								temp1 = beforeTrimStart;
								temp2 = beforeTrimEnd;
								//
								e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
								snp.start= snp.end=-1;
								removed4++; //seems ok
								continue;
							}
						if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
							{
								e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
								++removed2;
								snp.start= snp.end=-1;
								continue;
							}
						if( per_err > per_err_threshold ) //dropped due to pie
							{
								e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
								++removed1;
								snp.start= snp.end=-1;
								continue;
							}
						e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
						continue;
					}//end error1

				int snp1 = 0, snp2 = 0, hlength = 0;
				float noOfOppHom = 0;
				if( T_TRIM )
					{
						snp1 = e_obj.getNewSnp( temp1 );
						snp2 = e_obj.getNewSnp( temp2 );
						hlength = snp2 - snp1;
						if( hlength <= 0 )
							{
								hlength = 1;
							}
						noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
					}
				//update drop order 2/26/14
				if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
					{
						++removed2;
						snp.start= snp.end=-1;
						continue;
					}
				if( per_err > per_err_threshold )
					{
						++removed1;
						snp.start= snp.end=-1;
						continue;
					}
				if( option.compare("MovingAverages")==0 ) //make this ma2
					{
						if( T_TRIM)
							{
								e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
							}
						else
							{
								e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
							}
						continue;
					}
				if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
					{
						if( T_TRIM)
							{
								e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
							}
						else
							{
								e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
							}
						continue;
					}
				if ( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
					{
						++removed3;
						snp.start= snp.end=-1;
						continue;
					}
				if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
					{
						e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
					}


}//end of while loop

   std::vector<float>maxes;
   float cutoff = empirical_ma_threshold;
   if(empirical_ma_threshold < 0.0)
		{
			maxes = e_obj.getMaxAverages();
			std::sort(maxes.begin(),maxes.end());
			e_obj.setMaxAverage(maxes);
			cutoff = e_obj.getXthPercentile(ma_threshold);
		}
   e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH
	   //
	std::sort( consolidator.m_errors.begin(), consolidator.m_errors.end() );
	std::sort( consolidator.m_holdOutErrors.begin(), consolidator.m_holdOutErrors.end() );
	std::string str =  " \n No of elements in error check are: "
					  + NumberToString( consolidator.m_errors.size() );
	str  = str + " \n No of elements in hold  error check are: "
					  + NumberToString( consolidator.m_holdOutErrors.size() );
	e_obj.log( str );




	if( option.compare("weightedOutput") == 0 )
		{
			float snp_average_count = 0.0;
			int start_position = consolidator.find_genome_min();
			int end_position = consolidator.find_genome_max();
			int genome_length = (end_position - start_position);
			consolidator.genome_vector.resize(genome_length,0);
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
				consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
				}
			snp_average_count = consolidator.average_snp_count();
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
				consolidator.m_weighted_sh[i].snp_weight = consolidator.update_snp_weight(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
				}
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
				consolidator.m_weighted_sh[i].final_weight = ( snp_average_count / (consolidator.m_weighted_sh[i].snp_weight));
					 e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].final_weight);
				}
		}
         /*End weighted output*/
	str = " \n No of matches removed due to percentage error: " + NumberToString( removed1 );
	str = str+ " \n No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
	str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
	str = str+ " \n No of matches removed due to initial length: " + NumberToString( removed4 );
	e_obj.log( str );


}

else // for (	(situation_no == 2)	&&	(situation_no == 4)		)
{
	//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, T_ERROR, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
	ErrorCalculator &e_obj = eCalculator;
	float PIElength = PIELENGTH;
	//bool holdOut=T_ERROR;
	int window= WINDOW;
	float ma_threshold=MA_THRESHOLD;
	float empirical_ma_threshold=EMPIRICAL_MA_RESULT;
	std::string path=BMATCHFILE;
	int person_count=pers_count;
	int trueSNP=TRUESNP;
	float trueCM=TRUECM;
	//ErrorCalculator& e_obj=eCalculator;
	//int window=WINDOW;
	int ma_snp_ends=MA_SNP_END;
	//float ma_threshold=MA_THRESHOLD;
	int min_snp=MIN_SNP;
	float min_cm=MIN_CM;
	float per_err_threshold=PCT_ERR_THRESHOLD;
	std::string option=OPTION;
	float hThreshold=HO_THRESHOLD;
	//bool holdOut=false;
	float empirical_threshold=EMPIRICAL_MA_RESULT;
	//std::string path=BMATCHFILE;
	float empirical_pie_threshold=EMPIRICAL_PIE_RESULT;
	int removed1 =0, removed2 = 0, removed3 = 0, removed4 = 0;
	int not_removed = 0;
	int total_count = consolidator.global_initial; //for keeping track of total #SH
	//bool wrongOption = false;	//Not used anywhere in the actual version
	float per_err_threshold1;
	//begin streaming in bmatch data
	unsigned int pid[2];
	unsigned int sid[2];
	unsigned int dif,hom[2];
	std::ifstream file_bmatch(path.c_str(),std::ios::binary);
	int pers_count = person_count;
	SNP_lrf snp;
	if( !file_bmatch )
		{
			std::cerr<<"unable to open the bmatch file, exiting the program" << std::endl;
			exit( -1 );
		}
	//will need to replace all instances of m_matches in this function with pid[0], pid[1], etc.
	//will probably encapsulate all of the below code into a while !eof loop


	// Consolidator cons;	//glf object created for supporting m_errors and m_holdOutErrors
	if(empirical_pie_threshold >= 0.0)
		{
			per_err_threshold1 = empirical_pie_threshold;
		}
	else
		{
	    	per_err_threshold1 = consolidator.getPctErrThreshold( per_err_threshold );
		}
	str = " \n Corresponding empirical error threshold for " + NumberToString( per_err_threshold ) + " is: " + NumberToString( per_err_threshold1 )  + " \n";
	float hThreshold1 = 0;
	if( T_TRIM )//holdout performtrim
		{
			hThreshold1 = consolidator.getHoldOutThreshold( hThreshold );
			str = str + " \n corresponding holdout emperical error threshold for " + NumberToString( hThreshold ) + " is: " + NumberToString( hThreshold1 )  + " \n";
		}
	e_obj.log( str );
	per_err_threshold = per_err_threshold1;
	hThreshold = hThreshold1;

if (T_TRIM == false)
{

    if( OPTION.compare( "Error3") == 0 )
    {
  	  std::cerr<< " Error: You have provided option:Error3 " << std::endl
             <<" you can use Error3 only if you provided"
             <<" hold out ped and map file, program with not output anything" << std::endl;
        exit( -1 );
    }
}
else
{
	std::cerr<<"entering into HoldOut Mode Matching"<<std::endl;
}

	while ( !file_bmatch.eof())
		{
			//read in bmatch data
			pid[0] = -1;
			file_bmatch.read( (char*) &pid[0] , sizeof( unsigned int ) );
			if ( pid[0] == -1 ) continue;
			file_bmatch.read( (char*) &pid[1] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &sid[0] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &sid[1] , sizeof( unsigned int ) );
			file_bmatch.read( (char*) &dif , sizeof( int ) );
			file_bmatch.read( (char*) &hom[0] , sizeof( bool ) );
			file_bmatch.read( (char*) &hom[1] , sizeof( bool ) );
			if(pid[0]>=pers_count||pid[1]>=pers_count)
				{
					std::cerr<<"problem with bsid file, check it please/ problem with bsid file in findTrueSimplePctErrors, check it please"<<std::endl;
					return;
				}

			snp.start=sid[0];	//std::cout<<snp.start<<std::endl;//glf
			snp.end=sid[1];		//std::cout<<snp.end<<std::endl;//glf
			if(pid[0] > pid[1])
				{
					int temp = pid[0];
					pid[0] = pid[1];
					pid[1] = temp;
				}
			int i = pid[0];
			int j = pid[1];
			if(snp.end==-1)
				{
					continue;
				}


			//check to see if this SH is trulyIBD or not. do not currently need this, as germline input will be trueIBD
			if( ( e_obj.getCMDistance( sid[ 1 ] ) - e_obj.getCMDistance( sid[ 0 ] ) ) >= trueCM && ( sid[ 1 ] - sid[ 0 ] ) >= trueSNP &&  pid[0] != pid[1] )
				{
					//this indicates a trulyIBD segment, so we can go to work on it with the following code. Otherwise, just continue.
					//If you need to unroll these loops as well, you would add the stream here, and pluck out items that are
					//trulyIBD. This would result in streaming the data twice, as I don't think that there is a way to run non-IBD segments
					//without first having run all of the trulyIBD segments.
					//-------------------------------------------------------------------------------------------------
					int t1 = snp.start +( snp.end -snp.start ) * 0.25;
					int t2 = snp.end -( snp.end - snp.start ) * 0.25;

					std::vector<std::vector<int> > trueErrors=e_obj.checkErrors( i, j, t1, t2);
					std::vector<int>finalTrueErrors=e_obj.getFinalErrors( trueErrors );
					//This section handles finding the maximum moving averages amongst trulyIBD segments
					std::vector<float> av;
					float current_max;
					if(empirical_ma_threshold < 0.0)
						{ //this would indicate that the user did NOT supply an empirical-ma-threshold argument. So we need to calculate one.
							av = e_obj.getTrueMovingAverages(finalTrueErrors,t1,t2,window);
							current_max = av[0];
							for(int q = 1; q < av.size(); q++)
								{
								  if(av[q] > current_max)
									  {
										current_max = av[q];
									  }
								}
							e_obj.addMaxAverage(current_max);
						}
					//------------------------------------------------------------------------------
					int temp1 = snp.start;
					int temp2 = snp.end;
					float startCM = e_obj.getCMDistance( temp1 );
					float endCM = e_obj.getCMDistance( temp2 );
					float mid1CM = startCM + ( endCM - startCM ) / 2 - PIElength / 2;
					float mid2CM = startCM + ( endCM - startCM ) / 2 + PIElength / 2;
					while( e_obj.getCMDistance( temp1 ) <= mid1CM || e_obj.getCMDistance( temp2 ) >=mid2CM )
						{
							if( e_obj.getCMDistance( temp1 ) <= mid1CM )
								{
								  ++temp1;
								}
							if( e_obj.getCMDistance( temp2 ) >=mid2CM )
								{
								  --temp2;
								}
						}
					std::vector<std::vector<int> > errors=e_obj.checkErrors( i, j, temp1, temp2);
					std::vector<int>finalErrors=e_obj.getFinalErrors( errors );
					//                  float per_err = e_obj.getThreshold(finalErrors,temp1, temp2, 0 );
					float per_err = e_obj.getThreshold(finalErrors,temp1,temp2); //overload!
					consolidator.m_errors.push_back( per_err );
					if( T_ERROR  )//holdout in findTrueSimplePctErrors
						{
							float oppHom = ( e_obj.getOppHomThreshold( i, j, temp1, temp2 ) ) / ( temp2 -temp1 );
							consolidator.m_holdOutErrors.push_back( oppHom );
						}
				}
			else
				{//end check for trueIBD segment
					//continue; //this segment is not trueIBD, so go onto next SH
				}

			//TRIM
		    int temp1=snp.start;
		    int temp2=snp.end;
		    int pers1 = pid[0], pers2 = pid[1];
			if( option.compare( "ErrorRandom1" ) == 0 || option.compare( "ErrorRandom2" ) == 0 || option.compare( "ErrorRandom3" ) == 0 )
		    {
		      pers1 = std::rand() % e_obj.getNoOfPersons();
		      pers2 = std::rand() % e_obj.getNoOfPersons();
		      if( pers1 > pers2 )
		      {
		        pers1 = pers1 + pers2;
		        pers2 = pers1 - pers2;
		        pers1 = pers1 - pers2;
		      }
		    }

		    std::vector<std::vector<int> > errors=e_obj.checkErrors(pers1, pers2, temp1, temp2);
		    std::vector<int>finalErrors=e_obj.getFinalErrors(errors);
		    /*Inject implied error at start/end of SH here*/
		    std::vector<int>::iterator it;
		    it = finalErrors.begin(); //go to the start of the vector
		    if(finalErrors[0] != 1){
		      finalErrors.insert(it,1); //inject an error at position 1, if not already there
		    }
		    /*End inject implied error section*/
		    std::vector<int>trimPositions;
		    std::vector<float>movingAverages;
		    float threshold;


		    //length 32 should survive this, since the reduced arg is 32 (although the reduced cM may be the culprit)

		    if( (e_obj.isInitialCmDrop(temp1,temp2,min_cm)) || ((temp2-temp1) < min_snp) ){ //initial drop. Don't calculate MA
		      trimPositions.push_back(temp1);
		      trimPositions.push_back(temp2);
		      trimPositions.push_back(1);
		    }
		    else
		    {
		        movingAverages = e_obj.getMovingAverages(finalErrors,temp1,temp2,window);
		        if(empirical_threshold < 0.0){
		          threshold = e_obj.getCutoff(); //this is the empirical average threshold for moving averages
		        } else {
		          threshold = empirical_threshold;
		        }
		        trimPositions = e_obj.getTrimPositions(movingAverages,temp1,temp2,threshold,min_cm);
		    }
		    //-----------------
		    int beforeTrimStart = temp1;
		    int beforeTrimEnd = temp2;
		    snp.end = temp2 = temp1+trimPositions[1];
		    snp.start = temp1 = temp1+trimPositions[0];
		    int del0 = trimPositions[0];
		    int del1 = trimPositions[1];
		    //in order to update PIE calculation, may have to alter the finalErrors vector before passing it into getThreshold(). But the alteration
		    //would need to occur at the start and end of the trimmed SH, not just at the start and end points of the initial SH.
		    float per_err = e_obj.getThreshold(finalErrors,del0,del1,ma_snp_ends);
		    //add new weighted option
		    /*
		     For this new option, we only output SH that are not dropped. So, the output is finalOutput + weighted column.
		    */
		    if( (option.compare("weightedOutput") == 0) ){
		      int snp1 = 0, snp2 = 0, hlength = 0;
		      float noOfOppHom = 0;
		      /*if( holdOut )
		      {
		        snp1 = e_obj.getNewSnp( temp1 );
		        snp2 = e_obj.getNewSnp( temp2 );
		        hlength = snp2 - snp1;
		        if( hlength <= 0 )
		        {
		          hlength = 1;
		        }
		        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, m_matches[i][j][l].start, m_matches[i][j][l].end );
		      }*/
		      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
		        removed4++;
		        continue;
		      }
		      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){ //removed2 a tpos.size of 3 indicates trimming due ot cM
		        removed2++;
		        continue;
		      }
		      if( per_err > per_err_threshold){
		        removed1++;
		        continue;
		      }
		      if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
		        removed3++;
		        continue;
		      } //removed3

		      not_removed++;
		      consolidator.m_weighted_sh.push_back(Weighted_SH(temp1,temp2,i,j)); //build the vector of SH that passed
		      continue;
		    }//end weghtedOutput
		    if( (option.compare("FullPlusDropped") == 0) ){
		      int snp1 = 0, snp2 = 0, hlength = 0;
		      float noOfOppHom = 0;
		      if( T_TRIM )
		      {
		        snp1 = e_obj.getNewSnp( temp1 );
		        snp2 = e_obj.getNewSnp( temp2 );
		        hlength = snp2 - snp1;
		        if( hlength <= 0 )
		        {
		          hlength = 1;
		        }
		        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
		      }


		      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
		        std::vector<float>movingAverages;
		        temp1 = beforeTrimStart;
		        temp2 = beforeTrimEnd;
		        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,1);
		        snp.start= snp.end=-1;
		        removed4++;
		        continue;
		      }

		      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){
		        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,2);
		        removed2++;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      if( per_err > per_err_threshold){
		        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,3);
		        removed1++;
		        continue;
		      }
		      if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
		        e_obj.fullPlusDroppedOutput(i,j,temp1,temp2,min_snp,min_cm,finalErrors,per_err,4);
		        removed3++;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
		      continue;
		    } //end FullPlusDropped

		    //add finalOutput here
		    if(option.compare("finalOutput") == 0){
		    //only print if it passes all of the tests.
		      int snp1 = 0, snp2 = 0, hlength = 0;
		      float noOfOppHom = 0;
		      if( T_TRIM )
		      {
		        snp1 = e_obj.getNewSnp( temp1 );
		        snp2 = e_obj.getNewSnp( temp2 );
		        hlength = snp2 - snp1;
		        if( hlength <= 0 )
		        {
		          hlength = 1;
		        }
		        noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
		      }
		      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){
		        std::vector<float>movingAverages;
		        temp1 = beforeTrimStart;
		        temp2 = beforeTrimEnd;
		        snp.start= snp.end=-1;
		        removed4++;
		        continue;
		      }
		      if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) ){
		        removed2++;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      if( per_err > per_err_threshold){
		        removed1++;
		        continue;
		      }
		      if( T_TRIM && hThreshold < ( noOfOppHom ) / hlength ){
		        removed3++;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      e_obj.finalOutPut(i,j,temp1,temp2,min_cm);
		      continue;
		    }//end finaloutput

		    //Calculate Error1
		    if( (option.compare("Error1") == 0 ) || (option.compare("ErrorRandom1") == 0) || (option.compare("Error") == 0) ){

		      if( ( (beforeTrimEnd - beforeTrimStart) < min_snp) || ( (trimPositions.size() == 3) && (trimPositions[2] == 1) ) ){ //dropped before trimming
		        //don't bother printing out ma for this one. But go back and change it so that it doesn't actually calc it
		        std::vector<float>movingAverages;//null
		        //trying something special in this case. This can be removed once idrops aren't being trimmed
		        //test code
		        temp1 = beforeTrimStart;
		        temp2 = beforeTrimEnd;
		        //
		        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,1);
		        snp.start= snp.end=-1;
		        removed4++; //seems ok
		        continue;
		      }
		      if( (( temp2-temp1 ) < min_snp) || ((trimPositions.size() == 3) && (trimPositions[2] == 2) ) ) //dropped after trimming
		      {
		        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,2);
		        ++removed2;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      if( per_err > per_err_threshold ) //dropped due to pie
		      {
		        e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,3);
		        ++removed1;
		        snp.start= snp.end=-1;
		        continue;
		      }
		      e_obj.errorOutput(i,j,temp1,temp2,min_snp,min_cm,movingAverages,finalErrors,per_err,temp1,temp2,beforeTrimStart,beforeTrimEnd,0);//no drop
		      continue;
		    }//end error1

		    int snp1 = 0, snp2 = 0, hlength = 0;
		    float noOfOppHom = 0;
		    if( T_TRIM )
		    {
		      snp1 = e_obj.getNewSnp( temp1 );
		      snp2 = e_obj.getNewSnp( temp2 );
		      hlength = snp2 - snp1;
		      if( hlength <= 0 )
		      {
		        hlength = 1;
		      }
		      noOfOppHom = e_obj.getOppHomThreshold( pers1, pers2, snp.start, snp.end );
		    }
		    //update drop order 2/26/14
		    if( (( temp2-temp1 ) < min_snp) || (trimPositions.size() == 3) )
		    {
		      ++removed2;
		      snp.start= snp.end=-1;
		      continue;
		    }
		    if( per_err > per_err_threshold )
		    {
		      ++removed1;
		      snp.start= snp.end=-1;
		      continue;
		    }
		    if( option.compare("MovingAverages")==0 ) //make this ma2
		    {
		      if( T_TRIM)
		      {
		        e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp,min_cm,movingAverages,trimPositions,per_err, noOfOppHom, hlength );
		      }
		      else
		      {
		        e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm,movingAverages, trimPositions,per_err );
		      }
		      continue;
		    }
		    if(option.compare("Error2")==0 || option.compare( "ErrorRandom2" ) == 0)
		    {
		      if( T_TRIM)
		      {
		        e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
		      }
		      else
		      {
		        e_obj.middleOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err);
		      }
		      continue;
		    }
		    if ( T_TRIM && hThreshold < ( noOfOppHom ) / hlength )
		    {
		      ++removed3;
		      snp.start= snp.end=-1;
		      continue;
		    }
		    if( option.compare("Error3")==0 || option.compare( "ErrorRandom3" ) == 0  )
		    {
		      e_obj.middleHoldOutPut(i,j,temp1,temp2, min_snp, min_cm, finalErrors, trimPositions, per_err, noOfOppHom, hlength );
		    }





		}//end while looop --------------------------------(!eof)




	//this section actually handles the sorting of the max averages, and the setting of the user supplied percentile.
	std::vector<float>maxes;
	float cutoff = empirical_ma_threshold; //assume the user wanted to supply a value. This value will be overwritten shortly if they did not.
	if(empirical_ma_threshold < 0.0)
		{
			maxes = e_obj.getMaxAverages();
			std::sort(maxes.begin(),maxes.end());
			e_obj.setMaxAverage(maxes);
			cutoff = e_obj.getXthPercentile(ma_threshold); //<-make that an actual user input value
		}
	//std::cout<<maxes.size()<<"\t"<<m_errors.size()<<std::endl;//glf
	e_obj.setCutoff(cutoff);//set the actual threshold to be used when calculating MA in all other SH

/*	if(empirical_ma_threshold < 0.0)
		{
			std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
			outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
			e_obj.log(outt);
		}
	else
		{
			std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
			outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
			e_obj.log(outt);
		}*/
	std::sort( consolidator.m_errors.begin(), consolidator.m_errors.end() );
	std::sort( consolidator.m_holdOutErrors.begin(), consolidator.m_holdOutErrors.end() );
	std::string str =  " \n No of segments deemed to be IBD for finding empirical error threshold " + NumberToString( consolidator.m_errors.size() );
	str  = str + " \n No of segments deemed to be IBD for finding empirical error threshold in hold out are " + NumberToString( consolidator.m_holdOutErrors.size() );
	e_obj.log( str );

	if(empirical_ma_threshold < 0.0)
		{
			std::string outt = "\n User supplied ma-threshold is: " + NumberToString(ma_threshold);
			outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
			e_obj.log(outt);
		}
	else
		{
			std::string outt = "\n User supplied empirical-ma-threshold is: " + NumberToString(cutoff);
			outt = outt + "\n Moving Averages will be tested usign the empirical threshold: " + NumberToString(cutoff);
			e_obj.log(outt);
		}
	//-------------------------------TRIM CODE----------------------------------------------
	if( option.compare("weightedOutput") == 0 )
		{
			float snp_average_count = 0.0;
			int start_position = consolidator.find_genome_min();
			int end_position = consolidator.find_genome_max();
			int genome_length = (end_position - start_position);
			consolidator.genome_vector.resize(genome_length,0);
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
					consolidator.update_genome(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
				}
			snp_average_count = consolidator.average_snp_count();
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
					consolidator.m_weighted_sh[i].snp_weight = consolidator.update_snp_weight(consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2);
				}
			for(int i = 0; i < consolidator.m_weighted_sh.size(); i++)
				{
					consolidator.m_weighted_sh[i].final_weight = ( snp_average_count / (consolidator.m_weighted_sh[i].snp_weight));
					e_obj.weightedOutput(consolidator.m_weighted_sh[i].per1, consolidator.m_weighted_sh[i].per2, consolidator.m_weighted_sh[i].snp1, consolidator.m_weighted_sh[i].snp2, consolidator.m_weighted_sh[i].final_weight);
				}
		}
	/*End weighted output*/
	str = " \n No of matches removed due to percentage error: " + NumberToString( removed1 );
	str = str+ " \n No of matches removed due to length of trimming by moving averages: " + NumberToString( removed2 );
	str = str+ " \n No of matches removed due hold out ped file checking: "+ NumberToString( removed3 );
	str = str+ " \n No of matches removed due to initial length: " + NumberToString( removed4 );

	e_obj.log( str );


}//end of else
if (T_TRIM == true)
{
	std::cerr<<" Main Trim operation has completed "<< std::endl;
	std::cerr<< " Hold out trim has completed" <<std::endl;
	if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
	{
	   consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
	}
}

				//-----------------------------glf--------------------------------
time_t endTime;
time (&endTime);
std::string str2 = " The program ended at: " + std::string( ctime ( &endTime ) );
str2 = str2 +  "  Total time ( in seconds): " + NumberToString( ( endTime - startTime ) );
eCalculator.log( str1 );
eCalculator.log( str2 );


}


void ErrorFinderManager::initiateErrorFinder_PRE()
{
	eCalculator.readBmidFile(BMIDFILE);
	std::cerr<<"Reading bmid file completed"<<std::endl;
	eCalculator.readBsidFile(BSIDFILE);
	std::cerr<<"Reading bsid file completed"<<std::endl;
	eCalculator.readPedFile(PEDFILE, HO_MISSING);
	std::cerr<<"Reading ped file completed"<<std::endl;
/*	int pers_count=eCalculator.getNoOfPersons();*/
	consolidator.setPersonCount(eCalculator.getNoOfPersons());
}

void ErrorFinderManager::initiateErrorFinder_CURRENT(int pers_count)
{

	if(		( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	ISMOL)
		{
			SITUATION_NO	=	1;
			T_ERROR	=	true;
			//T_TRIM = true;
			eCalculator.changeMapFile( HMAPFILE );
			eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
			std::cerr<< " new map and ped File has been read" <<std::endl;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTruePctErrors( eCalculator, MA_SNP_END, T_VAL, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count, TRUESNP, TRUECM);//T_VAL = true
			CERRMESSAGE	=	" true percentage errors calculated ";

		}
	else if(	( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	!ISMOL)
		{
			SITUATION_NO	=	2;
			T_ERROR	=	true;
			//T_TRIM = true;
			eCalculator.changeMapFile( HMAPFILE );
			eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
			std::cerr<< " new map and ped File has been read" <<std::endl;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, T_VAL, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );//T_VAL = true
			CERRMESSAGE	=	" true percentage errors calculated ";
		}
	else if(	!( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	ISMOL)
		{
			SITUATION_NO	=	3;
			T_ERROR	=	false;
			//T_TRIM = false;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
			CERRMESSAGE	=	" true hold out percentage errors calculated ";


		}
	else if(	!( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ) &&	!ISMOL)
		{
			SITUATION_NO	=	4;
			T_ERROR	=	false;
			//T_TRIM = false;
			std::cerr<< " calculating true percentage errors" <<std::endl;
			//consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
			CERRMESSAGE	=	" true hold out percentage errors calculated ";
		}
	else
		{

		}

	if (( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 ))
		{
			T_TRIM = true;
		}
	else
		{
			T_TRIM = false;
		}



	if (!SITUATION_NO)
		{
			std::cerr<<"Wrong Situation:HPED HMAP, ISMOL Error "<<std::endl;
			exit(0);
		}


}





void ErrorFinderManager::displayError(std::string argv)
{
	std::cerr<<"these parameters are not allowed "<<wrongParam<<std::endl;
	std::cerr << "Usage: " << argv << " -bmatch [BMATCH FILE]  -bsid [BSID FILE] -bmid [BMID FILE] -reduced [min_snp] [min_cm] "
  <<" -ped-file [ped file] -window [window width to calculate moving averages] "
  <<" -gap [max gap to consolidate two matches]"
  <<" -pct-err-threshold [max percentage of errors in a match after the trim] OR -emp-pie-threshold" 
  <<" -ma-threshold [specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"
  <<" Note that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"
  <<"-output.type [ must provide any of these. it can be "
  << "MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 " 
  << "or ErrorRandom2 or Error3 or ErrorRandom3 or Full "
  <<  "look at the description about how these works in wiki ]"
  << "(optional) -holdout-ped [new ped file path] -holdout-map [new map file] "
  << "-holdout-threshold [threshold to drop a match with new ped file ]"
  << " -holdout-missing [missing value representation in new ped file] "
  << " -log.file [log file name]"
  << " -trueCM [ true match maximum cm length] " 
  << " - trueSNP [ true match SNP length]"
  << " -PIE.dist.length [ can be MOL or any cm distance length "
  << "please refer wiki for more details on how to use this option"
  << "-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"
  << std::endl;
}
/*void ErrorFinderManager::initiateErrorFinder()
{
        eCalculator.readBmidFile(BMIDFILE);
        std::cerr<<"Reading bmid file completed"<<std::endl;
        eCalculator.readBsidFile(BSIDFILE);
        std::cerr<<"Reading bsid file completed"<<std::endl;
        eCalculator.readPedFile(PEDFILE, HO_MISSING);
        std::cerr<<"Reading ped file completed"<<std::endl;
        int pers_count=eCalculator.getNoOfPersons();
        if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
        {
          eCalculator.changeMapFile( HMAPFILE );
          eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
          std::cerr<< " new map and ped File has been read" <<std::endl;
          std::cerr<< " calculating true percentage errors" <<std::endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count, TRUESNP, TRUECM);
          //    consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //  consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          std::cerr<< " true percentage errors calculated "<<std::endl;
       }
       else
       {
    	   std::cerr<< " calculating true percentage errors" <<std::endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          std::cerr<< " true hold out percentage errors calculated "<<std::endl;
      }   
  
}*/
