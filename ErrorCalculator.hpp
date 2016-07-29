
#ifndef ERRORCALCULATOR_HPP_
#define ERRORCALCULATOR_HPP_

//#include "Consolidator.hpp"
//#include "ErrorFinderManager.hpp"


#include <fstream>
#include <sstream>
#include <string>
#include <vector>
/*class SNPs;
class Marker;*/
class SNPs_lrf
{
       public:
       bool SNP1,SNP2;

};
class Marker
{
	public:
        short   chr;
        std::string  rsid;
        float   cm_distance;
        long    bp_distance;
};
class ErrorCalculator
{
public:

         ErrorCalculator():pers_count(0){}
         void createLogFile( std::string path  );
         void log( std::string& log );
         void readBmidFile(std::string path);
         void readBsidFile(std::string path);
         void readPedFile(std::string path, std::string missing);
         void finalOutPut( int pers1,int pers2,int snp1,int snp2,float min_cm );
         void weightedOutput( int pers1, int pers2, int snp1, int snp2, float weight);
         
         //Here 1
         template < class T >
         void fullPlusDroppedOutput( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,int reason );
         //Here 2
		 template < class T >
		 void middleOutPut( int pers1,int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions, float pct_err, int start, int end);
		 //Here 3
		 template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm, std::vector<T >positions, float pct_err );

		 //Here 4
		 template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm, std::vector< T >positions, std::vector< int > trims, float pct_err );

		 //Here 5
		 template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm,std::vector< T >positions, float pct_err, std::string reason);

		 //Here 6
		 template < class T >
         void middleOutPut( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm,std::vector< T >positions,std::vector< int > trims, float pct_err,int start, int end );

		 //Here 7
         template < class T >
         void middleHoldOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm, std::vector< T >positions,std::vector< int > trims,float pct_err, int nSNPS, int length );
         //this will be used for the new "erro" output option

         //Here 8
         template < class T >
		 void errorOutput(int p1, int p2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T >positions,std::vector< int > trims, float pct_err, int startTrim, int endTrim, int start, int end, int reason);

		 template < class T >
		 void errorOutput(int p1, int p2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T >positions,std::vector< int > trims, float pct_err, int startTrim, int endTrim, float pie, float maMAX, float randomPie, float randomMAmax);

         std::vector<int>getFinalErrors(std::vector<std::vector< int> > errors )const;

         std::vector<std::vector <int> >checkErrors(int pers1,int pers2,int snp1,int snp2)const;

         int getMax(int &a,int &b,int &c,int &d)const;

         std::vector<int>getTrimPositions(std::vector<float> averages,int snp1,int snp2,float threshold,float minLength );
         
         std::vector<int>getTrimPositions(std::vector<float> averages,int snp1,int snp2, float threshold_start, float threshold_end,float minLength,int ma_snp_ends);

         std::vector<float>getMovingAverages(std::vector<int> errors,int snp1,int snp2,int width);

         std::vector<float>getTrueMovingAverages(std::vector<int> errors,int snp1,int snp2,int width);
      
         std::vector<float>getTrueMovingAverages2(std::vector<int> errors,int snp1,int snp2, int width);
 
         float getThreshold(std::vector<int>trimPositions,std::vector<int>finalErrors );

		 float getThreshold(std::vector<int> finalErrors,int snp1,int snp2, int ma_err_ends );

		 float getThreshold(std::vector<int> finalErrors, int snp1, int snp2);

		 int getNoOfPersons()
			 {
				 return pers_count;
			 }

         virtual ~ErrorCalculator(){}; 

         float getCMDistance(int position);

	     void changeMapFile( std::string mapfile);
            
         void readHPedFile( std::string pedfile, std::string missing );

         float getOppHomThreshold( int pers1, int pers2, int snp1, int snp2 ) const;

         int getNewSnp ( int );
  
         void countGapErrors( bool );

         std::vector<float>getMaxAverages()
			 {
        	 	 return max_averages;
			 }//new

         void addMaxAverage(float av)
			 {
				 max_averages.push_back(av);
			 }//new

         void setMaxAverage(std::vector<float> av)
			 {
        	 	 max_averages = av;
			 }
         float getXthPercentile(float x);//new

         void setCutoff(float c)
			 {
				 cutoff = c;
			 }

		 float getCutoff()
			 {
				 return cutoff;
			 }
		 bool isInitialCmDrop(int snp1, int snp2, float minLength);

		 float newCmLength(int snp1, int snp2, float cut_length);

		 float adjustCmLength(int snp, float val);

		 int snpFinder(float cm_length,int start_snp, int end_snp, int direction);

private:
         int pers_count;
         std::vector<std::vector<SNPs_lrf > >ped_file;
         std::vector<std::vector<SNPs_lrf > >hped_file;
         std::vector< int > mapper;
         std::vector<Marker>marker_id;
         std::vector<Marker>hMarker_id;
         std::vector<std::string>sample_id;
         std::vector<float>max_averages; //new nate 2/3/2014
         std::ofstream m_logger;
         bool countGapError;
         float cutoff;
};    

template <class T>
void ErrorCalculator::fullPlusDroppedOutput( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,int reason ){
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err;
        std::cout << "\t"<< reason << std::endl;
}

template <class T>
void ErrorCalculator::middleOutPut(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector<T> positions, float pct_err, int start, int end){
	if( (snp1==-1) || (snp2==-1) ) return;
	if( (sample_id.size() <= (pers1*2)) || (sample_id.size() <= (pers2*2)) ){
		std::cerr << "Something went wrong with the bsid file while outputting data." << std::endl;
		return;
	}
	std::cout << sample_id[pers1*2] << "\t" << sample_id[pers2*2] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
	     << start << "/" << end << "\t" << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance - marker_id[snp1].cm_distance) << "\t" << "/";
	for(int i = 0; i < positions.size(); i++){
		std::cout << positions[i] << "/";
	}
	std::cout << "\t" << pct_err << std::endl;
}

template <class T>
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2, int min_snp, float min_cm,std::vector<T> positions, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
	//in theory, changing this to output trim-dropped data shouldn't affect anything. But if it does, looks here first,
	//and consider overloading this function(yet again).
        /*if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {

             return;
        }*/

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
		<<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err <<std::endl;

}

template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,int snp1,int snp2,int min_snp, float min_cm, std::vector<T> positions,float pct_err,std::string reason){
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
    <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err;
        std::cout << "\t"<< reason <<std::endl;
}


template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp,
                                  float min_cm,std::vector< T > positions,
                                  std::vector< int > trims, float pct_err )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
        	std::cout << "Testing a theory. " << std::endl;
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
        	std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err <<std::endl;

}

template < class T >
void ErrorCalculator::middleOutPut( int pers1,int pers2,
                                  int snp1,int snp2, int min_snp,
                                  float min_cm,std::vector< T > positions,
                                  std::vector< int > trims, float pct_err,int start,int end )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
/*        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
             cout << "Testing a theory. " << endl;
             return;
        }
*/
        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"
	       << pct_err << "\t" << start <<"/" << end << "/";
        for(int i=0;i<trims.size();++i)
        {
        	std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
//        cout<<"\t"<< pct_err <<endl;

}

template < class T >
void ErrorCalculator::middleHoldOutPut( int pers1,
                                    int pers2, int snp1, int snp2,
                                    int min_snp, float min_cm,
                                    std::vector< T > positions, std::vector< int > trims,
                                    float pct_err, int numberofOppHoms,
                                    int length )
{
        if(snp1==-1) return;
        if(snp2==-1)
        {
           return;
        }
        if((marker_id[snp2].cm_distance - marker_id[snp1].cm_distance)<min_cm || (snp2-snp1) < min_snp )
        {
             return;
        }

        if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
        	std::cerr<<"some thing went wrong with bsid file while outputing"<<std::endl;
               return;
        }

        std::cout<<sample_id[pers1*2]<<"\t"
             <<sample_id[pers2*2]<<"\t"
                <<marker_id[snp1].bp_distance<<"\t"
             <<marker_id[snp2].bp_distance<<"\t"
               <<(snp2-snp1)<<"\t"
               <<(marker_id[snp2].cm_distance-marker_id[snp1].cm_distance)<<"\t"<<"/";
        for(int i=0;i<trims.size();++i)
        {
        	std::cout<<trims[i]<<"/";
        }
        std::cout<<"\t";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout<<positions[i]<<"/";
        }
        std::cout<<"\t"<< pct_err << "\t"<< numberofOppHoms << "\t" << length <<std::endl;

}

template < class T > void ErrorCalculator::errorOutput(int pers1, int pers2, int snp1, int snp2, int min_snp, float min_cm, std::vector< T > positions,std::vector< int > errors, float pct_err, int startTrim,int endTrim, int start, int end, int reason)
{
	if( (snp1 == -1) || (snp2 == -1) ) return;
	if(sample_id.size()<=pers1*2||sample_id.size()<=pers2*2)
        {
		std::cerr<<"something went wrong with bsid file while outputing"<<std::endl;
               return;
        }
/*	cerr << "ABOUT TO COUT" << endl;//debug
	cerr << "Pers1 is: " << sample_id[pers1] << endl;
	cerr << "But the snps are: " << snp1 << ", " << snp2 << endl;
	cerr << "Therefore for snp1: " << marker_id[snp1].bp_distance << endl;
	cerr << "And snp2: " << marker_id[snp2].bp_distance << endl;*/
	std::cout << sample_id[pers1*2] << "\t" << sample_id[pers2*2] << "\t" << marker_id[snp1].bp_distance << "\t" << marker_id[snp2].bp_distance << "\t"
             << (snp2-snp1) << "\t" << (marker_id[snp2].cm_distance-marker_id[snp1].cm_distance) << "\t" << pct_err << "\t" << start << "/" << end << "\t" << startTrim <<"/" << endTrim << "\t" << reason
	     << "\t" << "/";
	for(int i=0;i<errors.size();++i)
        {
		std::cout << errors[i] << "/";
        }
	std::cout << "\t";
        for(int i=0;i<positions.size();++i)
        {
        	std::cout << positions[i] << "/";
        }
        std::cout << std::endl;
}
/*
class SNPs
{
       public:
       bool SNP1,SNP2;

};
class Marker
{
	public:
        short   chr;
        std::string  rsid;
        float   cm_distance;
        long    bp_distance;
};*/
#endif
