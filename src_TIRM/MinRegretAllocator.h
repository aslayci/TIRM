#ifndef _MIN_REGRET_ALLOC_H
#define _MIN_REGRET_ALLOC_H

#include <ctime>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <map>
#include <vector>
#include "Advertiser.h"
#include "anyoption.h"
#include "utils.h"
#include "TimGraph.h"

namespace _MinRegret {	
	
	typedef std::vector<Advertiser*> AdvertiserList;
	typedef std::vector<TimGraph*> TimGraphList;	
	
	class MinRegretAllocator {
		
	protected:
		
		AdvertiserList *advList;
		TimGraphList *timList;			
		
		multimap< float,int> allocQueue; // arranges allocation priority <reduction in regret, advertiserID>
		
		AnyOption *opt;
		
		int userAttentionConstraint; 
		int nrCompanies;	
		int nrTopics;
		float epsilon, lambda;
		int n;
		int m;
		
		string delim; 
		ofstream *outFileStreams;	
		ofstream outMasterStream;
		string outFolderName;
		string *outFileNames; // array of output file names for each advertiser
		string outMasterName;			
		
		time_t startTime;	// to keep the start-time of the algorithm
		// 		time_t *advStartTimes; //to keep the advertiser-specific start-times
		
		struct stat info;
		
		
	public:
		
		MinRegretAllocator(AnyOption* opt);
		~MinRegretAllocator();				
		
		void readKappasFile();
		void readItemDistsFile(); //z columns
		void readBudgetsFile(); // budget and cpc per node		
		void readTICGraph(); // keeps in transposed order
		void readNodeCTR();
		void allocateSeeds();
		
		void arrangeOutputFiles();
		void openOutputFiles(int advID); // adv-level output files
		void openOutputFiles(); // master-output file
		void writeInAdvertisersFile(int advID, int v, float mmg, float ctr, float budget, float cpe, float lambda);
		void writeInMasterOutputFile(int advertiserID, float totalCoverage, float budget, float totalRegret, int size, float duration, float memory, float lambda, float seedCosts);
        
        
		
	};	
}

#endif