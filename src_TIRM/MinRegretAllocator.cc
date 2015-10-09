#ifndef _MIN_REGRET_ALLOC_C
#define _MIN_REGRET_ALLOC_C

#include "MinRegretAllocator.h"
#include "memoryusage.h"

namespace _MinRegret {
	
	MinRegretAllocator::MinRegretAllocator(AnyOption *opt1) {
		
		opt = opt1;		
		nrTopics = strToInt(opt->getValue("nrTopics"));				
		nrCompanies = strToInt(opt->getValue("nrCompanies"));			
		userAttentionConstraint = strToInt(opt->getValue("attentionConstraint"));				
		outFolderName = opt->getValue("outputFolder"); 
		
		epsilon = strToFloat(opt->getValue("epsilon")); 
		lambda = strToFloat(opt->getValue("lambda"));

		n = strToInt(opt->getValue("n"));
		m = strToInt(opt->getValue("m"));
		
		delim = " \t";		
		
		for(int i = 0; i < n; i++) 
			graphT.push_back(vector<int>());
		
		attentionQuota = std::vector<int>(n,userAttentionConstraint);
		
		/* create Advertiser and TimGraph objects */		
		Advertiser *adv;
		TimGraph *tim;		
		advList = new AdvertiserList();
		timList = new TimGraphList();
		for(int i = 0; i < nrCompanies; i++) {
			adv = new Advertiser(i, nrTopics);
			advList->push_back(adv);						
			tim = new TimGraph(adv, epsilon, n, m, lambda);
			timList->push_back(tim);						
		} 			
		
		// read input item distributions, budgets, graph and kappas
		readNodeCTR();		
		readItemDistsFile();		
		readBudgetsFile();
		//readKappasFile();
		readTICGraph();
		
		for(int i = 0; i < nrCompanies; i++) {
			adv = advList->at(i);							
			adv->regret = adv->budget;
			adv->candidateRegret = adv->budget;
			adv->totalMonCoverage = 0.0;
		}
        
        /* kontrol theta
        for(int i = 0; i < nrCompanies; i++) {
            adv = advList->at(i);
            tim = timList->at(i);
            tim->doInitialGeneration();
        }*/
		
        arrangeOutputFiles();
        allocateSeeds();
		
	}
	
	void MinRegretAllocator::allocateSeeds() {
		
		allocQueue.clear(); 
		
		Advertiser *adv;
		TimGraph *timGraph;		
		int attentionTemp;
		infPair bests;
		
		cout << "MinRegret greedy-allocation with RR sets started for lambda = " << lambda << endl;
		time(&startTime);	// get it for the total algo running time		
		
		// initial RR sets computation and allocQueue update 
		for(int i = 0; i < nrCompanies; i++) {			
			adv = advList->at(i);	
			// cout << "doing adv " << adv->advertiserId << endl;
			timGraph = timList->at(i);
			timGraph->doInitialGeneration();
			timGraph->findBestNode(bests);
			adv->candidateRegret = (abs(adv->budget - (adv->totalMonCoverage + bests.second))) + (lambda * (float) (timGraph->seedSet.size() + 1)) ;
			allocQueue.insert(pair< float,int>((adv->regret - adv->candidateRegret), adv->advertiserId));
			// cout << getCurrentMemoryUsage() << endl;
		}
		
		while(!allocQueue.empty()) {	//stop allocation when no more advertiser is available for allocation			
			// give allocation priority to the advertiser who has the most reduction in regret wrt possible assignment			
			multimap< float, int>::iterator i = allocQueue.end();
			i--;		
			
			adv = advList->at(i->second); 
			timGraph = timList->at(i->second);
			allocQueue.erase(i); //remove from allocation queue, next best will be reinserted soon if allocation do not finish for this
			
			if(attentionQuota[timGraph->candidateNode] > 0) { // if the best candidate still has attention quota left for assignment		
				if(adv->regret > adv->candidateRegret) { //if this assignment is not increasing advertiser's regret
					timGraph->assignBestNode(); // inserts candidateNode into seed set and takes care of tim specific details		
					adv->seedSize = (int)timGraph->seedSet.size();
					adv->totalMonCoverage = adv->totalMonCoverage + timGraph->candidateMMG; 					
					adv->regret = adv->candidateRegret; 
					attentionTemp = attentionQuota[timGraph->candidateNode];  // decrease its attention quota
					attentionQuota[timGraph->candidateNode] = attentionTemp - 1;					
					//cout << "assignment done for advertiser " << adv->advertiserId << " with current regret " << adv->regret << endl;
					// now commenting this, will just take the tot adv results
					writeInAdvertisersFile(adv->advertiserId, timGraph->candidateNode, timGraph->candidateMMG, timGraph->nodeCTRs[timGraph->candidateNode], adv->budget, adv->cpe, lambda);
				
					// if the total coverage is already above budget then stop the allocation for this advertiser	
					if(adv->totalMonCoverage >= adv->budget) { // coz next seed will surely increase the regret
						//cout << "budget consumed, so stopping assignments to advertiser " << adv->advertiserId << endl;
						continue;
					}
					
					else { // calculate the next best node - // first check if an estimation update needed 
						if(timGraph->seedSet.size() == timGraph->kappa) { 						
							timGraph->updateEstimates();		
							//cout << "updating estimates for advertiser " << adv->advertiserId << endl;
						}						
						timGraph->findBestNode(bests);
						adv->candidateRegret = (abs(adv->budget - (adv->totalMonCoverage + bests.second))) + (lambda * (float) (timGraph->seedSet.size() + 1));
						allocQueue.insert(pair< float,int>((adv->regret - adv->candidateRegret), adv->advertiserId));
						continue;
					}					
				}//
				
				else { // this is the place if we want to continue to find the next best node that creates tighter regret
					// regret is increasing so no more assignment to this advertiser, already removed from allocQueue above
					//delete timGraph;
					continue;					
				}				
			}
			
			else { // if the best candidate node selected for that node is already out of the game due to attention quota, select another best
				//for the same iteration 
				timGraph->findBestNode(bests);
				adv->candidateRegret = (abs(adv->budget - (adv->totalMonCoverage + bests.second))) + (lambda * (float) (timGraph->seedSet.size() + 1));
				allocQueue.insert(pair< float,int>((adv->regret - adv->candidateRegret), adv->advertiserId));
			}
		}	
		
		float totalDuration = getRunningTime(startTime); // in seconds 
		float totalMemory = disp_mem_usage(); // in MB
		cout << "MinRegret allocation finished! " << endl;
        cout << "Results stored in " << outFolderName << " folder." << endl;
        
        // write master-level results to master-output
        float total_revenue = 0;
		float total_regret = 0;
        float total_budget = 0;
        float totSeedCost = 0;
        set<int> allocSeeds; // seeds of allocation -- no matter to which adv
   
        for(int i = 0; i < nrCompanies; i++) {
			adv = advList->at(i);
            timGraph = timList->at(i);
            float advSeedCost = 0;
            for(int j = 0; j < timGraph->seedSet.size(); j++) { // compute the money spent on seeds by ctr * cpe
                int seedNode = timGraph->seedSet[j];
                advSeedCost += (timGraph->nodeCTRs[seedNode] * adv->cpe);
                allocSeeds.insert(seedNode);
            }
            total_revenue += adv->totalMonCoverage;
			total_regret += adv->regret;
			total_budget += adv->budget;
            totSeedCost += advSeedCost;
            
            outMasterStream << "TIRM " << userAttentionConstraint << " " << lambda << " " << adv->advertiserId << " " << advSeedCost << " " << adv->seedSize << " " << totalDuration << " " << totalMemory << endl;
//
//            writeInMasterOutputFile(adv->advertiserId, adv->totalMonCoverage, adv->budget, adv->regret, adv->seedSize, totalDuration, totalMemory, lambda, advSeedCost);
            }
        
        int seedSizeTotal = (int) allocSeeds.size(); // total seeds of allocation
        outMasterStream << "TIRM " << userAttentionConstraint << " " << lambda << " " << "T" << " " << totSeedCost << " " << seedSizeTotal << " " << totalDuration << " " << totalMemory << endl;
        
//		writeInMasterOutputFile(-1, total_revenue, total_budget, total_regret, seedSizeTotal, totalDuration, totalMemory, lambda, totSeedCost);
        
    
    }
	
	MinRegretAllocator::~MinRegretAllocator() {
		
		string pid = intToStr(unsigned(getpid()));
		string outfile = "temp/tmp_" + pid + ".txt";
		
		//string command = string("rm -f ") + outfile ;
        string command = string("rm -r ") + "temp";
		system(command.c_str());		
	}
	
	void MinRegretAllocator::readTICGraph() {
		
		string probGraphFile = opt->getValue("probGraphFile");
		cout << "Reading file " << probGraphFile << endl;
		ifstream myfile (probGraphFile.c_str(), ios::in);		
		
		float *dists;		
		float p;
		Advertiser *advTemp;
		TimGraph *timGraphTemp;
		
		int nrEdges = 0;
		set<int> nodes; // kontrol amacli simdi
		
		if (myfile.is_open()) {
			while (! myfile.eof() )	{
				std::string line;
				getline (myfile,line);
				if (line.empty()) continue;
				nrEdges++;
				
				std::string::size_type pos = line.find_first_of(delim);
				int prevpos = 0;
				
				//first user
				string str = line.substr(prevpos, pos-prevpos);
				int u1 = strToInt(str);				
				
				//second user
				prevpos = line.find_first_not_of(delim, pos);
				pos = line.find_first_of(delim, prevpos);
				int u2 = strToInt(line.substr(prevpos, pos-prevpos));				
				
				if (u1 == u2) 
					continue;
				
				graphT[u2].push_back(u1); //insert to the transposed graph
				
				// kontrol amacli
				nodes.insert(u1);
				nodes.insert(u2);				
				
				prevpos = line.find_first_not_of(delim, pos);
				
				str = line.substr(prevpos);
				dists = new  float[nrTopics];
				stringTokenizer(str, dists, nrTopics, delim);
				
				for(int i = 0; i < nrCompanies; i++) {
					advTemp = advList->at(i);
					timGraphTemp = timList->at(i);					
					p = 0.0;
					for(int j = 0; j < nrTopics; j++)
                        p += (dists[j] * advTemp->gamma[j]);
					timGraphTemp->probT[u2].push_back(p);
				}
			}
			
			myfile.close();
		} 
		
		else 
			cout << "Can't open friendship graph file " << probGraphFile << endl;		
		
		cout << "Built transposed adj matrix from file " << endl;
		cout << "Number of nodes: " << nodes.size() << endl;
		cout << "Number of edges: " << nrEdges << endl;
		// 		cout << "size " << graphT.size() << endl;		
	}
	
	void MinRegretAllocator::readNodeCTR() {
		string ctrFile = opt->getValue("nodeCTRfile");
		cout << "Reading ctr file " << ctrFile << endl;
		
		int nodeIndex = 0;
		TimGraph *timGraphTemp;
		
		if(ctrFile.compare("NA") == 0)  { // set all ctr = 1
			cout << "assigning unit ctrs" << endl;
			while(nodeIndex < n) {				
				for(int i = 0; i < nrCompanies; i++) {
					timGraphTemp = timList->at(i);
					timGraphTemp->nodeCTRs[nodeIndex] = 0.01;
				}
				nodeIndex++;				
			}
		}
		
		else {
			
			ifstream myfile(ctrFile.c_str(), ios::in);			
			float *ctrs;							
			
			if(myfile.is_open()) {
				while(!myfile.eof()) {
					std::string line;				
					getline(myfile, line);
					if(line.empty())
						continue;	
					ctrs = new  float[nrCompanies]; //degisecek
					stringTokenizer(line, ctrs, nrCompanies, delim); //degisecek
					
					for(int i = 0; i < nrCompanies; i++) {
						timGraphTemp = timList->at(i);
						timGraphTemp->nodeCTRs[nodeIndex] = ctrs[i];					
					}
					nodeIndex++;
				}			
				
				myfile.close();
			}
			
			else {
				cout << "problem opening ctrs file, exiting... " <<  endl;
				exit(1);
			}
			
			
		}		
	}
	
	
	void MinRegretAllocator::readItemDistsFile() {	
		
		cout << "Reading item distributions file " << endl;
		string itemDistsFile = opt->getValue("itemDistsFile");
		ifstream myfile(itemDistsFile.c_str(), ios::in);
		float *tokens;
		
		int advIndex = 0;
		
		if(myfile.is_open()) {
			while(!myfile.eof()) {
				std::string line;				
				getline(myfile, line);
				if(line.empty())
					continue;	
				tokens = new  float[nrTopics];
				stringTokenizer(line, tokens, nrTopics, delim);		
				// 				for (int i = 0; i < nrTopics; i++)
				// 					cout << tokens[i] << " ";
				// 				cout << endl;
				advList->at(advIndex++)->setItemDist(tokens, nrTopics);
				if(advIndex >= nrCompanies)
					break; 
			}			
			
			myfile.close();
		}
		
		else {
			cout << "problem opening the item distributions file, exiting... " << itemDistsFile <<  endl;
			exit(1);
		}
		
	}

// Use default kappa = 1
//	void MinRegretAllocator::readKappasFile() {				
//		string kappasFile = opt->getValue("kappasFile");
//		ifstream myfile(kappasFile.c_str(), ios::in);		
//		TimGraph *tim;		
//		int counter = 0; 
//		
//		if(kappasFile.compare("NA") == 0) { // use default 1
//			//cout << "using default kappas " << endl;
//			while(counter < nrCompanies) {
//				timList->at(counter++)->kappa = 1;
//			}			
//		}
//		
//		else {			
//			if(myfile.is_open()) {
//				while(!myfile.eof()) {				
//					std::string line;				
//					getline(myfile, line);
//					if(line.empty())
//						continue;					
//					timList->at(counter++)->kappa = strToInt(line);				
//					if(counter >= nrCompanies)
//						break; 
//					
//				}			
//				
//				myfile.close();
//			}
//			
//			else {
//				cout << "problem opening the kappas file, exiting... " << kappasFile <<  endl;
//				exit(1);
//			}	
//			
//		}
//		
//	}
	
	void MinRegretAllocator::readBudgetsFile() {
		
		cout << "Reading advertisers' budgets and cpes " << endl;
		string budgetsFile = opt->getValue("budgetsFile");
		ifstream myfile(budgetsFile.c_str(), ios::in);
		float *bc;
		
		Advertiser *adv;
		int advIndex = 0;
		if(myfile.is_open()) {
			while(!myfile.eof()) {
				std::string line;
				
				getline(myfile, line);
				if(line.empty())
					continue;
				
				bc = new  float[2];				
				stringTokenizer(line, bc, 2, delim);
				
				adv = advList->at(advIndex++);
				adv->budget = bc[0];
				adv->cpe = bc[1];	
				
				if(advIndex >= nrCompanies)
					break; 
			}			
			
			myfile.close();
		}
		
		else {
			cout << "problem opening the budgets file, exiting... " << budgetsFile <<  endl;
			exit(1);
		}		
		
		/*		cout << "budget kontrol " << endl;
		 f o*r(*advIndex = 0; advIndex < nrCompanies; advIndex++) 
		 cout << advList->at(advIndex)->budget << endl;
		 
		 cout << "cpe kontrol " << endl;
		 for(advIndex = 0; advIndex < nrCompanies; advIndex++) 
			 cout << advList->at(advIndex)->cpe << endl;	*/	
	}
	
	void MinRegretAllocator::arrangeOutputFiles() {	
		
		string command = string("mkdir -p ") + outFolderName ;
		system(command.c_str());		
		
		// alti kapadim cunku sadece outputs folderina master output yazicam		
		/* check if the output folder requested already exits */		
		/*		if( access(outFolderName.c_str(), 0) == 0 ) { 
		 s t**at(outFolderName.c_st*r(), &info);
		 if (info.st_mode & S_IFDIR) { //if it is a directory
			 cout << "Output folder specified in the config file already exists, exiting ... " << endl;
			 exit(1);				
	} 					 
	}		
	else {			
		string command = string("mkdir -p ") + outFolderName;
		system(command.c_str());			
	}*/	
		
		/* arrange the output files - folders */	
		string masterFileName = "TIRM_h" +  intToStr(nrCompanies) + "_k" + intToStr(userAttentionConstraint) + "_L" + floatToStr(lambda) + ".txt";
		outMasterName = outFolderName + OS_SEP + masterFileName;
		openOutputFiles();
		
		//
		outFileStreams = new ofstream[nrCompanies];
		outFileNames = new string[nrCompanies];
		/* create output file names for each advertiser */
		for(int i =0; i < nrCompanies; i++) 
			outFileNames[i] = outFolderName + OS_SEP + "TIRM_adv" + intToStr(i) + "_k" + intToStr(userAttentionConstraint) + "_L" + floatToStr(lambda) + ".txt";
		
		for(int i =0; i < nrCompanies; i++) {
			openOutputFiles(advList->at(i)->advertiserId);								
		}
		//		
		// memory icin
		command = string("mkdir -p temp");
		system(command.c_str());
	}
	
	void MinRegretAllocator::writeInAdvertisersFile(int advID, int v, float mmg, float ctr, float budget, float cpe, float lambda) {
		// will write results per advertiser
		(outFileStreams[advID]) << v << " " << mmg << " " << ctr << " " << budget << " " << cpe << " " << lambda << endl;
	}
	
	void MinRegretAllocator::openOutputFiles(int advID) {
		
		// will create output files per advertiser			
		
		if((outFileStreams[advID]).is_open()) 
			(outFileStreams[advID]).close();
		
		(outFileStreams[advID]).open(outFileNames[advID].c_str());
		
		if ((outFileStreams[advID]).is_open() == false) {
			cout << "Can't open file " << outFileNames[advID]  << " for writing" << endl;
			exit(1);
		}
	}	
	
	void MinRegretAllocator::openOutputFiles() {		
		// will create master-output file			
		if(outMasterStream.is_open()) 
			outMasterStream.close();
		
		outMasterStream.open(outMasterName.c_str());
		
		if (outMasterStream.is_open() == false) {
			cout << "Can't open file " << outMasterName  << " for writing" << endl;
			exit(1);
		}
	}
	
    
	void MinRegretAllocator::writeInMasterOutputFile(int advertiserID, float totalCoverage, float budget, float totalRegret, int size, float duration, float memory, float lambda, float seedCost) {
		outMasterStream << advertiserID << " " << totalCoverage << " " << budget << " " << totalRegret << " " << size << " " << duration << " " << memory << " " << lambda << " " << seedCost << endl;
    }
	
	
	
}

#endif