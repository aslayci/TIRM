#ifndef _TIMGRAPH_H_
#define _TIMGRAPH_H_

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <vector>
#include <deque>
#include <utility>
#include "sfmt/SFMT.h"
#include "Advertiser.h"
#include "utils.h"

#define IF_TRACE(args) ;

namespace _MinRegret {		
	
	
	class TimGraph {
		
	public: 
		
		std::vector< std::vector< float> > probT; // keep the item specific influence probs in this pbject, in transposed graph style
		std::vector<float> nodeCTRs; 
		
		// maintain RR sets tru execution
		std::vector< std::vector<int> > hyperG_adv;
		std::vector< std::vector<int> > hyperGT_adv;
		
		// keep the temp versions to be used in the kpt estimations
		std::vector< std::vector<int> > hyperG_temp;
		std::vector< std::vector<int> > hyperGT_temp;		
		std::set<int> seedSetTemp; // temporary seed set created for Kpt estimation purposes	
		
		// acts local
		std::vector<bool> visit; 
		std::vector<int> visit_mark; 
		
		std::deque<int> q;
		sfmt_t sfmtSeed;	
		
		// needs to be actively updated
		std::vector<int> seedSet; // keep like this to be able to update in this order
		std::vector<int> nrRRCovered; //keep this with number of covered (not m.g.) to be able to update later
		
		std::vector<int> hyper_degree; // 0 to n, her node icin var - bu RR set yaratilirken doldurulmali 
		std::vector<bool> isCovered; // her RR_id icin - bu da RR set yapilirken de
		
		Advertiser *adv;
		int candidateNode;
		float candidateMMG; // candidate monetary marginal gain
		int candidateNrRR; // candidate nr of covered 
		int64 theta;
		int64 theta_old;
		int kappa;		
		float epsilon;
        float lambda;
		int n, m;			
		
		TimGraph(Advertiser *adv,  float eps, int nrNodes, int nrEdges, float lambda);
		~TimGraph(void);		
		
		void doInitialGeneration();
		void GenerateRRSets();
		void findBestNode(infPair &bestCandidate);		
		void assignBestNode();
		void updateEstimates();
		
		int BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT);
		
		void estimateTheta();
		void BuildHypergraphKPT(int64 R);
		void BuildSeedSetTemp();
		float MgT(int u);
		float logcnk(int n, int k);
		float sqr( float t);
		
		uint64 rdtsc(void)
		{
			unsigned a, d;
			//asm("cpuid");
			asm volatile("rdtsc" : "=a" (a), "=d" (d));
			return (((uint64)a) | (((uint64)d) << 32));
		}
		
	};
	
}


#endif