#ifndef _TIMGRAPH_C_
#define _TIMGRAPH_C_

#include "TimGraph.h"

namespace _MinRegret {	
	
	TimGraph::TimGraph(Advertiser *adv,  float eps, int nrNodes, int NrEdges, float lambda) {
		
		// 		cout << "........ Initializing TimGraph object ............" << endl;
		// initialize vars
		this->epsilon = eps;
		this->adv = adv;
		this->n = nrNodes;
		this->m = NrEdges;
        this->seedSet.clear();
        this->lambda = lambda;
		
		visit_mark = std::vector<int>(n,0);
		visit = std::vector<bool>(n,false);
		hyper_degree = std::vector<int>(n,0);	
		nodeCTRs = std::vector<float>(n,0);
		
		for(int i = 0; i < n; i++) 
			probT.push_back(std::vector< float>());
		
		for(int i = 0; i < 12; i++)
			sfmt_init_gen_rand(&sfmtSeed, i+1234);	
		
	}
	
	TimGraph::~TimGraph(void) {
		cout << "destructor called for timgraph" << endl;
		// 		delete[] hyperG_adv; 
		// 		delete[] hyperGT_adv;
		// 		delete[] hyperG_temp;
		// 		delete[] hyperGT_temp;
		// 		delete[] hyper_degree;
		// 		delete[] probT;		
		// 		delete[] nodeCTRs;
	}
	
	void TimGraph::doInitialGeneration() {		
		kappa = 1;
        theta = theta_old = 0;
		estimateTheta();		
		cout << "advertiser " << adv->advertiserId << " estimated theta " << theta << endl;
		GenerateRRSets();
	}
	
	void TimGraph::updateEstimates() { // the seed nodes selected "by default" covers the new rr sets in their seed-selection order	
		
		kappa = ceil(adv->regret / candidateMMG);
		// 		cout	<< "new kappa " << kappa << endl;
		theta_old = theta; 
		estimateTheta(); 
		// 		cout << "new theta " << theta << endl;
		
		if(theta <= theta_old) { // no updates, use the same set of RR sets
			theta = theta_old;
			// 			cout << "moving on with old theta : " << theta << endl;		
		}
		else { // // generate new RR sets and updates estimates
			GenerateRRSets();	
			// 			cout << "new RR sets generated for advertiser " << adv->advertiserId << endl;
			
			int oldNR;
			float totNR = 0; //ctr icin
			
			for (int i = 0; i < (int) seedSet.size(); i++) {				
				oldNR = nrRRCovered[i];
				nrRRCovered[i] = oldNR + hyper_degree[seedSet[i]];
				// 				cout << "old and new nrCovered for seed " << seedSet[i] << " is " << oldNR << " - " << nrRRCovered[i] << endl;
// 				totNR += nrRRCovered[i];
				//ctr integrate et buraya
				totNR += (float) nrRRCovered[i] * nodeCTRs[seedSet[i]];
				
				hyper_degree[seedSet[i]] = 0; // since we dont want to select this node as seed again						
							
				// oldCPP
				for(int j = 0; j < hyperG_adv[seedSet[i]].size(); j++){ //for each RR set covered by this node
					int t = hyperG_adv[seedSet[i]][j];
					if(!isCovered[t]){
						isCovered[t]=true; // thos old seed covers it						
						for(int z = 0; z < hyperGT_adv[t].size(); z++){ //decrease the degree of parents who are associated with this RR set
							int item = hyperGT_adv[t][z];
							hyper_degree[item]--;
							if(hyper_degree[item] < 0) { //myopic icin 					
								hyper_degree[item] = 0;	
							}
						}
					}
				}
				
				
	}
	
	// update advertisers tot coverage regret etc
	adv->totalMonCoverage = (( float) n * (( float) totNR / theta)) * adv->cpe;	//totnr a yukarida integrate ettim ctr i		
	adv->regret = abs(adv->budget - adv->totalMonCoverage) + (lambda * (float) seedSet.size()); // candidate regret in priority queue will updated after the new best node selection wrt this regret
}		
}

// returns best node and his "monetary" marginal gain
void TimGraph::findBestNode(infPair &bestCandidate) {					
	float mmg;
	//we can still do this since we set degree zero to the previously assigned seed nodes or to the nodes that are out of game		
	// new cpp
	// 		auto t = max_element(hyper_degree.begin(), hyper_degree.end()); 
	// 		int id = t - hyper_degree.begin();		
	// oldCPP
	
	// kendi max elementini yaz ctr ile carpmak icin
// 	std::vector<int>::iterator t = max_element(hyper_degree.begin(), hyper_degree.end()); 
// 	int id = std::distance(hyper_degree.begin(), t);			
	
// 	mmg = (( float) n * (( float) hyper_degree[id] / theta)) * adv->cpe;
// 	this->candidateNode = id;
// 	this->candidateNrRR = hyper_degree[id];
// 	this->candidateMMG = mmg;
// 	bestCandidate = std::make_pair(id,mmg);
// 	hyper_degree[id]=0;  // assign zero since either it will be allocated soon or will be out of the game due to attention bound
	
	int id = -1;
	float maxVal = 0.0;
	
	for(int i = 0; i < hyper_degree.size(); i++) { // i burada node id ye denk geliyor
		if((maxVal < hyper_degree[i] * nodeCTRs[i] ) &&  (attentionQuota[i] > 0) ) {
			maxVal = hyper_degree[i] * nodeCTRs[i];
			id = i;
		}
	}
	
	mmg = (( float) n * (( float) maxVal / theta)) * adv->cpe; //ctr integrated
	this->candidateNode = id;
	this->candidateNrRR = hyper_degree[id];
	this->candidateMMG = mmg;
	bestCandidate = std::make_pair(id,mmg);
	hyper_degree[id]=0;  // assign zero since either it will be allocated soon or will be out of the game due to attention bound
}

void TimGraph::assignBestNode() {
	
	seedSet.push_back(candidateNode);
	nrRRCovered.push_back(candidateNrRR); // needed for estimate updates		
	// update the degree of other parents
	// new cpp
	// 		for(int t : hyperG_adv[candidateNode]){ //for each RR set covered by this node
	// 			if(!isCovered[t]){ //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
	// 				isCovered[t]=true; // bestnode covers it
	// 				for(int item : hyperGT_adv[t]){ //decrease the degree of parents who are associated with this RR set
	// 					hyper_degree[item]--;
	// 				}
	// 			}
	// 		}	
	
	// oldCPP
	for(int j = 0; j < hyperG_adv[candidateNode].size(); j++){ //for each RR set covered by this node
		int t = hyperG_adv[candidateNode][j];
		if(!isCovered[t]){ //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
			isCovered[t]=true; // bestnode covers it
			for(int z = 0; z < hyperGT_adv[t].size(); z++){ //decrease the degree of parents who are associated with this RR set
				int item = hyperGT_adv[t][z];
				hyper_degree[item]--;
				if(hyper_degree[item] < 0) { //myopic icin 					
					hyper_degree[item] = 0;	
				}
			}
		}
	}	
	
	
	
	}	
	
	void TimGraph::GenerateRRSets() {
		
		while((int)hyperGT_adv.size() < theta) //resize it for new theta
			hyperGT_adv.push_back(std::vector<int>());
		
		for(int i = theta_old; i < theta; i++) {			
			isCovered.push_back(false); // resize isCovered for new theta with false value
			BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_adv);			
		}
		
		hyperG_adv.clear();
		for(int i = 0; i < n; i++)
			hyperG_adv.push_back(std::vector<int>());
		
		// burada -- kontrol et sonra da hyper_degree yi update et
		// new cpp
		// 		for(int i = 0; i < theta; i++){
		// 			if(!isCovered[i]) { // hyperG_adv boylece covered edilmislerin parentlarini tutmayacak
		// 				for(int t:hyperGT_adv[i]) { 			
		// 					hyperG_adv[t].push_back(i);
		// 				}				
		// 			}			
		// 		}
		
		// oldCPP
		for(int i = 0; i < theta; i++){
			if(!isCovered[i]) { // hyperG_adv boylece covered edilmislerin parentlarini tutmayacak
				for(int j = 0; j < hyperGT_adv[i].size(); j++) { 	
					int t = hyperGT_adv[i][j];
					hyperG_adv[t].push_back(i);
				}				
			}			
		}		
		
		
		// update hyper_degree, keep only number of uncovered rr sets!
		for(int i = 0; i < n; i++) {
			hyper_degree[i] = hyperG_adv[i].size();
		}		
		
	}
	
	void TimGraph::BuildHypergraphKPT(int64 R){
		
		// used for kpt estimation
		hyperG_temp.clear();
		for(int i = 0; i < n; i++)
			hyperG_temp.push_back(std::vector<int>());
		
		hyperGT_temp.clear();
		
		while((int)hyperGT_temp.size() <= R) // applies both to isTemp true and false
			hyperGT_temp.push_back(std::vector<int>());	
		
		for(int i = 0; i < R; i++)
			BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_temp); 	
		/*
		 / / n*ew cpp
		 // 		for(int i = 0; i < R; i++){
		 // 			for(int t:hyperGT_temp[i]) {
		 // 				hyperG_temp[t].push_back(i);
		 // 			}
		 // 		}
		 */
		
		// oldCPP
		for(int i = 0; i < R; i++){
			for(int j = 0; j < hyperGT_temp[i].size(); j++) {
				int t = hyperGT_temp[i][j];
				hyperG_temp[t].push_back(i);
			}
		}
		
	}	
	
	void TimGraph::estimateTheta() {
		
		// first estimateOPT via KptEstimation and refineKpt
		float ept, eps_prime; 
		
		float lb=1/2.0;	
		float c=0;
		int64 lastR=0;
		while(true){
			int loop= (6 * log(n)  +  6 * log(log(n)/ log(2)) )* 1/lb  ;			
			c=0;
			lastR=loop;
			// 			IF_TRACE(int64 now=rdtsc());
			float sumMgTu=0;
			for(int i=0; i<loop; i++){
				int u=rand()%n;
				float MgTu=MgT(u);
				float pu=MgTu/m;
				sumMgTu+=MgTu;				
				c+=1-pow((1-pu), kappa);			
			}		
			c/=loop;
			if(c>lb) break;                
			lb /= 2;
		}		
		BuildHypergraphKPT(lastR); // for refineKpt		
		
		ept = (c * n) / 2; //KptEstimation
		// 		cout << "first step ept " << ept << endl;
		// estimateOPT - refineKpt kismi burasi
		// select kappa seeds from the hypergraph builded above : coming from the last iteration of KptEstimation
		BuildSeedSetTemp();  
        eps_prime = 5.0 * std::pow(sqr(epsilon) / (kappa * 1.0), 1.0/3.0); //epsilon' used in RefineKPT
		int64 R = (8 + 2 * eps_prime) * ( n * log(n) +  n * log(2)  ) / ( eps_prime * eps_prime * ept)/4;
		BuildHypergraphKPT(R);
		
		std::set<int> s;
		// new cpp
		// 		for(auto t : seedSetTemp) {
		// 			for(auto tt : hyperG_temp[t]){
		// 				s.insert(tt);
		// 			}
		// 		}
		// oldCPP		
		for(set<int>::iterator t = seedSetTemp.begin(); t != seedSetTemp.end(); t++) {	
			int u = *t;
			for(int j = 0; j < hyperG_temp[u].size(); j++) {
				s.insert(hyperG_temp[u][j]);
			}
		}				
		
		ept = (( float)(n * s.size())) / R;		
		ept = ept / (1 + eps_prime);
		
// 		cout << "final ept " << ept << endl;
		
		// 		 float lambda = (8+2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, kappa) ) / ( epsilon * epsilon); 
		// 		cout << "lambda " << lambda << endl;
		
		theta = (8 + 2 * epsilon) * ( n * log(n) + n * log(2) +  n * logcnk(n, kappa) ) / ( epsilon * epsilon * ept);	
	}
	
	// will be used just to create the kappa seeds needed in RefineKPT 
	void TimGraph::BuildSeedSetTemp() {
		std::vector<int> degree;
		std::vector<int> visit_local(hyperGT_temp.size());		
		seedSetTemp.clear();
		
		for(int i = 0; i < n; i++) 	{
			degree.push_back(hyperG_temp[i].size());	
		}
		
		for(int i = 0; i < kappa; i++){
			// new cpp
			// 			auto t = max_element(degree.begin(), degree.end());
			// 			int id = t - degree.begin();
			// oldCPP
			
			// ctr buraya gelebilir
			std::vector<int>::iterator t = max_element(degree.begin(), degree.end()); 
			int id = std::distance(degree.begin(), t);			
			
			seedSetTemp.insert(id);
			degree[id]=0;			
			// new cpp
			// 			for(int t : hyperG_temp[id]){
			// 				if(!visit_local[t]){
			// 					visit_local[t]=true;
			// 					for(int item : hyperGT_temp[t]){
			// 						degree[item]--;
			// 					}
			// 				}
			// 			}
			
			// oldCPP
			for(int j = 0; j < hyperG_temp[id].size(); j++){				
				int t = hyperG_temp[id][j];
				if(!visit_local[t]){
					visit_local[t]=true;
					for(int z = 0; z < hyperGT_temp[t].size(); z++){
						int item = hyperGT_temp[t][z];
						degree[item]--;
					}
				}
			}
			
			
			
		}
	}
	
	float TimGraph::logcnk(int n, int k){
		float ans = 0;
		
		for(int i = n - k + 1; i <= n; i++){
			ans += log(i);
		}
		for(int i = 1; i <= k; i++){
			ans -= log(i);
		}
		return ans;
	}
	
	int TimGraph::BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT){
		int n_visit_edge=1;
		
		if(addHyperEdge) {
			// 			ASSERT((int)hyperGT.size() > hyperiiid);
			hyperGT[hyperiiid].push_back(uStart);
		}
		
		int n_visit_mark=0;
		
		q.clear();
		q.push_back(uStart);
		// 		ASSERT(n_visit_mark < n);
		
		// randomness kontrolu
		/*		if(addHyperEdge) {
		 i f(*hyperiiid < 10)
		 //cout << "hyperiiid : " << hyperiiid << endl;
		 cout << "generated random node id " << uStart << endl;		
		}*/		
		
		visit_mark[n_visit_mark++] = uStart;
		visit[uStart] = true;
		
		while(!q.empty()) {
			int expand = q.front();
			q.pop_front();
			
			int i = expand;
			for(int j = 0; j < (int)graphT[i].size(); j++){
				//int u=expand;
				int v = graphT[i][j]; //parent of u in the original graph G
				n_visit_edge++;
				float randFloat= float(sfmt_genrand_uint32(&sfmtSeed))/ float(RAND_MAX)/2;
				if(randFloat > probT[i][j])
					continue;
				if(visit[v])
					continue;
				if(!visit[v]) {
					// 					ASSERT(n_visit_mark < n);
					visit_mark[n_visit_mark++]=v;
					visit[v]=true;
				}
				q.push_back(v); 
				
				if(addHyperEdge) {				
					// 					ASSERT((int)hyperGT.size() > hyperiiid);
					hyperGT[hyperiiid].push_back(v);
				}
			}		
		}
		
		for(int i = 0; i < n_visit_mark; i++)
			visit[visit_mark[i]]=false; 		
		
		return n_visit_edge; // returns number of edges considered, i.e., width of the RR set created from uStart
		}
		
		float TimGraph::MgT(int u){            
			// 		ASSERT(u>=0);
			// 		ASSERT(u<n);
			return ( float)BuildHypergraphNode(u, 0, false, hyperGT_temp);
		}
		
		float TimGraph::sqr( float t) {
			return t*t;
		}
		
	}
	
	#endif