#ifndef _ADV_H_
#define _ADV_H_

namespace _MinRegret {
	
	class Advertiser {
		
	public: 
		
		int advertiserId;
		float *gamma; 
		float cpe;
		float budget; 
		int seedSize;
		
		float regret; 
		float candidateRegret; // assignment yaptiktan sonra ayarlaniyor
		float totalMonCoverage; // total monetary coverage 
		
		Advertiser(int id, int nrTopics) {
			this->advertiserId = id;		
			this->gamma = new  float[nrTopics];
            this->seedSize = 0;
		}
		
		~Advertiser() {}		
		
		void setItemDist( float *temp, int size) {
			for(int i = 0; i < size; i++) 
				this->gamma[i] = temp[i];
		}
		
	};	
	
}


#endif