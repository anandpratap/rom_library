#include "../src/main.hpp"
#include "../c_api/main.h"

     
int main(){
	auto solution = new GemsRom();
	solution->initialize();
	for(int i=0; i<20; i++){
		arma::vec r = arma::zeros(64000);
		for(int proc=0; proc<10; proc++){
			arma::vec r_p = arma::zeros(64000);
			
			int n;
			gemsrom_get_deim_n(solution, proc, n);
			arma::vec r_s = arma::zeros(n);
			r_s.load("r_s_"+std::to_string(i) + "_" +std::to_string(proc), arma::arma_ascii);
			assert(r_s.size() == n);
			if(n>0){
				gemsrom_calc_deim(solution, proc, n, r_s.memptr(), r_p.memptr());
			}
			r = r + r_p;
		}
		gemsrom_renormalize(solution, 64000, r.memptr(), r.memptr());
		r.save("deim_cxx_"+std::to_string(i), arma::arma_ascii);
	}
}
