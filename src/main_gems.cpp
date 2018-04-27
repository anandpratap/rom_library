#include "main.hpp"
int main(){
	auto m = GemsRom();
	m.initialize();
	m.m->calc_svd();
	

	int dim = 10;
	int ms = 10;

	m.m->calc_deim(dim);
	for(int i=0; i<2000; i++){
		if(i%100 == 0){
			m.m->use_deim(dim, i);
		}
	}
	for(int i=0; i<2000; i++){
		m.m->calc_adeim(dim, ms, i);
		if(i%100 == 0){
			m.m->use_adeim(dim, i);
		}
	}

}
