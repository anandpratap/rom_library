#include "../src/main.hpp"

int main(){
	// auto solution = new GemsRom();
	// solution->m = new MainRom();
	// arma::mat snapshots = solution->load_snapshots();
	// solution->m->set_snapshots(snapshots);
	// solution->m->calc_svd();
	// solution->m->save_modes();
	// delete solution->m;
	// delete solution;
	
	auto solution = new GemsRom();
	solution->m = new MainRom();
	auto snapshots = solution->load_snapshots("_deim");
	solution->m->set_snapshots(snapshots);
	solution->m->calc_svd();
	solution->m->save_modes("_deim");
	solution->m->load_modes("_deim");
	solution->m->calc_deim(100);
	for(int i=0; i<2000; i++){
		solution->m->use_deim(100, i);
	}
	delete solution->m;
	delete solution;
}
