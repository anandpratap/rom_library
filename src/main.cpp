#include "main.hpp"
int main(){
	auto m = Main();
	m.load_snapshots();
	//m.calc_svd();
	//m.reconstruct();
	//m.calc_deim();
	//m.use_deim();
	m.calc_adeim();
}
