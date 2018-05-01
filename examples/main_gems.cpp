#include "main.hpp"

int main(){
	auto m = GemsRom();
	m.load_partition_info();

	int gid = m.get_global_id(0, 100);
	println(gid);
	int lid = m.get_local_id(0, gid);
	println(lid);
	assert(lid == 100);

	int dim = 2000;
	int ms = 10;


	m.initialize();
	//m.m->calc_svd();
	m.m->load_modes();
	arma::vec u = m.m->snapshots.col(0);
	u = m.m->renormalize(u);
	
	arma::vec uhat = m.m->projection_on_basis(u, dim);
	arma::vec ut = m.m->projection_from_basis(uhat, dim);
	assert(u.size() == ut.size());
	println(u(1));
	println(uhat(1));
	println(ut(1));

	arma::uvec idx = {0, 2};
 	arma::uvec vidx = m.get_index(0, idx);
	vidx.print();

	arma::vec uhatp(dim, arma::fill::zeros);
	for(int i=0; i<m.num_processor; i++){
		idx = arma::find(m.partition_id == i);
		vidx = m.get_index(i, idx);
		//vidx.print();
		arma::vec up = u(vidx);
		uhatp += m.get_uhat(i, up, dim);
	}
	

	arma::vec utp = m.get_u(0, uhatp, dim);
	//assert(up.size() == utp.size());
	gid = m.get_global_id(0, 0);
	println(u(gid*8));
	println(uhatp(1));
	println(utp(0));
							 
	
	double *uptr = u.memptr();
	double *uhatptr = uhat.memptr();
	double *utptr = ut.memptr();
	//m.get_u
	
	// m.m->calc_deim(dim);
	// for(int i=0; i<2000; i++){
	// 	if(i%100 == 0){
	// 		m.m->use_deim(dim, i);
	// 	}
	// }
	// for(int i=0; i<2000; i++){
	// 	m.m->calc_adeim(dim, ms, i);
	// 	if(i%100 == 0){
	// 		m.m->use_adeim(dim, i);
	// 	}
	// }

}
