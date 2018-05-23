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
	solution->load_partition_info();
	solution->m = new MainRom();
	auto snapshots = solution->load_snapshots("_deim");
	solution->m->set_snapshots(snapshots, "_deim");
	solution->m->calc_svd();
	solution->m->save_modes("_deim");
	solution->m->load_modes("_deim");
	solution->m->calc_deim(100);
	//exit(1);
	arma::Col<arma::uword> deim_p;
	deim_p.load("deim_p.bin", arma::arma_binary);
	int ndim = deim_p.size();
	arma::Mat<arma::uword> tmp_mat(ndim, 6);
	for(int i=0; i<ndim; i++){
		int icv = deim_p(i)/8;
		int ivar = deim_p(i)%8;
		int pid = solution->get_partition_id(icv);
		int lid = solution->get_local_id(pid, icv);
		println("id = " + std::to_string(deim_p(i)) + " icv = "+std::to_string(icv) + " ivar = " + std::to_string(ivar) + " proc = " + std::to_string(pid) + " local id = " + std::to_string(lid));
		tmp_mat(i, 0) = deim_p(i);
		tmp_mat(i, 1) = icv;
		tmp_mat(i, 2) = ivar;
		tmp_mat(i, 3) = pid;
		tmp_mat(i, 4) = lid;
		tmp_mat(i, 5) = i;
	}

	int n_modes = ndim;
	arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
	arma::mat usub = solution->m->modes_spatial(arma::span::all, arma::span(0, n_modes-1));
	arma::mat PP = usub * arma::inv(usub.rows(psub));

	
	for(int i=0; i<10; i++){
		arma::uvec pids = tmp_mat.col(3);
		arma::uvec idx = arma::find(pids == i);
		arma::umat tmp_mat_s = tmp_mat.rows(idx);
		tmp_mat_s.print();
		tmp_mat_s.save("deim_p_"+std::to_string(i), arma::arma_ascii);
		arma::mat PP_s = PP.cols(tmp_mat_s.col(5));
		print_mat_shape(PP_s);
		PP_s.save("PP_p_"+std::to_string(i), arma::arma_ascii);
	}
	

	for(int i=0; i<20; i++){
		arma::vec deim_r = arma::zeros(64000);
		for(int proc=0; proc<10; proc++){
			arma::umat tmp_mat_s;
			tmp_mat_s.load("deim_p_"+std::to_string(proc));
			arma::vec r = solution->m->snapshots.col(i);
			arma::vec r_s = r(tmp_mat_s.col(0));
			deim_r = deim_r + solution->calc_deim(proc, r_s);
		}
		deim_r = solution->m->renormalize(deim_r);
		deim_r.save("deim_parallel_"+std::to_string(i), arma::arma_ascii);
	}

	
	for(int i=0; i<20; i++){
		for(int proc=0; proc<10; proc++){
			arma::umat tmp_mat_s;
			tmp_mat_s.load("deim_p_"+std::to_string(proc));
			arma::vec r = snapshots.col(i);
			arma::vec r_s = r(tmp_mat_s.col(0));
			r_s.save("r_s_"+std::to_string(i) + "_" +std::to_string(proc), arma::arma_ascii);
		}
	}

	//for(int i=0; i<2000; i++){
	//	solution->m->use_deim(100, i);
	//}
	delete solution->m;
	delete solution;
}
