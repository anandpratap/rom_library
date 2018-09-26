#include "../src/main.hpp"
#define ENABLE_CXXOPTS
#ifdef ENABLE_CXXOPTS
#include "../external/cxxopts/include/cxxopts.hpp"
#endif
int main(int argc, char **argv){
	bool calc_solution_modes = true, calc_residual_modes = false, calc_deim = false;
	bool use_solution_modes = false;
	bool flag_normalization = false;
	int normalization = 4;
	int n_deim_samples = 500;
#ifdef ENABLE_CXXOPTS
	cxxopts::Options options("MyProgram", "One line description of MyProgram");
	options.add_options()("calc_solution_modes", "Calculate solution modes");
	options.add_options()("calc_residual_modes", "Calculate residual modes");
	options.add_options()("use_solution_modes", "Calculate residual modes");
	options.add_options()("calc_deim", "Calculate DEIM sampling points and output values");
	options.add_options()("deim_mode", "Calculate DEIM mode", cxxopts::value<int>()->default_value("0"));
	options.add_options()("normalization", "Calculate DEIM mode", cxxopts::value<int>()->default_value("4"));
	options.add_options()("n_deim_samples", "Calculate DEIM mode", cxxopts::value<int>()->default_value("10"));
	options.add_options()("help", "Print help");

	auto parser = options.parse(argc, argv);
	if(parser["help"].as<bool>()){
		std::cout << options.help({"", "Group"}) << std::endl;
		exit(0);
	}

	calc_solution_modes = parser["calc_solution_modes"].as<bool>();
	calc_residual_modes = parser["calc_residual_modes"].as<bool>();
	use_solution_modes = parser["use_solution_modes"].as<bool>();

	calc_deim = parser["calc_deim"].as<bool>();
	n_deim_samples = parser["n_deim_samples"].as<int>();
	normalization = parser["normalization"].as<int>();
#endif

	if(calc_solution_modes){
		auto solution = new GemsRom();
		//solution->load_partition_info();
		solution->m = new MainRom();
		solution->m->read_config_file("");
		auto snapshots = solution->load_snapshots("");
		solution->m->set_snapshots(snapshots, "", normalization);
		//delete snapshots.memptr();
		solution->m->calc_svd();
		solution->m->save_modes("");
		delete solution->m;
		delete solution;
	}
	if(calc_residual_modes){
		auto solution = new GemsRom();
		int num_processor = solution->load_partition_info();
		solution->m = new MainRom();
		solution->m->read_config_file("");
		auto snapshots = solution->load_snapshots("_residual");
		solution->m->set_snapshots(snapshots, "_residual", 5);
		solution->m->calc_svd();
		solution->m->save_modes("_residual");

		delete solution->m;
		delete solution;

	}		
	if(calc_deim){
		auto solution = new GemsRom();
			solution->load_partition_info();
			solution->m = new MainRom();
			solution->m->read_config_file("");
			solution->m->isnormalize = 4;
			int ndim = n_deim_samples;
			if(use_solution_modes){
				solution->m->load_modes_full("");
			}
			else{
				solution->m->load_modes_full("_residual");
			}
			solution->m->calc_qdeim(ndim);
			
			arma::Col<arma::uword> deim_p;
			deim_p.load("residual_p.ascii", arma::arma_ascii);
			ndim = deim_p.size();
			arma::Mat<arma::uword> tmp_mat(ndim, 6);
			for(int i=0; i<ndim; i++){
				int icv = deim_p(i)/solution->m->nspecies;
				int ivar = deim_p(i)%solution->m->nspecies;
				int pid = solution->get_partition_id(icv);
				int lid = solution->get_local_id(pid, icv);
				
				tmp_mat(i, 0) = deim_p(i);
				tmp_mat(i, 1) = icv;
				tmp_mat(i, 2) = ivar;
				tmp_mat(i, 3) = pid;
				tmp_mat(i, 4) = lid;
				tmp_mat(i, 5) = i;
			}
			
			int n_modes = ndim;
			arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
			arma::Col<arma::uword> upsub = arma::find_unique(psub);
			std::cout<<"unque size: "<<upsub.size()<<std::endl;
			arma::mat usub = solution->m->modes_spatial(arma::span::all, arma::span(0, n_modes-1));

			// for(int i=0; i<usub.n_rows; i++){
			// 	for(int j=0; j<usub.n_cols; j++){
			// 		usub(i,j) = usub(i, j)*solution->m->snap_std_v(i%m->nspecies);
			// 	}
			// }
			arma::mat PP;
			try{
				PP = usub * arma::inv(usub.rows(psub));
			}
			catch (const std::runtime_error& error){
				println("xxxxxx Using pseudo inverse");
				PP = usub * arma::pinv(usub.rows(psub));
			}
			
			//std::cout<<"Condition number:"<<arma::cond(arma::inv(usub.rows(psub)))<<std::endl;
			int num_processor = solution->load_partition_info();
			for(int i=0; i<num_processor; i++){
				arma::uvec pids = tmp_mat.col(3);
				arma::uvec idx = arma::find(pids == i);
				arma::umat tmp_mat_s = tmp_mat.rows(idx);
				//tmp_mat_s.print();
				tmp_mat_s.save("residual_p_"+std::to_string(i), arma::arma_binary);
				arma::mat PP_s = PP.cols(tmp_mat_s.col(5));
				//print_mat_shape(PP_s);
				PP_s.save("PP_p_"+std::to_string(i), arma::arma_binary);
			}
			
	
		delete solution->m;
		delete solution;
	}
	return 0;
}
