#ifndef __MAIN_HPP
#define __MAIN_HPP
#include <cassert>
#include <iostream>
#include <armadillo>
template<typename T>
void print(T s){
	std::cout<<s;
}
template<typename T>
void println(T s){
	print(s);print("\n");
}

template<typename T>
void print_mat_shape(arma::Mat<T> mat, std::string header=""){
	println(header);
	print("n_rows = "); print(mat.n_rows); println("");
	print("n_cols = "); print(mat.n_cols); println("");
}



class MainRom{
private:
public:
	arma::mat snap_mean, snap_std, snap_norm;
	arma::vec snap_mean_v, snap_std_v;
	arma::mat snapshots;
	arma::mat modes_spatial;
	arma::mat modes_temporal;
	arma::vec singular_values;

	int isnormalize = 2;
	
	arma::Col<arma::uword> deim_p;

	void set_snapshots(arma::mat isnapshots, std::string suffix="");
	
	void save_modes(std::string suffix="");

	void load_modes(std::string suffix="");

	void calc_svd();

	void reconstruct(int n_mode=11);
	void calc_deim(int dim);

	void calc_adeim(int dim, int m, int tidx);
	arma::vec renormalize(arma::vec x);
	arma::vec normalize(arma::vec x);

	arma::vec renormalize(arma::vec x, arma::uvec var_idx);
	arma::vec normalize(arma::vec x, arma::uvec var_idx);

	void use_adeim(int dim, int t);
	void use_deim(int dim, int t);


	arma::vec projection_on_basis(arma::vec u, int n_mode);
	arma::vec projection_from_basis(arma::vec uhat, int n_mode);
	
	arma::vec projection_on_basis(arma::vec u, int n_mode, arma::uvec var_idx);
	arma::vec projection_from_basis(arma::vec uhat, int n_mode, arma::uvec var_idx);

	void renormalize(int isize, double *x, double *y);

};

class GemsRom{
private:
public:
	MainRom *m;
	int num_processor;

	arma::umat preload_tmp_idx;
	arma::mat preload_PP;

	arma::Col<int> local_id, partition_id,buffer;
	arma::mat load_snapshots(std::string suffix="");
	void initialize(int ipartition_id);
	void load_partition_info();
	int get_global_id(int ipartition_id, int ilocal_id);
	int get_local_id(int ipartition_id, int iglobal_id);
	int get_partition_id(int iglobal_id);
	arma::uvec get_index(int ipartition_id, arma::uvec idx);
	arma::vec calc_deim(int ipartition_id, arma::vec r_s);
	arma::vec get_uhat(int ipartition_id, arma::vec q, int n_mode);
	arma::vec get_u(int ipartition_id, arma::vec qhat, int n_mode);
	void get_uhat(int ipartition_id, double* q, int nq, double *qhat, int n_mode);
	void get_u(int ipartition_id, double* q, int nq, double *qhat, int n_mode);

	void get_uhat(int ipartition_id, int lq, double *q, int lqhat, double *qhat);
	void get_u(int ipartition_id, int lq, double *q, int lqhat, double *qhat);

	int get_deim_local_size(int ipartition_id);

	void get_deim_local_id_idx(int ipartition_id, int *ilocal_id, int *ivar);
	void calc_deim(int ipartition_id, double *r_s, double *deim_r);	
};



#endif
