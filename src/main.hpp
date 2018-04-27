#ifndef __MAIN_HPP
#define __MAIN_HPP

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

arma::vec test_function(arma::vec x, double t){
	double pi = arma::datum::pi;
	arma::vec y(x.size());
    y = (sin(t*x) + sin(t*pi/2*x) + sin(t*pi*x))*t;
	return y;
}

class Main{
private:
public:
	arma::mat snapshots_tmp;
	arma::mat snap_mean, snap_std, snap_norm;
	arma::mat snapshots;
	arma::mat modes_spatial;
	arma::mat modes_temporal;
	arma::vec singular_values;

	arma::Col<arma::uword> deim_p;
	
	void load_snapshots(){
		arma::Mat<int> shape;
		shape.load("shape.bin", arma::raw_binary);
		this->snapshots_tmp.load("snapshots.bin", arma::raw_binary);
		this->snapshots.set_size(shape[1]*shape[2], shape[0]);

		for(int i=0; i<shape[0]; i++){
			for(int j=0; j<shape[1]; j++){
				for(int k=0; k<shape[2]; k++){
					snapshots(j*shape[2] + k, i) = snapshots_tmp(i*shape[2]*shape[1] + j*shape[2] + k);
				}
			}
		}
		//snapshots = snapshots.submat(arma::span::all, arma::span(0, 20));

		for(int i=0; i< 20; i++){
			arma::vec snap_t = snapshots.col(i);
			snap_t.save("snap_t_"+std::to_string(i)+".dat", arma::arma_ascii);
		}


		snap_mean = arma::mean(snapshots, 1);
		snap_std = arma::stddev(snapshots, 0, 1);
		print_mat_shape(snap_mean, "Snapshots mean: ");
		print_mat_shape(snap_std, "Snapshots std: ");

		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean(i))/(1.0 + snap_std(i));
			}
		}

		//snapshots = snapshots.cols(0, 100);
		print_mat_shape(this->snapshots, "Snapshots: ");
		println(this->snapshots(0,0));
	}

	void calc_svd(){
		svd_econ(modes_spatial, singular_values, modes_temporal, this->snapshots);
		print_mat_shape(modes_spatial, "modes_spatial: ");
		print_mat_shape(modes_temporal, "modes_temporal: ");
		modes_spatial.save("modes_spatial.bin", arma::arma_binary);
		modes_temporal.save("modes_temporal.bin", arma::arma_binary);
		singular_values.save("singular_values.bin", arma::arma_binary);
		//singular_values.print();
	}

	void reconstruct(int n_mode=11){
		arma::mat sigma = arma::diagmat(singular_values.subvec(0, n_mode-1));
		arma::mat snapshots_recon = modes_spatial.cols(0, n_mode-1)*sigma*modes_temporal.cols(0, n_mode-1).t();
		arma::mat error = snapshots_recon - snapshots;
		print_mat_shape(error, "error: ");
		print("Error in reconstruction: "); println(arma::norm(error, 2)/error.size());
	}

	void calc_deim(){
		arma::Col<arma::uword> p(modes_spatial.n_cols);
		arma::vec u = arma::abs(modes_spatial.col(0));
		auto imax = u.index_max();
		p(0) = imax;
		for(int i=1; i < modes_spatial.n_cols; i++){
			println(i);
			arma::Col<arma::uword> p_sub = p.subvec(0, i-1);
			arma::mat usub = modes_spatial(arma::span::all, arma::span(0, i-1));
			arma::mat A = usub.rows(p_sub);
			//psub.print();
			arma::vec ui = modes_spatial.col(i);
			arma::vec b(i);
			for(arma::uword k=0; k<i; k++){
				b(k) = ui(p_sub(k));
			}
			arma::vec c = arma::solve(A, b);
			arma::vec r = ui - usub * c;
			auto rabs = arma::abs(r);
			arma::uword imax = rabs.index_max();
			p(i) = imax;
		}
		p.print();
		//P.save("P.bin", arma::arma_binary);
		p.save("p_idx.bin", arma::arma_binary);

		//deim_p = p;
		//deim_P = P;
		//arma::vec snap_t = snapshots.col(1500);
		//arma::vec snap_t_recon = modes_spatial * arma::inv(P.t() * modes_spatial) * (P.t() * snap_t);

		//arma::vec snap_t_recon = PP * snap_t;
		//auto error = arma::norm(snap_t - snap_t_recon, 2)/snap_t.size();
		//println("Error DEIM = ");
		//println(error);
	}


	void calc_adeim(){

		int t = 1500;
		
		modes_spatial.load("modes_spatial.bin", arma::arma_binary);
		modes_temporal.load("modes_temporal.bin", arma::arma_binary);
		singular_values.load("singular_values.bin", arma::arma_binary);

		deim_p.load("p_idx.bin", arma::arma_binary);
		int n_modes = 100;
		arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
		arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		int m = 10;
		auto ms = arma::randi<arma::Col<arma::uword>>(m, arma::distr_param(0,usub.n_rows));
		arma::Col<arma::uword> psub_new(psub.size() + m);
		psub_new.subvec(0, psub.size()-1) = psub;
		psub_new.subvec(psub.size(), psub.size()+m-1) = ms;

		psub_new = arma::unique(psub_new);

		arma::vec F(psub_new.size());
		auto ut = snapshots.col(t);
		for(int k=0; k<F.size(); k++){
			F(k) = ut(psub_new(k));
		}
		println(psub_new.size());
		//psub_new.print();

		double ranktol = 1e-14;
		arma::mat modes_new(usub);

		arma::mat A = usub.rows(psub_new);
		print_mat_shape(A, "A");
		arma::vec b = F;
		arma::vec c = arma::solve(A, b);
		arma::vec res = A*c - F;
		arma::mat Q, R;
		arma::qr(Q, R, c);
		arma::vec Rabs = arma::abs(R);
		arma::vec num = arma::mean(Rabs, 1);
		double den = arma::norm(R, "fro");
		arma::vec I = num/(den*den);
		arma::uvec I_idx = arma::find(I > ranktol);
		I_idx.print();
		arma::mat R_tilde = R.rows(I_idx);
		arma::mat RR = R_tilde * R_tilde.t();
		//arma::mat E = arma::eye(arma::size(Q));
		//./print_mat_shape(R, "R");
		//print_mat_shape(Q, "Q");
		//print_mat_shape(R_tilde, "R_tilde");

		//print_mat_shape(E, "E");
		//println(res.size());
		//println(c.size());
		
		arma::mat E_1 = R_tilde*(res).t()*(res)*R_tilde.t();
		arma::cx_vec eigval;
		arma::cx_mat eigvec;
		arma::eig_pair(eigval, eigvec, E_1, RR);
		eigval.print();
		println(E_1.size());
		arma::vec eigval_r = arma::real(eigval);
		arma::mat eigvec_r = arma::real(eigvec);
		auto idx_max = eigval_r.index_max();
		arma::vec col_max = eigvec_r.col(idx_max);
		arma::mat bnorm_square =  col_max.t()*RR*col_max;
		// assert bnorm_square is size 1x1
		print_mat_shape(bnorm_square, "bnorm");

		arma::mat beta = Q.cols(0, R_tilde.n_rows-1) * col_max;
		print_mat_shape(beta, "beta");

		arma::mat alpha = -1.0/bnorm_square(0,0) * res * c.t() * beta;
		print_mat_shape(alpha, "alpha");

		arma::mat delta = alpha*beta.t();
		modes_new.rows(psub_new) = usub.rows(psub_new) + alpha*beta.t();
		delta.print();
		println(arma::max(arma::max(delta)));
		println(arma::min(arma::min(delta)));
	}


	arma::vec renormalize(arma::vec x){
		for(arma::uword k=0; k<x.size(); k++){
			x(k) = x(k)*(1.0 + snap_std(k)) + snap_mean(k);
		}
		return x;
	}
	
	void use_deim(){
		deim_p.load("p_idx.bin", arma::arma_binary);
		int n_modes = 100;
		arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
		arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		
		arma::mat sigma = arma::diagmat(singular_values.subvec(0, n_modes-1));
		arma::mat snapshots_recon = modes_spatial.cols(0, n_modes-1)*sigma*modes_temporal.cols(0, n_modes-1).t();
		singular_values.save("singular_values.data", arma::arma_ascii);
		deim_p.save("p_idx.dat", arma::arma_ascii);
		arma::mat PP = usub * arma::inv(usub.rows(psub));

		usub.save("spatial_modes.dat", arma::arma_ascii);
		for(int i=0; i< 20; i++){
			arma::vec snap_t = snapshots.col(100*i);
			arma::vec snap_t_pod = snapshots_recon.col(100*i);

			arma::vec b(n_modes);
			for(arma::uword k=0; k<n_modes; k++){
				b(k) = snap_t(psub(k));
			}

			arma::vec snap_t_deim = PP * b;
			snap_t = renormalize(snap_t);
			snap_t_pod = renormalize(snap_t_pod);
			snap_t_deim = renormalize(snap_t_deim);
			snap_t.save("snap_t_"+std::to_string(i)+".dat", arma::arma_ascii);
			snap_t_pod.save("snap_t_pod_"+std::to_string(i)+".dat", arma::arma_ascii);
			snap_t_deim.save("snap_t_deim_"+std::to_string(i)+".dat", arma::arma_ascii);
			psub.save("snap_iblank_"+std::to_string(i)+".dat", arma::arma_ascii);
			//arma::vec snap_t_recon = PP * snap_t;
			auto error_pod_deim = arma::norm(snap_t - snap_t_deim, 2)/snap_t.size();
			auto error_pod = arma::norm(snap_t - snap_t_pod, 2)/snap_t.size();
			auto error_deim = arma::norm(snap_t_pod - snap_t_deim, 2)/snap_t.size();
			println(i);
			println("Error POD DEIM = ");
			println(error_pod_deim);
			println("Error POD = ");
			println(error_pod);
			println("Error DEIM = ");
			println(error_deim);
			println("-----------");
		}
	};
};

#endif
