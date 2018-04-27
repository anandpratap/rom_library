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
	//y = 10*x%x*t*t;
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

	void set_snapshots(arma::mat isnapshots){
		snapshots = isnapshots;
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
	
	void load_snapshots(){
		arma::Mat<int> shape;
		arma::mat isnapshots;
		shape.load("shape.bin", arma::raw_binary);
		this->snapshots_tmp.load("snapshots.bin", arma::raw_binary);
		isnapshots.set_size(shape[1]*shape[2], shape[0]);

		for(int i=0; i<shape[0]; i++){
			for(int j=0; j<shape[1]; j++){
				for(int k=0; k<shape[2]; k++){
					isnapshots(j*shape[2] + k, i) = snapshots_tmp(i*shape[2]*shape[1] + j*shape[2] + k);
				}
			}
		}
		//snapshots = snapshots.submat(arma::span::all, arma::span(0, 20));
		set_snapshots(isnapshots);
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

	void calc_deim(int dim){
		arma::Col<arma::uword> p(dim);
		arma::vec u = arma::abs(modes_spatial.col(0));
		//u.print();
		auto imax = u.index_max();
		println("imax");println(imax);
		p(0) = imax;
		for(int i=1; i < dim; i++){
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


	void calc_adeim(int dim, int m, int tidx){



		
		int t = tidx;
		
		modes_spatial.load("modes_spatial.bin", arma::arma_binary);
		modes_temporal.load("modes_temporal.bin", arma::arma_binary);
		singular_values.load("singular_values.bin", arma::arma_binary);

		deim_p.load("p_idx.bin", arma::arma_binary);
		int n_modes = dim;
		arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
		arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		
		auto ms = arma::randi<arma::Col<arma::uword>>(m, arma::distr_param(0,usub.n_rows-1));
		arma::Col<arma::uword> psub_new(psub.size() + m);
		psub_new.subvec(0, psub.size()-1) = psub;
		psub_new.subvec(psub.size(), psub.size()+m-1) = ms;

		psub_new = arma::unique(psub_new);

		arma::vec F(psub_new.size());
		auto ut = snapshots.col(t);
		for(int k=0; k<F.size(); k++){
			F(k) = ut(psub_new(k));
		}
		//psub_new.print();
		//F.print();

// rankTol = 1e-14;

// % compute DEIM coefficients matrix and residual at sampling points
// size(Uold(S, :))
// size(F)
// C = Uold(S, :)\F; % DEIM coefficients
// size(C)
// %C
// res = Uold(S, :)*C - F; % residual at sampling points
// % reveal rank of C matrix following Lemma 3.4
// [Q, R, E] = qr(C);
// %Q, R
// I = mean(abs(R), 2)/norm(R, 'fro')^2 > rankTol;
// RTilde = R(I, :);
// RR = RTilde*RTilde';


		
		//println(psub_new.size());
		//psub_new.print();

		double ranktol = 1e-14;
		arma::mat modes_new(usub);

		arma::mat A = usub.rows(psub_new);
		//print_mat_shape(A, "A");
		arma::mat b = F;
		//print_mat_shape(b, "b");
		arma::mat c = arma::solve(A, b);
		arma::mat res = A*c - b;
		//res.print();
		arma::mat Q, R;
		arma::qr(Q, R, c);
		arma::mat Rabs = arma::abs(R);
		arma::mat num = arma::mean(Rabs, 1);
		double den = arma::norm(R, "fro");
		arma::mat I = num/(den*den);
		arma::uvec I_idx = arma::find(I > ranktol);
		
		//I.print();
		//I_idx.print();
		arma::mat R_tilde = R.rows(I_idx);
		arma::mat RR = R_tilde * R_tilde.t();
		//arma::mat E = arma::eye(arma::size(Q));
		//./print_mat_shape(R, "R");
		//print_mat_shape(Q, "Q");
		//print_mat_shape(R_tilde, "R_tilde");

		//print_mat_shape(E, "E");
		//println(res.size());
		//println(c.size());

// % now use RTilde instead of C to solve eigenvalue problem (Lemma 3.5)
// % compute update vectors alpha and beta
// [rMat, eMat] = eig(RTilde*(res*E)'*(res*E)*RTilde', RR);
// [~, maxI] = max(diag(eMat));
// normBSquare = rMat(:, maxI)'*RR*rMat(:, maxI);
// 'Q'
// size( Q(:, 1:size(RTilde, 1)))
// 'rmat'
// size(rMat(:, maxI))
// beta = Q(:, 1:size(RTilde, 1))*rMat(:, maxI);
// alpha = -1/normBSquare*res*C'*beta;

// % apply update to basis
// Unew = Uold;
// Unew(S, :) = Unew(S, :) + alpha*beta';

		
		arma::mat E_1 = R_tilde*res.t()*res*R_tilde.t();
		arma::cx_vec eigval;
		arma::cx_mat eigvec;
		arma::eig_pair(eigval, eigvec, E_1, RR);
		//eigval.print();
		//println(E_1.size());
		arma::vec eigval_r = arma::real(eigval);
		arma::mat eigvec_r = arma::real(eigvec);
		arma::uword idx_max = eigval_r.index_max();
		arma::cx_vec col_max = eigvec.col(idx_max);
		arma::mat bnorm_square =  arma::real(col_max.t()*RR*col_max);
		// assert bnorm_square is size 1x1
		//print_mat_shape(bnorm_square, "bnorm");

		arma::mat beta = Q.cols(0, R_tilde.n_rows-1) * arma::real(col_max);
		//print_mat_shape(beta, "beta");

		arma::mat alpha = -1.0/bnorm_square(0,0) * res * c.t() * beta;
		//print_mat_shape(alpha, "alpha");

		arma::mat delta = alpha*beta.t();
		
		modes_new.rows(psub_new) = usub.rows(psub_new) + alpha*beta.t();
		
		modes_spatial(arma::span::all, arma::span(0, n_modes-1)) = modes_new;

		modes_spatial.save("modes_spatial.bin", arma::arma_binary);

		//delta.print();
		//println(arma::max(arma::max(delta)));
		//println(arma::min(arma::min(delta)));

		arma::mat Unew = modes_new;
		arma::mat Uold = usub;
		arma::mat Uno = Unew.t()*Uold;
		arma::vec Uno_diag = Uno.diag();

		arma::mat Un = Unew.t()*Unew;
		arma::vec Un_diag = Un.diag();
		
		arma::mat Uo = Uold.t()*Uold;
		arma::vec Uo_diag = Uo.diag();

// [~, minI] = min(abs(diag(Unew'*Uold))./(sqrt(diag(Unew'*Unew)).*sqrt(diag(Uold'*Uold))));

// % remove maxI-th vector from Unew
// Utmp = Unew;
// Utmp(:, minI) = [];
// % remove maxI-th point from P
// Ptmp = P;
// Ptmp(minI) = [];

// % find new point following standard DEIM
// c = Utmp(Ptmp, :)\Unew(Ptmp, minI);
// r = abs(Unew(:, minI) - Utmp*c); % residual
// [~, I] = sort(r, 'descend');
// curI = 1;
// while(~isempty(find(Ptmp == I(curI), 1))) % check if point is already interp
//     curI = curI + 1;
// end
// P(minI) = I(curI); % add first point that is not already an interp point

		
		int idx_min = arma::index_min(arma::abs(Uno_diag)/(arma::sqrt(Un_diag)%arma::sqrt(Uo_diag)));

		arma::mat Utmp(Unew);
		Utmp.shed_col(idx_min);
		//print_mat_shape(Utmp, "Utmp");
		arma::uvec ptmp(psub);
		ptmp.shed_row(idx_min);
		
		arma::mat AA = Utmp.rows(ptmp);
		arma::mat BBt = Unew.rows(ptmp); 
		arma::mat BB = BBt.col(idx_min);
		arma::mat cc = arma::solve(AA, BB);
		//println("c");println(c.size());
		arma::mat r = Unew.col(idx_min) - Utmp*cc;
		r = arma::abs(r);
		arma::uvec ridx = sort_index(r, "descend");

		int curl = 0;
		
		while (arma::any(ptmp == ridx(curl))){
			curl += 1;
		}
		//psub(idx_min) = ridx(curl);
		//psub.print();
		deim_p.save("p_idx.bin", arma::arma_binary);
		

		//modes_spatial
// % find new point following standard DEIM
// c = Utmp(Ptmp, :)\Unew(Ptmp, minI);
// r = abs(Unew(:, minI) - Utmp*c); % residual
// [~, I] = sort(r, 'descend');
// curI = 1;
// while(~isempty(find(Ptmp == I(curI), 1))) % check if point is already interp
//     curI = curI + 1;
// end
// P(minI) = I(curI); % add first point that is not already an interp point
		

	}


	arma::vec renormalize(arma::vec x){
		for(arma::uword k=0; k<x.size(); k++){
			x(k) = x(k)*(1.0 + snap_std(k)) + snap_mean(k);
		}
		return x;
	}

	void use_adeim(int dim, int t){
		deim_p.load("p_idx.bin", arma::arma_binary);
		int n_modes = dim;
		arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
		arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		arma::mat PP = usub * arma::inv(usub.rows(psub));
		arma::vec snap_t = snapshots.col(t);
		arma::vec b(n_modes);
		for(arma::uword k=0; k<n_modes; k++){
			b(k) = snap_t(psub(k));
		}
		arma::vec snap_t_adeim = PP * b;
		snap_t_adeim = renormalize(snap_t_adeim);
		snap_t_adeim.save("snap_t_adeim_"+std::to_string(t)+".dat", arma::arma_ascii);
		
	}
	
	void use_deim(int dim, int t){
		deim_p.load("p_idx.bin", arma::arma_binary);
		int n_modes = dim;
		arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
		arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		
		arma::mat sigma = arma::diagmat(singular_values.subvec(0, n_modes-1));
		arma::mat snapshots_recon = modes_spatial.cols(0, n_modes-1)*sigma*modes_temporal.cols(0, n_modes-1).t();
		singular_values.save("singular_values.data", arma::arma_ascii);
		deim_p.save("p_idx.dat", arma::arma_ascii);
		arma::mat PP = usub * arma::inv(usub.rows(psub));

		usub.save("spatial_modes.dat", arma::arma_ascii);
		//		for(int i=1; i < 11; i++){
			arma::vec snap_t = snapshots.col(t);
			arma::vec snap_t_pod = snapshots_recon.col(t);

			arma::vec b(n_modes);
			for(arma::uword k=0; k<n_modes; k++){
				b(k) = snap_t(psub(k));
			}

			arma::vec snap_t_deim = PP * b;
			snap_t = renormalize(snap_t);
			snap_t_pod = renormalize(snap_t_pod);
			snap_t_deim = renormalize(snap_t_deim);
			snap_t.save("snap_t_"+std::to_string(t)+".dat", arma::arma_ascii);
			snap_t_pod.save("snap_t_pod_"+std::to_string(t)+".dat", arma::arma_ascii);
			snap_t_deim.save("snap_t_deim_"+std::to_string(t)+".dat", arma::arma_ascii);
			psub.save("snap_iblank_"+std::to_string(t)+".dat", arma::arma_ascii);
			//arma::vec snap_t_recon = PP * snap_t;
			auto error_pod_deim = arma::norm(snap_t - snap_t_deim, 2)/snap_t.size();
			auto error_pod = arma::norm(snap_t - snap_t_pod, 2)/snap_t.size();
			auto error_deim = arma::norm(snap_t_pod - snap_t_deim, 2)/snap_t.size();
			println(t);
			println("Error POD DEIM = ");
			println(error_pod_deim);
			println("Error POD = ");
			println(error_pod);
			println("Error DEIM = ");
			println(error_deim);
			println("-----------");
		
	};
};

#endif
