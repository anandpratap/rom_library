#include <iostream>
#include <fstream>
#include <armadillo>
#include <cassert>
#include "main.hpp"
#include "mkl_lapacke.h"
#define USE_LAPACK
void MainRom::read_config_file(std::string directory){
	std::fstream config_file(directory + "rom.config", std::ios_base::in);
	config_file >> ndimension;
	config_file >> ncv;
	config_file >> nspecies;
	config_file.close();
	ndof = ncv*nspecies;
	std::cout<<ndimension<<std::endl;
	std::cout<<ncv<<std::endl;
	std::cout<<nspecies<<std::endl;
	
}
arma::mat load_arma_binary_partial(std::string filename, int n_cols){
	std::ifstream f;
	f.open(filename.c_str(), std::fstream::binary);
	std::string f_header;
	arma::uword f_n_rows;
	arma::uword f_n_cols;
	arma::mat x;
	f >> f_header;
	f >> f_n_rows;
	f >> f_n_cols;
	assert(n_cols <= f_n_cols);
	arma::uword skipcols = f_n_cols - n_cols;
	
	f.get();
	x.set_size(f_n_rows, n_cols);
	arma::uword n_elem = f_n_rows*n_cols;
	f.read(reinterpret_cast<char *>(x.memptr()), std::streamsize(n_elem*sizeof(double)));
	return x;
}

using MatrixXd = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;

MatrixXd arma_matrix_to_eigen(arma::mat m){
	Eigen::MatrixXd me = Eigen::Map<Eigen::MatrixXd>(m.memptr(), m.n_rows, m.n_cols);
	return me;
}

arma::mat eigen_matrix_to_arma(MatrixXd me){
	arma::mat m = arma::mat(me.data(), me.rows(), me.cols(), false, false);
	return m;
}

void MainRom::set_snapshots(int idim, int insamples, double *isnapshots, int normalization){
	assert(isnapshots != nullptr);
	arma::mat tmp_snapshots = arma::mat(isnapshots, idim, insamples, false, true);
	set_snapshots(tmp_snapshots, "", normalization);
}

arma::Col<int> calc_qr_lapack(arma::mat A){
	int m = A.n_rows;
	int n = A.n_cols;
	int lda = m;
	arma::Col<int> jpvt = arma::zeros<arma::Col<int>>(n);
	arma::Col<double> tau = arma::zeros<arma::Col<double>>(m);
	LAPACKE_dgeqpf(LAPACK_COL_MAJOR, m, n, A.memptr(), lda, jpvt.memptr(), tau.memptr());
	// this outputs indexing starting from 1
	jpvt = jpvt - 1;
	//jpvt.print();
	return jpvt;
}


void MainRom::set_snapshots(arma::mat isnapshots, std::string suffix, int normalization){
	isnormalize = normalization;
	snapshots = isnapshots;
	if(isnormalize==1){
		snap_mean = arma::mean(snapshots, 1);
		snap_std = arma::stddev(snapshots, 0, 1);
		print_mat_shape(snap_mean, "Snapshots mean: ");
		print_mat_shape(snap_std, "Snapshots std: ");
		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean(i))/(1.0 + snap_std(i));
			}
		}
		snap_mean.save("snap_mean.bin"+suffix, arma::arma_binary);
		snap_std.save("snap_std.bin"+suffix, arma::arma_binary);
	}
	else if(isnormalize == 2){
		snap_mean_v = arma::vec(nspecies);
		snap_std_v = arma::vec(nspecies);
		int ncv = snapshots.n_rows/nspecies;
		int nt = snapshots.n_cols;
		arma::vec Var(ncv*nt);
		for(int i=0; i<nspecies; i++){
			int k=0;
			for(int icv=0; icv<ncv; icv++){
				for(int t=0; t<nt; t++){
					Var(k) = snapshots(icv*nspecies + i, t);
					k+=1;
				}
			}
			assert(k == ncv*nt);
			snap_mean_v(i) = arma::mean(Var);
			snap_std_v(i) = arma::stddev(Var, 0);
		}
		snap_mean_v(2) = snap_mean_v(1);
		snap_std_v(2) = snap_std_v(1);

		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean_v(i%nspecies))/snap_std_v(i%nspecies);
			}
		}
		
		snap_mean_v.save("snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
	
	}
	
	else if(isnormalize == 3){
		snap_mean_v = arma::vec(nspecies);
		snap_std_v = arma::vec(nspecies);
		int ncv = snapshots.n_rows/nspecies;
		int nt = snapshots.n_cols;
		arma::vec Var(ncv*nt);
		for(int i=0; i<nspecies; i++){
			int k=0;
			for(int icv=0; icv<ncv; icv++){
				for(int t=0; t<nt; t++){
					Var(k) = snapshots(icv*nspecies + i, t);
					k+=1;
				}
			}
			assert(k == ncv*nt);
			snap_mean_v(i) = arma::mean(Var);
			snap_std_v(i) = arma::abs(Var).max();

			if(i>3){
				snap_std_v(i) = 1.0;
			}

		}
		//snap_mean_v(2) = snap_mean_v(1);
		snap_std_v(2) = snap_std_v(1);
		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean_v(i%nspecies))/snap_std_v(i%nspecies);
			}
		}
		
		snap_mean_v.save("snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
		
	}
	else if(isnormalize == 4 || isnormalize == 5){
		snap_mean = arma::mean(snapshots, 1);
		snap_std_v = arma::vec(nspecies);
		int ncv = snapshots.n_rows/nspecies;
		int nt = snapshots.n_cols;
		arma::vec Var(ncv*nt);
		for(int i=0; i<nspecies; i++){
			int k=0;
			for(int icv=0; icv<ncv; icv++){
				for(int t=0; t<nt; t++){
					Var(k) = snapshots(icv*nspecies + i, t) - snap_mean(icv*nspecies + i);
					k+=1;
				}
			}
			assert(k == ncv*nt);
			snap_std_v(i) = arma::abs(Var).max();
			if(i>= 4 + (ndimension-2)){
				snap_std_v(i) = 1.0;
			}
			
		}
		//snap_mean_v(2) = snap_mean_v(1);
		double umag = 0;
		for(int i=0; i<ndimension; i++){
			umag += snap_std_v(1+i)*snap_std_v(1+i);
		}
		umag = std::sqrt(umag);

		for(int i=0; i<ndimension; i++){
			snap_std_v(1+i) = umag;
		}
		
		if(isnormalize == 5){
			snap_std_v.load("snap_std_v.bin", arma::arma_binary);
		}
		
		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean(i))/snap_std_v(i%nspecies);
			}
		}
		snap_mean.save("snap_mean.bin"+suffix, arma::arma_binary);
		if(isnormalize == 5){
			snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
			snap_std_v.save("snap_std_v.ascii"+suffix, arma::arma_ascii);
		}
		else{
			snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
			snap_std_v.save("snap_std_v.ascii"+suffix, arma::arma_ascii);
		}
	}
	else{
		throw std::invalid_argument("isnormalize not valid.");
	}


	//snapshots = snapshots.cols(0, 100);
	//print_mat_shape(this->snapshots, "Snapshots: ");
	//println(this->snapshots(0,0));
}

arma::vec MainRom::projection_on_basis(arma::vec u, int n_mode){
	assert(u.size() == modes_spatial.n_rows);
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_mode-1));
	return usub.t()*normalize(u);
}

arma::vec MainRom::projection_on_basis(arma::vec u, int n_mode, arma::uvec var_idx){
	assert(u.size() == var_idx.size());
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_mode-1));
	arma::mat uusub = usub.rows(var_idx);
	u = normalize(u, var_idx);
	return uusub.t()*u;
}


arma::vec MainRom::projection_from_basis(arma::vec uhat, int n_mode){
	assert(uhat.size() == n_mode);
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_mode-1));
	arma::vec u= usub*uhat;
	return renormalize(u);
}


arma::vec MainRom::projection_from_basis(arma::vec uhat, int n_mode, arma::uvec var_idx){
	assert(uhat.size() == n_mode);
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_mode-1));
	arma::mat uusub = usub.rows(var_idx);
	arma::vec u = uusub*uhat;
	u = renormalize(u, var_idx);
	return u;
}


void MainRom::save_modes(std::string suffix){
	modes_spatial.save("modes_spatial.bin"+suffix, arma::arma_binary);
	modes_temporal.save("modes_temporal.bin"+suffix, arma::arma_binary);
	singular_values.save("singular_values.bin"+suffix, arma::arma_binary);
	singular_values.save("singular_values.ascii"+suffix, arma::arma_ascii);
	
}

void MainRom::load_modes(std::string suffix, std::string directory){
	//isnormalize = 3;
	
	if(isnormalize == 1){
		snap_mean.load(directory + "snap_mean.bin"+suffix, arma::arma_binary);
		snap_std.load(directory + "snap_std.bin"+suffix, arma::arma_binary);
		
	}
	else if(isnormalize == 2 || isnormalize == 3){
		snap_mean_v.load(directory + "snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.load(directory + "snap_std_v.bin"+suffix, arma::arma_binary);
	}
	else if(isnormalize == 4 || isnormalize == 5){
		snap_mean.load(directory+"snap_mean.bin"+suffix, arma::arma_binary);
		//std::cout<<"Loading "<<"snap_std_v.bin"+suffix<<std::endl;
		snap_std_v.load(directory+"snap_std_v.bin"+suffix, arma::arma_binary);
		//std::cout<<"tempfix Loading "<<"snap_std_v.bin"<<std::endl;

		//* tempfix
		//* for now using solution std for normalization
		//snap_std_v.load(directory+"snap_std_v.bin", arma::arma_binary);
		
	}
	
	// if(n_cols > 0){
	// 	modes_spatial = load_arma_binary_partial("modes_spatial.bin"+suffix, n_cols);
	// }
	// else{
	// 	modes_spatial.load("modes_spatial.bin"+suffix, arma::arma_binary);
	// }
	//	modes_spatial.save("modes_spatial_raw.bin"+suffix, arma::raw_binary);
	
	//modes_temporal.load("modes_temporal.bin"+suffix, arma::arma_binary);
	//singular_values.load("singular_values.bin"+suffix, arma::arma_binary);
}

void MainRom::load_modes_full(std::string suffix, std::string directory){
	//isnormalize = 3;
	
	load_modes(suffix, directory);
	// if(n_cols > 0){
	// 	modes_spatial = load_arma_binary_partial("modes_spatial.bin"+suffix, n_cols);
	// }
	// else{
	modes_spatial.load(directory + "modes_spatial.bin"+suffix, arma::arma_binary);
	// }
	//	modes_spatial.save("modes_spatial_raw.bin"+suffix, arma::raw_binary);
	
	//modes_temporal.load("modes_temporal.bin"+suffix, arma::arma_binary);
	//singular_values.load("singular_values.bin"+suffix, arma::arma_binary);
}


void MainRom::calc_svd(){
	svd_econ(modes_spatial, singular_values, modes_temporal, this->snapshots);
	print_mat_shape(modes_spatial, "modes_spatial: ");
	print_mat_shape(modes_temporal, "modes_temporal: ");
	
	//singular_values.print();
}

void MainRom::reconstruct(int n_mode){
	arma::mat sigma = arma::diagmat(singular_values.subvec(0, n_mode-1));
	arma::mat snapshots_recon = modes_spatial.cols(0, n_mode-1)*sigma*modes_temporal.cols(0, n_mode-1).t();
	arma::mat error = snapshots_recon - snapshots;
	print_mat_shape(error, "error: ");
	print("Error in reconstruction: "); println(arma::norm(error, 2)/error.size());
}

arma::mat MainRom::get_var_modes(int ivar_idx){
	arma::uvec var_indices = arma::regspace<arma::uvec>(ivar_idx, nspecies, ndof-1);
	var_indices.print();
	return modes_spatial.rows(var_indices);
}

void MainRom::calc_qdeim(int dim){
	arma::mat sub_basis;
	if(deim_mode == DEIM_MODE_VECTOR){
		sub_basis = modes_spatial(arma::span::all, arma::span(0, dim-1));
	}
	else if(deim_mode == DEIM_MODE_CELL){
		arma::mat tmp_modes = get_var_modes(3);
		sub_basis = tmp_modes(arma::span::all, arma::span(0, dim-1));
	}

	arma::mat sub_basis_t = arma::trans(sub_basis);

	std::cout<<" Calculating QR"<<std::endl;

#if defined(USE_LAPACK)
	arma::Col<int> p_ind = calc_qr_lapack(sub_basis_t);
#else
	MatrixXd eigen_sub_basis = arma_matrix_to_eigen(sub_basis);
	// do qr
	Eigen::ColPivHouseholderQR<MatrixXd> qr(eigen_sub_basis.transpose());
	//qr.compute(eigen_sub_basis.transpose());
	std::cout<<" DONE QR"<<std::endl;
	
	int m = sub_basis.n_cols;
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType P(qr.colsPermutation());
	auto p_ind = P.indices();
	std::cout<<P.rows()<<" "<<P.cols();
	std::cout<<p_ind.rows()<<" "<<p_ind.cols();
#endif
	
	arma::uvec pp(dim);
	for(int i=0; i<dim; i++){
		pp(i) = p_ind(i);
		std::cout<<p_ind(i)<<std::endl;
	}
	//std::cout<<P;

	if(deim_mode == DEIM_MODE_VECTOR){
		assert(pp.min() >= 0);
		assert(pp.max() < ndof);
		pp.save("residual_p.bin", arma::arma_binary);
		pp.save("residual_p.ascii", arma::arma_ascii);

	}
	else if(deim_mode == DEIM_MODE_CELL){
		assert(pp.min() >= 0);
		assert(pp.max() < 8000);
		arma::uvec ppp = arma::uvec(pp.size()*nspecies);
		for(int i=0; i<pp.size(); i++){
			for(int j=0; j<nspecies; j++){
				ppp(i*nspecies + j) = pp(i)*nspecies + j;
			}
		}
		ppp.save("residual_p.bin", arma::arma_binary);
		ppp.save("residual_p.ascii", arma::arma_ascii);
	}
}

	

void MainRom::calc_deim(int dim){
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
	arma::sort(p).print();
	deim_p = p;
	
	//P.save("P.bin", arma::arma_binary);
	//p.save("p_idx.bin", arma::arma_binary);
	deim_p.save("residual_p.bin", arma::arma_binary);
	deim_p.save("residual_p.ascii", arma::arma_ascii);
	//deim_p = p;
	//deim_P = P;
	//arma::vec snap_t = snapshots.col(1500);
	//arma::vec snap_t_recon = modes_spatial * arma::inv(P.t() * modes_spatial) * (P.t() * snap_t);

	//arma::vec snap_t_recon = PP * snap_t;
	//auto error = arma::norm(snap_t - snap_t_recon, 2)/snap_t.size();
	//println("Error DEIM = ");
	//println(error);
}


void MainRom::calc_adeim(int dim, int m, int tidx){



		
	int t = tidx;
		


	int n_modes = dim;
	arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));

	arma::Col<arma::uword> ms(m);
	if(0){
		ms = arma::randi<arma::Col<arma::uword>>(m, arma::distr_param(0,usub.n_rows-1));
	}
	else{
		ms.load("../../matlab/adeim/rand_"+std::to_string(t+1) + ".mat", arma::raw_ascii);
		ms = ms - 1;
	}
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

	//modes_spatial.save("modes_spatial.bin", arma::arma_binary);

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
		
	psub(idx_min) = ridx(curl);
	deim_p = psub;
	//psub.print();
	//deim_p.save("p_idx.bin", arma::arma_binary);
		

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


arma::vec MainRom::renormalize(arma::vec x){
	arma::vec y(x);
	if(isnormalize==1){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*(1.0 + snap_std(k)) + snap_mean(k);
		}
	}
	else if(isnormalize==2 || isnormalize==3){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*snap_std_v(k%nspecies) + snap_mean_v(k%nspecies);
		}
	}
	else if(isnormalize==4 || isnormalize == 5){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*snap_std_v(k%nspecies) + snap_mean(k);
		}
	}

	else{
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k);
		}
	}
	return y;
}

arma::vec MainRom::renormalize(arma::vec x, arma::uvec var_idx){
	arma::vec y(x);
	if(isnormalize==1){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*(1.0 + snap_std(var_idx(k))) + snap_mean(var_idx(k));
		}
	}
	else if(isnormalize==2 || isnormalize ==3){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*snap_std_v(var_idx(k)%nspecies) + snap_mean_v(var_idx(k)%nspecies);
		}
	}
	else if(isnormalize==4 || isnormalize == 5){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k)*snap_std_v(var_idx(k)%nspecies) + snap_mean(var_idx(k));
		}
	}

	else{
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k);
		}
	}
	return y;
}



void MainRom::renormalize(int isize, double *x, double *y){
	if(isnormalize==1){
		for(arma::uword k=0; k<isize; k++){
			y[k] = x[k]*(1.0 + snap_std(k)) + snap_mean(k);
		}
	}
	else if(isnormalize==2 || isnormalize==3){
		for(arma::uword k=0; k<isize; k++){
			y[k] = x[k]*snap_std_v(k%nspecies) + snap_mean_v(k%nspecies);
		}
	}
	else if(isnormalize==4  || isnormalize == 5){
		//std::cout<<"calling renormalize"<<std::endl;
		for(arma::uword k=0; k<isize; k++){
			y[k] = x[k]*snap_std_v(k%nspecies);// +  snap_mean(k);
		}
	}

	else{
		for(arma::uword k=0; k<isize; k++){
			y[k] = x[k];
		}
	}
}


arma::vec MainRom::normalize(arma::vec x){
	arma::vec y(x);
	if(isnormalize==1){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean(k))/(1.0 + snap_std(k));
		}
	}
	else if(isnormalize==2 || isnormalize==3){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean_v(k%nspecies))/snap_std_v(k%nspecies);
		}
	}
	else if(isnormalize==4  || isnormalize == 5){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean(k))/snap_std_v(k%nspecies);
		}
	}

	else{
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k);
		}
	}
	return y;
}


arma::vec MainRom::normalize(arma::vec x, arma::uvec var_idx){
	arma::vec y(x);
	if(isnormalize==1){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean(var_idx(k)))/(1.0 + snap_std(var_idx(k)));
		}
	}
	else if(isnormalize==2 || isnormalize==3){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean_v(var_idx(k)%nspecies))/snap_std_v(var_idx(k)%nspecies);
		}
	}
	else if(isnormalize==4  || isnormalize == 5){
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = (x(k) - snap_mean(var_idx(k)))/snap_std_v(var_idx(k)%nspecies);
		}
	}

	else{
		for(arma::uword k=0; k<x.size(); k++){
			y(k) = x(k);
		}
	}
	return y;
}


void MainRom::use_adeim(int dim, int t){
		
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
	psub.save("iblank_adeim_"+std::to_string(t)+".dat", arma::arma_ascii);
}
	
void MainRom::use_deim(int dim, int t){
	
	int n_modes = dim;
	arma::Col<arma::uword> psub = deim_p.subvec(0, n_modes-1);
	arma::mat usub = modes_spatial(arma::span::all, arma::span(0, n_modes-1));
		
	arma::mat sigma = arma::diagmat(singular_values.subvec(0, n_modes-1));
	arma::mat snapshots_recon = modes_spatial.cols(0, n_modes-1)*sigma*modes_temporal.cols(0, n_modes-1).t();
	//singular_values.save("singular_values.data", arma::arma_ascii);
	//deim_p.save("p_idx.dat", arma::arma_ascii);
	arma::mat PP = usub * arma::inv(usub.rows(psub));

	//		usub.save("spatial_modes.dat", arma::arma_ascii);
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
	//snap_t_pod.save("snap_t_pod_"+std::to_string(t)+".dat", arma::arma_ascii);
	snap_t_deim.save("snap_t_deim_"+std::to_string(t)+".dat", arma::arma_ascii);
	//	psub.save("snap_iblank_"+std::to_string(t)+".dat", arma::arma_ascii);
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
		
}


arma::vec GemsRom::calc_deim(int ipartition_id, arma::vec r_s){
	int n_cols = preload_PP.n_cols;
	assert(n_cols == r_s.size());
	return preload_PP*r_s;
}



int GemsRom::get_deim_local_size(int ipartition_id){
	int n_cols = preload_PP.n_cols;
	return n_cols;
}

void GemsRom::get_deim_local_id_idx(int ipartition_id, int *ilocal_id, int *ivar){
	arma::umat tmp_s;
	tmp_s.load(directory + "residual_p_"+std::to_string(ipartition_id));
	int n_rows = tmp_s.n_rows;
	assert(ilocal_id != nullptr);
	assert(ivar != nullptr);
	for(int i=0; i< n_rows; i++){
		ilocal_id[i] = tmp_s(i, 4);
		ivar[i] = tmp_s(i, 2);
	}
	std::cout<<"get_deim_local_id_idx"<<std::endl;
}

void GemsRom::calc_deim(int ipartition_id, double *r_s, double *deim_r){
	assert(r_s != nullptr);
	assert(deim_r != nullptr);
	arma::vec r_s_v = arma::vec(r_s, get_deim_local_size(ipartition_id), false, true);
	double mean = 0.0;
	double stddev = 1.0;
	//	arma::umat tmp_s;
	//tmp_s.load("deim_p_"+std::to_string(ipartition_id));
	//m->snap_mean_v.print();
	//std::cout<<"HERE before "<<deim_r_v.size()<<" "<<r_s_v.size()<<std::endl;	
	for(int i=0; i<r_s_v.size(); i++){
		int ivar = preload_tmp_idx(i, 2);
		if(m->isnormalize < 4){
			//mean = m->snap_mean_v(ivar);
		}
		else{
			int gid = preload_tmp_idx(i, 0);
			//mean = m->snap_mean(gid);
		}
		stddev = m->snap_std_v(ivar);
		r_s_v(i) = (r_s_v(i) - mean)/stddev;
	}


	arma::vec deim_r_v = calc_deim(ipartition_id, r_s_v);
	
	for(int i=0; i<m->ndof; i++){
		deim_r[i] = deim_r_v(i);
	}

	//std::cout<<"HERE after "<<deim_r_v.size()<<" "<<r_s_v.size()<<std::endl;
	// this is a temporary fix to avoid memory leak, ideally the object
	// should delete itself
	try{
		//delete[] deim_r_v.memptr();
	}
	catch(...){
		
	}
}

arma::mat GemsRom::load_snapshots(std::string suffix){
	arma::mat snapshots;
	std::string filename = directory + "snapshots_solution.bin" + suffix;
	println("Reading file: " + filename);
	snapshots.load(filename, arma::arma_binary);
	return snapshots;
}
void GemsRom::initialize(int ipartition_id){
	try{
		m = new MainRom();
	}
	catch(std::bad_alloc &e){
		std::cout<<e.what()<<std::endl;
	}
	m->isnormalize = 4;
	directory = "RomFiles/";
	m->read_config_file(directory);
	//arma::mat snapshots = load_snapshots("_deim");
	//m->set_snapshots(snapshots);
	int num_processor = load_partition_info();
	m->load_modes("", directory);
	
	preload_tmp_idx.load(directory + "residual_p_"+std::to_string(ipartition_id), arma::arma_binary);
	preload_PP.load(directory + "PP_p_"+std::to_string(ipartition_id), arma::arma_binary);
}

int GemsRom::load_partition_info(){
	buffer.load(directory + "partition.bin", arma::arma_binary);
	local_id = buffer.col(0);
	partition_id = buffer.col(1);
	num_processor = arma::max(partition_id) + 1; 
	println("Number of processors = "); println(num_processor);
	return num_processor;
}

int GemsRom::get_local_id(int ipartition_id, int iglobal_id){
	assert(iglobal_id >= 0 && ipartition_id >= 0);
	assert(partition_id[iglobal_id] == ipartition_id);
	return local_id[iglobal_id];
}

int GemsRom::get_partition_id(int iglobal_id){
	return partition_id[iglobal_id];
}

int GemsRom::get_global_id(int ipartition_id, int ilocal_id){
	assert(ilocal_id >= 0 && ipartition_id >= 0);
	arma::uvec idx = arma::find(local_id == ilocal_id);
	for(arma::uword i=0; i<idx.size(); i++){
		int global_id = idx[i];
		int partition_idx = partition_id[idx[i]];
		if(partition_idx == ipartition_id){
			return idx[i];
		}
	}
	return -1;
}

arma::uvec GemsRom::get_index(int ipartition_id, arma::uvec idx){
	arma::uvec var_idx(m->nspecies*idx.size());
	for(int i=0; i<idx.size(); i++){
		//		int gid = get_global_id(ipartition_id, i);
		for(int k=0; k<m->nspecies; k++){
			var_idx(i*m->nspecies + k) = idx(i)*m->nspecies + k;
		}
	}
	return var_idx;
}

arma::vec GemsRom::get_uhat(int ipartition_id, arma::vec q, int n_mode){
	arma::uvec idx = arma::find(partition_id == ipartition_id);
	assert(q.size() == m->nspecies*idx.size());
	arma::uvec var_idx = get_index(ipartition_id, idx);
	arma::vec uhat = m->projection_on_basis(q, n_mode, var_idx);
	return uhat;
}

arma::vec GemsRom::get_u(int ipartition_id, arma::vec qhat, int n_mode){
	assert(qhat.size() == n_mode);
	arma::uvec idx = arma::find(partition_id == ipartition_id);
	arma::uvec var_idx = get_index(ipartition_id, idx);
	arma::vec u = m->projection_from_basis(qhat, n_mode, var_idx);
	assert(u.size() == idx.size()*m->nspecies);
	return u;
}


void GemsRom::get_uhat(int ipartition_id, double* q, int nq, double *qhat, int n_mode){
	assert(q != nullptr);
	assert(qhat != nullptr);
	arma::vec qvec = arma::vec(q, nq);
	arma::vec qhatvec = get_uhat(ipartition_id, qvec, n_mode);
	assert(qhatvec.size() == n_mode);

	for(int i=0; i<qhatvec.size(); i++){
		qhat[i] = qhatvec(i);
	}

}

void GemsRom::get_u(int ipartition_id, double* q, int nq, double *qhat, int n_mode){
	assert(q != nullptr);
	assert(qhat != nullptr);
	arma::vec qhatvec = arma::vec(qhat, n_mode);
	arma::vec qvec = get_u(ipartition_id, qhatvec, n_mode);
	assert(qvec.size() == nq);
	for(int i=0; i<qvec.size(); i++){
		q[i] = qvec(i);
	}
}


void GemsRom::get_uhat(int ipartition_id, int lq, double* q, int lqhat, double *qhat){
	get_uhat(ipartition_id, q, lq, qhat, lqhat);
}

void GemsRom::get_u(int ipartition_id, int lq, double* q, int lqhat, double *qhat){
	get_uhat(ipartition_id, q, lq, qhat, lqhat);
}
