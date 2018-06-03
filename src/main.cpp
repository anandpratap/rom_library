#include <iostream>
#include <fstream>
#include <armadillo>
#include <cassert>
#include "main.hpp"

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

void MainRom::set_snapshots(int idim, int insamples, double *isnapshots){
	assert(isnapshots != nullptr);
	arma::mat tmp_snapshots = arma::mat(isnapshots, idim, insamples, false, true);
	set_snapshots(tmp_snapshots, "");
}


void MainRom::set_snapshots(arma::mat isnapshots, std::string suffix){
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
		snap_mean_v = arma::vec(8);
		snap_std_v = arma::vec(8);
		int ncv = snapshots.n_rows/8;
		int nt = snapshots.n_cols;
		arma::vec Var(ncv*nt);
		for(int i=0; i<8; i++){
			int k=0;
			for(int icv=0; icv<ncv; icv++){
				for(int t=0; t<nt; t++){
					Var(k) = snapshots(icv*8 + i, t);
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
				snapshots(i, j) = (snapshots(i,j) - snap_mean_v(i%8))/snap_std_v(i%8);
			}
		}
		
		snap_mean_v.save("snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
	
	}
	
	else if(isnormalize == 3){
		snap_mean_v = arma::vec(8);
		snap_std_v = arma::vec(8);
		int ncv = snapshots.n_rows/8;
		int nt = snapshots.n_cols;
		arma::vec Var(ncv*nt);
		for(int i=0; i<8; i++){
			int k=0;
			for(int icv=0; icv<ncv; icv++){
				for(int t=0; t<nt; t++){
					Var(k) = snapshots(icv*8 + i, t);
					k+=1;
				}
			}
			assert(k == ncv*nt);
			snap_mean_v(i) = arma::mean(Var);
			snap_std_v(i) = arma::abs(Var).max();
		}
		//snap_mean_v(2) = snap_mean_v(1);
		snap_std_v(2) = snap_std_v(1);

		for(arma::uword j=0; j<snapshots.n_cols; j++){
			for(arma::uword i=0; i<snapshots.n_rows; i++){
				snapshots(i, j) = (snapshots(i,j) - snap_mean_v(i%8))/snap_std_v(i%8);
			}
		}
		
		snap_mean_v.save("snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.save("snap_std_v.bin"+suffix, arma::arma_binary);
	
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

void MainRom::load_modes(std::string suffix){
	isnormalize = 3;
	
	if(isnormalize == 1){
		snap_mean.load("snap_mean.bin"+suffix, arma::arma_binary);
		snap_std.load("snap_std.bin"+suffix, arma::arma_binary);
		
	}
	else if(isnormalize == 2 || isnormalize == 3){
		snap_mean_v.load("snap_mean_v.bin"+suffix, arma::arma_binary);
		snap_std_v.load("snap_std_v.bin"+suffix, arma::arma_binary);
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
	arma::uvec var_indices = arma::regspace<arma::uvec>(ivar_idx, 8, 64000-1);
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
	MatrixXd eigen_sub_basis = arma_matrix_to_eigen(sub_basis);
	std::cout<<" Calculating QR"<<std::endl;
	// do qr
	Eigen::ColPivHouseholderQR<MatrixXd> qr(eigen_sub_basis.transpose());
	//qr.compute(eigen_sub_basis.transpose());
	std::cout<<" DONE QR"<<std::endl;
	
	int m = sub_basis.n_cols;
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType P(qr.colsPermutation());
	auto p_ind = P.indices();
	std::cout<<P.rows()<<" "<<P.cols();
	std::cout<<p_ind.rows()<<" "<<p_ind.cols();
	arma::uvec pp(dim);
	for(int i=0; i<dim; i++){
		pp(i) = p_ind(i);
		std::cout<<p_ind(i)<<std::endl;
	}
	//std::cout<<P;

	if(deim_mode == DEIM_MODE_VECTOR){
		assert(pp.min() >= 0);
		assert(pp.max() < 64000);
		pp.save("deim_p.bin", arma::arma_binary);
		pp.save("deim_p.ascii", arma::arma_ascii);

	}
	else if(deim_mode == DEIM_MODE_CELL){
		assert(pp.min() >= 0);
		assert(pp.max() < 8000);
		arma::uvec ppp = arma::uvec(pp.size()*8);
		for(int i=0; i<pp.size(); i++){
			for(int j=0; j<8; j++){
				ppp(i*8 + j) = pp(i)*8 + j;
			}
		}
		ppp.save("deim_p.bin", arma::arma_binary);
		ppp.save("deim_p.ascii", arma::arma_ascii);
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
	deim_p.save("deim_p.bin", arma::arma_binary);
	deim_p.save("deim_p.ascii", arma::arma_ascii);
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
			y(k) = x(k)*snap_std_v(k%8) + snap_mean_v(k%8);
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
			y(k) = x(k)*snap_std_v(var_idx(k)%8) + snap_mean_v(var_idx(k)%8);
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
			y[k] = x[k]*snap_std_v(k%8) + snap_mean_v(k%8);
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
			y(k) = (x(k) - snap_mean_v(k%8))/snap_std_v(k%8);
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
			y(k) = (x(k) - snap_mean_v(var_idx(k)%8))/snap_std_v(var_idx(k)%8);
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
	tmp_s.load("deim_p_"+std::to_string(ipartition_id));
	int n_rows = tmp_s.n_rows;
	assert(ilocal_id != nullptr);
	assert(ivar != nullptr);
	for(int i=0; i< n_rows; i++){
		ilocal_id[i] = tmp_s(i, 4);
		ivar[i] = tmp_s(i, 2);
	}
}

void GemsRom::calc_deim(int ipartition_id, double *r_s, double *deim_r){
	assert(r_s != nullptr);
	assert(deim_r != nullptr);
	arma::vec r_s_v = arma::vec(r_s, get_deim_local_size(ipartition_id), false, true);
	double mean, stddev;
	//	arma::umat tmp_s;
	//tmp_s.load("deim_p_"+std::to_string(ipartition_id));
	//m->snap_mean_v.print();
	for(int i=0; i<r_s_v.size(); i++){
		int ivar = preload_tmp_idx(i, 2);
		mean = m->snap_mean_v(ivar);
		stddev = m->snap_std_v(ivar);
		r_s_v(i) = (r_s_v(i) - mean)/stddev;
	}
	arma::vec deim_r_v = calc_deim(ipartition_id, r_s_v);
	for(int i=0; i<64000; i++){
		deim_r[i] = deim_r_v(i);
	}
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
	std::string filename = "snapshots_solution.bin"+suffix;
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
	//arma::mat snapshots = load_snapshots("_deim");
	//m->set_snapshots(snapshots);
	load_partition_info();
	m->load_modes("_residual");
	
	preload_tmp_idx.load("deim_p_"+std::to_string(ipartition_id), arma::arma_binary);
	preload_PP.load("PP_p_"+std::to_string(ipartition_id), arma::arma_binary);
}

void GemsRom::load_partition_info(){
	arma::Mat<int> shape;
	shape.load("snapshots_shape.bin", arma::raw_binary);
	//	arma::Col<int> buffer(shape[1]*2);
	println(shape[1]);
	buffer.load("partition.bin", arma::raw_binary);

	//local_id(shape[1]);
	//partition_id(shape[1]);

	local_id = buffer.subvec(0, shape[1]-1);
	partition_id = buffer.subvec(shape[1], shape[1]*2-1);


	// for(int i=0; i<local_id.size(); i++){
	// 	if(i>=4000){
	// 		partition_id(i) = 1;
	// 		local_id(i) = i-4000;
	// 	}
	// 	else{
	// 		partition_id(i) = 0;
	// 		local_id(i) = i;
	// 	}
	// }
	num_processor = arma::max(partition_id) + 1; 
	println("Number of processors = "); println(num_processor);
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
	arma::uvec var_idx(8*idx.size());
	for(int i=0; i<idx.size(); i++){
		//		int gid = get_global_id(ipartition_id, i);
		for(int k=0; k<8; k++){
			var_idx(i*8 + k) = idx(i)*8 + k;
		}
	}
	return var_idx;
}

arma::vec GemsRom::get_uhat(int ipartition_id, arma::vec q, int n_mode){
	arma::uvec idx = arma::find(partition_id == ipartition_id);
	assert(q.size() == 8*idx.size());
	arma::uvec var_idx = get_index(ipartition_id, idx);
	arma::vec uhat = m->projection_on_basis(q, n_mode, var_idx);
	return uhat;
}

arma::vec GemsRom::get_u(int ipartition_id, arma::vec qhat, int n_mode){
	assert(qhat.size() == n_mode);
	arma::uvec idx = arma::find(partition_id == ipartition_id);
	arma::uvec var_idx = get_index(ipartition_id, idx);
	arma::vec u = m->projection_from_basis(qhat, n_mode, var_idx);
	assert(u.size() == idx.size()*8);
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
