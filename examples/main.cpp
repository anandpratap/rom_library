#include "main.hpp"
arma::vec test_function(arma::vec x, double t){
	double pi = arma::datum::pi;
	arma::vec y(x.size());
    y = (sin(t*x) + sin(t*pi/2*x) + sin(t*pi*x))*t;
	//y = 10*x%x*t*t;
	return y;
}

int main(){
	int N = 64; // number of grid points in spatial domain
	int Nt = 10000; // number of time steps
	int Tstart = 1; // start time
	int Tend = 3; // end time
	int dim = 2; // dimension of (a)DEIM space
	int w = 1; // window size
	int ms = 5; // number of additional(!) sampling points

	arma::vec t = arma::linspace<arma::vec>(Tstart, Tend, Nt);
	arma::vec x = arma::linspace<arma::vec>(0.0, 2.0*arma::datum::pi, N);
	arma::mat snapshots(N, Nt);
	for(arma::uword i=0; i<Nt; i++){
		snapshots.col(i) = test_function(x, t(i));
	}
	
	auto m = Main();
	m.set_snapshots(snapshots);
	m.calc_svd();
	m.calc_deim(dim);
	for(int i=0; i<Nt; i++){
		if(i%1 == 0){
			m.use_deim(dim, i);
		}
	}

	for(int i=0; i<Nt; i++){
		m.calc_adeim(dim, ms, i);
		if(i%1 == 0){
			m.use_adeim(dim, i);
		}
	}
	//m.load_snapshots();
	//m.calc_svd();
	//m.reconstruct();
	//m.calc_deim();
	
	//m.calc_adeim();
}
