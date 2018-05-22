#include "main.h"
#include "../src/main.hpp"
#include <iostream>

gemsrom *create_gemsrom(void){
	return new GemsRom();
}

void delete_gemsrom(gemsrom* gr){
	delete gr;
};

void gemsrom_initialize(gemsrom* gr){
	gr->initialize();
};

void gemsrom_get_uhat(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode){
	gr->get_uhat(partition_id, q, nq, qhat, n_mode);
};

void gemsrom_get_u(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode){
	gr->get_u(partition_id, q, nq, qhat, n_mode);
};


int gemsrom_get_deim_n(gemsrom *gr, int partition_id){
	gr->get_deim_local_size(partition_id);
};
void gemsrom_get_deim_idx(gemsrom *gr, int partition_id, int *ilocal_id, int *ivar){
	gr->get_deim_local_id_idx(partition_id, ilocal_id, ivar);
};
void gemsrom_calc_deim(gemsrom *gr, int partition_id, double *r_s, double *deim_r){
	gr->calc_deim(partition_id, r_s, deim_r);
};
