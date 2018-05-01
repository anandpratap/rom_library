#include "main.h"
#include "main.hpp"
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
