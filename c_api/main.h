#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
	class GemsRom;
	typedef GemsRom gemsrom;
#else
	typedef struct gemsrom gemsrom;
#endif
	gemsrom* create_gemsrom(void);
	void delete_gemsrom(gemsrom* gr);
	void gemsrom_initialize(gemsrom *gr);
	void gemsrom_get_uhat(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode);
	void gemsrom_get_u(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode);
	int gemsrom_get_deim_n(gemsrom *gr, int partition_id);
	void gemsrom_get_deim_idx(gemsrom *gr, int partition_id, int *ilocal_id, int *ivar);
	void gemsrom_calc_deim(gemsrom *gr, int partition_id, double *r_s, double *deim_r);
#ifdef __cplusplus
}
#endif



#endif
