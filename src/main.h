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
#ifdef __cplusplus
}
#endif



#endif
