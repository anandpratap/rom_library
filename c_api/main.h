#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
	class GemsRom;
	typedef GemsRom gemsrom;
	class MainRom;
	typedef MainRom mainrom;
#else
	typedef struct GemsRom gemsrom;
	typedef struct MainRom mainrom;
#endif

	mainrom* create_mainrom(void);
	void mainrom_save_modes(mainrom *mr, char* suffix="");
	void mainrom_load_modes(mainrom *mr, char* suffix="");
	void mainrom_calc_svd(mainrom *mr);

	
	gemsrom* create_gemsrom(void);
	void delete_gemsrom(gemsrom* gr);
	void gemsrom_initialize(gemsrom *gr, int ipartition_id);
	void gemsrom_get_uhat(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode);
	void gemsrom_get_u(gemsrom *gr, int partition_id, double *q, int nq, double *qhat, int n_mode);
	void gemsrom_get_deim_n(gemsrom *gr, int partition_id, int &in);
	void gemsrom_get_deim_idx(gemsrom *gr, int partition_id, int isize, int *ilocal_id, int *ivar);
	void gemsrom_calc_deim(gemsrom *gr, int partition_id, int isize, double *r_s, double *deim_r);
	void gemsrom_renormalize(gemsrom *gr, int isize, double *x, double *y);
	void gemsrom_get_global_id(gemsrom *gr, int ipartition_id, int ilocal_id, int &iglobal_id);

#ifdef __cplusplus
}
#endif



#endif
