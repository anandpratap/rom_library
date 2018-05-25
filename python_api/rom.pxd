cdef extern from "../c_api/main.h":
    ctypedef struct gemsrom:
        pass
    ctypedef struct mainrom:
        pass

    mainrom* create_mainrom();
    void mainrom_save_modes(mainrom *mr, char* suffix);
    void mainrom_load_modes(mainrom *mr, char* suffix);
    void mainrom_calc_svd(mainrom *mr);
    
    gemsrom* create_gemsrom();
    void delete_gemsrom(gemsrom* gr);
    void gemsrom_initialize(gemsrom *gr, int ipartition_id);

    

