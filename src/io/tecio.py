import numpy as np
import scipy as sp
import scipy.io as io
import tecplot
from pylab import *
n = 20


var_list = ["P", "U", "V", "T", "Y_CH4", "Y_O2", "Y_H2O", "Y_CO2"]



for i in range(0, 2000, 100):
    dataset = tecplot.data.load_tecplot("../../data_from_flux/UnsteadyFieldResults/gemsma_cmb_%i.plt"%(i+1), read_data_option=2)
    zone = dataset.zone(0)

    snap_t = loadtxt("snap_t_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")
    #snap_t_pod = loadtxt("snap_t_pod_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")
    snap_t_deim = loadtxt("snap_t_deim_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")
    snap_t_adeim = loadtxt("snap_t_adeim_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")

    idxx = loadtxt("iblank_adeim_%i.dat"%i, skiprows=2).astype(int)
    p = zeros(64000)
    p[idxx] = 1.0
    p = p.reshape([8000, 8], order="C")
    
    
    for idx, var in enumerate(var_list):
        varname = "%s_iblank"%(var)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = p[:,idx]

    
    for idx, var in enumerate(var_list):
        varname = "%s_snap_t"%(var)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = snap_t[:,idx]

 #       varname = "%s_snap_t_pod_%i"%(var, i)
 #       dataset.add_variable(varname, dtypes=2, locations=0)
 #       zone.values(varname)[:] = snap_t_pod[:,idx]

        varname = "%s_snap_t_deim"%(var)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = snap_t_deim[:,idx]

        varname = "%s_snap_t_adeim"%(var)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = snap_t_adeim[:,idx]
    tecplot.data.save_tecplot_plt("all_%i.plt"%i)
