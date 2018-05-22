import numpy as np
import scipy as sp
import scipy.io as io
import tecplot
from pylab import *
n = 20

ymid = 0.001
x = np.linspace(0.0009, 0.00999, 1000)
y = np.ones_like(x)*ymid

var_list = ["P", "U", "V", "T", "Y_CH4", "Y_O2", "Y_H2O", "Y_CO2"]


for i in range(200):
    print i
    dataset = tecplot.data.load_tecplot("../../../gems_2/rom_demo/UnsteadyFieldResults/gemsma_cmb_%i.plt"%(i+1), read_data_option=2)
    zone = dataset.zone(0)
    idxx = loadtxt("../../examples/deim_p.ascii", skiprows=2).astype(int)
    p = zeros(64000)
    p[idxx] = 1.0
    p = p.reshape([8000, 8], order="C")

    snap_t_pod = loadtxt("../../examples/snap_t_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")
    snap_t_deim = loadtxt("../../examples/snap_t_deim_%i.dat"%i, skiprows=2).reshape([8000, 8], order="C")

#    snap_t = loadtxt("../../examples/true_%i"%i, skiprows=2).reshape([8000, 8], order="C")

    
    for idx, var in enumerate(var_list):
        varname = "res_%i_iblank"%(idx+1)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = p[:,idx]
    for idx, var in enumerate(var_list):
        varname = "res_%i_snap"%(idx+1)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = snap_t_pod[:,idx]

        varname = "res_%i_deim"%(idx+1)
        dataset.add_variable(varname, dtypes=2, locations=0)
        zone.values(varname)[:] = snap_t_deim[:,idx]

    # idx = 0
    # data = np.load("snap_deim.npz")
    # figure()
    # plot(zone.values("Eq_1_residual").as_numpy_array().astype(np.float64), label="true")
    # plot(snap_t_pod[:,idx], label="should be true")
    # plot(snap_t_deim[:,idx], label="deim")
    # plot(snap_t[:,idx], label="true s")
    # plot(data["data"][i,:,idx], label="true ss")

    # legend()
    # show()
        
    #x = zone.values("x")[:]
    #y = zone.values("y")[:]
    #ind = np.lexsort((x, y))
    #print len(x)
    #print ind

    
    tecplot.data.save_tecplot_plt("all_%i.plt"%i)
