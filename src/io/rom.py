import numpy as np
import scipy as sp
import scipy.io as io
import tecplot
from pylab import *
def get_filename(fidx):
    filename = "../../../gems_2/rom_demo/UnsteadyFieldResults/gemsma_cmb_%i.plt"%fidx
    return filename

class TecplotImporter(object):
    def __init__(self):
        self.variables_of_interest = ["Static_Pressure","U","V","Temperature","CH4_mf","O2_mf","H2O_mf","CO2_mf"]
        self.variables_of_interest_deim = ["Eq_%i_residual"%i for i in range(1, len(self.variables_of_interest)+1)]

        
    def read_parition_file(self, filename="../../../data_from_flux/dfd_cell.plt"):
        print filename
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        num_elements = zone.num_elements
        global_id = np.arange(1, num_elements+1).astype(np.int32)
        local_id = zone.values("id").as_numpy_array().astype(np.int32) - 1
        processor_id = zone.values("ipart").as_numpy_array().astype(np.int32) - 1
        data_array = np.concatenate([local_id, processor_id])
        data_array.tofile("partition.bin")
        print data_array.shape
        
    def read_data(self, filename, deim):
        print filename
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        num_variables = zone.num_variables
        num_points = zone.num_points
        num_elements = zone.num_elements
        time = zone.solution_time
        variables = list(zone.dataset.variables())
        data = {}
        for v in variables:
            if deim:
                if v.name in self.variables_of_interest_deim:
                    data[v.name] = zone.values(v.name).as_numpy_array().astype(np.float64)
            else:
                if v.name in self.variables_of_interest:
                    data[v.name] = zone.values(v.name).as_numpy_array().astype(np.float64)
        return data

    def read_plts(self, start=1, end=100, step=1, deim=False):
        data = self.read_data(get_filename(start), deim)
        if deim:
            ncv = data[self.variables_of_interest_deim[0]].size
        else:
            ncv = data[self.variables_of_interest[0]].size
        nvar = len(self.variables_of_interest)
        nt = (end + 1 - start)/step
        data_array = np.zeros([nt, ncv, nvar], dtype=np.float64)
        for fidx, time in enumerate(range(start, end+1, step)):
            print "Loading plts %i of %i"%(fidx+1, end+1-start) 
            data = self.read_data(get_filename(time), deim)
            if deim:
                for vidx, var in enumerate(self.variables_of_interest_deim):
                    data_array[fidx, :, vidx] = data[var]
            else:
                for vidx, var in enumerate(self.variables_of_interest):
                    data_array[fidx, :, vidx] = data[var]
            # figure()
            # plot(data["Eq_1_residual"])
            # plot(data_array[fidx,:,0])
            # print "fidx-1", fidx
            # show()
        return data_array

    def save_plts(self):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        

if __name__ == "__main__":
    t = TecplotImporter()
    t.read_parition_file()
    end = 2000
    data_array = t.read_plts(start=1, end=end, step=1)
    data_array_deim = t.read_plts(start=1, end=end, step=1, deim=True)
    nt, ncv, nvar = data_array.shape
    sizes = np.array([nt, ncv, nvar], dtype=np.int32)
    sizes.tofile("snapshots_shape.bin")
    data_array.tofile("snapshots_solution.bin")
    data_array_deim.tofile("snapshots_solution.bin"+"_deim", format="%.16f")
    np.savez("snap_deim.npz", data=data_array_deim)
