import numpy as np
import scipy as sp
import scipy.io as io
import tecplot

def get_filename(fidx):
    filename = "../../../data_from_flux/UnsteadyFieldResults/gemsma_cmb_%i.plt"%fidx
    return filename

class TecplotImporter(object):
    def __init__(self):
        self.variables_of_interest = ["Static_Pressure","U","V","Temperature","CH4_mf","O2_mf","H2O_mf","CO2_mf"]

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
        
    def read_data(self, filename):
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
            if v.name in self.variables_of_interest:
                data[v.name] = zone.values(v.name).as_numpy_array().astype(np.float64)
        return data

    def read_plts(self, start=1, end=100, step=1):
        data = self.read_data(get_filename(start))
        ncv = data['U'].size
        nvar = len(self.variables_of_interest)
        nt = (end + 1 - start)/step
        data_array = np.zeros([nt, ncv, nvar], dtype=np.float64)
        for fidx, time in enumerate(range(start, end+1, step)):
            print "Loading plts %i of %i"%(fidx-start+1, end+1-start) 
            data = self.read_data(get_filename(time))
            for vidx, var in enumerate(self.variables_of_interest):
                data_array[fidx-1, :, vidx] = data[var]
        return data_array

    def save_plts(self):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        

if __name__ == "__main__":
    t = TecplotImporter()
    t.read_parition_file()
    data_array = t.read_plts(start=1, end=2000, step=1)
    nt, ncv, nvar = data_array.shape
    sizes = np.array([nt, ncv, nvar], dtype=np.int32)
    sizes.tofile("snapshots_shape.bin")
    data_array.tofile("snapshots_solution.bin")
    import sys
    sys.exit()
