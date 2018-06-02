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
        data = self.read_vars_to_dict(dataset, self.variables_of_interest, num_elements)
        if deim:
            data_deim = self.read_vars_to_dict(dataset, self.variables_of_interest_deim, num_elements)

        output = [data]
        if deim:
            output.append(data_deim)
        return output

    def read_vars_to_dict(self, dataset, variables_list, num_elements):
        zone = dataset.zone(0)
        variables = list(zone.dataset.variables())
        nvar = len(variables_list)
        data = np.zeros([num_elements, nvar], dtype=np.float64)
        for v in variables:
            if v.name in variables_list:
                idx = variables_list.index(v.name)
                data[:,idx] = zone.values(v.name).as_numpy_array().astype(np.float64)
        return data

    def query_ncv(self, filename):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        num_variables = zone.num_variables
        num_points = zone.num_points
        num_elements = zone.num_elements
        return num_elements
    
    def read_plts(self, start=1, end=100, step=1, deim=False):
        ncv = self.query_ncv(get_filename(start))
        nvar = len(self.variables_of_interest)
        nt = (end + 1 - start)/step

        data_array = np.zeros([nt, ncv, nvar], dtype=np.float64)
        if deim:
            data_array_deim = np.zeros([nt, ncv, nvar], dtype=np.float64)
        for fidx, time in enumerate(range(start, end+1, step)):
            print "Loading plts %i of %i"%(fidx+1, end+1-start) 
            output_data = self.read_data(get_filename(time), deim)
            if deim:
                    data_array_deim[fidx, :, :] = output_data[1][:,:]
            else:
                    data_array[fidx, :, :] = output_data[0][:,:]
            # figure()
            # plot(data["Eq_1_residual"])
            # plot(data_array[fidx,:,0])
            # print "fidx-1", fidx
            # show()
        output = [data_array]
        if deim:
            output.append(data_array_deim)
        return output

    def save_plts(self):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        

if __name__ == "__main__":
    t = TecplotImporter()
    t.read_parition_file()
    start = 5000
    end = 11000
    data_array, data_array_deim = t.read_plts(start=start, end=end, step=1, deim=True)
    nt, ncv, nvar = data_array.shape
    sizes = np.array([nt, ncv, nvar], dtype=np.int32)
    sizes.tofile("snapshots_shape.bin")
    data_array.tofile("snapshots_solution.bin")
    data_array_deim.tofile("snapshots_solution.bin"+"_deim", format="%.16f")
    np.savez("snap_deim.npz", data=data_array_deim)
