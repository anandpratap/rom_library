import os
import numpy as np
import scipy as sp
import scipy.io as io
import argparse
import tecplot
import arma
from pylab import *

def get_filename(directory, fidx):
    filename = os.path.join(directory, "gemsma_cmb_%i.plt"%fidx)
    return filename

class TecplotImporter(object):
    def __init__(self):
        self.variables_of_interest = ["Static_Pressure","U","V","Temperature","CH4_mf","O2_mf","H2O_mf","CO2_mf"]
        self.variables_of_interest_residual = ["Eq_%i_residual"%i for i in range(1, len(self.variables_of_interest)+1)]

        
    def read_parition_file(self, filename):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        num_elements = zone.num_elements
        global_id = np.arange(1, num_elements+1).astype(np.int32)
        local_id = zone.values("id").as_numpy_array().astype(np.int32) - 1
        processor_id = zone.values("ipart").as_numpy_array().astype(np.int32) - 1
        data_array = np.concatenate([local_id, processor_id])
        data_array.tofile("partition.bin")
        print data_array.shape
        
    def read_data(self, filename, residual):
        print filename
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        num_variables = zone.num_variables
        num_points = zone.num_points
        num_elements = zone.num_elements
        time = zone.solution_time
        variables = list(zone.dataset.variables())
        data = self.read_vars_to_dict(dataset, self.variables_of_interest, num_elements)
        if residual:
            data_residual = self.read_vars_to_dict(dataset, self.variables_of_interest_residual, num_elements)

        output = [data]
        if residual:
            output.append(data_residual)
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
    
    def read_plts(self, start=1, end=100, step=1, residual=False, directory=None):
        ncv = self.query_ncv(get_filename(directory, start))
        nvar = len(self.variables_of_interest)
        nt = (end + 1 - start)/step

        data_array = np.zeros([nt, ncv, nvar], dtype=np.float64)
        if residual:
            data_array_residual = np.zeros([nt, ncv, nvar], dtype=np.float64)
        for fidx, time in enumerate(range(start, end+1, step)):
            print "Loading plts %i of %i"%(fidx+1, end+1-start) 
            output_data = self.read_data(get_filename(directory, time), residual)
            data_array[fidx, :, :] = output_data[0][:,:]
            if residual:
                data_array_residual[fidx, :, :] = output_data[1][:,:]
            
            
            
            # figure()
            # plot(data["Eq_1_residual"])
            # plot(data_array[fidx,:,0])
            # print "fidx-1", fidx
            # show()
        output = [data_array]
        if residual:
            output.append(data_array_residual)
        return output

    def save_plts(self):
        dataset = tecplot.data.load_tecplot(filename, read_data_option=2)
        zone = dataset.zone(0)
        

class TecplotExporter(object):
    def __init__(self):
        pass
    
    def save(self, nmodes, directory=os.getcwd(), residual=False, template="template.plt"):
        dataset = tecplot.data.load_tecplot(template, read_data_option=2)
        if not residual:
            modes = arma.load_mat(os.path.join(directory, "modes_spatial.bin"))
        else:
            modes = arma.load_mat(os.path.join(directory, "modes_spatial.bin_residual"))
        modes_vars_list = ["P_pod_mode", "U_pod_mode", "V_pod_mode", "T_pod_mode",
                               "Y_CH4_pod_mode", "Y_O2_pod_mode", "Y_H2O_pod_mode", "Y_CO2_pod_mode"]
        for v in modes_vars_list:
            dataset.add_variable(v, locations=tecplot.constant.ValueLocation.CellCentered)

        for i in range(nmodes):
            zone = dataset.zone(0)
            mode = modes[:,i]
            zone.solution_time = i
            for idx, v in enumerate(modes_vars_list):
                zone.values(v)[:] = mode[idx::len(modes_vars_list)]
            vars_to_write = ["x", "y"]
            vars_to_write.extend(modes_vars_list)
            variables = [dataset.variable(v) for v in vars_to_write]
            if not residual:
                file_format = "spatial_mode%i.dat"
            else:
                file_format = "spatial_mode_residual%i.dat"
            print "Write mode "+file_format%(i+1)
            tecplot.data.save_tecplot_ascii(file_format%(i+1), variables=variables)
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Utility for ROM framework')
    parser.add_argument('--iterations', nargs=2, type=int, help="Start and end (is included) iteration number for parsing and storing snapshots.")
    parser.add_argument("--process-plts",  action='store_true', help="Read plt files for solution (and residual) snapshots")
    parser.add_argument("--residual", action="store_true", help="Also process residual snapshots")
    parser.add_argument("--data-directory", type=str, help="Data directory", default=os.getcwd())
    
    parser.add_argument("--save-modes-plts", action="store_true", help="Save modes to tecplot format")
    parser.add_argument("--nmodes", type=int, help="Number of modes to save")
    parser.add_argument("--template", type=str, help="path to plt template")
    
    args = parser.parse_args()
    if args.process_plts:
        t = TecplotImporter()
        t.read_parition_file(os.path.join(args.data_directory, "../dfd_cell.plt"))
        start = args.iterations[0]
        end = args.iterations[1]
       
        if not args.residual:
            data_array, = t.read_plts(start=start, end=end, step=1, residual=False, directory=args.data_directory)
            
            nt, ncv, nvar = data_array.shape
            sizes = np.array([nt, ncv, nvar], dtype=np.int32)
            sizes.tofile("snapshots_shape.bin")

            data_array = np.reshape(data_array, [nt, ncv*nvar]).T
            arma.save_mat(data_array, "snapshots_solution.bin")
        else:
            data_array, data_array_residual = t.read_plts(start=start, end=end, step=1, residual=True, directory=args.data_directory)
            
            nt, ncv, nvar = data_array.shape
            sizes = np.array([nt, ncv, nvar], dtype=np.int32)
            sizes.tofile("snapshots_shape.bin")
            
            data_array = np.reshape(data_array, [nt, ncv*nvar]).T
            arma.save_mat(data_array, "snapshots_solution.bin")
            
            data_array_residual = np.reshape(data_array_residual, [nt, ncv*nvar]).T
            arma.save_mat(data_array_residual, "snapshots_solution.bin_residual")
    if args.save_modes_plts:
        t = TecplotExporter()
        t.save(args.nmodes, args.data_directory, False, args.template)
        if args.residual:
            t.save(args.nmodes, args.data_directory, True, args.template)
