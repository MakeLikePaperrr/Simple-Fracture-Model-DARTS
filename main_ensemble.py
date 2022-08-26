import numpy as np
from sim_single_model import simulate_single_realization
import os


# Specify some other time-related properties (NOTE: all time parameters are in [days])
max_time_step_size = 20  # Adjust the maximum time-step as desired (this is overwriting the max_ts specified in model.py)
size_report_step = 80  # Size of the reporting step (when output is writen to .vtk format)
end_time = 1500  # End time of the simulation

NR_REAL = 2
DIR_INPUT = os.path.join('..', '2D-Fracture-Generation-Tool', 'ensemble_1_meshes')
# DIR_INPUT = os.path.join('..', '2D-Fracture-Generation-Tool', 'ensemble_1_clean_meshes')
DIR_OUTPUT = 'ensemble_1_meshes_output_deadoil'
# DIR_OUTPUT = 'ensemble_1_clean_meshes_output_deadoil'

BASE_FILENAME_MESH = lambda ith_real: os.path.join(DIR_INPUT, f'ensemble_1_mesh_{ith_real}.msh')
# BASE_FILENAME = lambda ith_real: f'ensemble_1_mesh_{ith_real}_mergefac_0.86_clean_lc_16.msh'
BASE_FILENAME_LOG = lambda ith_real: os.path.join(DIR_OUTPUT, f'log_{ith_real}.log')
BASE_FILENAME_OUTPUT = lambda ith_real: os.path.join(DIR_OUTPUT, f'prod_data_real_{ith_real}.xlsx')

inj_well_coords = [[100, 100, 25], [900, 900, 25]]
prod_well_coords = [[100, 900, 25], [900, 100, 25]]

physics_type = 'dead_oil'  # geothermal, dead_oil
bound_cond = 'wells_in_nearest_cell'  # wells_in_frac, wells_in_mat, wells_in_nearest_cell
const_perm = 10
frac_aper = 1e-3
poro = 0.2

if not os.path.exists(DIR_OUTPUT):
    os.makedirs(DIR_OUTPUT)

for i in range(1, NR_REAL + 1):
    filename_mesh = BASE_FILENAME_MESH(i)
    filename_log = BASE_FILENAME_LOG(i)
    filename_simdata = BASE_FILENAME_OUTPUT(i)
    OUT_DIR = os.path.join(DIR_OUTPUT, f'real_{i}')
    simulate_single_realization(filename_mesh, physics_type, bound_cond, const_perm,
                                frac_aper, poro, filename_log, filename_simdata,
                                inj_well_coords, prod_well_coords, OUT_DIR, max_time_step_size,
                                size_report_step, end_time)
