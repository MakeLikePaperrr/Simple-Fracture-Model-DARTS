import numpy as np
from sim_single_model import simulate_single_realization
import os


# Specify some other time-related properties (NOTE: all time parameters are in [days])
max_time_step_size = 20  # Adjust the maximum time-step as desired (this is overwriting the max_ts specified in model.py)
size_report_step = 80  # Size of the reporting step (when output is writen to .vtk format)
end_time = 1500  # End time of the simulation

NR_REAL = 10
BASE_DIR_NAME = 'ensemble_2'

inj_well_coords = [[100, 100, 25], [900, 900, 25]]
prod_well_coords = [[100, 900, 25], [900, 100, 25]]

physics_type = 'dead_oil'  # geothermal, dead_oil
bound_cond = 'wells_in_nearest_cell'  # wells_in_frac, wells_in_mat, wells_in_nearest_cell
const_perm = 2
frac_aper = 5e-4
poro = 0.2

use_clean = True
char_len = 20
merg_fac = 0.86

if not use_clean:
    DIR_INPUT = os.path.join('..', '2D-Fracture-Generation-Tool', f'{BASE_DIR_NAME}_meshes')
    DIR_OUTPUT = f'{BASE_DIR_NAME}_meshes_output_{physics_type}'
    BASE_FILENAME_MESH = lambda ith_real: os.path.join(DIR_INPUT, f'{BASE_DIR_NAME}_mesh_{ith_real}.msh')
else:
    DIR_INPUT = os.path.join('..', '2D-Fracture-Generation-Tool', f'{BASE_DIR_NAME}_clean_meshes')
    DIR_OUTPUT = f'{BASE_DIR_NAME}_clean_meshes_output_{physics_type}'
    BASE_FILENAME_MESH = lambda ith_real: os.path.join(DIR_INPUT,
                                                       f'{BASE_DIR_NAME}_mesh_{ith_real}_mergefac_{merg_fac}_clean_lc_{char_len}.msh')

BASE_FILENAME_LOG = lambda ith_real: os.path.join(DIR_OUTPUT, f'log_{ith_real}.log')
BASE_FILENAME_OUTPUT = lambda ith_real: os.path.join(DIR_OUTPUT, f'prod_data_real_{ith_real}.xlsx')

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

# todo: write file for plotting all sim.data (reading outputs from excel files)
#   and test a few different configurations etc. to make sure it works!
