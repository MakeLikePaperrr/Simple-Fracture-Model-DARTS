# Section of the Python code where we import all dependencies on third party Python modules/libaries or our own
# libraries (exposed C++ code to Python, i.e. darts.engines && darts.physics)
from model import Model
import numpy as np
import meshio
from darts.engines import *

"""
Some general comments on the code:
    - This is the unstructured reservoir example in DARTS
    - Most code is using basic Object Oriented Programming principles
        * Please see:   https://www.programiz.com/python-programming/object-oriented-programming &&
                        https://www.tutorialspoint.com/python/python_classes_objects.htm &&
                        https://python.swaroopch.com/oop.html
    - Unstructured reservoir module is based on: https://doi.org/10.2118/79699-MS
    - Physics are dead-oil, but should be able to replace with anything else!
    - Have fun! 
"""
# Typical structure of the main.py file (which is the actual file that is being run in PyCharm) is the following:
# 1) Create model object by calling the Model() constructor from the file model.py
#   --> This model class contains everything related to the model which is run in DARTS
#   --> From permeability, to porosity, to the physics used in the simulator, as well as the simulation parameters
m = Model()

# After constructing the model, the simulator needs to be initialized. The init() class method is called, which is
# inherited (https://www.python-course.eu/python3_inheritance.php) from the parent class DartsModel (found in
# darts/models/darts_model.py (NOTE: This is not the same as the__init__(self, **) method which each class (should)
# have).
m.init()
# redirect_darts_output('')

# Specify some other time-related properties (NOTE: all time parameters are in [days])
m.params.max_ts = 20  # Adjust the maximum time-step as desired (this is overwriting the max_ts specified in model.py)
size_report_step = 80  # Size of the reporting step (when output is writen to .vtk format)
num_report_steps = 25  # Number of reporting steps (see above)
start_time = 0  # Starting time of the simulation
end_time = size_report_step * num_report_steps  # End time of the simulation

# Before starting the simulation, store initial condition also in .vtk format:
ith_step = 0  # Store initial conditions as ../solution0.vtk
pressure_field = m.reservoir.mesh.pressure  # Extract initial pressure field
saturation_field = m.reservoir.mesh.composition  # Extract initial saturation field

# Read mesh using meshio (for writing output to VTK):
Mesh = meshio.read(m.reservoir.file_path)

# Properties for writing to vtk format:
# output_directory = 'trial_dir'  # Specify output directory here
output_directory = 'sol_{:s}_{:s}_{:s}'.format(m.mesh_type, m.bound_cond, m.physics_type)
num_wells_tot = len(m.reservoir.well_perf_loc[0]) + len(m.reservoir.well_perf_loc[1])  # Specify here how much wells are being used
# Specify here the number of properties you want to extract (properties are based on selected physics, see model):
tot_properties = 2

# Calculate the size of the properties vector:
tot_unknws = m.reservoir.unstr_discr.fracture_cell_count + m.reservoir.unstr_discr.matrix_cell_count + num_wells_tot*2

# Allocate and store the properties in an array:
property_array = np.empty((tot_unknws, tot_properties))
property_array[:, 0] = pressure_field
property_array[:, 1] = saturation_field

# Write to vtk using class methods of unstructured discretizer (uses within meshio write to vtk function):
m.reservoir.unstr_discr.write_to_vtk(output_directory, property_array, m.cell_property, ith_step)

# Run over all reporting time-steps:
for ith_step in range(num_report_steps):
    # Run engine for reporting_step [days]:
    # print('\n---------------------------SELF-PRINT---------------------------')
    # print('Current simulation time: {:f}'.format((ith_step+1)*size_report_step))
    # print('---------------------------SELF-PRINT---------------------------\n')
    m.run_python(size_report_step)

    # Export some results to VTK:
    # Store engine result:
    pressure_field = m.physics.engine.X[:-1:2]
    saturation_field = m.physics.engine.X[1::2]

    # Allocate and store the properties in an array:
    property_array = np.empty((tot_unknws, tot_properties))
    property_array[:, 0] = pressure_field
    property_array[:, 1] = saturation_field

    # Write to vtk using class methods of unstructured discretizer (uses within meshio write to vtk function):
    m.reservoir.unstr_discr.write_to_vtk(output_directory, property_array, m.cell_property, ith_step+1)

# After the simulation, print some of the simulation timers and statistics,
# newton iters, etc., how much time spent where:
m.print_timers()
m.print_stat()

import pandas as pd
time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
writer = pd.ExcelWriter('time_data.xlsx')
time_data.to_excel(writer, 'Sheet1')
writer.save()

from darts.tools.plot_darts import *
w = m.reservoir.wells[1]
ax2 = plot_temp_darts(w.name, time_data)

plt.show()
