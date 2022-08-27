from model_ensemble import Model
import numpy as np
import meshio
from darts.engines import *
import pandas as pd


def simulate_single_realization(filename_mesh, physics_type, bound_cond, const_perm,
                                frac_aper, poro, filename_log, filename_simdata, inj_well_coords, prod_well_coords,
                                OUT_DIR, max_time_step_size, size_report_step, end_time):
    redirect_darts_output(filename_log)
    m = Model(filename_mesh, physics_type, bound_cond, const_perm, frac_aper, poro, inj_well_coords, prod_well_coords)
    m.init()
    m.params.max_ts = max_time_step_size
    start_time = 0

    num_report_steps = np.ceil((end_time - start_time) / size_report_step).astype('int')  # Number of reporting steps (see above)
    sim_steps = np.zeros((num_report_steps,))
    sim_steps[:np.floor((end_time - start_time) / size_report_step).astype('int')] = size_report_step
    sim_steps[-1] = end_time - np.floor((end_time - start_time) / size_report_step).astype('int') * size_report_step
    max_time_steps = np.ones((num_report_steps,)) * m.params.max_ts
    # max_time_steps[0] = 5
    # max_time_steps[1] = 10

    # Before starting the simulation, store initial condition also in .vtk format:
    ith_step = 0  # Store initial conditions as ../solution0.vtk
    pressure_field = m.reservoir.mesh.pressure  # Extract initial pressure field
    saturation_field = m.reservoir.mesh.composition  # Extract initial saturation field

    # Read mesh using meshio (for writing output to VTK):
    Mesh = meshio.read(m.reservoir.file_path)

    # Properties for writing to vtk format:
    # output_directory = 'trial_dir'  # Specify output directory here
    output_directory = OUT_DIR
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
        m.run_python(sim_steps[ith_step])

        property_array = np.empty((tot_unknws, tot_properties))
        property_array[:, 0] = m.physics.engine.X[:-1:2]
        property_array[:, 1] = m.physics.engine.X[1::2]
        m.reservoir.unstr_discr.write_to_vtk(output_directory, property_array, m.cell_property, ith_step+1)

    m.print_timers()
    m.print_stat()

    time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
    writer = pd.ExcelWriter(filename_simdata)
    time_data.to_excel(writer, 'Sheet1')
    writer.save()
