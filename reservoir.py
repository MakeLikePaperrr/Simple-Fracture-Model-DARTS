from darts.engines import conn_mesh, ms_well, ms_well_vector, index_vector, value_vector
import numpy as np
from math import inf, pi
from darts.mesh.unstruct_discretizer import UnstructDiscretizer
from itertools import compress


# Definitions for the unstructured reservoir class:
class UnstructReservoir:
    def __init__(self, permx, permy, permz, frac_aper, mesh_file, poro, bound_cond, physics_type):
        """
        Class constructor for UnstructReservoir class
        :param permx: Matrix permeability in the x-direction (scalar or vector)
        :param permy: Matrix permeability in the y-direction (scalar or vector)
        :param permz: Matrix permeability in the z-direction (scalar or vector)
        :param frac_aper: Aperture of the fracture (scalar or vector)
        :param mesh_file: Name and relative path to the mesh-file (string)
        :param poro: Matrix (and fracture?) porosity (scalar or vector)
        :param bound_cond: switch which determines the type of boundary conditions used (string)
        """
        # Create mesh object (C++ object used by DARTS for all mesh related quantities):
        self.mesh = conn_mesh()

        # Specify well index and store matrix geometry:
        self.file_path = mesh_file

        # Construct instance of Unstructured Discretization class:
        self.unstr_discr = UnstructDiscretizer(permx=permx, permy=permy, permz=permz, frac_aper=frac_aper,
                                               mesh_file=mesh_file)

        # Use class method load_mesh to load the GMSH file specified above:
        self.unstr_discr.load_mesh()

        # Calculate cell information of each geometric element in the .msh file:
        self.unstr_discr.calc_cell_information()

        # Store volumes and depth to single numpy arrays:
        self.unstr_discr.store_volume_all_cells()
        self.unstr_discr.store_depth_all_cells()
        self.unstr_discr.store_centroid_all_cells()

        # Perform discretization:
        cell_m, cell_p, tran, tran_thermal = self.unstr_discr.calc_connections_all_cells()

        # Write to files (in case someone needs this for Eclipse or other simulator):
        self.unstr_discr.write_conn2p_to_file(cell_m, cell_p, tran, file_name='conn2p.dat')
        self.unstr_discr.write_conn2p_therm_to_file(cell_m, cell_p, tran, tran_thermal, file_name='conn2p.dat.connsn')
        self.unstr_discr.write_volume_to_file(file_name='vol.dat')
        self.unstr_discr.write_depth_to_file(file_name='depth.dat')

        # Initialize mesh using built connection list
        if physics_type == 'dead_oil':
            # Initialize mesh with just three parameters (cell_m, cell_p, trans):
            self.mesh.init(index_vector(cell_m), index_vector(cell_p), value_vector(tran))

        elif physics_type == 'geothermal':
            # Initialize mesh with all four parameters (cell_m, cell_p, trans, trans_D):
            self.mesh.init(index_vector(cell_m), index_vector(cell_p), value_vector(tran), value_vector(tran_thermal))
            # self.mesh.init('conn2p.dat.connsn')

        # Store number of control volumes (NOTE: in case of fractures, this includes both matrix and fractures):
        self.nb = self.unstr_discr.volume_all_cells.size
        self.num_frac = self.unstr_discr.fracture_cell_count
        self.num_mat = self.unstr_discr.matrix_cell_count

        # Create numpy arrays wrapped around mesh data (no copying, this will severely slow down the process!)
        self.poro = np.array(self.mesh.poro, copy=False)
        self.depth = np.array(self.mesh.depth, copy=False)
        self.volume = np.array(self.mesh.volume, copy=False)

        # rock thermal properties
        self.hcap = np.array(self.mesh.heat_capacity, copy=False)
        self.conduction = np.array(self.mesh.rock_cond, copy=False)

        # Since we use copy==False above, we have to store the values by using the Python slicing option, if we don't
        # do this we will overwrite the variable, e.g. self.poro = poro --> overwrite self.poro with the variable poro
        # instead of storing the variable poro in self.mesh.poro (therefore "numpy array wrapped around mesh data!!!):
        self.poro[:] = poro
        self.depth[:] = self.unstr_discr.depth_all_cells
        self.volume[:] = self.unstr_discr.volume_all_cells

        # Calculate well_index (very primitive way....):
        self.well_index = np.mean(tran) * 1

        # Store type of boundary condition:
        self.bound_cond = bound_cond

        if bound_cond == 'const_pres_rate':
            # Set-up dictionary with data for boundary cells:
            boundary_data = {}  # Dictionary containing boundary condition data (coordinate and value of boundary):
            boundary_data['first_boundary_dir'] = 'X'  # Indicates the boundary is located at constant X (in this case!)
            # Constant X-coordinate value at which the boundary is located (used to be 3.40885):
            boundary_data['first_boundary_val'] = np.min(self.unstr_discr.mesh_data.points[:, 0])

            # Same as above but for the second boundary condition!
            boundary_data['second_boundary_dir'] = 'X'
            # Constant X-coordinate value at which the boundary is located (used to be 13.0014):
            boundary_data['second_boundary_val'] = np.max(self.unstr_discr.mesh_data.points[:, 0])

            # Create empty containers for cells which lay on the boundary:
            self.left_boundary_cells = np.array([])
            self.right_boundary_cells = np.array([])

            # Calculate boundary cells using the calc_boundary_cells method:
            self.calc_boundary_cells(boundary_data)

            # Calc maximum size of well cells (used to have more homogeneous injection conditions by scaling the WI):
            dummy_vol = np.array(self.volume, copy=True)
            self.max_well_vol = np.max([np.max(dummy_vol[self.left_boundary_cells]),
                                        np.max(dummy_vol[self.right_boundary_cells])])
        elif bound_cond == 'wells_in_frac':
            if 0:
                # Find fractures in the bottom-left and top-right of reservoir:
                bot_left_frac = np.array([inf, inf, inf])
                top_right_frac = np.array([0, 0, 0])
                bot_left_id = 0
                top_right_id = self.unstr_discr.fracture_cell_count - 1

                for ith_frac in self.unstr_discr.frac_cell_info_dict:
                    current_frac = self.unstr_discr.frac_cell_info_dict[ith_frac].centroid

                    # Find top right corner fracture:
                    if np.sqrt(current_frac[0]**2 + current_frac[1]**2) > \
                            np.sqrt(top_right_frac[0] ** 2 + top_right_frac[1] ** 2):
                        # Store new larger centroid of fracture:
                        top_right_frac = current_frac
                        top_right_id = ith_frac

                    # Find bottom left corner fracture:
                    if np.sqrt(current_frac[0] ** 2 + current_frac[1] ** 2) < \
                            np.sqrt(bot_left_frac[0] ** 2 + bot_left_frac[1] ** 2):
                        # Store new larger centroid of fracture:
                        bot_left_frac = current_frac
                        bot_left_id = ith_frac

                self.well_perf_loc = np.array([bot_left_id, top_right_id])
            else:
                # Find closest control volume to dummy_well point:
                self.injection_wells = []
                dummy_well_inj = [[50, 160, 25]]

                self.store_dist_to_well_inj = np.zeros((len(dummy_well_inj),))
                self.store_coord_well_inj = np.zeros((len(dummy_well_inj), 3))
                ii = 0
                for ith_inj in dummy_well_inj:
                    dist_to_well_point = np.linalg.norm(self.unstr_discr.centroid_all_cells[:self.num_frac] - ith_inj,
                                                        axis=1)
                    cell_id = np.argmin(dist_to_well_point)
                    self.injection_wells.append(cell_id)

                    self.store_coord_well_inj[ii, :] = self.unstr_discr.centroid_all_cells[cell_id]
                    self.store_dist_to_well_inj[ii] = np.min(dist_to_well_point)
                    ii += 1

                self.production_wells = []
                dummy_well_prod = [[925, 960, 25]]

                self.store_dist_to_well_prod = np.zeros((len(dummy_well_prod),))
                self.store_coord_well_prod = np.zeros((len(dummy_well_prod), 3))
                ii = 0
                for ith_prod in dummy_well_prod:
                    dist_to_well_point = np.linalg.norm(self.unstr_discr.centroid_all_cells[:self.num_frac] - ith_prod,
                                                        axis=1)
                    cell_id = np.argmin(dist_to_well_point)
                    self.production_wells.append(cell_id)

                    self.store_coord_well_prod[ii, :] = self.unstr_discr.centroid_all_cells[cell_id]
                    self.store_dist_to_well_prod[ii] = np.min(dist_to_well_point)
                    ii += 1

                self.well_perf_loc = np.array([self.injection_wells, self.production_wells])

        else:
            print("--------ERROR SPECIFY CORRECT PHYSICS NAME--------")

        # Create empty list of wells:
        self.wells = []

    def calc_boundary_cells(self, boundary_data):
        """
        Class method which calculates constant boundary values at a specif constant x,y,z-coordinate
        :param boundary_data: dictionary with the boundary location (X,Y,Z, and location)
        :return:
        """
        # Specify boundary cells, simply set specify the single coordinate which is not-changing and its value:
        # First boundary:
        index = []  # Dynamic list containing indices of the nodes (points) which lay on the boundary:
        if boundary_data['first_boundary_dir'] == 'X':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 0] == boundary_data['first_boundary_val']
        elif boundary_data['first_boundary_dir'] == 'Y':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 1] == boundary_data['first_boundary_val']
        elif boundary_data['first_boundary_dir'] == 'Z':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 2] == boundary_data['first_boundary_val']

        # Convert dynamic list to numpy array:
        left_boundary_points = np.array(list(compress(range(len(index)), index)))

        # Second boundary (same as above):
        index = []
        if boundary_data['second_boundary_dir'] == 'X':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 0] == boundary_data['second_boundary_val']
        elif boundary_data['second_boundary_dir'] == 'Y':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 1] == boundary_data['second_boundary_val']
        elif boundary_data['second_boundary_dir'] == 'Z':
            # Check if first coordinate of points is on the boundary:
            index = self.unstr_discr.mesh_data.points[:, 2] == boundary_data['second_boundary_val']

        right_boundary_points = np.array(list(compress(range(len(index)), index)))

        # Find cells containing boundary cells, for wedges or hexahedrons, the boundary cells must contain,
        # on the X or Y boundary four nodes exactly!
        #     0------0          0
        #    /     / |         /  \
        #  0------0  0        0----0
        #  |      | /         |    |
        #  0------0           0----0
        # Hexahedron       Wedge (prism)
        # Create loop over all matrix cells which are of the geometry 'matrix_cell_type'
        left_count = 0  # Counter for number of left matrix cells on the boundary
        left_boundary_cells = {}  # Dictionary with matrix cells on the left boundary
        for geometry in self.unstr_discr.geometries_in_mesh_file:
            if geometry in self.unstr_discr.available_matrix_geometries:
                # Matrix geometry found, check if any matrix control volume has exactly 4 nodes which intersect with
                # the left_boundary_points list:
                for ith_cell, ith_row in enumerate(
                        self.unstr_discr.mesh_data.cells_dict[geometry]):

                    if len(set.intersection(set(ith_row), set(left_boundary_points))) == 4:
                        # Store cell since it is on the left boundary:
                        left_boundary_cells[left_count] = ith_cell
                        left_count += 1

        right_count = 0
        right_boundary_cells = {}
        for geometry in self.unstr_discr.geometries_in_mesh_file:
            if geometry in self.unstr_discr.available_matrix_geometries:
                # Matrix geometry found, check if any matrix control volume has exactly 4 nodes which intersect with
                # the right_boundary_points list:
                for ith_cell, ith_row in enumerate(
                        self.unstr_discr.mesh_data.cells_dict[geometry]):
                    if len(set.intersection(set(ith_row), set(right_boundary_points))) == 4:
                        # Store cell since it is on the left boundary:
                        right_boundary_cells[right_count] = ith_cell
                        right_count += 1

        self.left_boundary_cells = np.array(list(left_boundary_cells.values()), dtype=int) + \
                                   self.unstr_discr.fracture_cell_count
        self.right_boundary_cells = np.array(list(right_boundary_cells.values()), dtype=int) + \
                                    self.unstr_discr.fracture_cell_count
        return 0

    def add_well(self, name, depth):
        """
        Class method which adds wells heads to the reservoir (Note: well head is not equal to a perforation!)
        :param name:
        :param depth:
        :return:
        """
        well = ms_well()
        well.name = name
        well.segment_volume = 0.0785 * 40  # 2.5 * pi * 0.15**2 / 4
        well.well_head_depth = depth
        well.well_body_depth = depth
        well.segment_transmissibility = 1e5
        well.segment_depth_increment = 1
        self.wells.append(well)
        return 0

    def add_perforation(self, well, res_block, well_index):
        """
        Class method which ads perforation to each (existing!) well
        :param well: data object which contains data of the particular well
        :param res_block: reservoir block in which the well has a perforation
        :param well_index: well index (productivity index)
        :return:
        """
        well_block = 0
        well.perforations = well.perforations + [(well_block, res_block, well_index)]
        return 0

    def init_wells(self):
        """
        Class method which initializes the wells (adding wells and their perforations to the reservoir)
        :return:
        """
        # Add injection well:
        self.add_well("I1", 0.5)
        if self.bound_cond == 'const_pres_rate':
            # Perforate all boundary cells:
            for nth_perf in range(len(self.left_boundary_cells)):
                well_index = self.mesh.volume[self.left_boundary_cells[nth_perf]] / self.max_well_vol * self.well_index
                self.add_perforation(well=self.wells[-1], res_block=self.left_boundary_cells[nth_perf],
                                     well_index=well_index)

        elif self.bound_cond == 'wells_in_frac':
            # Only perforating the single fracture/matrix block
            self.add_perforation(self.wells[-1], res_block=self.well_perf_loc[0], well_index=self.well_index)

        # Add production well:
        self.add_well("P1", 0.5)
        if self.bound_cond == 'const_pres_rate':
            # Perforate all boundary cells:
            for nth_perf in range(len(self.right_boundary_cells)):
                well_index = self.mesh.volume[self.right_boundary_cells[nth_perf]] / self.max_well_vol * self.well_index
                self.add_perforation(self.wells[-1], res_block=self.right_boundary_cells[nth_perf],
                                     well_index=well_index)

        elif self.bound_cond == 'wells_in_frac':
            # Only perforating the single fracture/matrix block
            self.add_perforation(self.wells[-1], res_block=self.well_perf_loc[1], well_index=self.well_index)

        # Add wells to the DARTS mesh object and sort connection (DARTS related):
        self.mesh.add_wells(ms_well_vector(self.wells))
        self.mesh.reverse_and_sort()
        self.mesh.init_grav_coef()
        return 0
