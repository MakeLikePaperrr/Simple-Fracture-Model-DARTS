# Section of the Python code where we import all dependencies on third party Python modules/libaries or our own
# libraries (exposed C++ code to Python, i.e. darts.engines && darts.physics)
from darts.engines import value_vector, sim_params
from darts.models.physics.dead_oil import DeadOil
from darts.models.physics.geothermal import Geothermal
from darts.models.darts_model import DartsModel
from reservoir import UnstructReservoir
import os


# Here the Model class is defined (child-class from DartsModel) in which most of the data and properties for the
# simulation are defined, e.g. for the reservoir/physics/sim_parameters/etc.
class Model(DartsModel):
    def __init__(self, filename_mesh, physics_type, bound_cond, const_perm, frac_aper, poro,
                 inj_well_coords, prod_well_coords, n_points=64):
        """
        Class constructor of Model class
        :param n_points: number of discretization points for the parameter space
        """
        # Call base class constructor (see darts/models/darts_model.py for more info as well as the links in main.py
        # on OOP)
        super().__init__()

        # Measure time spend on reading/initialization
        self.timer.node["initialization"].start()
        self.bound_cond = bound_cond
        self.physics_type = physics_type

        # Some permeability input data for the simulation
        permx = const_perm  # Matrix permeability in the x-direction [mD]
        permy = const_perm  # Matrix permeability in the y-direction [mD]
        permz = const_perm  # Matrix permeability in the z-direction [mD]
        mesh_file = filename_mesh

        # Instance of unstructured reservoir class from reservoir.py file.
        # When calling this class constructor, the def __init__(self, arg**) is executed which created the instance of
        # the class and constructs the object. In the process, the mesh is loaded, mesh information is calculated and
        # the discretization is executed. Besides that, also the boundary conditions of the simulations are
        # defined in this class --> in this case constant pressure/rate at the left (x==x_min) and right (x==x_max) side
        self.reservoir = UnstructReservoir(permx=permx, permy=permy, permz=permz, frac_aper=frac_aper,
                                           mesh_file=mesh_file, poro=poro, bound_cond=self.bound_cond,
                                           physics_type=self.physics_type, inj_well_coords=inj_well_coords,
                                           prod_well_coords=prod_well_coords)

        if self.physics_type == 'dead_oil':
            # Create physics using DeadOil constructor (imported from darts.models.physics.dead_oil)
            self.cell_property = ['pressure', 'saturation']
            self.physics = DeadOil(timer=self.timer, physics_filename='physics.in',
                                   n_points=n_points, min_p=0, max_p=800, min_z=1e-12)

        elif self.physics_type == 'geothermal':
            self.cell_property = ['pressure', 'enthalpy']
            self.physics = Geothermal(timer=self.timer, n_points=n_points, min_p=0, max_p=500, min_e=0, max_e=50000)

            # Fill some additional values for geothermal runs:
            self.reservoir.hcap.fill(2200)
            self.reservoir.conduction.fill(181.44)

        else:
            print("--------ERROR SPECIFY CORRECT PHYSICS NAME--------")

        # Some tuning parameters:
        self.params.first_ts = 1e-6  # Size of the first time-step [days]
        self.params.mult_ts = 8  # Time-step multiplier if newton is converged (i.e. dt_new = dt_old * mult_ts)
        self.params.max_ts = 20  # Max size of the time-step [days]
        self.params.tolerance_newton = 1e-4  # Tolerance of newton residual norm ||residual||<tol_newt
        self.params.tolerance_linear = 1e-6  # Tolerance for linear solver ||Ax - b||<tol_linslv
        self.params.newton_type = sim_params.newton_local_chop  # Type of newton method (related to chopping strategy?)
        self.params.newton_params = value_vector([0.2])  # Probably chop-criteria(?)

        self.runtime = 2000  # Total simulations time [days], this parameters is overwritten in main.py!

        # End timer for model initialization:
        self.timer.node["initialization"].stop()

    def set_initial_conditions(self):
        """
        Class method called in the init() class method of parents class
        :return:
        """
        if self.physics_type == 'dead_oil':
            # Takes care of uniform initialization of pressure and saturation (composition) in this example:
            self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=350,
                                                        uniform_composition=[0.2357])
        elif self.physics_type == 'geothermal':
            # Takes care of uniform initialization of pressure and temperature in this example:
            self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=350,
                                                        uniform_temperature=348.15)
        return 0

    def set_boundary_conditions(self):
        """
        Class method called in the init() class method of parents class
        :return:
        """
        # Takes care of well controls, argument of the function is (in case of bhp) the bhp pressure and (in case of
        # rate) water/oil rate:
        for i, w in enumerate(self.reservoir.wells):
            if 'I' in w.name or 'INJ' in w.name:
                # Add controls for injection well:
                if self.physics_type == 'dead_oil':
                    # For BHP control in injection well we usually specify pressure and composition (upstream) but here
                    # the method is wrapped such  that we only need to specify bhp pressure (see lambda for more info)
                    w.control = self.physics.new_bhp_water_inj(375)

                elif self.physics_type == 'geothermal':
                    # Specify both pressure and temperature (since it's upstream for injection well)
                    w.control = self.physics.new_bhp_water_inj(375, 308.15)
                    # w.control = self.physics.new_rate_water_inj(4800, 308.15)

            else:
                # Add controls for production well:
                if self.physics_type == 'dead_oil':
                    # Specify bhp for particular production well:
                    w.control = self.physics.new_bhp_prod(325)

                elif self.physics_type == 'geothermal':
                    # Specify bhp for particular production well:
                    w.control = self.physics.new_bhp_prod(325)
                    # w.control = self.physics.new_rate_water_prod(4800)

        return 0
