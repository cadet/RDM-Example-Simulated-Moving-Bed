# %% [markdown]
# (Example_SMB)=
# # Simulated Moving Bed

# %% [markdown]
# The following example is a reproduction of part of the research results published in "Efficient numerical simulation of simulated moving bed chromatography with a single-column solver" (Qiao-Le He, Samuel Leweke, Eric von Lieres, Computers & Chemical Engineering 2018;111:183-198. doi:10.1016/j.compchemeng.2017.12.022.) <br>
# https://www.sciencedirect.com/science/article/pii/S0098135417304520
#

# %% [markdown]
# The first case study depicted in the paper evaluates the separation of fructose `A` and glucose `B` in a four-zone simulated moving bed (SMB) with eight columns. The binding behavior follows a linear isotherm.  
#
#
# SMB:
# contains 4 different external units, connected to the column by valves:
# 1. Feed: mixture A+B
# 2. Raffinate: faster/more easily eluting component (elutes before feed)
# 3. Extract: slower elution component (elutes after feed input), interacts more strongly with the column solid phase
# 4. Desorbant: Solution to eluate Extract from column before Raffinate plug flow enters the column again. -> want to have an empty column to prevent Mixing of the separated components
#
# `w_e` defines the volume flow percentile with 

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, LumpedRateModelWithPores

# Component System
component_system = ComponentSystem(['A', 'B'])

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = True
binding_model.adsorption_rate = [2, 3]
binding_model.desorption_rate = [1, 1]

# Column
column = LumpedRateModelWithPores(component_system, name='column')
column.binding_model = binding_model

column.length = 0.536  # L [m]
column.diameter = 2.6e-2  # d [m]
column.bed_porosity = 0.38  # ε_c [-]

column.particle_porosity = 1.0e-5  # ε_p [-] 
column.particle_radius = 1.63e-3  # r_p [m]
column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]
column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
column.axial_dispersion = 3.81e-6  # D_ax [m² / s]

eluent = Inlet(component_system, name='eluent')
eluent.c = [0, 0from CADETProcess.modelBuilder import SMBBuilder

smb_builder = SMBBuilder(feed, eluent, column)
smb_builder.switch_time = 1552  # t_s [s]

w_e = 
smb_builder.set_output_state('zone_I', [w_e, 1-w_e])

w_r = 
smb_builder.set_output_state('zone_III', [w_r, 1-w_r])


eluent.flow_rate = 

feed = Inlet(component_system, name='feed')
feed.c = [, ]
feed.flow_rate = 

# %%
from CADETProcess.modelBuilder import SMBBuilder

smb_builder = SMBBuilder(feed, eluent, column)
smb_builder.switch_time = 100

w_e = 0.14
smb_builder.set_output_state('zone_I', [w_e, 1-w_e])

w_r = 0.13
smb_builder.set_output_state('zone_III', [w_r, 1-w_r])process = smb_builder.build_process()

# %%
process = smb_builder.build_process()

from CADETProcess.simulator import Cadet
process_simulator = Cadet()
process_simulator.n_cycles = 3

simulation_results = process_simulator.simulate(process)

_ = simulation_results.solution.raffinate.inlet.plot()
_ = simulation_results.solution.extract.inlet.plot()
