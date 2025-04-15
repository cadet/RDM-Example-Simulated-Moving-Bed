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
# `w_e` defines the volume flow percentile with which the fluid leaves a particular `zone` and enters another.

# %% [markdown]
# As seen in the  Four-zone SMB schemes with eight (right) columns indicating positions of the associated hold-up volumes (Fig.5, He et al.) Hold up-volumes between all columns and external valves exist and should may be evaluated in thier effect on the SMB elution.
#
# Linear Isotherm:     {\displaystyle q=K_{\text{H}}\cdot C_{\text{eq}}}
#
#     q – Beladung des Sorbents (Masse Sorbat bezogen auf Masse Sorbent)
#     KH – Henry-Koeffizient
#     Ceq – Konzentration des Sorbats in Lösung

# %% [markdown]
# ```{figure} ./figures/case_study1_practical_setup.jpg
# ## Four-zone SMB schemes with eight (right) columns indicating positions of the associated hold-up volumes (Fig.5, He et al.)](src/figures/case_study1_practical_setup.jpg)

# %% [markdown]
# `eluent`= desorbant 

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, LumpedRateModelWithPores

# Component System
component_system = ComponentSystem(['A', 'B'])

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.54, 0.28]  # Henry_1 = 	0.54; Henry_2 = 0.28
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
#column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
column.axial_dispersion = 3.81e-6  # D_ax [m² / s]

eluent = Inlet(component_system, name='eluent')
eluent.c = [0, 0]  # c_in_D [mol / m^3]
eluent.flow_rate = 4.14e-8  # Q_D [m^3 / s]

feed = Inlet(component_system, name='feed')
feed.c = [2.78e3, 2.78e3]  # c_in [mol / m^3]
feed.flow_rate = 2.0e-8  # Q_F [m^3 / s]

# %% [markdown]
# Henry coefficient can be assumed to equal the equilibrium constant under ideal, linear conditions. The percentile of the volume flow that leaves zone I for [extract, zone II] can be deducted by examining the differences in the volumetric flow rate as all columns have the same physical properties. 
# -> A_1 * v_1 = A_2 * v_2
# von zone I in Extract (A? v = Q_E) und zone II (A = column.cross..., v = Q_II)
# ```
# 1.4e-7 * 5.31e-4 = (3.48e-8 * A_extract + 1.05e-7 * 5.31e-4 ) 
# 1.4e-7 = (3.48e-8 * (A_extract / A) + 1.05e-7)
# A_extract = (1.4e-7  - 1.05e-7 ) / (3.48e-8 / A) = 5.34e-4
# 1.4e-7 = (3.48e-8 * 1.006 + 1.05e-7)
# 100 * (3.48e-8 * 1.006) / (1.4e-7) = 25% für Extrakt
#
# zoneIII = Raff + zone IV
# 1.25e-7 * 5.31e-4 = (2.66e-8 * A_raff + 9.81e-8 * 5.31e-4 ) 
# 1.25e-7 = (2.66e-8 * (A_raff / A) + 9.81e-8)
# A_raff = (1.25e-7  -  9.81e-8) / (2.66e-8 / 5.31e-4) = 5.37e-4
# 1.25e-7 = (2.66e-8 * 1.011 +  9.81e-8)
# 100 * (2.66e-8 * 1.011) / (9.81e-8) = 27.4% für Raffinate
# ```
# -> kann einfach Q vergleichen, der Rest der Gleichung wird durch die andere Fläche des Extrakts gestellt 
# interstitial velocities berechnet CADET selbst aus Q / A * porosity

# %%
(1.25e-7  -  9.81e-8) / (2.66e-8 / 5.31e-4) / 5.31e-4

100 * (2.66e-8 * 1.011) / (9.81e-8)


# %%
from CADETProcess.modelBuilder import SMBBuilder

smb_builder = SMBBuilder(feed, eluent, column)
smb_builder.switch_time = 1552

#builder.add_connection(zone_I, extract)
#builder.add_connection(zone_I, zone_II)
w_e = 0.25
smb_builder.set_output_state('zone_I', [w_e, 1-w_e])

#builder.add_connection(zone_III, raffinate)
#builder.add_connection(zone_III, zone_IV)
w_r = 0.274
smb_builder.set_output_state('zone_III', [w_r, 1-w_r])
process = smb_builder.build_process()

# %%
process = smb_builder.build_process()

from CADETProcess.simulator import Cadet
process_simulator = Cadet()
process_simulator.n_cycles = 1

simulation_results = process_simulator.simulate(process)

_ = simulation_results.solution.raffinate.inlet.plot()
_ = simulation_results.solution.extract.inlet.plot()

# %%
