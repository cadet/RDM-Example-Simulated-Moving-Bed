# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# (Example_3_SMB)=
# # Ternary Separation

# %% [markdown]
# Case study III examines the separation of the three components 2’-deoxycytidine `A`, 2’-deoxyguanosine `B` and 2’-deoxyadenosine `C`. A standard five-zone SMB system with two extract ports and without partial-feeding or partial-closing is assumed.
#
# This experiment was originally simulated and published by Mun et al. in "Improving performance of a five-zone simulated moving bed chromatography for ternary separation by simultaneous use of partial-feeding and partial-closing of the product port in charge of collecting the intermediate-affinity solute molecules" (Sungyong Mun, Journal of Chromatography A 2011;  1218(44):8060-8074) <br> https://doi.org/10.1016/j.chroma.2011.09.015.
#

# %% [markdown]
# ```{figure} ./figures/case_study3.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 4, He et al.) Integrated five-zone SMB scheme with two extract ports
# <div>

# %% [markdown]
# (Ternary separation with a five zone system is particularly effective when the target component is present in much higher abundance than the more strongly retained component.)

# %% [markdown]
# ## Differences in Parameters He's paper / He's Matlab code / Mun's paper
#    
#     1.) Feed concentrations: 
#     He paper text and He table: cF = 1.0 g/cm^3 = 1.0 g/mL
#     original Mun paper: Mun Table 1 cF = 1.0 g/L => 1.0e3 g/cm^3
#
#     Matlab code:     
#     concentrationFeed 	= [1.0, 1.0, 1.0];    % g/cm^3 [concentration_compA, concentration_compB]  => wrong unit 
#     opt.molMass         = [227.217, 267.24, 251.24192]; % The molar mass of each components g/mol
#     
#     He's paper: feed.c = [4.41e3, 3.75e3, 3.98e3]  # c_in [mol / m^3]  
#     Muns paper c + Matlab code M = [4.401079144606257, 3.74195479718605, 3.9802275034357324]  [mol / m^3]
#
#     2.) Particle Porosity
#     He Paper: column.particle_porosity = 1.0e-5  # ε_p [-] 
#     Matlab code: opt.porosityParticle    = 0.00000001;   % e_p very small to ensure e_t = e_c  => would be 1.0e-8
#
#     3.) Film Diffusion:
#     He Paper: column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]  
#     Matlab code: opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];  % K_f
#
#     4.) Pore Diffusion: 
#     He Paper: column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
#     Matlab code: opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];  % D_p
#
# # volumetric flow rates, time checked
#
#     % Purities archieved at operation point a for components A,B,C: 99.62  71.16	96.78
#     flowRate.recycle    = 2.9230e-7;      % m^3/s  == QI, not QIV or QV!
#
# # Parameters missing
#     opt.tolIter         = 1e-4;  % tolerance of the SMB stopping criterion
#     opt.nMaxIter        = 1000;  % the maximum iteration step in SMB
#     opt.nThreads        = 8;     % threads of CPU, up to your computer
#
#     opt.timePoints      = 1000;  % the observed time-points
#     -> For Feed concentration setup
#     Feed.time = linspace(0, opt.switch, opt.timePoints);  # 1000 evaluated points during 1 switching time
#     Feed.concentration = zeros(length(Feed.time), opt.nComponents);
#
#     
#     opt.Purity_extract1_limit   = 0.95;  % used for constructing constraints
#     opt.Purity_extract2_limit   = 0.65;  % used for constructing constraints
#     opt.Purity_raffinate_limit  = 0.99;  % used for constructing constraints
#     opt.Penalty_factor          = 10;    % penalty factor in penalty function
#
#     opt.yLim            = max(concentrationFeed ./ opt.molMass) * 1.1; % the magnitude for plotting
#
#
#  #   Capable of placing a CSTR or DPFR apparatues before and after the calculated column
#
#     Continuous Stirred Tank Reactor
#     opt.enable_CSTR = false;
#     opt.CSTR_length = 0.01;
#
#     Dispersive Plug Flow Reactor
#     opt.enable_DPFR = false;
#
#     opt.DPFR_length = 0.0066;
#     opt.DPFR_nCells = 50;
#
#     opt.DPFR_velocity   = 0.00315;
#     opt.DPFR_dispersion = 2.5e-20;
#

# %% [markdown]
# The following process parameters are taken from Table 1 (Mun et al.). The feed concentrations `feed.c` (mol/m^3) of all three components are calculated based on 1.0 g/L. The axial dispersion is assumed to be a typical value (Table 1, He et al.). The mass transfer follows a linear binding model, in which the `adsorption_rate` is given as the product of the `mass-transfer coefficient` and the `Henry constant` for each component. The `desorption_rate` is equal to the mass-transfer coefficient of each component (Equation 9b, 9j, Mun et al.). 
#

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, LumpedRateModelWithoutPores
import numpy as np

# Component System
component_system = ComponentSystem(['A', 'B', 'C'])

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = True
binding_model.adsorption_rate = [3.15, 7.40*0.5, 23.0*0.1]  # k_a = apkm * H [1/s]
binding_model.desorption_rate = [1, 0.5, 0.1]  # kd = apkm [1/s]


# Column
#class CADETProcess.processModel.(total_porosity, _q, length, diameter, axial_dispersion,
 #                                                           flow_direction, c, name)
column = LumpedRateModelWithoutPores(component_system, name='column')
column.binding_model = binding_model
column.length = 0.150 # L [m]
column.diameter = 1.0e-2  # d [m]
column.total_porosity = 0.80  # ε_c [-]
column.axial_dispersion = 3.81e-5  # D_ax [m² / s  #Matlab code: 3.8148e-10;
#column.discretization.npar = 1  # N_r
column.discretization.ncol = 40  # N_z
column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent')
eluent.c = [0, 0, 0]  # c_in_D [mol / m^3]
eluent.flow_rate = 2.34e-7  # Q_D [m^3 / s]  Matlab code: flowRate.desorbent  = 2.3412e-7;      % m^3/s



feed = Inlet(component_system, name='feed')
#feed.c = [4.41e3, 3.75e3, 3.98e3]  # c_in [mol / m^3]  
feed.c = [4.40, 3.74, 3.98]  # c_in [mol / m^3] 
# Muns paper c + Matlab code M = [4.401079144606257, 3.74195479718605, 3.9802275034357324]
feed.flow_rate = 1.67e-8  # Q_F [m^3 / s]  Matlab code: flowRate.feed       = 1.6667e-8;      % m^3/s

# %% [markdown]
# All zones are connected to each other in series `SerialZone`. The flow rates are taken from Table 2 (Mun at al.) and the fractions of the flow going into extract port I `w_e1`, extract port 2 `w_e2` and the raffinate port `w_r` are calculated. 
# ```
# zone_I -> extract_1 + zone_II
# Q_I = Q_E1 + Q_II 
#
# zone_II -> extract_2 + zone_III
# Q_II = Q_E2 + Q_III 
#
# zoneIV -> raffinate + zone_V
# Q_IV = Q_R + Q_V 
#
# w_e1 = Q_E1 / Q_I = 0.6418
# w_e2 = Q_E2 / Q_II = 0.443
# w_r = Q_R / Q_IV = 0.224

# %%
Q_R = 1.009  # Operating point a, Table 3 Mun et al.
Q_E1 = 11.256  # Operating point a, Table 3
Q_E2 = 2.782  
Q_I = 17.538
Q_II = 6.282
Q_IV = 4.500

Q_E1 / Q_I 
Q_E2 / Q_II 
Q_R / Q_IV 

# %%
extract_1 = Outlet(component_system, name = 'extract_1')
extract_2 = Outlet(component_system, name = 'extract_2')
raffinate = Outlet(component_system, name = 'raffinate')
from CADETProcess.modelBuilder import SerialZone

zone_I = SerialZone(component_system, 'zone_I', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_II = SerialZone(component_system, 'zone_II', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_III = SerialZone(component_system, 'zone_III', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_IV = SerialZone(component_system, 'zone_IV', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_V = SerialZone(component_system, 'zone_V', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})

from CADETProcess.modelBuilder import CarouselBuilder

builder = CarouselBuilder(component_system, 'smb')
builder.valve_dead_volume = 1e-9
builder.column = column
builder.add_unit(eluent)
builder.add_unit(feed)
builder.add_unit(extract_1)
builder.add_unit(extract_2)
builder.add_unit(raffinate)

builder.add_unit(zone_I)
builder.add_unit(zone_II)
builder.add_unit(zone_III)
builder.add_unit(zone_IV)
builder.add_unit(zone_V)

builder.add_connection(eluent, zone_I)

builder.add_connection(zone_I, extract_1)
builder.add_connection(zone_I, zone_II)
w_e1 = 0.642
builder.set_output_state(zone_I, [w_e1, 1 - w_e1])

builder.add_connection(zone_II, extract_2)
builder.add_connection(zone_II, zone_III)
w_e2 = 0.443
builder.set_output_state(zone_II, [w_e2, 1 - w_e2])

builder.add_connection(zone_III, zone_IV)

builder.add_connection(feed, zone_IV)
builder.add_connection(zone_IV, raffinate)
builder.add_connection(zone_IV, zone_V)
w_r = 0.224
builder.set_output_state(zone_IV, [w_r, 1 - w_r])

builder.add_connection(zone_V, zone_I)

builder.switch_time = 264  

process = builder.build_process()

# %%
from CADETProcess.simulator import Cadet
process_simulator = Cadet()
#process_simulator.evaluate_stationarity = True
process_simulator.n_cycles = 41  #200.99 steps = 200.99 switch times -> 40.198
process_simulator.use_dll = True

process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code!, not in Klatt paper 
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

#process_simulator.time_resolution = builder.switch_time / 1000  # default value is 1 second.  Matlab code: Feed.time = linspace(0, opt.switch, opt.timePoints);

simulation_results = process_simulator.simulate(process)

# %%
import numpy as np
raff_mM = simulation_results.solution.raffinate.inlet.solution
ext1_mM = simulation_results.solution.extract_1.inlet.solution
ext2_mM = simulation_results.solution.extract_2.inlet.solution
t = simulation_results.time_complete

# Transformation from mM to g/L
molar_mass = [227.217, 267.24, 251.24192]  # molar mass of each component
raff = np.multiply(raff_mM, molar_mass) * 1e-3
ext_1 = np.multiply(ext1_mM, molar_mass) * 1e-3
ext_2 = np.multiply(ext2_mM, molar_mass) * 1e-3

# %% [markdown]
# To compare the simulation results to those of Mun et al., the concentration of every component is averaged over one switching period. This results in a new average every 264s. As there are 5 columns that switch a total of 41 times, this results in 205 total average concentrations for every component. Dividing the total simulation time by the switch time yields the same number of 205 steps for `n_averages`. This is done for the raffinate, extract 1 and extract 2 ports.
#
#     np.shape() of concentration arrays during reshaping:
#     (54121,3)
#     (205,264,3)
#     (205,264,3)

# %%
n_averages = int(len(raff) // builder.switch_time)  

raff_average = (
    raff[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)
ext1_average = (
    ext_1[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)
ext2_average = (
    ext_2[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)

# %%
import matplotlib.pyplot as plt

n_steps = range(0,n_averages)
#n_steps = int(len(t)/builder.switch_time)
#n_steps = n_columns * n_cycles

fig, axs = plt.subplots(2, 2, figsize=(20, 17))
ax1 = axs[0, 0]  # Raff
ax2 = axs[0, 1]  # Ext2
ax3 = axs[1, 0]  # Ext1

ax1.plot(n_steps, raff_average)
ax1.set_title("Raffinate")
ax1.set_xlabel("Step number")
ax1.set_ylabel("c [g / L]")
ax1.set_ylim(0, 0.8)

ax2.plot(n_steps, ext2_average)
ax2.set_title("Extract 2")
ax2.set_xlabel("Step number")
ax2.set_ylabel("c [g / L]")
ax2.set_ylim(0, 0.8)

ax3.plot(n_steps, ext1_average)
ax3.set_title("Extract 1")
ax3.set_xlabel("Step number")
ax3.set_ylabel("c [g / L]")
ax3.set_ylim(0, 0.8)


# %%
# Axial concentrations in mM
from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
axial_conc = CarouselSolutionBulk(builder, simulation_results)
before_switch = axial_conc.plot_at_time(t=200.99 * builder.switch_time)
after_switch = axial_conc.plot_at_time(t=200.01 * builder.switch_time)

# Conversion of axial concentration plot from mM to g/L
t = 200.01 * builder.switch_time
n_cols = axial_conc.builder.n_columns

fig, axs = plt.subplots(
ncols=n_cols,
figsize=(n_cols*4, 6),
gridspec_kw=dict(wspace=0.0, hspace=0.0),
sharey='row')
t_i = np.where(t <= axial_conc.time)[0][0]

x = axial_conc.axial_coordinates

y_min_data = 0
y_max_data = 0
zone_counter = 0
column_counter = 0
_lines = []

for position, ax in enumerate(axs):
    col_index = axial_conc.builder.column_indices_at_time(t, position) 
    y_data = axial_conc.solution[f'column_{col_index}'].bulk.solution[t_i, :]
    y = np.multiply(y_data, molar_mass) * 1e-3
    y_min_data = min(y_min_data, min(0, np.min(y)))
    y_max_data = max(y_max_data, 1.1*np.max(y))

    l = ax.plot(x, y)
  
    _lines.append(l)

    zone = axial_conc.builder.zones[zone_counter]
    if zone.n_columns > 1:
        ax.set_title(f'{zone.name}, position {column_counter}')
    else:
        ax.set_title(f'{zone.name}')

    if column_counter < (zone.n_columns - 1):
        column_counter += 1
    else:
        zone_counter += 1
        column_counter = 0
ax.set_ylabel("c(g/L)")        
np.shape(y)

# %% [markdown]
# ```{figure} ./figures/ternary.png
# :width: 800px
# <div style="text-align: center">
# (Fig. 8, Mun et al.) internal concentration profiles of the five-zone SMBs, (a) Standard (at 200.01 steps), (b) Standard (at 200.99 steps), Blue line: component A, red line: component B, green line: component C.; numerical simulations where all the mass-transfer effects were considered in accordance with the information in Table 1. 
# <div>

# %% [markdown]
# ```{figure} ./figures/ternary_separation_Mun.png
# :width: 800px
# <div style="text-align: center">
# (Fig. 8, Mun et al.) internal concentration profiles of the five-zone SMBs, Operation mode: Standard mode (at 200.99 steps), numerical simulations where the mass-transfer effects were minimized to approach an equilibrium (or ideal) state.
# <div>

# %%
