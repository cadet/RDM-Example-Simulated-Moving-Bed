# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# (Example_SMB)=
# # Simulated Moving Bed

# %% [markdown]
# The following example is a reproduction of part of the research results published in "Efficient numerical simulation of simulated moving bed chromatography with a single-column solver" (Qiao-Le He, Samuel Leweke, Eric von Lieres, Computers & Chemical Engineering 2018; 111:183-198. doi:10.1016/j.compchemeng.2017.12.022.) <br>
# https://www.sciencedirect.com/science/article/pii/S0098135417304520
#

# %% [markdown]
# The first case study depicted in the paper evaluates the separation of the two components fructose `A` and glucose `B` in a four-zone simulated moving bed (SMB) system with eight columns. The binding behavior follows a linear isotherm.  
#
# In addition to the columns a SMB system contains four different external units connected to the columns:
# 1. Feed: component mixture
# 2. Raffinate: faster eluting component (elutes before feed plug flow)
# 3. Extract: slower eluting component (elutes after feed plug flow), interacts more strongly with the column solid phase
# 4. Desorbant (Eluent): Solution to elute extract from column before raffinate plug flow enters the column again to prevent mixing of the separated components
#
#

# %% [markdown]
# A continuous process can be approximated by periodically switching the inlet and outlet valves of the columns in the system ("column switching"). The position of each column relative to these external units determine their specific "Zone" in the SMB system. Based on thier external units, every zone has a different function in the seperation process and a different flow rate within the columns. 
#

# %% [markdown]
# ```{figure} ./figures/case_study1_practical_setup.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 5, He et al.) Four-zone SMB schemes with eight columns indicating positions of the associated hold-up volumes
# <div>

# %% [markdown]
# As seen in Fig. 5, hold up-volumes between all columns and external units exist and should ideally be considered in thier effect on the SMB elution. (triangle theory) They generally increase retention time and dispersion. (The hold-up volume on either side of a column, i.e., tubing and frits, can be described by a CSTR that is moved through the network together with that column. )(four categories: (1) tubing between multi-port valve and column inlet plus frit before packed bed, (2) frit after packed bed plus tubing between column outlet and multi-port valve, (3) tubing between injection point and multi-port valve, and (4) tubing between multi-port valve and detector. Each of these categories can be modeled as one or more PFR, CSTR and/or DPFR in series. )
# -> CADET-SMB allows to consider hold-up volumes in the column network. This is demonstrated by introducing CSTR models, Eq. (10), in case study I as illustrated by Fig. 5. The residence time, τCSTR, is varied between 0s, 5s and 10s. Fig. 14 shows the impact of these hold-up volumes on the column states in CSS.

# %% [markdown]
# To simulate the SMB process, first the physical properties of the columns and the Inlet are defined. The mass transfer within the column is characterized by the equilibrium-dispersive model (EDM) which can be derived from the `GeneralRateModel` by defining the spatial discretization `column.discretization.npar` as 1. In the finite volume method, only one radial cell is assumed. The axial column dimension `column.discretization.ncol` is set to 40 axial cells. All numerical values are taken from Table 1.(4. Case Studies). The Henry coefficient can be assumed to equal the equilibrium constant under ideal, linear conditions. 

# %% [markdown]
# # Differences in He's Matlab code / He's paper / original case study in [Klatt's paper](https://www.sciencedirect.com/science/article/pii/S0959152401000051?ref=pdf_download&fr=RR-2&rr=94b07706292368ec#TBL1):
#
# ## Parameters from Matlab code not in CarouselBuilder Example:
#
#         % The parameter setting for simulator
#         opt.tolIter         = 1e-3;
#
#         % The parameter settting for the SMB
#         opt.Purity_limit    = [0.99, 0.99];
#         opt.Penalty_factor  = 10;
#
#         % opt.compTargID = [2, 1];
#         opt.structID    = [2 2 2 2];
#         opt.diffusionParticleSurface  = [0.0 0.0];  # probably automatically 0 if not defined
#         opt.enableDebug = true;
#
#         % The parameter setting for simulator
#         opt.nMaxIter        = 1000;
#         opt.nThreads        = 4;
#         opt.timePoints      = 1000;
#
#         opt.yLim            = max(concentrationFeed ./ opt.molMass);
#         
#         (Viscosity in Klatt paper, not in He paper)
#
# ## Interstitial velocities missing (Calculation already checked for all flow rates Q): 
#         Interstitial velocities calculated explicitly in Matlab code, done automatically by CarouselBuilder?
#         
#         % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
#         interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
#         interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
#         interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
#         interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
#         interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
#
#
#         process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code, not in Klatt paper, but in He paper
#         

# %% [markdown]
# ##  Matlab code capable of placing a CSTR or DPFR apparatues before and after the calculated column:
#
# %   Continuous Stirred Tank Reactor
#     opt.enable_CSTR = false;
#     opt.CSTR_length = 0.01;
#
# %   Dispersive Plug Flow Reactor
#     opt.enable_DPFR = false;
#
#     opt.DPFR_length = 0.0066;
#     opt.DPFR_nCells = 50;
#
#     opt.DPFR_velocity   = 0.00315;
#     opt.DPFR_dispersion = 2.5e-20;
#
#
# ## Differences in parameters:
#     
#     1.)  He paper: "The inlet concentrations are converted from 0.55 g/cm^3 assuming that fructose and glucose have the same molar mass of 180 g/mol." 
#     => Would be: (0.55/180)*1e6 = 3055.555 = 3.06e3 mol/m^3
#     
#     Text does not match table in paper
#     
#     # He paper Table 1: 2.78e3 
#     # original Klatt paper:  "cF = 0.5 g/cm3" => would also be 2.78e3
#
#     Matlab code: 
#         concentrationFeed 	= [0.5, 0.5];   % g/m^3 [concentration_compA, concentration_compB]   => wrong unit, is actually g/cm^3
#         opt.molMass         = [180.16, 180.16];
#         opt.yLim            = max(concentrationFeed ./ opt.molMass);
#
#         Feed.time = linspace(0, opt.switch, opt.timePoints);
#         Feed.concentration = zeros(length(Feed.time), opt.nComponents);  
#     
#     2.) Film diffusion and particle diffusion exchanged in He's paper and Matlab code
#     Matlab code: opt.filmDiffusion             = [5e-5 5e-5];
#     Matlab code: opt.diffusionParticle         = [1.6e4 1.6e4];
#
#
#
#
#     % Example in Matlab code is not the exact same example as Klatt et al.: different values for one flow rate and component adsorption, but otherwise every value taken from Klatt et al. 
#     for Recycle(= Zone IV) flow rate:
#     flowRate.recycle    = 0.1395e-6;      % m^3/s = MATLAB code; 
#     Klatt/He paper = Recycle flow rate QIV = 0.0981 cm 3 /s = 9.81e-8 m^3/s 
#     
#     Different value for Henry coefficients:
#     He, Klatt paper: Henry_1 = 	0.54; Henry_2 = 0.28
#     Matlab code: opt.KA = [0.28, 0.61]; % [comp_A, comp_B], A for raffinate, B for extract

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, GeneralRateModel

# Component System
component_system = ComponentSystem(['A', 'B'])

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.54, 0.28]  # Henry_1 = 	0.54; Henry_2 = 0.28
# Matlab code: opt.KA = [0.28, 0.61]; % [comp_A, comp_B], A for raffinate, B for extract
binding_model.desorption_rate = [1, 1]
 

# Column
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model

column.length = 0.536  # L [m]
column.diameter = 2.6e-2  # d [m]
column.bed_porosity = 0.38  # ε_c [-] 
# Matlab: porosityColumn

column.particle_porosity = 1.0e-5  # ε_p [-] 
column.particle_radius = 1.63e-3  # r_p [m]
#Matlab code:  opt.particleRadius      = 0.325e-2 /2 = 0,001625; 
#column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]
column.film_diffusion = component_system.n_comp * [5e-5]
#Matlab code: opt.filmDiffusion             = [5e-5 5e-5];

#column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
column.pore_diffusion = component_system.n_comp * [1.6e4]
#Matlab code: opt.diffusionParticle         = [1.6e4 1.6e4];

column.axial_dispersion = 3.81e-6  # D_ax [m² / s]
#Matlab code: opt.dispersionColumn          = ones(1,opt.nZone) .* 3.8148e-6;
column.discretization.npar = 1  # N_r
column.discretization.ncol = 40  # N_z

column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent') #Name in paper = "desorbent"
eluent.c = [0, 0]  # c_in_D [mol / m^3]
#Matlab code: Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);
eluent.flow_rate = 4.14e-8  # Q_D [m^3 / s] 
#

feed = Inlet(component_system, name='feed')
feed.c = [2.78e3, 2.78e3]  # c_in [mol / m^3] => He Matlab, Klatt, NOT HE PAPER
#He paper: "The inlet concentrations are converted from 0.55 g/m^3 assuming that fructose and glucose have the same molar mass of 180 g/mol." => [3052.84, 3052.84])

#Matlab code: Feed.time = linspace(0, opt.switch, opt.timePoints)
#Matlab code: Feed.concentration = zeros(length(Feed.time), opt.nComponents)
feed.flow_rate = 2.0e-8  # Q_F [m^3 / s]

# %%
import numpy as np
concentrationFeed 	= np.array([0.5, 0.5])   # g/cm^3 [concentration_compA, concentration_compB]
opt_molMass         = np.array([180, 180]) #g/mol
concentrationFeed/opt_molMass*1e6

# %% [markdown]
# The unit system of the SMB is implemented using the `CarouselBuilder` from CADETProcess. Four zones with two columns in each zone and the connections to their respective external units are implemented as seen in Fig. 5. For more information please refer: [here](https://cadet-process.readthedocs.io/en/stable/user_guide/tools/carousel_builder.html#). (Not using SMB builder because `n_columns` = 2) <br>
# The percentile of the volume flow that leaves `zone_I` for `extract`(`w_e`) or `zone_II` (`1-w_e`) can be deducted by examining the differences in the volumetric flow rate as all columns have the same cross section `A`(Table 1.). The same can be done for the percentile of the volume flow that leaves `zone_III` for `raffinate`(`w_r`) and `zone_IV` (`1-w_r`) <br>
# The continuity equation for laminar flow is assumed:
# `A_1 * v_1 = A_2 * v_2`
# ```
# zone_I -> extract + zone_II
# w_e = Q_E / Q_I 
# w_e = (3.48e-8) / (1.4e-7) = 0.249  (= 0.24857142857142855)
#
# zoneIII -> raffinate + zone_IV
# w_r = Q_R / Q_III 
# w_r = (2.66e-8 ) / (1.25e-7) = 0.213 (= 0.2128000000000000)
# ``` 

# %%
column.volume_liquid
(3.48e-8) / (1.4e-7)
(2.66e-8 ) / (1.25e-7)
0.325e-2 /2

# %%
extract = Outlet(component_system, name='extract')
raffinate = Outlet(component_system, name='raffinate')
from CADETProcess.modelBuilder import SerialZone, ParallelZone

zone_I = SerialZone(component_system, 'zone_I', n_columns = 2, valve_dead_volume=1e-9)
zone_II = SerialZone(component_system, 'zone_II', n_columns = 2, valve_dead_volume=1e-9)
zone_III = SerialZone(component_system, 'zone_III', n_columns = 2, valve_dead_volume=1e-9)
zone_IV = SerialZone(component_system, 'zone_IV', n_columns = 2, valve_dead_volume=1e-9)

from CADETProcess.modelBuilder import CarouselBuilder

builder = CarouselBuilder(component_system, 'smb')
builder.valve_dead_volume = 1e-9
builder.column = column
builder.add_unit(eluent)
builder.add_unit(feed)

builder.add_unit(extract)
builder.add_unit(raffinate)

builder.add_unit(zone_I)
builder.add_unit(zone_II)
builder.add_unit(zone_III)
builder.add_unit(zone_IV)

builder.add_connection(eluent, zone_I)

builder.add_connection(zone_I, extract)
builder.add_connection(zone_I, zone_II)
w_e = 0.249  
builder.set_output_state(zone_I, [w_e, 1-w_e])

builder.add_connection(zone_II, zone_III)

builder.add_connection(feed, zone_III)

builder.add_connection(zone_III, raffinate)
builder.add_connection(zone_III, zone_IV)
w_r = 0.213
builder.set_output_state(zone_III, [w_r, 1-w_r])

builder.add_connection(zone_IV, zone_I)

builder.switch_time = 1552

process = builder.build_process()

# %% [markdown]
# ((CSS - cyclic steady state: dynamic trajectory is repeated after every switch => all columns have same state after specific time laps ))

# %%
#process = smb_builder.build_process()

from CADETProcess.simulator import Cadet
process_simulator = Cadet()
#process_simulator.evaluate_stationarity = True
process_simulator.n_cycles = 6
process_simulator.use_dll = True
#process_simulator.timeout = 15*60
#simulate first 8 switch times (1 iteration), conc bis 1mol

process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code!, not in Klatt paper 
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

simulation_results = process_simulator.simulate(process)
cycle = 6
_ = simulation_results.solution.raffinate.inlet.plot(start = cycle * 0, end = (cycle) * 8 * builder.switch_time)
_ = simulation_results.solution.extract.inlet.plot(start = 0, end = 8 * builder.switch_time)
_ = simulation_results.solution.extract.inlet.plot(start = cycle * 0, end = (cycle) * 8 * builder.switch_time) 

# %%
raff = simulation_results.solution.raffinate.inlet.solution
ext = simulation_results.solution.extract.inlet.solution
t = simulation_results.time_complete
t

# %%
import matplotlib.pyplot as plt
import numpy as np

plt.plot(t, ext*np.pi*0.536*(2.6e-2/2)**2)
plt.ylim(0,1)

# %%
plt.plot(t, raff*np.pi*0.536*(2.6e-2/2)**2)
plt.ylim(0,1)


# %%
#class CarouselSolutionBulk(SolutionBase): 
from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
axial_conc = CarouselSolutionBulk(builder, simulation_results)
axial_conc.component_system
axial_conc.solution
axial_conc.axial_coordinates
axial_conc.time
simulation_results.solution
axial_conc.plot_at_time(t = 74496.0)

# t = 48*switchtime = 6 cycles 
#for t = switching time -> only 1 column switch, have to at least switch once for every column
#8x switching time = 1 cycle
#needs a few cycles to get to CSS => columns are filled completely 

# looks like graph is shifted by 1 switch time 

# %%
