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
# The first of the five case studies depicted in the paper evaluates the separation of the two components fructose `A` and glucose `B` in a four-zone simulated moving bed (SMB) system with eight columns. The binding behavior follows a linear isotherm. 
# This experiment was originally simulated and published by Klatt et al. in "Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose" (Karsten-Ulrich Klatt, Felix Hanisch, Guido Dünnebier, Journal of Process Control 2002; 12(2):203-219. https://doi.org/10.1016/S0959-1524(01)00005-1.) <br>
# https://www.sciencedirect.com/science/article/abs/pii/S0959152401000051
#
#
# Continuous processes generally have many benefits over batch processes, like higher efficiency and throughput. However, a truly continuous chromatography process is not practically feasable. This so called true moving bed chromatography would entail the solid phase of the chromatography column (the bed) moving in the opposite direction of the mobile phase. This would induce the local separation of components in the feed solution based on their column binding behaviour. Those components could then by retrieved by different outlet streams located upstream and downstream of the feed inlet. 
#
# ```{figure} ./figures/true_MB.png
# :width: 400px
# <div style="text-align: center">
#
#  [Link to Youtube Video](https://www.youtube.com/watch?v=xhhJxb48tgc)
# <div>
#
#

# %% [markdown]
# Simulated Moving Bed Chromatography is a way to approach such a continuous process in practice. This is realized by having multiple chromatography columns connected to each other in a caroussel. By periodically switching the inlet and outlet valves connected to the columns ("column switching") the movement of the solid phase can be mimicked. 
#
# A general SMB system contains four different external units connected to the columns:
# 1. Feed (Inlet): Component mixture
# 2. Raffinate (Outlet): Faster eluting component (elutes before feed plug flow)
# 3. Extract (Outlet): Slower eluting component (elutes after feed plug flow), interacts more strongly with the column solid phase<br>
# There can be multiple extract or raffinate outlets in a SMB system depending on the number of components to be separated.
# 4. Desorbant / Eluent / Solvent (Inlet): Solution to elute extract from column before raffinate plug flow enters the column again to prevent mixing of the separated components
#
# The position of each column relative to these external units determine their specific “Zone” in the SMB system. Based on thier external units, every zone has a different function in the seperation process and a different flow rate within the columns. The most basic SMB system for binary separation is made up of four zones with one chromatography column in each zone. 
#

# %% [markdown]
# ```{figure} ./figures/four_zone.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 1, He et al.) Schematic of four-zone SMB chromatography for binary separations. Column positions are periodically switched in opposite direction of liquid flow.

# %% [markdown]
# ```{figure} ./figures/case_study1_practical_setup.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 5, He et al.) Four-zone SMB schemes with eight columns indicating positions of the associated hold-up volumes
# <div>

# %% [markdown]
# As seen in Fig. 5, hold up-volumes between all columns and external units exist and should ideally be considered in thier effect on the SMB elution. (triangle theory) They generally increase retention time and dispersion. (The hold-up volume on either side of a column, i.e., tubing and frits, can be described by a CSTR that is moved through the network together with that column. )(four categories: (1) tubing between multi-port valve and column inlet plus frit before packed bed, (2) frit after packed bed plus tubing between column outlet and multi-port valve, (3) tubing between injection point and multi-port valve, and (4) tubing between multi-port valve and detector. Each of these categories can be modeled as one or more PFR, CSTR and/or DPFR in series. )
# -> CADET-SMB allows to consider hold-up volumes in the column network. This is demonstrated by introducing CSTR models, Eq. (10), in case study I as illustrated by Fig. 5. The residence time, τ
# CSTR, is varied between 0s, 5s and 10s. Fig. 14 shows the impact of these hold-up volumes on the column states in CSS.

# %% [markdown]
# To simulate a SMB process, first the physical properties of the columns and the Inlet are defined. The mass transfer within the column is characterized by the equilibrium-dispersive model (EDM) which can be derived from the `GeneralRateModel` by defining the spatial discretization `column.discretization.npar` as 1. In the finite volume method, only one radial cell is assumed. The axial column dimension `column.discretization.ncol` is set to 40 axial cells. All numerical values are taken from Table 1.(4. Case Studies). The Henry coefficient can be assumed to equal the equilibrium constant under ideal, linear conditions. 

# %% [markdown]
# # Differences in He's Matlab code / He's paper / original case study in [Klatt's paper](https://www.sciencedirect.com/science/article/pii/S0959152401000051?ref=pdf_download&fr=RR-2&rr=94b07706292368ec#TBL1):
#
# ## Parameters from Matlab code (getParameters_binary_case2.m) not in CarouselBuilder Example:
#
#         % The parameter setting for simulator
#         opt.tolIter         = 1e-4;
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
#         opt.nThreads        = 4;(/8)
#         opt.timePoints      = 1000         !!
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
#     3.) Recycle flow rate nomenclature:
#     Matlab code: Recycle flow rate QI: flowRate.recycle    = 0.1395e-6;      % m^3/s 
#     Klatt/He paper: Recycle flow rate QIV = 0.0981 cm 3 /s = 9.81e-8 m^3/s 
#     
#     4.) Order of Henry coefficients exchanged:
#     # Matlab code: opt.KA = [0.28 0.54]; % [comp_A, comp_B], A for raffinate, B for extract
#     opt.comp_raf_ID = 1;  % the target component withdrawn from the raffinate ports
#     opt.comp_ext_ID = 2;  % the target component withdrawn from the extract ports    
