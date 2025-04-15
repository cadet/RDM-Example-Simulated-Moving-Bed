# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# (lwe_example_concentration)=
# # Multi-State Steric Mass Action
#

# %% [markdown]
# The following example is a reproduction of part of the research results published in "Multi-state steric mass action model and case study on complex high loading behavior of mAb on ion exchange tentacle resin" (Diedrich J, Heymann W, Leweke S, et al., J Chromatogr A. 2017;1525:60-70. doi:10.1016/j.chroma.2017.09.039). <br>
# https://pubmed.ncbi.nlm.nih.gov/29055527/

# %% [markdown]
# In their study, Diedrich et al. examined the binding behavior of a therapeutic monoclonal antibody (mAb) from a CHO culture on Fractogel EMD SO₃⁻. This tentacle resin is used for cation exchange chromatography (CEX) with a high selectivity. Under four different loading experiments, the elution profile of the protein was examined and found to exhibit a "shoulder under overloaded conditions". This is a result of complex binding behaviour between mAb and tentacle resin. The experimental data can be replicated with CADET-Process using the Multi-State Steric Mass Action model. The mAb is assumed to have two distinct binding states with the stationary phase in the column. Components are able to alternate between the bound states with given conversion rates. They exhibit different binding behavior based on their binding state such as different sorption rates, characteristic charges and steric factors. The protein is treated with a Load-Wash-Elute process with a linear salt gradient for elution from the column.
#
# In the example covered in the paper, the mAb in the mobile phase is able to bind to the stationary phase in two different states. This is implemented in the binding model, where the salt component is assigned one, and the mAb component "A" is assigned two `bound_states`. `is_kinetic`is set to `False`to simulate the establishment of a rapid equilibrium between bound and unbound particles in the column. This is sensible, because of the high values of the adsorption and desorption rates (Table 3).
# All numerical values from `adsoption_rate` to `conversion_rate` are taken from Table 3. The maximal salt concentration in the mobile phase during elution was used for `reference_liquid_phase_conc` and the column capacity for `reference_solid_phase_conc` (Table A1).  <br>
# The parameters of the `conversion_rate` are listed in a component-row-major ordering for all `bound_states`: <br>
#
# ```
# [
#     comp0fromBnd0toBnd0,
#     comp1fromBnd0toBnd0, comp1fromBnd0toBnd1,
#     comp1fromBnd1toBnd0, comp1fromBnd1toBnd1
# ]
# ```
# For more detail please refer [here](https://forum.cadet-web.de/t/reference-simulation-for-multi-state-sma/818/8). The conversion rates within the same bound state are set to 0.0. 

# %%
import numpy as np

from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import MultistateStericMassAction
from CADETProcess.processModel import Inlet, GeneralRateModel, Outlet
from CADETProcess.processModel import FlowSheet
from CADETProcess.processModel import Process

# Component System
component_system = ComponentSystem()
component_system.add_component('Salt')
component_system.add_component('A')

# Binding Model
binding_model = MultistateStericMassAction(component_system, name='MultistateSMA')
binding_model.bound_states = [1, 2]
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.0, 1.1e31, 7.7e26 ]  # k_a [m_MP³ / (m_SP³ * s)]
binding_model.desorption_rate = [0.0, 5.9e31, 2.0e36]  # k_d [1 / s]
binding_model.characteristic_charge = [0.0, 9.6, 24.7]  # ν [-]
binding_model.steric_factor = [0.0, 47.8, 65.9]  # σ [-]
binding_model.conversion_rate = [0.0, 0.0, 9.4e39, 9.5, 0.0]  # k [1 / s]
binding_model.capacity = 223.55  # Λ [mM]
binding_model.reference_liquid_phase_conc = 520.0  # c_ref [mM]
binding_model.reference_solid_phase_conc = 223.55  # q_ref [mM]

# %% [markdown]
# ```{figure} ./figures/flow_sheet_concentration.svg
# Flow sheet for load-wash-elute process using a single inlet.
# ```

# %% [markdown]
# The unit operation model for this process is the General Rate Model. 
# The `diameter` of the column is inferred from the interstitial velocity u (Table A1, 0.0011438 m/s) which is given by division of the volumetric `flow_rate` (3.1 Experimental) by the product of the `cross_section_area` and the `bed_porosity`. 
# The inital concentrations of salt and protein (3.1 Experimental) are given by `column.c` for the mobile phase and `column.q` for the bound states. 
# All other numerical values are taken from Table A1.

# %%
#Unit Operations
inlet = Inlet(component_system, name='inlet')
inlet.flow_rate = 4.0333e-8  # 2.42 mL / min -> 4.03e−8 m³ / s 

#Transport Model
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model
column.length = 0.215  # L [m]
column.diameter = 0.0115  # L [m]
column.bed_porosity = 0.34  # ε_c [-]
column.particle_radius = 3.25e-5  # r_p [m]
column.particle_porosity = 0.39  # ε_p [-] 
column.axial_dispersion = 10.0e-7  # D_ax [m² / s]
column.film_diffusion = column.n_comp * [2.0e-5]  # k_f [m / s]
column.pore_diffusion = column.n_comp * [9.0e-12]  # D_p [m² / s]
column.surface_diffusion = column.n_bound_states * [0.0]  # [m² / s]
column.c = [69.97, 0.0]  # [mM]
column.q = [binding_model.capacity, 0.0, 0.0]  # [mM]

outlet = Outlet(component_system, name='outlet')

# Flow Sheet
flow_sheet = FlowSheet(component_system)

flow_sheet.add_unit(inlet)
flow_sheet.add_unit(column)
flow_sheet.add_unit(outlet, product_outlet=True)

flow_sheet.add_connection(inlet, column)
flow_sheet.add_connection(column, outlet)


# %% [markdown]
# The following process simulates the load-wash-elute (LWE) CEX under overloaded conditions with a mAb feed concentration of 118.2 g/L = 0.106 mM (3.2. Model calibration). The salt and protein concentrations of the inlet during every step of the LWE are specified using `events`. A linear salt gradient is implemented for elution. The process protocol is taken from Table A1.
#
# The plot of the `simulation_results` shows a "characteristic 'knive blade' shape" (4.1. Standart SMA model) of the large elution peak of the protein. This is the result of complex binding behaviour between the mAb and the tentacle resin (5. Conclusions and outlook). The Multi-State Steric Mass Action Model is able to "quantitatively reproduce" the experimental data (Fig.&nbsp;7d, 4.3. Discussion).

# %%
process = Process(flow_sheet, 'lwe')

load_duration = 74.15 * 60  # [s]
wash_duration = 27.7 * 60  # [s]
t_gradient_start = load_duration + wash_duration
elute_duration = 88.1 * 60  # [s]
elution_slope = 0.053  # [mM / s]

process.cycle_time = t_gradient_start + elute_duration

c_load = np.array([69.97, 0.106])  
c_wash = np.array([69.97, 0.0])
c_elute = np.array([69.97, elution_slope])  

process.add_event('load', 'flow_sheet.inlet.c', c_load)
process.add_event('wash', 'flow_sheet.inlet.c',  c_wash, load_duration)
process.add_event('elute','flow_sheet.inlet.c', list(c_elute), t_gradient_start, indices =[(0,0), (0,1)])

# %%
if __name__ == '__main__':
    from CADETProcess.simulator import Cadet
    process_simulator = Cadet()

    simulation_results = process_simulator.simulate(process)

    from CADETProcess.plotting import SecondaryAxis
    sec = SecondaryAxis()
    sec.components = ['Salt']
    sec.y_label = '$c_{salt}$'

    simulation_results.solution.column.outlet.plot(secondary_axis=sec)
