```julia
include("Network_workplace_model/src/param_sweeps.jl")
```

# Overview

This library accompanies the paper "Modelling the impact of non-pharmaceutical interventions on workplace transmission of SARS-CoV-2 in the home-delivery sector", currently available in preprint form here (final peer-reviewed version will be linked when published).

This julia-based model simulates the stochastic transmission of SARS-CoV-2 on a dynamic contact network representing staff in a delivery/logistics warehouse. The model uses viral load trajectories from data in the literature, as well as test-sensitivity data, to simulate the impact of various non-pharmaceutical interventions including regular PCR and lateral-flow testing. 

Below is a brief guide to how the model works, the underlying code is designed to be easily repurposed for simulating outbreaks in other workplaces and networks of moderate size (< 1000 nodes).

# Pre-requisites

This library also requires files in the Github repository <a href="https://github.com/CarlWhitfield/Viral_load_testing_COV19_model">Viral_load_testing_COV19_model</a>. Paths to this repository are defined relative to the currrent directory, which assume the two repositories be stored in the same directory. Therefore, for ease of setup we recommend downloading these two repositories into the same directory.

# Code structure

The key files in this directory are organised hierarchically, as follows:

## Base level: `src/transmission_model_framework.jl`

This file calls `Viral_load_testing_COV19_model/src/viral_load_infectivity_testpos.jl` from the Github repository <a href="https://github.com/CarlWhitfield/Viral_load_testing_COV19_model">Viral_load_testing_COV19_model</a>. Note that the path to this file in the include statement at the top of this file is relative, so may need to be altered.

This file contains functions to run the generic contact and transmission model underlying the dynamics of the workplace model. This includes the following key functions:


```julia
@doc(init_transmission_model)
```


```julia
@doc(setup_transmission_model!)
```


```julia
@doc(get_infectivities)
```


```julia
@doc(get_network_infections)
```


```julia
@doc(do_infections_randomly!)
```


```julia
@doc(update_in_work!)
```


```julia
@doc(update_all_statuses!)
```


```julia
@doc(update_contact_network!)
```


```julia
@doc(update_testing_state!)
```

## Top level: `src/network_transmission_workplace.jl`

This file calls `transmission_model_framework.jl` and applies the components specific to the delivery workplaces. The main function is described below



```julia
@doc(run_sim_delivery_wp)
```

# Examples

## Example: Reproducing results in pre-print publication

All of the results in the pre-print publication can be reproduced by calling any of the following functions from `Network_workplace_model/src/param_sweeps.jl`. These generate .csv files, which have been converted to .pkl files and stored <a href="">here</a>. The plots can be reproduced using the python notebooks contained <a href="">here</a>. Note that the simulations are stochastic, so outcomes will vary between each run. 

- `run_testing_sweep_outbreak_parcel`: Generates data used to plot Figure 3
- `run_testing_sweep_outbreak_pairs` and `run_testing_sweep_outbreak_pairs_alt`: Generates data used to plot Figure 4
- `run_house_share_sweep_parcel` and `run_house_share_sweep_pairs`: Generates data used to plot Figure 5

The following functions require 3 inputs `Prev::Array{Float64,1}`, `Inc::Array{Float64,1}`, and `Demand::Array{Float64,1}`, which are the daily prevalence, incidence and demand used in the simulations. Data used for each workplace can be found <a href="">here</a>.

- `run_all_interventions_variableprev_scenario_parcel` and `run_all_interventions_variableprev_scenario_parcel_isolfirst`: Generates data used to plot Figure 6
- `run_all_interventions_variableprev_scenario_pairs` and `run_all_interventions_variableprev_scenario_parcel_isolfirst`: Generates data used to plot Figure 7
- `run_all_interventions_separately_scenario_parcel` and `run_all_interventions_separately_scenario_pairs`: Generates data used to plot Figure 8


The following functions are used to reproduce the supplementary figures
- `run_param_sweep_outbreak_parcel`: Generates data used to plot Supplementary Fig S6
- `run_param_sweep_outbreak_pairs`: Generates data used to plot Supplementary Fig S7
- `run_presenteeism_param_sweep_outbreak_parcel`: Generates data used to plot Supplementary Fig S8
- `run_presenteeism_param_sweep_outbreak_pairs`: Generates data used to plot Supplementary Fig S9 and S10
- `run_param_sweep_outbreak_transmod_parcel` and `run_param_sweep_outbreak_transmod_pairs`: Generates data used to plot Supplementary Fig S11
- `run_param_sweep_outbreak_fomite_parcel` and `run_param_sweep_outbreak_fomite_pairs`: Generates data used to plot Supplementary Fig S12
- `run_param_sweep_outbreak_wpsize_parcel` and `run_param_sweep_outbreak_wpsize_pairs`: Generates data used to plot Supplementary Fig S13
- `run_contact_sweeps_outbreak_parcel` and `run_contacts_sweep_outbreak_pairs`: Generates data used to plot Supplementary Fig S14


## Example: Running a stand-alone delivery workplace simulation

The following demonstrates how to run 10,000 simulations with the same parameters (corresponding to the Baseline parameters for the SPDD setting) and write to "baseline_params.csv". Outputs are also stored in the returned dataframe.


```julia
OccPattern = repeat(ParcelOccPattern,NweeksDefault)
PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
NPvec = Int64.(round.(NPparcel*PkgPattern))
PP = Array{Dict{Any,Any},1}(undef,0)
PkgP = Array{Dict{Any,Any},1}(undef,0)
Infinit = [1,2,3]
for ii in Infinit
    PPh = copy(BasicParcelParams)
    PPh["InfInit"] = ii
    push!(PP,PPh)
    push!(PkgP,copy(BasicPkgParams))
end
df = run_many_sims(PP, 10000, OccPattern; NPPerDay = NPvec,
                    PkgParams = PkgP, filename="baseline_params.csv")
```

## Example: Creating a simulation with user-generated contact networks

The script `simple_workplace_model.jl` provides an example of a pared-back workplace model. 

The function `generate_random_contact_network!(sim::Dict, i_day::Int)` updates the contact network `sim["rand_contact_network"]` with random contacts between job roles, set by the matrix `sim["contact_prob_mat"]`. Only contacts of infectious individuals are included, for efficiency. This is overwritten each day.

The function `generate_graph!(sim::Dict, degree_logmean::Float64, degree_logstd::Float64)` generates a random contact network with lognormal degree distribution, which remains the same throughout the simulation. This is stored in `sim["contact_graph"]`.

The function `collate_networks(sim::Dict)` then combined these two graphs into the combined contact network returned by this function, which is then passed to `update_contact_network!`. 


