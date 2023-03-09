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




### Description

`init_transmission_model(N_per_role::Array{Int64,1}, Pisol::Float64, Psusc::Float64)`

Initialises the transmission model.

### Arguments

`N_per_role::Array{Int64,1}` = Array where each entry is the number of employees in a different job role. 

`Pisol::Float64` = Probability an individual will isolate due to symptom onset

`Psusc::Float64` = Probability an individual is susceptible at simulation start

`Inc::Array{Float64,1}` = Community incidence values for each day of the simulation                            (vector length defines the maximum length of the simulation)

`Prev::Array{Float64,1}` = Community prevalence values for each day of the simulation                             (length must match length of `Inc`) 

### Returns:

`sim::Dict()` = Framework for storing simulation data, to be passed to other functions





```julia
@doc(setup_transmission_model!)
```




### Description

`setup_transmission_model!(sim::Dict, Params::Dict, TestParams::Dict, NDays::Int)`

Setup model based on initialisation

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`Params::Dict` = Dict containing model parameters. This function requires:

  * `Params["NContactTypes"]::Int` = Number of contact types
  * `Params["SimType"]::Int` = Outbreak*sim (1) or Scenario*sim (2)
  * `Params["InfInit"]::Int` = If Outbreak sim, group that index case is in (0 is random)

`TestParams::Dict` = Dict containting parameters for testing model. This function uses the following args:

  * `TestParams["is_testing"]::Bool` = whether testing or not
  * `TestParams["protocol"]::Int` = `PCR_mass_protocol` (1) or `LFD_mass_protocol` (2) or `LFD_pattern` (3)
  * `TestParams["tperiod"]::Int` = number of days between tests (if protocol is 1 or 2)
  * `TestParams["Pattern"]::Vector{Bool}` = LFD testing pattern (if protocol is 3)
  * `TestParams["testing_enforced"]::Bool` = Whether testing is mandated (so adherence is overidden)
  * `TestParams["policy_adherence"]::Float` = Probability that an individual adheres to (or ignores) testing protocol
  * `TestParams["test_miss_prob"]::Float` = Probability of missing a scheduled test at random (optional, assumed 1 otherwise)
  * `TestParams["sens_rel"]::Float` = Scaling factor for maximum sensitivity (optional, assumed 1 otherwise)

`Ndays::Int` = maximum number of simulation days

### Returns

`sim_summary::Dict` = Dictionary containing keys for simulation summary outputs

`i_day::Int` = Day number for first day of simulation 

`Anyinf::Bool` = Whether any nodes are infected on first day





```julia
@doc(get_infectivities)
```




### Description

`get_infectivities(sim::Dict, i_day::Int)`

function to list infectious nodes and their current infectivity

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

### Returns:

`inf::Array{Int}` = List of infectious nodes

`inf_scales::Array{Float64}` = Infectivity values for infectious nodes

### See also

[`init_transmission_model`](@ref).





```julia
@doc(get_network_infections)
```




### Description

`get_network_infections(sim::Dict, i_day::Int)`

Generates all infection events on current day based on contact network and node infectivities/susceptibilities

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`i_day::Int` = current day number

### Return

`ipairs::Array{Int64,2}` = 3 x N array listing all successful infection events, potentially including repeat infections. 

  * Row 1 is infector node index
  * Row 2 is infectee node index
  * Row 3 is the mode of infection (defined in the contact network)

### See also

[`init_transmission_model`](@ref), [`get_infectivities`](@ref).





```julia
@doc(do_infections_randomly!)
```




### Description

`do_infections_randomly!(infpairs::Array{Int,2}, sim::Dict, i_day::Int)`

Takes the output of `get_network_infections` and filters to remove repeat infections by selecting  the successful event at random in the event of repeat infectees. Infections are then applied by calling `infect_node!`.

### Arguments

`ipairs::Array{Int64,2}` == 3 x N array listing all successful infection events, returned by `get_network_infections`. 

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`i_day::Int` = current day number

### Return

`infpairs_kept::Array{Int64,2}` = 3 x N array listing final successful infection events

  * Row 1 is infector node index
  * Row 2 is infectee node index
  * Row 3 is the mode of infection (defined in the contact network)

### See also

[`infect_node!`](@ref), [ `init_transmission_model`](@ref), [`get_network_infections`](@ref).





```julia
@doc(update_in_work!)
```




### Description

`update_in_work!(sim::Dict, in_work::Array{Bool,1})`

Updates which nodes are in work on the current day

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`in_work::Array{Bool,1}` =  Boolean vector where true/false entries indicate whether     the node with that index is in work that day. This is stored in sim["at_work"]. This may be used     is assmebling the contact network but note that the contact network can include contacts      between nodes not in work.

### Return

Null





```julia
@doc(update_all_statuses!)
```




### Description

`update_all_statuses!(sim::Dict, i_day::Int)`

Called after day increment to update node infectivity and isolation statuses

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`i_dat::Int` = index of current day

### Return

Null





```julia
@doc(update_contact_network!)
```




Description

`update_contact_network!(sim::Dict, new_network::MetaGraphs.MetaGraph)`

Updates the contact network sim["contact_network"]

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`new_network::MetaGraphs.MetaGraph` =  Network object, in MetaGraphs format, must  contain the following edge metadata:

  * `:weights` -> a `Vector{Float64,1}` of infection rates associated with each contact between   the two nodes `nin` and `nout`. <br>Set using:   `set_prop(sim["contact_network"]),nin,nout,:weights,w)`   where `w` is the vector of weights, and access using:   `get_prop(sim["contact_network"],nin,nout,:weights)`<br>
  * `:types` -> a `Vector{Int64,1}` of indices denoting contact routes for each contact listed in   `:weights`. Can be set and accessed in the same way. Note that the lengths of the contact    weights and types MUST BE EQUAL.

### Return

Null





```julia
@doc(update_testing_state!)
```




### Description

`update_testing_state!(sim::Dict, i_day::Int)`

Called after day increment to update testing status (called if testing is being simulated)

### Arguments

`sim::Dict` = simulation framework (returned by `init_transmission_model`)

`i_dat::Int` = index of current day




## Top level: `src/network_transmission_workplace.jl`

This file calls `transmission_model_framework.jl` and applies the components specific to the delivery workplaces. The main function is described below



```julia
@doc(run_sim_delivery_wp)
```




### Description

`run_sim_delivery_wp(Params::Dict, OccPerDay::Array{Float64,1}, NPPerDay::Array{Int64,1}; PkgParams::Dict=DefaultPkgParams, PairParams::Dict=DefaultPairParams, TestParams::Dict=DefaultTestParams, Incidence::Array{Float64,1} = zeros(length(OccPerDay)), Prevalence::Array{Float64,1} = zeros(length(OccPerDay)))`

Function to run the simulation

### Arguments

`Params::Dict` = Parameter dictionary, as well as parameters required for `init_transmission_model`  and `setup_transmission_model` this uses the following entries:

  * `Params["is_cohorts"]::Bool` = Whether the employees are sorted into cohorts/teams
  * `Params["NDteams"]::Int` = Number of Driver teams/cohorts used in this simulation
  * `Params["NLteams"]::Int` = Number of Loader teams/cohorts used in this simulation
  * `Params["NOteams"]::Int` = Number of Office teams/cohorts used in this simulatio
  * `Params["CohortChangeRate"]::Float64` = Probability per day of an employee changing cohort/team
  * `Params["TeamTimes"]::Vector{Float64}` = F2F contact time between team members in hours used to calculate infection risk
  * `Params["TeamsOutside"]::Vector{Float64}` = To what extent F2F contacts between team members are outside, used to calculate infection risk
  * `Params["TeamDistances"]::Vector{Float64}` = F2F contact distance between team members in m used to calculate infection risk
  * `Params["outdoor_contact_frac"]`::Float64 [Optional] = The fraction of random contacts that occur outside, used to calculate infection risk
  * `Params["distance"]`::Float64 [Optional] = F2F contact distance in m for random contacts used to calculate infection risk
  * `Params["PairIsolation"]::Bool` [Optional] = Whether close-contact pairs isolate when one isolates
  * `Params["CohortIsolation"]::Bool` [Optional] = Whether whole cohort isolates when one member isolates
  * `Params["CarShareIsolation]::Bool` [Optional] = Whether employees who share transport all isolate when one member isolates
  * `Params["HouseShareIsolation"]::Bool` [Optional] = Whether employees who share a house all isolate when one member isolates
  * `Params["Office_WFH"]::Bool` = Whether office staff are working from home in this simulation
  * `Params["AbsenceRate"]::Float64` = Probability of each individual being absent from work for non-COVID reasons per day
  * `Params["F2F_mod"]::Float64` [Optional] = Multiplier to modify F2F tranmission rate
  * `Params["Aerosol_mod"]::Float64` [Optional] = Multiplier to modify aerosol tranmission rate

`OccPerDay::Array{Float64,1}` = Proportion of staff scheduled to be in work for each day of the      simulation up to maximum simulation length

`NPPerDay::Array{Int64,1}` = Number of packages/items to be delivered in each day of the      simulation up to maximum simulation length

`PkgParams::Dict=DefaultPkgParams` = Parameter dictionary defining the interactions of staff with packages. Includes the following elements:

  * `PkgParams["p_fomite_trans"]::Float64` = Probability that handling a package immediately after

an infected individual results in a transmission event

  * `PkgParams["Ltime"]::Float64` = Time window for all loading to occur in hours (only required

if PkgParams["p*fomite*trans"] > 0)

  * `PkgParams["Dtime"]::Float64` = Time window for all deliveries to occur in hours (only required

if PkgParams["p*fomite*trans"] > 0)

  * `PkgParams["PkgHlife"]::Float64` = Half life of viable virus on package surface (only required

if PkgParams["p*fomite*trans"] > 0)

`PairParams::Dict=DefaultPairParams` = Parameter dictionary defining whether close-contact pair work     is simulated. Includes the following elements:

  * `PairParams["is_driver_pairs"]::Bool` = Whether drivers work in close-contact pairs to perform deliveries
  * `PairParams["is_loader_pairs"]::Bool` = Whether loaders work in close-contact pairs to load packages
  * `PairParams["is_window_open"]::Bool` [Optional] = Whether shared delivery vehicles have the "window-open" modifier (Bool)
  * `PairParams["fixed_driver_pairs"]::Bool` = Whether driver pairings are always the same where possible (only required

if `PairParams["is_driver_pairs"]` == true`

  * `PairParams["fixed_loader_pairs"]::Bool` = Whether loader pairings are always the same where possible (only required

if `PairParams["is_driver_pairs"]` == true`

`TestParams::Dict=DefaultTestParams` = Parameter dictionary defining testing intervention

  * `TestParams["is_testing"]::Bool` = Whether testing intervention is in place
  * `TestParams["protocol"]::Int` = '1' or '2' to indicate PCR or LFD repeat testing interventions
  * `TestParams["tperiod"]::Int` = Number of days between scheduled tests
  * `TestParams["specificity"]::Float64` = Probability that a positive test is a true positive
  * `TestParams["delay"]::Float64` = Mean time in days between day of test and result
  * `TestParams["test_miss_prob"]::Float64` = Probability that an individual fails to perform a scheduled test
  * `TestParams["test_pause"]::Float64` = Time in days that a person stops testing after a positive result
  * `TestParams["testing_enforced"]::Float64` = Whether testing intervention is mandated and forced

`Incidence::Array{Float64,1} = zeros(length(OccPerDay))` = Array of community incidence values used to simulate rate of workplace introductions

`Prevalence::Array{Float64,1} = zeros(length(OccPerDay))` =  = Array of community incidence values used to simulate risk of introduction via customer contact




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

The following demonstrates how to run 500 simulations with the same parameters (corresponding to the Baseline parameters for the SPDD setting) for each index case possibility. Outputs are written to "baseline_params.csv" and also stored in the returned dataframe.


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
df = run_many_sims(PP, 500, OccPattern; NPPerDay = NPvec,
                    PkgParams = PkgP, filename="baseline_params.csv")
```

    1/3
    2/3
    3/3





<div class="data-frame"><p>6,000 rows × 48 columns (omitted printing of 43 columns)</p><table class="data-frame"><thead><tr><th></th><th>TPAsympIsolatorsFrac</th><th>IndexCaseViralLoad</th><th>CarShareFactor</th><th>PairInfectionFrac</th><th>IsolatorsFrac</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>0.0</td><td>7.05396</td><td>0.05</td><td>0.0</td><td>0.0263158</td></tr><tr><th>2</th><td>0.0</td><td>7.05396</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>3</th><td>0.0</td><td>7.05396</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>4</th><td>0.0</td><td>7.05396</td><td>0.05</td><td>0.0</td><td>0.0144928</td></tr><tr><th>5</th><td>0.0</td><td>7.91163</td><td>0.05</td><td>0.0</td><td>0.0263158</td></tr><tr><th>6</th><td>0.0</td><td>7.91163</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>7</th><td>0.0</td><td>7.91163</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>8</th><td>0.0</td><td>7.91163</td><td>0.05</td><td>0.0</td><td>0.0144928</td></tr><tr><th>9</th><td>0.0</td><td>7.24192</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>10</th><td>0.0</td><td>7.24192</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>11</th><td>0.0</td><td>7.24192</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>12</th><td>0.0</td><td>7.24192</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>13</th><td>0.0</td><td>7.7016</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>14</th><td>0.0</td><td>7.7016</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>15</th><td>0.0</td><td>7.7016</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>16</th><td>0.0</td><td>7.7016</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>17</th><td>0.0</td><td>8.43963</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>18</th><td>0.0</td><td>8.43963</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>19</th><td>0.0</td><td>8.43963</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>20</th><td>0.0</td><td>8.43963</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>21</th><td>0.0</td><td>7.7977</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>22</th><td>0.0</td><td>7.7977</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>23</th><td>0.0</td><td>7.7977</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>24</th><td>0.0</td><td>7.7977</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>25</th><td>0.0</td><td>7.28999</td><td>0.05</td><td>0.0</td><td>0.0263158</td></tr><tr><th>26</th><td>0.0</td><td>7.28999</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>27</th><td>0.0</td><td>7.28999</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>28</th><td>0.0</td><td>7.28999</td><td>0.05</td><td>0.0</td><td>0.0144928</td></tr><tr><th>29</th><td>0.0</td><td>7.31074</td><td>0.05</td><td>0.0</td><td>0.0263158</td></tr><tr><th>30</th><td>0.0</td><td>7.31074</td><td>0.05</td><td>0.0</td><td>0.0</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



## Example: Creating a simulation with user-generated contact networks

The script `simple_workplace_model.jl` provides an example of a pared-back workplace model. 

The function `generate_random_contact_network!(sim::Dict, i_day::Int)` updates the contact network `sim["rand_contact_network"]` with random contacts between job roles, set by the matrix `sim["contact_prob_mat"]`. Only contacts of infectious individuals are included, for efficiency. This is overwritten each day.

The function `generate_graph!(sim::Dict, degree_logmean::Float64, degree_logstd::Float64)` generates a random contact network with lognormal degree distribution, which remains the same throughout the simulation. This is stored in `sim["contact_graph"]`.

The function `collate_networks(sim::Dict)` then combined these two graphs into the combined contact network returned by this function, which is then passed to `update_contact_network!`. 


