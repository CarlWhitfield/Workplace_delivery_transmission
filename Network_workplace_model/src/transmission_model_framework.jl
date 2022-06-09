"""How to use the functions in this code




"""

#include("../../../Viral_load_testing_COV19_model/src/viral_load_infectivity_testpos_v2.jl")
include("../../../Viral_load_testing_COV19_model/src/viral_load_infectivity_testpos.jl")

using Graphs
using MetaGraphs

#sim type definitions
const Outbreak_sim = 1
const Scenario_sim = 2


const isol_time = 10

#Status definitions
const Susc = 0
const Expd = 1
const Symp = 2
const Recd = 3

#Isolation reasons
const Suspected_Isol = 1
const Confirmed_Isol = 2
const SuspContact_Isol = 3
const ConfContact_Isol = 4
const NIsolStatuses = 4
const IsolStatusNames = ["Symptomatic","Positive","SymptomContact","PositiveContact"]

#Infection rates
infection_rate_F2F = 0.15   #infectivity of F2F interactions per hour = microCOVID approx
const outside_factor = 0.2
const no_talking_factor = 0.2
const mask_factor_infector = 0.25
const mask_factor_infectee = 0.5
const distance_factor_per_m = 0.5

#P_inf = 1 - exp(-w(t,x)) -- this function returns w
function return_infection_weight(distance::Float64, duration::Float64,
                                 outdoor::Bool, talking::Bool)
    w = infection_rate_F2F * duration * (distance_factor_per_m^(distance-1.0))
    if outdoor
        w = w * outside_factor
    end
    if talking == false
        w = w * no_talking_factor
    end

    return w
end

#beta = w * mask1_factor * mask2_factor
#returns p_inf
function convert_weight_to_prob(weight::Float64, inf::Float64,
                     susc::Float64; mask_infector::Bool = false,
                     mask_infectee::Bool = false)
    beta = weight*inf*susc
    if mask_infector
        beta *= mask_factor_infector
    end
    if mask_infectee
        beta *= mask_factor_infectee
    end

    return 1 - exp(-beta)
end

"""
    init_transmission_model(N_per_role::Array{Int64,1}, Pisol::Float64, Psusc::Float64)

Initialises the transmission model.

## Arguments

`N_per_role` = Array where each entry is the number of employees in a different job role.

`contact_graph` = A graph of contacts (can be updated each day), with the following metadata on the edges:

        `:weights`, a float array of weights listing the duration*transmission rate for each type of contact that occured

        `:types`, an integer array (the same size as weights) indicating the type of contact associated with each

`Pisol` = Probability an individual will isolate

`Psusc` = Probability an individual is susceptible at simulation start

## Returns:

`sim::Dict()` = Framework for storing simulation data, to be passed to other functions
"""
function init_transmission_model(N_per_role::Array{Int64,1}, Pisol::Float64, Psusc::Float64, Inc::Array{Float64,1},
    Prev::Array{Float64,1})
    Ntot = sum(N_per_role)
    Nj = length(N_per_role)
    Ndays = length(Inc)
    sim = Dict("Ntot"=>Ntot, "N"=>N_per_role, "job"=>zeros(Int8, Ntot),
               "Njobs"=>Nj, "inf_time"=>(zeros(Int64, Ntot) .- 1),
                "symp_day"=>-ones(Int64, Ntot),
                "asymptomatic"=>zeros(Bool, Ntot),
                "would_isolate"=>zeros(Bool, Ntot),
                "isolation_status"=>zeros(Int64, Ntot),
                "isolation_start_day"=>zeros(Int64,(NIsolStatuses,Ntot)),
                "isolation_end_day"=>zeros(Int64,(NIsolStatuses,Ntot)),
        
#                 "isolation_contactor"=>zeros(Int64, Ntot),
#                 "isolation_contacts"=>Array{Array{Int64,1},1}(undef, Ntot),
#                 "isolation_action"=>Array{Array{Int64,1},1}(undef, Ntot),
#                 "isolation_action_day"=>Array{Array{Int64,1},1}(undef, Ntot),
#                 "isolation_action_counter"=>ones(Int64,Ntot),
                "new_isolators"=>Array{Int64,1}(undef, 0),
#                 "isolation_status_true_test"=>zeros(Bool, Ntot),
#                 "isolation_status_false_test"=>zeros(Bool, Ntot),
#                 "isolation_status_contact_true_test"=>zeros(Bool, Ntot),
#                 "isolation_status_contact_false_test"=>zeros(Bool, Ntot),
#                 "isolating_due_to_true_test"=>zeros(Bool, Ntot),
#                 "isolating_due_to_false_test"=>zeros(Bool, Ntot),
#                 "isolating_due_to_contact_true_test"=>zeros(Bool, Ntot),
#                 "isolating_due_to_contact_false_test"=>zeros(Bool, Ntot),
                "susceptibility"=>ones(Float64,Ntot),
                "VL_mag"=>zeros(Float64, Ntot),
                "inf_mag"=>zeros(Float64, Ntot),
                "infection_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
                "VL_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
                "at_work"=>zeros(Bool, Ntot),
                "infection_status"=>zeros(Int8,Ntot),
                "days_infectious" => zeros(Int64,Ntot),
                "contact_network" => MetaGraphs.MetaGraph(SimpleGraph(Ntot)),
                "job_sorted_nodes"=>Array{Array{Int64,1},1}(undef, Nj),
                "Incidence"=>Inc, "Prevalence"=>Prev)
    istart = 0
    for (i,n) in enumerate(N_per_role)
        nrh = (istart+1):(istart + n)
        sim["job"][nrh] .= i   #assign job/role labels
        sim["job_sorted_nodes"][i] = nrh
        #null initialise these
        sim["infection_profiles"][nrh] .= fill(zeros(0),n)
        sim["VL_profiles"][nrh] .= fill(zeros(0),n)
        istart += n
    end
    nr = 1:Ntot
    sim["would_isolate"][generate_isolations(sim["Ntot"], Pisol)] .= true   #generate propensity to isolate
    sim["non_isolators"] = nr[(sim["would_isolate"] .== false)]
    sim["susceptibility"][randsubseq(1:sim["Ntot"], 1.0 - Psusc)] .= 0.0    #generate non-susceptibles

    return sim
end

function generate_symp_day_from_time(stime::Float64)
    return Int64(floor(rand(0:1) + stime))
end

function infect_node!(sim::Dict, i::Int, time::Int)
    if haskey(sim,"isolation_network")
        contact_tracing_network = sim["isolation_network"]
    else
        contact_tracing_network = SimpleGraph(sim["Ntot"])
    end
    
    sim["infection_status"][i] = Expd
    sim["susceptibility"][i] = 0.0
    sim["inf_time"][i] = time

    #Generate viral load
    build_viral_load_distribution!(sim, i)
    if sim["asymptomatic"][i]
        will_isolate = false
    else
        will_isolate = sim["would_isolate"][i]
    end
    sim["days_infectious"][i] = length(sim["infection_profiles"][i])

    if haskey(sim,"test_protocol")
       sim["test_pos_profiles"][i] = get_pos_profile(sim, i, sim["test_protocol"])
    end

    #Sort out symptomatic isolation and contact isolation
    #Decide what happens with multiple statuses
    #Assume close contact isolation is guaranteed to isolate if infected node isolates (i.e. enforced)
    if will_isolate
        isol_length = get_symp_isolation_days(sim, i)    #depends on confimatory PCR policy
        if sim[isolation_status[i]] != Confirmed_Isol   
            #this should almost always be the case, but if somebody is already isolated due to a (false)
            #positive test, then they ignore their symptoms and continue with their isolation as before
            ic = copy(sim["isolation_action_counter"][i])
            ignore_symp_isol_start = false
            ignore_symp_isol_end = false
            if length(sim["isolation_action"][i] >= ic)   #find upcoming isolation actions 
                while sim["isolation_action_day"][i][ic] < time + sim["symp_day"][i]
                    #if sim["isolation_action"] is a confirmed case then ignore symp_isol
                    if sim["isolation_action"][i] == Confirmed_Isol
                        ignore_symp_isol_start == true
                        ignore_symp_isol_end == true
                    end
                    ic += 1
                end
                while sim["isolation_action_day"][i][ic] < time + sim["symp_day"][i] + isol_length
                    #if sim["isolation_action"] is a confirmed case then ignore symp_isol
                    if sim["isolation_action"][i] == Confirmed_Isol
                        ignore_symp_isol_end == true
                    end
                    ic += 1
                end
            end
            if ignore_symp_isol_start == false
                push!(sim["isolation_action_day"][i], time + sim["symp_day"][i])
                push!(sim["isolation_action"][i], Suspected_Isol)
            end
            if ignore_symp_isol_end == false
                push!(sim["isolation_action_day"][i], time + sim["symp_day"][i]) + isol_length
                push!(sim["isolation_action"][i], 0)
            end
        
        
        
        
        
        if sim["isolation_time"][i] > time + sim["symp_day"][i] ||
           sim["isolation_time"][i] == 0 #might already be isolating
               sim["isolation_time"][i] = time + sim["symp_day"][i]
               sim["isolation_reason"][i] = Suspected_Isol
        end
#         print("Symp isol: ", i, ' ', time, ' ', sim["isolation_time"][i], '\n')
        traced = neighbors(contact_tracing_network, i)
        if length(traced) > 0
            later_isol = (sim["isolation_time"][traced] .>  (time + sim["symp_day"][i]))
            no_isol = (sim["isolation_time"][traced] .== 0)
            to_update = no_isol .| later_isol
            #if not isolating already or isolating later than this
            if sum(to_update) > 0
                sim["isolation_time"][traced[to_update]] .=  time + sim["symp_day"][i]
                sim["isolation_reason"][traced[to_update]] .= SuspContact_Isol
                sim["isolation_contactor"][traced[to_update]] .= i
                sim["isolation_contacts"][i] = traced[to_update]   
            end
        end
    end
end

#Need to work out a better way to deal with isolations
#Isolation status: True or False
#Isolation time: If 0: never isolated, otherwise, last isolation start day (can isolate multiple times)
#Isolation reason: Should exist for anybody with isolation time > 0. Options
                #Suspected (symptoms, no test) -> Either isolate for 10-days or until conf PCR depending on policy
                #Confirmed (positive test) -> If previously: no-status/close-contact status -> 10-day isolation from test-date
                #                          -> If previously: suspected -> unchanged
                #Close-Contact (policy defined) of (suspected + isolating) -> Either isolate or regularly test depending on policy
                #Close-contact of confirmed case -> Either isolate or regularly test, depending on policy
#Leaving isolation:
                #Suspected (symptoms, no test) -> If PCR conf negative leave with some prob. (may still be ill)
                #Confirmed -> After 10-days automatically, could have testing policy options
                #Close-contact of suspected -> If PCR conf of index case negative
                #Close-contact of confirmed -> 10-days or X-days + return to work PCR?
#Need an option for individual testing

function infect_random!(sim::Dict, InfInit::Int, i_day::Int)
    #infect from non infectious
    if InfInit > 0
        nrj = sim["job_sorted_nodes"][InfInit]   #all nodes with job no. InfInit
        nrsub = nrj[sim["infection_status"][nrj] .== Susc]
    else  #infinit 0 means any random person
        nr = 1:sim["Ntot"]
        nrsub = nr[sim["infection_status"] .== Susc]
    end
    inf = rand(nrsub,1)[1]
#     partner = 0
#     if length(pairs) > 0   #for pair isolation
#         partner = find_partner(pairs, inf)
#     end
    infect_node!(sim, inf, i_day)

    return inf
end


function get_infected_bool(sim::Dict)
    return (sim["infection_status"] .!= Susc) .* (sim["infection_status"] .!= Recd)
end


"""
    get_infectivities(sim::Dict, i_day::Int)

function to list infectious nodes and their current infectivity

## Arguments

`sim` = simulation framework (returned by `init_transmission_model`)

## Returns:

`inf` = List of infectious nodes

`inf_scales` = Infectivity values for infectious nodes

## See also

`init_transmission_model`
"""
function get_infectivities(sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    infected = get_infected_bool(sim)
    inf = nr[infected .* sim["at_work"]]
    inf_scales = zeros(length(inf))
    j = 1
    for i in inf
        inf_scales[j] = get_infectivity(sim, i_day, i)
        j += 1
    end
    return inf, inf_scales
end

function get_infectivity(sim::Dict, i_day::Int, index::Int)
    return sim["infection_profiles"][index][i_day - sim["inf_time"][index] + 1]
end

function increment_infectivity!(sim::Dict, i_day::Int)
    #based on dates, stop infectivity
    nr = 1:sim["Ntot"]
    inf = nr[get_infected_bool(sim)]
    #people becoming Symptomatic (or asymptomatic)
    incend = inf[((i_day .- sim["inf_time"][inf]) .== sim["symp_day"][inf])]
    sim["infection_status"][incend] .= Symp
    #people recovered
    infend = inf[((i_day .- sim["inf_time"][inf]) .== sim["days_infectious"][inf])]
    sim["infection_status"][infend] .= Recd
end

function increment_isolations!(sim::Dict, i_day::Int)
    nr = collect(1:sim["Ntot"])
    sim["new_isolators"] = Array{Int64,1}(undef,0)
    for i in nr[length.(sim["isolation_action_days"]) .<= sim["isolation_action_counter"]]  #loop over people with tests remaining
        if i_day == sim["isolation_action_days"][i][sim["isolation_action_counter"][i]]
            action = sim["isolation_action"][i][sim["isolation_action_counter"][i]]
            #make a record of people entering isolation
            if sim["isolation_status"][i] == 0 && action > 0
                push!(sim["new_isolators"],i)
            end
            #update isolation status
            sim["isolation_status"][i] = action
            sim["isolation_action_counter"][i] += 1
        end
    end
    
    print("People in isolation: ", nr[sim["isolation_status"] .> 0],'\n')
end

"""

"""
function get_introductions(sim::Dict, i_day::Int, introduction_ID::Int; 
                           rel_rates::Array{Float64,1} = ones(sim["Njobs"]))
    #doesn't matter if individual is in work or not
    iprob = zeros(sim["Ntot"])
    for j in 1:sim["Njobs"]
        iprob[sim["job_sorted_nodes"][j]] .= rel_rates[j]*sim["Incidence"][i_day]
    end
    ind = collect(1:sim["Ntot"])
    newinfs = ind[rand(sim["Ntot"]) .< (iprob)]
    NI = length(newinfs)
    infpairs = Array{Int64,2}(undef,3,0)
    if NI > 0
        infpairs = [transpose(zeros(Int64,NI)); transpose(newinfs);
                       introduction_ID*transpose(ones(Int64,NI))]
    end
    
    
    return infpairs
end


"""

"""
function get_network_infections(sim::Dict, i_day::Int)
    inf, inf_scales = get_infectivities(sim, i_day)
    ipairs = Array{Int64,2}(undef,3,0)
    #collect all unique edges
    #all have different infection rates
    #get indices of all infectious edges
    #print(inf,'\n')
    if length(inf) > 0
        #get all neighbours of infectious nodes
        nbrs = neighbors.(Ref(sim["contact_network"]), inf)
        #create array of infectious nodes for each contact edge
        nin = vcat(fill.(inf,length.(nbrs))...)
        #flatten the other array to match
        nout = vcat(nbrs...)
        #get infectious node infectivity for each edge
        j_inf = vcat(fill.(inf_scales,length.(nbrs))...)
        #get target node susceptibility for each edge
        j_susc = sim["susceptibility"][nout]

        #get the edge weights (total transmission rate * contact time, can be multiple per edge)
        w = get_prop.(Ref(sim["contact_network"]),nin,nout,:weights)
        t = get_prop.(Ref(sim["contact_network"]),nin,nout,:types)
    
        #fill and flatten all arrays so there is one entry for each entry
        nin_all = vcat(fill.(nin,length.(w))...)
        nout_all = vcat(fill.(nout,length.(w))...)
        j_inf_all = vcat(fill.(j_inf,length.(w))...)
        j_susc_all = vcat(fill.(j_susc,length.(w))...)
        beta = vcat(w...)
        etype = vcat(t...)

        #get probability of infection on each edge
        eprob = convert_weight_to_prob.(beta, j_inf_all, j_susc_all)

        #draw which will be infectious (duplicates don't matter)
        #print(ecprob,"\n\n")
        eind = collect(1:length(eprob))
        einf = eind[rand(length(eprob)) .< (eprob)]

        if length(einf) > 0
            ipairs = hcat(ipairs, [transpose(nin_all[einf]);
                    transpose(nout_all[einf]); transpose(etype[einf])])
        end
    end

    return ipairs
end

"""

"""
function do_infections_randomly!(infpairs::Array{Int,2}, sim::Dict, i_day::Int)
    N = size(infpairs,2)
    nkeep = Array{Int64,1}()
    if N > 0
        ind = randperm(N)
        for i in ind
            k = infpairs[2,i]   #node to be infected
            if sim["infection_status"][k] == Susc  #if susceptible
                infect_node!(sim, k, i_day)
                push!(nkeep, i)
            end
        end
    end
    infpairs_kept = infpairs[:,nkeep]
    return infpairs_kept
end

"""


"""
function update_in_work!(sim::Dict, in_work::Array{Bool,1})
    sim["at_work"] = in_work
end

"""


"""
function update_all_statuses!(sim::Dict, i_day::Int)
    increment_infectivity!(sim, i_day)
    increment_isolations!(sim, i_day)
end


"""update contact network with new version



"""
function update_contact_network!(sim::Dict, new_network::MetaGraphs.MetaGraph)
    sim["contact_network"] = new_network
end


"""

"""
function do_testing!(sim::Dict, testing_params::Dict, i_day::Int,
          isolation_network::Graph)

    #do testing if test day
    nr = collect(1:sim["Ntot"])
    i_testers = Array{Int64,1}(undef,0)
    for i in nr[length.(sim["test_days"]) .<= sim["test_day_counter"]]  #loop over people with tests remaining
        if i_day == sim["test_days"][i][sim["test_day_counter"][i]]
            push!(i_testers,i)   #add people who will test today to this list
        end
    end
    Ntesters = length(i_testers)
    if Ntesters > 0        
        print("Day: ", i_day, '\n')
        bool_exposed = (sim["infection_status"][i_testers] .!= Susc)
        exposed = i_testers[bool_exposed]
        should_be_positive = exposed[length.(sim["test_pos_profiles"][exposed]) .>
                                       i_day .- sim["inf_time"][exposed]]
        print("Exposed: ", exposed,'\n')
        print("Should be positive: ", should_be_positive,'\n')
        pos_tests = zeros(Bool,Ntesters)
        #false positives due to specificity
        FPprob = (1-testing_params["specificity"])
        false_pos = randsubseq(i_testers, FPprob)
        # print("False pos: ",false_pos,'\n')
        pos_tests[false_pos] .= true
        LP = length(should_be_positive)
        if LP > 0
            pos_probs_exposed = zeros(LP)
            for (i,c) in enumerate(should_be_positive)
                pos_probs_exposed[i] = sim["test_pos_profiles"][c][1 + i_day - sim["inf_time"][c]]
            end
            #print(pos_probs_exposed,'\n')
            # print("Pos prob: ", pos_probs_exposed, '\n')
            true_pos = should_be_positive[rand(LP) .< pos_probs_exposed]
            pos_tests[true_pos] .= true
        end
        print("Positive: ", i_testers[pos_tests],'\n')
        # print("Pos tests: ",nr[pos_tests],'\n')

        isolators = apply_positive_tests!(sim, i_testers[pos_tests], testing_params, i_day, isolation_network)
        
        print("Isolators: ", isolators, '\n')
        
        sim["test_day_counter"][i_testers] .= sim["test_day_counter"][i_testers] .+ 1

        return isolators
    else
        return Array{Int64,1}(undef,0)
    end
end

"""

"""
function apply_positive_tests!(sim::Dict, detected::Array{Int64,1}, testing_params::Dict, i_day::Int,
                               isolation_network::Graph)
    will_isolate = detected[(sim["will_isolate_with_test"][detected])]
    isol_day = Int64.(round.(i_day + testing_params["delay"] .+ rand([0,0.5],length(will_isolate))))
    bool_update_isol_time = Bool.(((sim["isolation_time"][will_isolate] .== 0) .+
                                   (sim["isolation_time"][will_isolate] .> isol_day)))
    to_update = will_isolate[bool_update_isol_time]

    if length(will_isolate) > 0
        if length(to_update) > 0
            sim["isolation_time"][to_update] .= isol_day[bool_update_isol_time]

            bool_false_test = (sim["infection_status"][to_update] .== Susc) .|
                             (sim["infection_status"][to_update] .== Recd)
            sim["isolating_due_to_true_test"][to_update[.!bool_false_test]] .= true
            sim["isolating_due_to_false_test"][to_update[bool_false_test]] .= true
        end

        #locate contacts for those who will isolate
        nbrs = neighbors.(Ref(isolation_network), will_isolate)
        srcs = vcat(fill.(will_isolate, length.(nbrs))...)
        dests = vcat(nbrs...)

        for (i,ip) in enumerate(dests)
            #all partners of people who report do isolate
            #if not isolating or contact's isolation time is before theirs
            if sim["isolation_time"][ip] == 0 || sim["isolation_time"][ip] > sim["isolation_time"][srcs[i]]
                sim["isolation_time"][ip] = sim["isolation_time"][srcs[i]]
                if sim["isolating_due_to_true_test"][srcs[i]]
                    sim["isolating_due_to_contact_true_test"][ip] = true
                elseif sim["isolating_due_to_false_test"][srcs[i]]
                    sim["isolating_due_to_contact_false_test"][ip] = true
                end
            end
        end
        will_isolate = unique(vcat(will_isolate, dests))
        #remove all scheduled tests in pause window
        for i in will_isolate 
            bool_elements = (sim["test_days"][i] .> i_day) & 
                            (sim["test_days"][i] .<= i_day + testing_params["test_pause"])
            i_elements = collect(1:length(sim["test_days"][i]))
            deleteat!(sim["test_days"][i],i_elements[bool_elements])
        end
    end
    return will_isolate
end

"""

"""
function basic_sim_setup(sim::Dict, i_day::Int64, Ndays::Int64)
    Nj = sim["Njobs"]
    summary = Dict{Any,Any}("time"=>collect(1:Ndays),
                "Susceptible"=>zeros(Int64,(Nj,Ndays)),
                "Infectious"=>zeros(Int64,(Nj,Ndays)),
                "Exposed"=>zeros(Int64,(Nj,Ndays)),
                "Isolated"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToSymptoms"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToTestAsymp"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToTestSymp"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToFalsePos"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToContactTruePos"=>zeros(Int64,(Nj,Ndays)),
                "IsolatedDueToContactFalsePos"=>zeros(Int64,(Nj,Ndays)),
                "NewIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewSympIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewTestAsympIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewTestSympIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewFalseIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewContactTruePosIsolators"=>zeros(Int64,(Nj,Ndays)),
                "NewContactFalsePosIsolators"=>zeros(Int64,(Nj,Ndays)),
                "Presenting"=>zeros(Int64,(Nj,Ndays)),
                "Asymptomatic"=>zeros(Int64,(Nj,Ndays)),
                "Recovered"=>zeros(Int64,(Nj,Ndays)),
                "InfsByType"=>Array{Array{Int64,2},1}(undef,sim["NTypes"]),
                "ExternalIntroductions"=>zeros(Int64,(Nj,Ndays)))
    for i in 1:sim["NTypes"]
        summary["InfsByType"][i] = zeros(Int64,(Nj,Ndays))
    end

    return summary
end

"""

"""
function sim_setup!(sim::Dict, InfInit::Int64, i_day::Int64, Ndays::Int64)
    index_case = infect_random!(sim,InfInit,i_day)
    sim_summary = basic_sim_setup(sim,i_day,Ndays)
    sim_summary["IndexCase"] = index_case
    sim_summary["IndexCaseInfections"] = 0
    sim_summary["IndexCaseInfectivity"] = sim["inf_mag"][index_case]
    sim_summary["IndexCasePeakVL"] = sim["VL_mag"][index_case]
    for j in 1:sim["Njobs"]
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, Array{Int64,2}(undef,3,0), i_day)
    return sim_summary
end

"""

"""
function scenario_sim_setup!(sim::Dict, i_day::Int64, Ndays::Int64)

    sim_summary = basic_sim_setup(sim, i_day, Ndays)
    sim_summary["Incidence"] = sim["Incidence"]
    sim_summary["Prevalence"] = sim["Prevalence"]
    for j in 1:sim["Njobs"]
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, Array{Int64,2}(undef,3,0), i_day)
    return sim_summary
end

"""

"""
function update_sim_summary!(summary::Dict, sim::Dict, inf_pairs::Array{Int64,2}, i_day::Int)
    nr = 1:sim["Ntot"]
    new_isolator_bool =  zeros(Bool,nr)
    new_isolator_bool[sim["new_isolators"]] .= true
    at_work_bool = (sim["at_work"] .== true)  
    after_onset_bool = (sim["infection_status"] .== Symp)  #after symptom onset time (which they have even if asymptomatic)
    for j in 1:sim["Njobs"]
        jt = nr[sim["job"] .== j]  #update this

        summary["Susceptible"][j,i_day] = sum(sim["infection_status"][jt] .== Susc)
        summary["Recovered"][j,i_day] = sum(sim["infection_status"][jt] .== Recd)
        summary["Infectious"][j,i_day] = length(jt) - summary["Susceptible"][j,i_day]
                                                    - summary["Recovered"][j,i_day])
        summary["Exposed"][j,i_day] = sum(sim["infection_status"][jt] .== Expd)
        summary["Isolated"][j,i_day] = sum(sim["isolation_status"][jt] .> 0)
        summary["NewIsolators"][j,i_day] = sum(new_isolator_bool[jt])
        
        for n in 1:NIsolStatuses
            name = "Isolated" + IsolStatusNames[n]
            summary[name][j,i_day] = sum(sim["isolation_status"][jt] .== n)
            name2 = "New" + IsolStatusNames[n] + "Isolators"
            summary[name2][j,i_day] = sum((sim["isolation_status"][jt] .== n).*new_isolator_bool[jt])
        end
        summary["Presenting"][j,i_day] = sum(jt .* after_onset_bool .* (sim["asymptomatic"].==false) .* at_work_bool)
        summary["Asymptomatic"][j,i_day] = sum(jt .* after_onset_bool .* sim["asymptomatic"])
    end

    for k in 1:length(inf_pairs[1,:])
        j = sim["job"][inf_pairs[2,k]]  #job of infectee
        summary["InfsByType"][inf_pairs[3,k]][j,i_day] += 1
    end
    
    if haskey(summary,"IndexCase")
        summary["IndexCaseInfections"] += sum(inf_pairs[1,:] .== summary["IndexCase"])
    end
end



"""

"""
function trim_sim_summary!(sim_summary, i_day, Ndays)
    if i_day < Ndays
        for (key,value) in sim_summary
            if isa(value, Array)
                if key == "InfsByType"
                    for n in 1:length(sim_summary[key])
                        sim_summary[key][n] = sim_summary[key][n][:,1:i_day]
                    end
                elseif length(size(value)) == 1
                    sim_summary[key] = sim_summary[key][1:i_day]
                elseif length(size(value)) == 2
                    sim_summary[key] = sim_summary[key][:,1:i_day]
                end
            end
        end
    end
end

"""

"""
# function simulate_day!(sim::Dict, i_day::Int)

#     #add parts of previous function that go here


#     infpairs = get_network_infections(sim::Dict, i_day::Int)
#     infpairs_reduced = do_infections_randomly!(sim, infpairs)

# end
function setup_transmission_model!(sim::Dict, Params::Dict, TestParams::Dict,
                                   NDays::Int)
    i_day = rand(1:7)
    if any(TestParams["is_testing"])
        if (TestParams["protocol"] == PCR_mass_protocol
         || TestParams["protocol"] == LFD_mass_protocol)
            sim["non_testers"] = zeros(Int,sim["Ntot"])
            for j in 1:sim["Njobs"]
                if TestParams["is_testing"][j] == false
                    sim["non_testers"][sim["job_sorted_nodes"][j]] .= true
                end
            end
            sim["test_days"], sim["test_day_counter"] =
                     init_testing!(sim, TestParams, i_day, NDays; fill_pos_profiles=false)
        end  #add options for other protocols here
    else
        sim["test_day_counter"] = 1
        sim["test_days"] = [0]
    end
    Anyinf = true
    if Params["SimType"] == Outbreak_sim
        sim_summary = sim_setup!(sim, Params["InfInit"], i_day, NDays)
        Anyinf = any(((sim["infection_status"] .== Susc)
                      .+ (sim["infection_status"] .== Recd)) .== 0)
    elseif Params["SimType"] == Scenario_sim
        sim_summary = scenario_sim_setup!(sim, i_day, NDays)
    end

    return sim_summary, i_day, Anyinf
end


function example_sim_loop!(sim::Dict, i_day::Int,  testing::Bool, TestParams::Dict,
                           isolation_network::Graph, Prev::Float64)

#     #do_testing
#     if testing
#         new_isolators = do_testing!(sim, TestParams, i_day, isolation_network)
#     end

#     #update infectivity and isolation status
#     update_all_statuses!(sim, i_day)

#     #update_in_work
#     #simple thing
#     update_in_work!(sim, at_work)

#     infpairs = Array{Int64,2}(undef,3,0)

#     #introductions
#     infpairs = get_introductions(infpairs, sim, i_day)

#     #get all contacts
#     g = #gen network

#     update_contact_network!(sim, g)
#     infpairs = get_network_infections(infpairs, sim, i_day)

#     infpairs_final = do_infections_randomly!(infpairs, sim, i_day)

#     return infpairs_final
end


function run_example_sim(Params::Dict; visualise::Bool = false, testing::Bool=false,
           TestParams::Dict=Dict(), Incidence::Array{Float64,1} = zeros(length(OccPerDay)),
           Prevalence::Array{Float64,1} = zeros(length(OccPerDay)))

#     #     if visualise
# #         sim, node_x, node_y = initialise_withvis(Params,
# #           PairParams, degree_logmean, degree_logstd, is_network, is_pairs, Incidence, Prevalence)
# #     else
# #         sim = initialise_novis(Params,
# #           PairParams, degree_logmean, degree_logstd, is_network, is_pairs, Incidence, Prevalence)
# #     end

#     i_day, Go = setup_transmission_model(sim, testing, TestParams)
#     while Go && (i_day <= length(OccPerDay))
#         infpairs = sim_loop_delivery_wp!(sim, sim_summary, i_day, OccPerDay[i_day],
#                   NPPerDay[i_day], Params, PkgParams, PairParams, TestParams,
#                   is_pairs, is_network, testing, test_days[test_day_counter])

# #         fname = Printf.@sprintf("infection_network_%03d.png",i_day)
# #         if visualise
# #             print_infection_network(sim, fname, infpairs, node_x, node_y, pairs)
# #         end

# #         if testing && (i_day == test_days[test_day_counter])
# #             test_day_counter += 1
# #         end

#         if Params["SimType"] == Outbreak_sim
#             Go = any(((sim["infection_status"] .== Susc)
#                   .+ (sim["infection_status"] .== Recd)) .== 0)
#         end

#         i_day += 1
#     end
#     trim_sim_summary!(sim_summary, i_day-1, length(OccPerDay))

#     return sim_summary
end
