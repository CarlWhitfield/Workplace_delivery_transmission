"""How to use the functions in this code




"""




include("../../../Viral_load_testing_COV19_model/src/viral_load_infectivity_testpos.jl")

using LightGraphs
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

#Infection rates
const infection_rate_F2F = 0.06   #infectivity of F2F interactions per hour = microCOVID approx
const outside_factor = 0.2
const mask_factor_infector = 0.25
const mask_factor_infectee = 0.5
const distance_factor_per_m = 0.5

function return_infection_weight(distance::Float64, duration::Float64,
                                 outdoor::Bool)
    w = infection_rate_F2F * duration * (distance_factor_per_m^distance)
    if outdoor
        w = w * outside_factor
    end
    return w
end

function convert_weight_to_prob(weight::Float64, inf::Float64,
                     susc::Float64; mask_infector::Bool = false,
                     mask_infectee::Bool = true)
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
    sg = LightGraphs.SimpleGraph(Ntot)
    Nj = length(N_per_role)
    sim = Dict("Ntot"=>Ntot, "N"=>N_per_role, "job"=>zeros(Int8, Ntot),
               "Njobs"=>Nj, "inf_time"=>(zeros(Int64, Ntot) .- 1),
                "symp_day"=>-ones(Int64, Ntot),
                "asymptomatic"=>zeros(Bool, Ntot),
                "would_isolate"=>zeros(Bool, Ntot),
                "isolation_time"=>zeros(Int64, Ntot),
                "isolation_status"=>zeros(Bool, Ntot),
                "isolation_status_true_test"=>zeros(Bool, Ntot),
                "isolation_status_false_test"=>zeros(Bool, Ntot),
                "isolation_status_contact_true_test"=>zeros(Bool, Ntot),
                "isolation_status_contact_false_test"=>zeros(Bool, Ntot),
                "isolating_due_to_true_test"=>zeros(Bool, Ntot),
                "isolating_due_to_false_test"=>zeros(Bool, Ntot),
                "isolating_due_to_contact_true_test"=>zeros(Bool, Ntot),
                "isolating_due_to_contact_false_test"=>zeros(Bool, Ntot),
                "susceptibility"=>ones(Float64,Ntot),
                "VL_mag"=>zeros(Float64, Ntot),
                "inf_mag"=>zeros(Float64, Ntot),
                "infection_profiles"=>Array{Array{Float64,1},1}(undef,0),
                "VL_profiles"=>Array{Array{Float64,1},1}(undef,0),
                "at_work"=>zeros(Bool, Ntot),
                "infection_status"=>zeros(Int8,Ntot),
                "days_infectious" => zeros(Int64,Ntot),
                "contact_network" => MetaGraphs.MetaGraph(sg),
                "job_sorted_nodes"=>Array{Array{Int64,1},1}(undef, Nj),
                "Incidence"=>Inc, "Prevalence"=>Prev)
    istart = 0
    for (i,n) in enumerate(N_per_role)
        nrh = (istart+1):(istart + n)
        sim["job"][nrh] .= i   #assign job/role labels
        sim["job_sorted_nodes"][i] = nrh
        #null initialise these
        push!.(Ref(sim["infection_profiles"]), fill(zeros(0),n))
        push!.(Ref(sim["VL_profiles"]), fill(zeros(0),n))
        istart += n
    end

    sim["would_isolate"][generate_isolations(sim["Ntot"], Pisol)] .= true   #generate propensity to isolate
    sim["susceptibility"][randsubseq(1:sim["Ntot"], 1.0 - Psusc)] .= 0.0    #generate non-susceptibles

    return sim
end

function generate_symp_day_from_time(stime::Float64)
    return Int64(floor(rand(0:1) + stime))
end

function infect_node!(sim::Dict, i::Int, time::Int,
                      contact_tracing_network::MetaGraphs.MetaGraph
                       = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"])))

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

    #if testing generate test positivity -- TODO
    if haskey(sim,"test_protocol")
       sim["test_pos_profiles"][i] = get_pos_profile(sim, i, sim["test_protocol"])
    end

    #Sort out isolation
    #Assume close contact isolation is guaranteed to isolate if infected node isolates (i.e. enforced)
    if will_isolate
        if sim["isolation_time"][i] > time + sim["symp_day"][i] ||
           sim["isolation_time"][i] == 0 #might already be isolating
               sim["isolation_time"][i] = time + sim["symp_day"][i]
        end
        #       print("Symp isol: ", i, ' ', time, ' ', sim["isolation_time"][i], '\n')
        traced = neighbors(contact_tracing_network, i)
        if length(traced) > 0
            to_update = ((sim["isolation_time"][traced] .>  (time + sim["symp_day"][i])) .+
                        (sim["isolation_time"][traced] .== 0))
            #if not isolating already or isolating later than this
            if sum(to_update) > 0
                sim["isolation_time"][traced[to_update]] .=  time + sim["symp_day"][i]
            end
        end

    end
end

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
    infect_node!(sim, inf, i_day)#, partner)

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
    # if length(inf) > 0
    #     print(i_day, ' ', inf, '\n')
    # end
    #people becoming Symptomatic (or asymptomatic)
    incend = inf[((i_day .- sim["inf_time"][inf]) .== sim["symp_day"][inf])]
    sim["infection_status"][incend] .= Symp
    #people recovered
    infend = inf[((i_day .- sim["inf_time"][inf]) .== sim["days_infectious"][inf])]
    sim["infection_status"][infend] .= Recd
end

function increment_isolations!(sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    #isolators
    cond0 = (sim["isolation_time"] .> 0)
    isolators = nr[cond0]
    #isolating today (over-cautious definition)
    cond1 = (sim["isolation_time"][isolators] .<= i_day) #isolation time today or before
    cond2 = (sim["isolation_status"][isolators] .== false)   #not currently isolating
    infisol = isolators[cond1 .* cond2]

    sim["isolation_status"][infisol] .= true

    is_true_test = sim["isolating_due_to_true_test"][infisol]
    sim["isolation_status_true_test"][infisol[is_true_test]] .= true
    sim["isolating_due_to_true_test"][infisol[is_true_test]] .= false #reset

    is_true_contact_test = sim["isolating_due_to_contact_true_test"][infisol]
    sim["isolation_status_contact_true_test"][infisol[is_true_contact_test]] .= true
    sim["isolating_due_to_contact_true_test"][infisol[is_true_contact_test]] .= false #reset

    is_false_test = sim["isolating_due_to_false_test"][infisol]
    sim["isolation_status_false_test"][infisol[is_false_test]] .= true
    sim["isolating_due_to_false_test"][infisol[is_false_test]] .= false #reset

    is_false_contact_test = sim["isolating_due_to_contact_false_test"][infisol]
    sim["isolation_status_contact_false_test"][infisol[is_false_contact_test]] .= true
    sim["isolating_due_to_contact_false_test"][infisol[is_false_contact_test]] .= false #reset

    #people leaving isolation
    isols = nr[sim["isolation_status"] .== true]
    isol_leavers = isols[(i_day  .- sim["isolation_time"][isols]) .== isol_time ]
    sim["isolation_status"][isol_leavers] .= false
    #reset flags
    sim["isolation_time"][isol_leavers] .= 0   #can isolate again, even though recovered
    sim["isolation_status_true_test"][isol_leavers] .= false
    sim["isolation_status_false_test"][isol_leavers] .= false
    sim["isolation_status_contact_true_test"][isol_leavers] .= false
    sim["isolation_status_contact_false_test"][isol_leavers] .= false
end

"""

"""
function get_introductions(sim::Dict, i_day::Int)
    #doesn't matter if individual is in work or not
    newinfs = randsubseq(1:sim["Ntot"],sim["Incidence"][i_day])
    NI = length(newinfs)
    infpairs = Array{Int64,2}(undef,3,0)
    if NI > 0
        infpairs = [transpose(zeros(Int64,NI)); transpose(newinfs);
                       introduction*transpose(ones(Int64,NI))]
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
function update_testing_state!(sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    testing_paused = nr[sim["testing_paused"] .== true]
    resume_testing = testing_paused[i_day .> sim["resume_testing"][testing_paused]]
    sim["testing_paused"][resume_testing] .= false
end

"""

"""
function do_testing!(sim::Dict, testing_params::Dict, i_day::Int,
          isolation_network::MetaGraphs.MetaGraph)

    #do testing if test day
    if i_day == sim["test_days"][sim["test_day_counter"]]
        #get all non susceptibles
        nr = 1:sim["Ntot"]
        update_testing_state!(sim, i_day)
        bool_can_be_tested = (sim["testing_paused"] .== false)

        bool_exposed = (sim["infection_status"] .!= Susc)
        exposed = nr[bool_exposed]
        should_be_positive = exposed[length.(sim["test_pos_profiles"][exposed]) .>
                                       i_day .- sim["inf_time"][exposed]]
        pos_tests = zeros(Bool,sim["Ntot"])
        #false positives due to specificity
        false_pos = randsubseq(nr[bool_can_be_tested],1-testing_params["specificity"])
        # print("Should test pos: ", should_be_positive,'\n')
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
        # print("Pos tests: ",nr[pos_tests],'\n')

        isolators = apply_positive_tests!(sim, nr[pos_tests], testing_params, i_day, isolation_network)

        sim["test_day_counter"] += 1

        return isolators
    else
        return Array{Int64,1}(undef,0)
    end
end

"""

"""
function apply_positive_tests!(sim::Dict, detected::Array{Int64,1}, testing_params::Dict, i_day::Int,
                               isolation_network::MetaGraphs.MetaGraph)
    will_isolate = detected[(sim["will_isolate_with_test"][detected])]
    isol_day = Int64.(round.(i_day + testing_params["delay"] .+ rand([0,0.5],length(will_isolate))))
    bool_update_isol_time = Bool.(((sim["isolation_time"][will_isolate] .== 0) .+
                                   (sim["isolation_time"][will_isolate] .> isol_day)))
    to_update = will_isolate[bool_update_isol_time]
    sim["isolation_time"][to_update] .= isol_day[bool_update_isol_time]

    bool_true_test = Bool.(sim["infection_status"][to_update] .== Susc .+
                           sim["infection_status"][to_update] .== Recd)
    sim["isolating_due_to_true_test"][bool_true_test] .= true
    sim["isolating_due_to_false_test"][.!bool_true_test] .= true

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

    sim["testing_paused"][will_isolate] .= true
    sim["resume_testing"][will_isolate] .= i_day + testing_params["test_pause"]

    return will_isolate
end

"""

"""
# function print_infection_network(sim::Dict, fname::String, infpairs::Array{Int64,2},
#                                  x_pos::Array{Float64,1}, y_pos::Array{Float64,1},
#                                  pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))

#     lg = LightGraphs.SimpleGraph(sim["Ntot"])
#     if haskey(sim, "social_graph")
#         ges = collect(Graphs.edges(sim["social_graph"]))
#         LightGraphs.add_edge!.(Ref(lg), Graphs.source.(ges), Graphs.target.(ges))
#     end
#     # if length(pairs) > 0
#     #     for i in 1:size(pairs,2)
#     #         if !LightGraphs.has_edge(lg,pairs[1,i],pairs[2,i])
#     #             LightGraphs.add_edge!(lg,pairs[1,i],pairs[2,i])
#     #         end
#     #     end
#     # end
#     for i in 1:size(infpairs,2)
#         if !LightGraphs.has_edge(lg,infpairs[1,i],infpairs[2,i])
#             LightGraphs.add_edge!(lg,infpairs[1,i],infpairs[2,i])
#         end
#     end

#     mg = MetaGraphs.MetaGraph(lg)
#     if haskey(sim, "social_graph")
#         for e in ges
#             set_prop!(mg, LightGraphs.Edge(Graphs.source(e), Graphs.target(e)), :color, length(edge_colours))
#         end
#     end
#     # if length(pairs) > 0
#     #     for i in 1:size(pairs,2)
#     #         set_prop!(mg, LightGraphs.Edge(pairs[1,i], pairs[2,i]), :color, 2)
#     #     end
#     # end
#     for i in 1:size(infpairs,2)
#         if infpairs[1,i] > 0 &&  infpairs[2,i] > 0
#             set_prop!(mg, LightGraphs.Edge(infpairs[1,i], infpairs[2,i]), :color, infpairs[3,i])
#         end
#     end
#     ecolors = zeros(Int8,ne(lg))
#     for (j,e) in enumerate(LightGraphs.edges(lg))
#         ecolors[j] = get_prop(mg,e,:color)
#     end
#     #inflabels = 2 .* ones(Int8,sim["Ntot"])
#     #inflabels[sim["infection_status"] .== Susc] .= 1
#     #inflabels[sim["infection_status"] .== Recd] .= 3
#     #here
#     inflabels = inf_ref .* ones(Int8,sim["Ntot"])
#     inflabels[sim["infection_status"] .== Susc] .= susceptible_ref
#     inflabels[sim["infection_status"] .== Recd] .= recovered_ref
#     inflabels[.!sim["at_work"]] .= not_at_work_ref
#     inflabels[sim["isolation_status"]] .= isolating_ref


#     draw(PNG(fname,19cm,19cm,dpi=150),gplot(lg, x_pos, y_pos,
#                 nodefillc=inf_colours[inflabels],
#                 #nodestrokec=inf_colours[infstrokes],
#                 nodestrokelw=1.0,
#                 nodelabel=job_labels[sim["job"]],
#                 NODESIZE = 0.3/sqrt(sim["Ntot"]),
#                 NODELABELSIZE = 3.0,
#                 #nodelabel=1:sim["Ntot"],
#                 #edgelabel=1:ne(lg),
#                 edgestrokec=edge_colours[ecolors]))
# end

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
    sim_summary["IndexCaseInfections"] = zeros(Int64,Ndays)
    sim_summary["IndexCase"] = index_case
    sim_summary["IndexCaseInfections"] = 0
    sim_summary["IndexCaseInfectivity"] = sim["inf_mag"][index_case]
    sim_summary["IndexCasePeakVL"] = sim["VL_mag"][index_case]
    for j in 1:sim["Njobs"]
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, i_day)
    return sim_summary
end

"""

"""
function scenario_sim_setup!(sim::Dict, inc::Array{Float64,1}, prev::Array{Float64,1}, i_day::Int64, Ndays::Int64)

    sim_summary = basic_sim_setup(sim, i_day, Ndays)
    sim_summary["Incidence"] = inc
    sim_summary["Prevalence"] = prev
    for j in 1:sim["Njobs"]
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, i_day)
    return sim_summary
end

"""

"""
function update_sim_summary!(summary::Dict, sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    new_isolator_bool = sim["isolation_status"] .* (sim["isolation_time"] .== i_day)
    isol_no_test = (sim["isolation_status"] .- sim["isolation_status_true_test"] .-
                   sim["isolation_status_false_test"] .- sim["isolation_status_contact_true_test"] .-
                   sim["isolation_status_contact_false_test"])
    after_onset_bool = (sim["infection_status"] .== Symp)
    at_work_bool = (sim["at_work"] .== true)
    for j in 1:sim["Njobs"]
        jt = (sim["job"] .== j)

        summary["Susceptible"][j,i_day] = sum(jt .* (sim["infection_status"] .== Susc))
        summary["Recovered"][j,i_day] = sum(jt .* (sim["infection_status"] .== Recd))
        summary["Infectious"][j,i_day] = (sum(jt) - summary["Susceptible"][j,i_day]
                                         - summary["Recovered"][j,i_day])
        summary["Exposed"][j,i_day] = sum(jt .* (sim["infection_status"] .== Expd))
        summary["Isolated"][j,i_day] = sum(jt .* sim["isolation_status"])
        summary["NewIsolators"][j,i_day] = sum(jt .* new_isolator_bool)


        summary["IsolatedDueToSymptoms"][j,i_day] = sum(jt .* isol_no_test)
        summary["IsolatedDueToTestAsymp"][j,i_day] = sum(jt .* sim["isolation_status_true_test"] .* sim["asymptomatic"])
        summary["IsolatedDueToTestSymp"][j,i_day] = sum(jt .* sim["isolation_status_true_test"] .* .!(sim["asymptomatic"]))
        summary["IsolatedDueToFalsePos"][j,i_day] = sum(jt .* sim["isolation_status_false_test"])
        summary["IsolatedDueToContactTruePos"][j,i_day] = sum(jt .* sim["isolation_status_contact_true_test"])
        summary["IsolatedDueToContactFalsePos"][j,i_day] = sum(jt .* sim["isolation_status_contact_false_test"])

        summary["NewSympIsolators"][j,i_day] = sum(jt .* isol_no_test .* new_isolator_bool)
        summary["NewTestAsympIsolators"][j,i_day] = sum(jt .* sim["isolation_status_true_test"] .*
                                                new_isolator_bool .* sim["asymptomatic"])
        summary["NewTestSympIsolators"][j,i_day] = sum(jt .* sim["isolation_status_true_test"] .*
                                                new_isolator_bool .* .!(sim["asymptomatic"]))
        summary["NewFalseIsolators"][j,i_day] = sum(jt .* sim["isolation_status_false_test"] .* new_isolator_bool)
        summary["NewContactTruePosIsolators"][j,i_day] = sum(jt .* sim["isolation_status_contact_true_test"] .*
                                                         new_isolator_bool)
        summary["NewContactFalsePosIsolators"][j,i_day] = sum(jt .* sim["isolation_status_contact_false_test"] .*
                                                         new_isolator_bool)

        summary["Presenting"][j,i_day] = sum(jt .* after_onset_bool .* at_work_bool)
        summary["Asymptomatic"][j,i_day] = sum(jt .* after_onset_bool .* sim["asymptomatic"])
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
                                   OccPerDay::Array{Float64,1})
    i_day = rand(1:7)
    if TestParams["is_testing"]
        if (TestParams["protocol"] == PCR_mass_protocol
         || TestParams["protocol"] == LFD_mass_protocol)
             sim["test_days"], sim["test_day_counter"] =
                     init_testing!(sim, TestParams, i_day, length(OccPerDay), false)
        end  #add options for other protocols here
    else
        sim["test_day_counter"] = 1
        sim["test_days"] = [0]
    end
    Anyinf = true
    if Params["SimType"] == Outbreak_sim
        sim_summary = sim_setup!(sim, Params["InfInit"], i_day, length(OccPerDay))
        Anyinf = any(((sim["infection_status"] .== Susc)
                      .+ (sim["infection_status"] .== Recd)) .== 0)
    elseif Params["SimType"] == Scenario_sim
        sim_summary = scenario_sim_setup!(sim, Incidence, Prevalence, i_day, length(OccPerDay))
    end

    return sim_summary, i_day, Anyinf
end


function example_sim_loop!(sim::Dict, i_day::Int,  testing::Bool, TestParams::Dict,
                           isolation_network::MetaGraphs.MetaGraph, Prev::Float64)

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
