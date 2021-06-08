"""How to use the functions in this code




"""




include("viral_load_infectivity_testpos.jl")

using LightGraphs
using MetaGraphs





"""
Initialises: the model.

Args:
N_per_role::Array{Int64,1}: Array where each entry is the number of employees in a different job role.

contact_graph::MetaGraphs.MetaGraph: A graph of contacts (can be updated each day), with the following metadata on the edges:
        :weights, a float array of weights listing the duration*transmission rate for each type of contact that occured
        :types, an integer array (the same size as weights) indicating the type of contact associated with each

Pisol::Float64: Probability an individual will isolate

Psusc::Float64: Probability an individual is susceptible at simulation start

Returns: 

sim::Dict(): Framework for storing simulation data, to be passed to other functions
"""
function init_transmission_model(N_per_role::Array{Int64,1}, contact_graph::MetaGraphs.MetaGraph,
                                 Pisol::Float64, Psusc::Float64)
    sim = Dict("Ntot"=>sum(N_per_role), "job"=>zeros(Int8, Ntot), "inf_time"=>(zeros(Int64, Ntot) .- 1),
               "symp_time"=>-ones(Int64, Ntot), 
                "asymptomatic"=>zeros(Bool, Ntot),
                "would_isolate"=>zeros(Bool, Ntot),
                "isolation_time"=>zeros(Int64, Ntot),
                "isolation_status"=>zeros(Bool, Ntot),
                "susceptibility"=>ones(Float64,Ntot),
                "inf_mag"=>zeros(Float64, Ntot),
                "infection_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
                "VL_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
                "at_work"=>zeros(Bool, Ntot),
                "infection_status"=>zeros(Int8,Ntot),
                "days_infectious" => zeros(Int64,Ntot),
                "fixed_contact_network" => contact_graph)
    istart = 0
    for (i,n) in enumerate(N_per_role)
        sim["job"][(istart+1):(istart + n)] .= i   #assign job/role labels
    end
    
    
    sim["would_isolate"][generate_isolations(sim["Ntot"], Pisol)] .= true   #generate propensity to isolate 
    sim["susceptibility"][randsubseq(1:sim["Ntot"], 1.0 - Psusc)] .= 0.0    #generate non-susceptibles
#     
#     build_viral_load_distributions!(sim)
#     sim["symp_day"] = Int64.(floor.(rand(0:1,sim["Ntot"]) .+ sim["symp_time"]))    #day isolation would start
#     nr = 1:sim["Ntot"]
#     sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]  #those without propensity to isolate
#     sim["will_isolate"][sim["asymptomatic"]] .= false   #asymptomatics don't self isolate because they lack symptoms, 
#                                                         #but are not "non isolators"
    return sim
end

function infect_node!(sim::Dict, i::Int, time::Int)
    sim["infection_status"][i] = Expd
    sim["susceptibility"][i] = 0.0
    sim["inf_time"][i] = time
    
    #Generate viral load
    sim["symp_time"][i], sim["VL_mag"][i], sim["asymptomatic"][i] .= 
            build_viral_load_distribution!(sim["VL_profiles"][i], sim["infection_profiles"][i])
    
    #if testing generate test positivity -- TODO
    
    
    #Sort out isolation
    #Assume close contact isolation is guaranteed to isolate if infected node isolates (i.e. enforced)
    if sim["will_isolate"][i]
        if sim["isolation_time"][i] > time + sim["symp_day"][i] ||
           sim["isolation_time"][i] == 0 #might already be isolating
               sim["isolation_time"][i] = time + sim["symp_day"][i]
        end
#       print("Symp isol: ", i, ' ', time, ' ', sim["isolation_time"][i], '\n')
        if length(sim["close_contacts_identified"][i]) > 0   #pair isolation
            to_update = ((sim["isolation_time"][sim["close_contacts_identified"][i]] .> 
                        (time + sim["symp_day"][i])) .+
                        (sim["isolation_time"][sim["close_contacts_identified"][i]] .== 0))   
            #if not isolating already or isolating later than this
            if sum(to_update) > 0
                sim["isolation_time"][sim["close_contacts_identified"][i][to_update]] .= 
                           time + sim["symp_day"][i]
            end
        end
    end
end

"""

"""
function infect_random!(sim::Dict, InfInit::Int, i_day::Int,
                        pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))
    #infect from non infectious
    nr = 1:sim["Ntot"]
    if InfInit > 0
        nrsub = nr[(sim["infection_status"] .== Susc) .* (sim["job"] .== InfInit)]
        else  #infinit 0 means any random person
        nrsub = nr[(sim["infection_status"] .== Susc)]
    end
    inf = rand(nrsub,1)[1]
    partner = 0
    if length(pairs) > 0   #for pair isolation
        partner = find_partner(pairs, inf)
    end
    infect_node!(sim, inf, i_day, partner)

    return inf
end


"""function to list infectious nodes and their current infectivity

Args:

sim::Dict simulation framework (returned by init_transmission_model)

Returns: 

inf::Array{Int64,1}: List of infectious nodes

inf_scales::Array{Float64,1}: Infectivity values for infectious nodes
"""
function get_infectivities(sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    infected = (sim["infection_status"] .!= Susc) .* (sim["infection_status"] .!= Recd)
    inf = nr[infected .* sim["at_work"]]
    inf_scales = zeros(length(inf))
    j = 1
    for i in inf
        inf_scales[j] = sim["infection_profiles"][i][i_day - sim["inf_time"][i] + 1]
        j += 1
    end
    return inf, inf_scales
end

function increment_infectivity!(sim::Dict, i_day::Int)
    #based on dates, stop infectivity
    nr = 1:sim["Ntot"]
    inf = nr[(sim["infection_status"] .!= Susc) .* (sim["infection_status"] .!= Recd)]
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

    is_true_partner_test = sim["isolating_due_to_partner_true_test"][infisol]
    sim["isolation_status_partner_true_test"][infisol[is_true_partner_test]] .= true
    sim["isolating_due_to_partner_true_test"][infisol[is_true_partner_test]] .= false #reset

    is_false_test = sim["isolating_due_to_false_test"][infisol]
    sim["isolation_status_false_test"][infisol[is_false_test]] .= true
    sim["isolating_due_to_false_test"][infisol[is_false_test]] .= false #reset

    is_false_partner_test = sim["isolating_due_to_partner_false_test"][infisol]
    sim["isolation_status_partner_false_test"][infisol[is_false_partner_test]] .= true
    sim["isolating_due_to_partner_false_test"][infisol[is_false_partner_test]] .= false #reset

    #people leaving isolation
    isols = nr[sim["isolation_status"] .== true]
    isol_leavers = isols[(i_day  .- sim["isolation_time"][isols]) .== isol_time ]
    sim["isolation_status"][isol_leavers] .= false
    #reset flags
    sim["isolation_time"][isol_leavers] .= 0   #can isolate again, even though recovered
    sim["isolation_status_true_test"][isol_leavers] .= false
    sim["isolation_status_false_test"][isol_leavers] .= false
    sim["isolation_status_partner_true_test"][isol_leavers] .= false
    sim["isolation_status_partner_false_test"][isol_leavers] .= false
end

"""

"""
function get_introductions(infpairs::Array{Int64,2}, sim::Dict, i_day::Int)
    #doesn't matter if individual is in work or not

    newinfs = randsubseq(1:sim["Ntot"],sim["Incidence"][i_day])
    NI = length(newinfs)
    if NI > 0
        infpairs = hcat(infpairs,[transpose(zeros(Int64,NI)); transpose(newinfs);
                       introduction*transpose(ones(Int64,NI))])
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
        eprob = 1 .- exp.(-j_inf_all .* beta .* j_susc_all)

        #draw which will be infectious (duplicates don't matter)
        #print(ecprob,"\n\n")
        eind = collect(1:length(eprob))
        einf = eind[rand(length(eprob)) .< (eprob)]

        if length(einf) > 0
            ipairs = hcat(ipairs, [nin_all[einf], nout_all[einf], etype[einf]])
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
            if sim["infection_status"][k] .== Susc  #if susceptible
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
                     workpairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))
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

    isolators = apply_positive_tests!(sim, nr[pos_tests], testing_params, i_day, workpairs)

    return isolators
end

"""

"""
function apply_positive_tests!(sim::Dict, detected::Array{Int64,1}, testing_params::Dict, i_day::Int,
                               pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))
    will_isolate = detected[(sim["will_isolate_with_test"][detected])]
    isol_day = Int64.(round.(i_day + testing_params["delay"] .+ rand([0,0.5],length(will_isolate))))
    bool_update_isol_time = Bool.(((sim["isolation_time"][will_isolate] .== 0) .+
                       (sim["isolation_time"][will_isolate] .> isol_day)))
    to_update = will_isolate[bool_update_isol_time]
    sim["isolation_time"][to_update] .= isol_day[bool_update_isol_time]


    #locate partners for those who will isolate
    partners = zeros(Int64, length(to_update))
    if length(pairs) > 0
        partners .= find_partner.(Ref(pairs),to_update)
    end

    # print(isol_day,'\n')
    # print(sim["isolation_time"][will_isolate],'\n')



    for (i,ip) in enumerate(to_update)
        #all partners of isolators do isolate
        if partners[i] > 0
            if sim["isolation_time"][partners[i]] == 0 ||
               sim["isolation_time"][partners[i]] > sim["isolation_time"][ip]
                sim["isolation_time"][partners[i]] = sim["isolation_time"][ip]
                if sim["infection_status"][ip] == Susc || sim["infection_status"][ip] == Recd
            sim["isolating_due_to_partner_false_test"][partners[i]] = true
        else
            sim["isolating_due_to_partner_true_test"][partners[i]] = true
        end
            end
        end
        if sim["infection_status"][ip] == Susc || sim["infection_status"][ip] == Recd
            sim["isolating_due_to_false_test"][ip] = true
        else
            sim["isolating_due_to_true_test"][ip] = true
        end
    end

#     if length(to_update) > 0
#         print(i_day, '\n')
#         print("Isolators: ", to_update, '\n')
#         print("Partners: ", partners, '\n')
#         print("Isol times: ", sim["isolation_time"][to_update], '\n')
#         print("Isol times p: ", sim["isolation_time"][partners[partners .> 0]], '\n')
#     end

    sim["testing_paused"][will_isolate] .= true
    sim["resume_testing"][will_isolate] .= i_day + testing_params["delay"] +
                                           testing_params["test_pause"]

    return will_isolate
end

"""

"""
function print_infection_network(sim::Dict, fname::String, infpairs::Array{Int64,2},
                                 x_pos::Array{Float64,1}, y_pos::Array{Float64,1},
                                 pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))

    lg = LightGraphs.SimpleGraph(sim["Ntot"])
    if haskey(sim, "social_graph")
        ges = collect(Graphs.edges(sim["social_graph"]))
        LightGraphs.add_edge!.(Ref(lg), Graphs.source.(ges), Graphs.target.(ges))
    end
    # if length(pairs) > 0
    #     for i in 1:size(pairs,2)
    #         if !LightGraphs.has_edge(lg,pairs[1,i],pairs[2,i])
    #             LightGraphs.add_edge!(lg,pairs[1,i],pairs[2,i])
    #         end
    #     end
    # end
    for i in 1:size(infpairs,2)
        if !LightGraphs.has_edge(lg,infpairs[1,i],infpairs[2,i])
            LightGraphs.add_edge!(lg,infpairs[1,i],infpairs[2,i])
        end
    end

    mg = MetaGraphs.MetaGraph(lg)
    if haskey(sim, "social_graph")
        for e in ges
            set_prop!(mg, LightGraphs.Edge(Graphs.source(e), Graphs.target(e)), :color, length(edge_colours))
        end
    end
    # if length(pairs) > 0
    #     for i in 1:size(pairs,2)
    #         set_prop!(mg, LightGraphs.Edge(pairs[1,i], pairs[2,i]), :color, 2)
    #     end
    # end
    for i in 1:size(infpairs,2)
        if infpairs[1,i] > 0 &&  infpairs[2,i] > 0
            set_prop!(mg, LightGraphs.Edge(infpairs[1,i], infpairs[2,i]), :color, infpairs[3,i])
        end
    end
    ecolors = zeros(Int8,ne(lg))
    for (j,e) in enumerate(LightGraphs.edges(lg))
        ecolors[j] = get_prop(mg,e,:color)
    end
    #inflabels = 2 .* ones(Int8,sim["Ntot"])
    #inflabels[sim["infection_status"] .== Susc] .= 1
    #inflabels[sim["infection_status"] .== Recd] .= 3
    #here
    inflabels = inf_ref .* ones(Int8,sim["Ntot"])
    inflabels[sim["infection_status"] .== Susc] .= susceptible_ref
    inflabels[sim["infection_status"] .== Recd] .= recovered_ref
    inflabels[.!sim["at_work"]] .= not_at_work_ref
    inflabels[sim["isolation_status"]] .= isolating_ref


    draw(PNG(fname,19cm,19cm,dpi=150),gplot(lg, x_pos, y_pos,
                nodefillc=inf_colours[inflabels],
                #nodestrokec=inf_colours[infstrokes],
                nodestrokelw=1.0,
                nodelabel=job_labels[sim["job"]],
                NODESIZE = 0.3/sqrt(sim["Ntot"]),
                NODELABELSIZE = 3.0,
                #nodelabel=1:sim["Ntot"],
                #edgelabel=1:ne(lg),
                edgestrokec=edge_colours[ecolors]))
end

"""

"""
function basic_sim_setup(sim::Dict, i_day::Int64, Ndays::Int64)
    return Dict{Any,Any}("time"=>collect(1:Ndays),
                "Susceptible"=>zeros(Int64,(3,Ndays)),
                "Infectious"=>zeros(Int64,(3,Ndays)),
                "Exposed"=>zeros(Int64,(3,Ndays)),
                "Isolated"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToSymptoms"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToTestAsymp"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToTestSymp"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToFalsePos"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToPartnerTruePos"=>zeros(Int64,(3,Ndays)),
                "IsolatedDueToPartnerFalsePos"=>zeros(Int64,(3,Ndays)),
                "NewIsolators"=>zeros(Int64,(3,Ndays)),
                "NewSympIsolators"=>zeros(Int64,(3,Ndays)),
                "NewTestAsympIsolators"=>zeros(Int64,(3,Ndays)),
                "NewTestSympIsolators"=>zeros(Int64,(3,Ndays)),
                "NewFalseIsolators"=>zeros(Int64,(3,Ndays)),
                "NewPartnerTruePosIsolators"=>zeros(Int64,(3,Ndays)),
                "NewPartnerFalsePosIsolators"=>zeros(Int64,(3,Ndays)),
                "Presenting"=>zeros(Int64,(3,Ndays)),
                "Asymptomatic"=>zeros(Int64,(3,Ndays)),
                "Recovered"=>zeros(Int64,(3,Ndays)),
                "PackagesInfectiousOnDelivery"=>zeros(Int64,Ndays),
                "FomiteInfs"=>zeros(Int64,(3,Ndays)),
                "CustomersInfectedByPkgs"=>zeros(Int64,Ndays),
                "CustomersInfectedByDrivers"=>zeros(Int64,Ndays),
                "NetworkInfs"=>zeros(Int64,(3,Ndays)),
                "ContactInfs"=>zeros(Int64,(3,Ndays)),
                "RoomInfs"=>zeros(Int64,(3,Ndays)),
                "PairInfs"=>zeros(Int64,(3,Ndays)),
                "CustomerIntroductions"=>zeros(Int64,(3,Ndays)),
                "ExternalIntroductions"=>zeros(Int64,(3,Ndays)))
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
    for j in 1:3
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, i_day)
    return sim_summary
end

"""

"""
function scenario_sim_setup!(sim::Dict, inc::Array{Float64,1}, prev::Array{Float64,1},
                             i_day::Int64, Ndays::Int64)

    sim_summary = basic_sim_setup(sim, i_day, Ndays)
    sim_summary["Incidence"] = inc
    sim_summary["Prevalence"] = prev
    for j in 1:3
        sim_summary["Susceptible"][j, 1:(i_day-1)] .= sim["N"][j]
    end
    update_sim_summary!(sim_summary, sim, i_day)
    return sim_summary
end

"""

"""
function trim_sim_summary!(sim_summary, i_day, Ndays)
    if i_day < Ndays
        for (key,value) in sim_summary
            if isa(value, Array)
                if length(size(value)) == 1
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
function simulate_day!(sim::Dict, i_day::Int)
    
    #add parts of previous function that go here
    
    
    infpairs = get_network_infections(sim::Dict, i_day::Int)
    infpairs_reduced = do_infections_randomly!(sim, infpairs)
    
end
