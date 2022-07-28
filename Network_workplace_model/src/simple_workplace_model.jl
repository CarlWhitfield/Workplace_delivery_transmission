include("transmission_model_framework.jl")

using Plots
using Cairo
using Compose
using Colors

using DataFrames
using CSV

using Distributions
using Random
using StatsBase

using Graphs
using MetaGraphs
using GraphPlot

using Printf
using Base.Threads

#infection type definitions
const Ninftypes = 3
const network_contact = 1
const non_network_contact = 2
const external_contact = 3
const edge_colours = [colorant"red",colorant"blue"]

#Contact duration parameters
const t_F2F_random = 1.0/4.0
const x_random = 1.5
# const job_labels = ["W"]
# const job_names = ["Workers"]
# const job_colours = [colorant"blue"]
const stat_labels = ["S","I","R"]
const stat_colours = [colorant"white", colorant"red",colorant"blue"]

const susceptible_ref = 1
const inf_ref = 2
const recovered_ref = 3
const not_at_work_ref = 4
const isolating_ref = 5
const inf_colours = [colorant"lightgrey",colorant"red",colorant"blue",
                     colorant"black",colorant"purple"]

DefaultTestParams = Dict("is_testing"=>false)

function generate_graph!(sim::Dict, degree_logmean::Float64,
                         degree_logstd::Float64)
    #generate graph with lognormal node distribution and weighted probability for node connections
    #generate degree distribution
    Tlognorm_pdf = truncated(LogNormal(degree_logmean, degree_logstd), 0, sim["Ntot"] - 1)
    node_degrees = Int64.(round.(rand(Tlognorm_pdf,sim["Ntot"])))

    #lognormal truncated at N-1
    #connect all loose edges
    G = SimpleGraph(sim["Ntot"])
    sim["contact_graph"] =  MetaGraphs.MetaGraph(G)
    nr = 1:sim["Ntot"]
    iorder = sortperm(node_degrees, rev = true)
    w = return_infection_weight(x_random, t_F2F_random, false, true)
    for i in iorder
        ijob = sim["job"][i]
        NDweights = copy(node_degrees)
        NDweights[i] = 0
        jvals = neighbors(sim["contact_graph"], i)
        NDweights[jvals] .= 0
        while node_degrees[i] > 0 && sum(NDweights .> 0) > 1
            j = sample(nr, Weights(NDweights))
            add_edge!(sim["contact_graph"],i,j)
            NDweights[j] = 0
            set_prop!(sim["contact_graph"],i,j,:weight,w)
            set_prop!(sim["contact_graph"],i,j,:type,network_contact)
            node_degrees[i] -= 1
            node_degrees[j] -= 1
        end
    end
end

#repeats not accepted for these contacts
function add_to_random_contact_network!(sim::Dict, source::Int64,
                               dest::Int64, weight::Float64)
    if has_edge(sim["rand_contact_network"],source,dest) == false
        add_edge!(sim["rand_contact_network"], source, dest)
        set_prop!(sim["rand_contact_network"], source, dest, :weight, weight)
    end
end

function generate_random_contact_network!(sim::Dict, i_day::Int)
    inf, inf_scales = get_infectivities(sim, i_day)
    #for those who are infectious and in work, generate contacts randomly with others
    #in work
    if length(inf) > 0
        nr = 1:sim["Ntot"]
        #nw = nr[sim["at_work"]]
        #nw in work
        #j0 = sim["job"][nw]          #j0 is jobs of in work
        w_rand = return_infection_weight(x_random, t_F2F_random, false, true)
        for (k, i) in enumerate(inf)
            j = sim["job"][i]
            contacts = Array{Int64,1}(undef,0)
            for j0 in 1:sim["Njobs"]
                p1 = sim["contact_prob_mat"][j,j0]
                nwh = sim["job_sorted_nodes"][j0]
                if j == j0
                    keep = ones(Bool,length(nwh))
                    keep[nwh .== i] .= false
                    nwj0 = nwh[keep]
                else
                    nwj0 = nwh[sim["at_work"][nwh]]
                end
                new_cs = randsubseq(nwj0,p1)
                contacts = push!(contacts,new_cs...)
            end

            add_to_random_contact_network!.(Ref(sim), Ref(i), contacts, Ref(w_rand))
        end
    end
end

function reset_daily_contact_networks!(sim::Dict)
    sim["rand_contact_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
end

function generate_random_absences(available::Array{Int64,1}, AbsRate::Float64)
    #random absences
    away = randsubseq(available, AbsRate)
    #return availables with aways filtered out
    return available[.!in.(available,Ref(away))]
end

function init(Params::Dict, Inc::Array{Float64,1},
              Prev::Array{Float64,1}; F2F_mod::Float64=1.0)
    sim = init_transmission_model(Params["N"], Params["Pisol"], Params["Psusc"], Inc, Prev)
    if haskey(Params, "contact_mat")
        sim["contact_prob_mat"] = Params["contact_mat"]
    else
        sim["contact_prob_mat"] = zeros(sim["Njobs"],sim["Njobs"])
    end
    
    return sim
end

function initialise(Params::Dict, Incidence::Array{Float64,1},
                    Prevalence::Array{Float64,1}, TransModifiers::Dict)
    sim = init(Params, Incidence, Prevalence; F2F_mod=TransModifiers["F2F_mod"])
    sim["Infections"] = zeros(Int,(3,0))
    sim["isolation_network"] = SimpleGraph(sim["Ntot"])
    return sim
end

function setup_wp_model!(sim::Dict, Params::Dict, TestParams::Dict, NDays::Int, 
                         ShiftPattern::Array{Bool})
    Params["NContactTypes"] = Ninftypes
    #assign shift pattern
    Npattern = size(ShiftPattern,1)
    if Params["random_start"]
        sim["start_day"] = rand(1:Npattern,sim["Ntot"])
    else
        sim["start_day"] = ones(Int,sim["Ntot"])
    end
    
    sim["shift_schedule"] = Array{Array{Bool,1},1}(undef,sim["Ntot"])
    if sim["Njobs"] > 1
        for j in 1:sim["Njobs"]
            shift_repeat = repeat(ShiftPattern[:,j],Int64(ceil(NDays/Npattern)+1))
            for i in sim["job_sorted_nodes"][j]
                sim["shift_schedule"][i] = shift_repeat[
                    sim["start_day"][i]:(sim["start_day"][i] + NDays - 1)]
            end
        end
    else
        shift_repeat = repeat(ShiftPattern,Int64(ceil(NDays/Npattern)+1))
        for i in 1:sim["Ntot"]
            sim["shift_schedule"][i] = shift_repeat[
                    sim["start_day"][i]:(sim["start_day"][i] + NDays - 1)]
        end
    end
    
    summary, i_day, Anyinf = setup_transmission_model!(sim, Params, TestParams, NDays)
    at_work = ones(sim["Ntot"])
    if haskey(Params, "degree_logmean") && haskey(Params, "degree_logstd")
        generate_graph!(sim, Params["degree_logmean"], Params["degree_logstd"])
        update_contact_network!(sim, sim["contact_graph"])
    #else
        #print("No fixed contact network created (parameters", 
        #      "\"degree_logmean\" and \"degree_logstd\" not given)\n")
    end
   
    
    return summary, i_day, Anyinf
end

function get_at_work(sim::Dict, AbsRate::Float64, i_day::Int)
    bool_sheduled_at_work = zeros(Bool,sim["Ntot"])
    nr = collect(1:sim["Ntot"])
    for i in nr
        if sim["shift_schedule"][i][i_day]
            bool_sheduled_at_work[i] = true
        end
    end
    at_work = generate_random_absences(nr[bool_sheduled_at_work], AbsRate)
    bool_at_work = zeros(Bool,sim["Ntot"])
    bool_at_work[at_work] .= true
    bool_at_work[sim["isolation_status"]] .= false
    
    return bool_at_work
end

function update_sim_summary_wp!(summary::Dict, sim::Dict, i_day::Int,
        inf_pairs::Array{Int64,2})
    update_sim_summary!(summary, sim, inf_pairs, i_day)
end

function add_edge_to_composite_network!(g::MetaGraphs.MetaGraph, snode::Int64, dnode::Int64,
                                        w::Float64, typ::Int64)
    if has_edge(g, snode, dnode)
        wvec = get_prop(g, snode, dnode, :weights)
        tvec = get_prop(g, snode, dnode, :types)
        push!(wvec, w)
        push!(tvec, typ)
    else
        wvec = fill(w,1)
        tvec = fill(typ,1)
        add_edge!(g, snode, dnode)
    end
    set_prop!(g, snode, dnode, :weights, wvec)
    set_prop!(g, snode, dnode, :types, tvec)
end

function collate_networks(sim::Dict)
    composite_graph = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    #add edges from each graph, with vector of weights and types
    nr = 1:sim["Ntot"]
    b_infected = get_infected_bool(sim)
    inf = nr[b_infected]
    nets = ["rand_contact_network", "contact_graph"]
    work_only = [true, true]
    indices = [non_network_contact, network_contact]

    for i in 1:length(nets)
        if haskey(sim, nets[i])
            #for e in edges(sim[nets[i]])
            #only add infectious node edges to save time
            nbrs = neighbors.(Ref(sim[nets[i]]), inf)
            #create array of infectious nodes for each contact edge
            n_in = vcat(fill.(inf,length.(nbrs))...)
            #flatten the other array to match
            n_out = vcat(nbrs...)
            w = get_prop.(Ref(sim[nets[i]]), n_in, n_out, Ref(:weight))
            if work_only[i]
                b1 = sim["at_work"][n_in]
                b2 = sim["at_work"][n_out]
                bcond = b1 .* b2
                add_edge_to_composite_network!.(Ref(composite_graph), n_in[bcond], n_out[bcond],
                    w[bcond], Ref(indices[i]))
            else
                add_edge_to_composite_network!.(Ref(composite_graph), n_in, n_out, w, Ref(indices[i]))
            end
        end
    end
    
    return composite_graph
end

function sim_loop_wp!(sim::Dict, sim_summary::Dict, i_day::Int, Params::Dict, TestParams::Dict, TransModifiers::Dict)
    
    reset_daily_contact_networks!(sim)
    
    #do_testing
    if any(TestParams["is_testing"])
        new_isolators = do_testing!(sim, TestParams, i_day, sim["isolation_network"])
    end

    #update infectivity and isolation status
    update_all_statuses!(sim, i_day)
    if haskey(Params, "absence_rate")
        at_work = get_at_work(sim, Params["absence_rate"], i_day)
    else
        at_work = get_at_work(sim, 0.0, i_day)
    end
    update_in_work!(sim, at_work)
    
    generate_random_contact_network!(sim, i_day)
    
    infpairs = Array{Int64,2}(undef,3,0)

    #introductions
    if haskey(Params,"intro_rate")
        intro_pairs = get_introductions(sim, i_day, external_contact; rel_rates = Params["intro_rate"])
    else
        intro_pairs = get_introductions(sim, i_day, external_contact)
    end
    g = collate_networks(sim)
    
    update_contact_network!(sim, g)
    infpairs = get_network_infections(sim, i_day)

    all_infpairs = hcat(intro_pairs, infpairs)

    infpairs_final = do_infections_randomly!(all_infpairs, sim, i_day)
    update_sim_summary_wp!(sim_summary, sim, i_day, infpairs_final)

    return infpairs_final
end

#run a either single outbreak with 1 initial case, or a sim with no cases but introductions
function run_sim_wp(Params::Dict, ShiftPattern::Array{Bool}, NDays::Int; 
        TestParams::Dict=DefaultTestParams,
        Incidence::Array{Float64,1} = zeros(NDays), 
        Prevalence::Array{Float64,1} = zeros(NDays))
    
    if length(TestParams) == 0
        TestParams = DefaultTestParams
    end

    #Transmission modifiers
    TransModifiers = Dict("F2F_mod"=>1.0)
    #modifies relative rate of F2F and aerosol transmission
    for key in keys(TransModifiers)  #if given in params, change them
        if haskey(Params,key)==true
            TransModifiers[key] = Params[key]   #add to Params Dict if not already there
        end
    end
    
    sim = initialise(Params, Incidence, Prevalence, TransModifiers)
    sim_summary, i_day, Go = setup_wp_model!(sim, Params, TestParams, NDays, ShiftPattern)
    
    
    #edit these based on params
    
    while Go && (i_day <= NDays)
        infpairs = sim_loop_wp!(sim, sim_summary, i_day,
                                Params, TestParams, TransModifiers)

        if Params["SimType"] == Outbreak_sim
            Go = any(get_infected_bool(sim))
        end
        
        Ninfs = size(infpairs,2)
        internal_infpairs = infpairs[:,infpairs[3,:].!=external_contact]
        if size(internal_infpairs,2) > 0
            days_since_inf = transpose(i_day .- sim["inf_time"][internal_infpairs[1,:]])
            sim["Infections"] = hcat(sim["Infections"],vcat(internal_infpairs[1:2,:], days_since_inf))
        end
        
        i_day += 1
    end
    trim_sim_summary!(sim_summary, i_day-1, NDays)
    
    return sim_summary, sim["Infections"]
end