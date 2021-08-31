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

using LightGraphs
using MetaGraphs
using GraphPlot

using Printf
using Base.Threads


#infection type definitions
const Ninftypes = 9   #should match the number of types listed below
const network_contact = 1
const non_network_contact = 2
const room_transmission = 3
const package_contact = 4
const customer_contact = 5
const introduction = 6
const pair_contact = 7
const car_share = 8
const house_share = 9
const edge_colours = [colorant"red", colorant"red", colorant"green", colorant"orange",
                      colorant"white",colorant"white",colorant"blue",colorant"lightgrey"]


#const per_person_office_volume = 20   #m^3
#const office_turnover = 4*24    #ventilation rate per day
#const cabin_volume = 5           #m^3
#const cabin_turnover_closed = 5*24    #ventilation rate per day
#const cabin_turnover_open = 50*24    #ventilation rate per day
#const breathing_rate = 0.72*24  #lung turnovers per day m^3
#const shift_pattern = 7     #days to draw shifts -- to be included

#Contact duration parameters
const t_F2F = 1.0          #Total duration of F2F contacts
const t_F2F_random = 1.0/6.0
const t_office = 7.0
const t_lunch = 0.75

#const t_breaks = 1/24       #Time in break/lunch areas (1 hr)
#around 10 hrs per day: 85 deliveries per day,
const ParcelTimesPerDelivery = Dict("t_cabin"=>1/(12),"t_doorstep"=>1/(120), "t_handling"=>1/(60),
                                "t_picking"=>1/(60))
const BulkTimesPerDelivery = Dict("t_cabin"=>1/(12),"t_doorstep"=>1/(60), "t_handling"=>1/(12), "t_picking"=>1/(12))
const job_labels = ["D","L","O"]
const job_names = ["Drivers","Pickers","Other"]
const job_colours = [colorant"blue",colorant"orange",colorant"green"]
const stat_labels = ["S","I","R"]
const stat_colours = [colorant"white", colorant"red",colorant"blue"]

const susceptible_ref = 1
const inf_ref = 2
const recovered_ref = 3
const not_at_work_ref = 4
const isolating_ref = 5
const inf_colours = [colorant"lightgrey",colorant"red",colorant"blue",colorant"black",colorant"purple"]

DefaultPkgParams = Dict("p_fomite_trans"=>0.0)
DefaultPairParams = Dict("is_driver_pairs"=>false, "is_loader_pairs"=>false)
DefaultTestParams = Dict("is_testing"=>false)

"""
    add_team_to_graph!(G::MetaGraphs.MetaGraph, team::Array{Int64,1}, w::Float64)

Add a set of edges that fully connects a subset of nodes to the MetaGraph G

## Arguments

`G` = MetaGraph to be edited

`team` = List of nodes in G to connect

`w` = value for all edge weights

## See also

`generate_cohort_graph!`
"""

function add_cohort_to_graph!(G::MetaGraphs.MetaGraph, team::Array{Int64,1}, w::Float64)
    for i in 1:length(team)
        jadd = (i+1):length(team)
        add_edge!.(Ref(G), Ref(team[i]), team[jadd])
        set_prop!.(Ref(G), Ref(team[i]), team[jadd], Ref(:weight), Ref(w))
    end  
end


function get_team_edge_weight(team::Array{Int64,1}, modifiers::Dict)
    #assume contacts are one-to-one, total time is preserved
    edge_weight = return_infection_weight(modifiers["distance"],
       modifiers["rel_time"] * t_F2F / (length(team) - 1),
        modifiers["outside"])
    
    return edge_weight
end

#split drivers into teams who all have pre-shift meeting with one office person? Or one warehouse person?
#If pairs delivery, pairs are assigned within teams. 
#Could do same with warehouse people and office people.
#NPteams, NOteams
#Regular network is just team interactions
#Random interactions are just that, but with t_d and phi factors
#Could reassign daily by regenerating this graph, or keep fixed
"""
    generate_cohort_graph!(sim::Dict, Nteams::Array{Int64,1}, TeamF2FTime::Array{Float64,1},
                               outside::Array{Bool,1}, distance::Array{Float64})

## Arguments

`sim` = Simulation dictionary

`Nteams` = Number of teams for each job role

`TeamF2FTime` = Time (for each job role) that a team spends in contact (in total)

`outside` = Bool (for each job role) indicating whether contacts occur outside

`distance` = Distance (for each job role) between contacts
"""
function generate_cohort_graph!(sim::Dict, Nteams::Array{Int64,1}, TeamF2FTime::Array{Float64,1},
                               outside::Array{Bool,1}, distance::Array{Float64,1})
    team_assign = Array{Array{Int64,1},1}(undef, length(Nteams))
    teams = Array{Array{Int64,1},1}(undef, sum(Nteams))
    Nstart = 0
    Ntstart = 0
    for j in 1:length(sim["N"])
        team_assign[j] = rand(1:Nteams[j], sim["N"][j])
        nr = Nstart .+ collect(1:sim["N"][j])
        
        #sort out fixed pairings, if they exist
        if haskey(sim,"fixed_job_pairings")
            fj_pairs_bool = (sim["fixed_job_pairings"][1,:] .> Nstart) .* 
                   (sim["fixed_job_pairings"][1,:] .<= Nstart + sim["N"][j])
            fj_pairs = sim["fixed_job_pairings"][:,fj_pairs_bool]
            if length(fj_pairs) > 0  #fixed job pairings exist
                for k in 1:size(fj_pairs,2) #make sure pairs go in the same cohort
                    pair_teams = zeros(Int64,2)
                    for m in 1:2    #get team no.s for each pair
                        pair_teams[m] = team_assign[j][fj_pairs[m,k] - Nstart]
                    end
                    team_no = rand(pair_teams)   #pick one at random
                    for m in 1:2    #assign to both
                        team_assign[j][fj_pairs[m,k] - Nstart] = team_no
                    end

                end
            end
        end
        
        #assign teams
        for n in 1:Nteams[j]      
            teams[Ntstart + n] = nr[team_assign[j] .== n]   #each team has a new entry
        end 
        
        #assign a member of warehouse staff to each driver team
        if j == 1
            nw = sim["N"][j] .+ collect(1:sim["N"][j+1])
            nwm = sample(nw,Nteams[j],replace=false) 
            for n in 1:Nteams[j]
                append!(teams[n],nwm[n])
            end
        end
        
        Nstart += sim["N"][j]
        Ntstart += Nteams[j]
    end
    
    graph = LightGraphs.SimpleGraph(sim["Ntot"])
    team_graph = MetaGraphs.MetaGraph(graph)
    Ntstart = 0
    for j in 1:length(sim["N"])
        for n in 1:Nteams[j]
            w = get_team_edge_weight(teams[Ntstart+n], Dict("outside"=>outside[j],
                           "distance"=>distance[j],"rel_time"=>TeamF2FTime[j]))
            add_cohort_to_graph!(team_graph, teams[Ntstart+n], w)
        end
        Ntstart += Nteams[j]
    end
    
    sim["cohort_network"] = team_graph
end


"""

"""
function generate_car_share_and_house_share_graphs!(sim::Dict, HomeParams::Dict)
    sim["car_share_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    sim["house_share_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    w_house_share = 0.0
    w_car_share = 0.0
    
    nr = collect(1:sim["Ntot"])
    nr_rand = shuffle(nr)
    pph = ones(sim["Ntot"])
    cumul_pph = cumsum(pph)
    #HomeParams["HouseShareFactor"]   #1 - mean number of employees per house
    if (HomeParams["HouseShareFactor"] > 0)
        Nhouses = Int64(ceil(sim["Ntot"]/(1 + HomeParams["HouseShareFactor"])))
        pph = 1 .+ rand(Multinomial(sim["Ntot"] - Nhouses, ones(Nhouses)./Nhouses))
        cumul_pph = cumsum(pph)
        for h in 1:Nhouses
            if (pph[h] > 1)
                if(h>1)
                    iph = (cumul_pph[h-1] + 1):cumul_pph[h]
                else
                    iph = 1:cumul_pph[h]
                end
                add_cohort_to_graph!(sim["house_share_network"], nr_rand[iph], w_house_share)
            end 
        end
    end
    
    #fill vector of hhs
    hhs = Array{Array{Int64,1},1}(undef,0)
    for h in 1:length(pph)
        if(h>1)
            iph = (cumul_pph[h-1] + 1):cumul_pph[h]
        else 
            iph = 1:cumul_pph[h]
        end
        push!(hhs, nr_rand[iph])
    end
    
    #HomeParams["CarShareFactor"]   #no. of houses / no. of cars (assume all house shares share vehicle)
    NH = length(pph)
    nh = collect(1:NH)
    nh_rand = shuffle(nh)
    hpc = ones(NH)
    cumul_hpc = cumsum(hpc)
    cars = Array{Array{Int64,1},1}(undef,0)
    if (HomeParams["CarShareFactor"] > 0)
        nh = collect(1:NH)
        Ncars = Int64(ceil(NH/(1 + HomeParams["CarShareFactor"])))
        hpc = 1 .+ rand(Multinomial(NH - Ncars, ones(Ncars) ./ Ncars))
        cumul_hpc = cumsum(hpc)
        for c in 1:Ncars
            car_share = Array{Int64,1}(undef,0)
            if(c>1)
                ipc = (cumul_hpc[c-1] + 1):cumul_hpc[c]
            else
                ipc = 1:cumul_hpc[c]
            end
            for h in nh_rand[ipc]
                push!.(Ref(car_share), hhs[h])
            end
            push!(cars, car_share)
            if length(car_share) > 1
                add_cohort_to_graph!(sim["car_share_network"], car_share, w_car_share)
            end 
        end
    end
    
    return hhs, cars
end


"""
    apply_contact_mixing_params!(sim::Dict, Params::Dict)


"""
function apply_contact_mixing_params!(sim::Dict, Params::Dict)
    for i = 1:3, j=1:3
        if i == 1 || j == 1
            sim["contact_prob_mat"][i,j] = Params["tD"]*sim["contact_prob_mat"][i,j]
        end
        if i != j
            sim["contact_prob_mat"][i,j] = Params["phi"]*sim["contact_prob_mat"][i,j]
        end
    end
    norm = Params["tD"]*sim["N"][1]*((sim["N"][1]-1)/2 + Params["phi"]*
           (sim["N"][2] + sim["N"][3])) + sim["N"][2]*(sim["N"][2]-1)/2 +
           sim["N"][3]*(sim["N"][3]-1)/2 + Params["phi"]*sim["N"][2]*sim["N"][3]
    coeff = 0.5*sim["Ntot"]*(sim["Ntot"]-1)/norm
    sim["contact_prob_mat"] .*= coeff
end

"""

"""
function init(Params::Dict, Inc::Array{Float64,1},
        Prev::Array{Float64,1})
    sim = init_transmission_model([Params["ND"],Params["NL"],Params["NO"]],
            Params["Pisol"], Params["Psusc"], Inc, Prev)
    sim["contact_prob_mat"] = Params["p_contact"]*ones(3,3)
    HSparams = Dict("HouseShareFactor"=>0.0, "CarShareFactor"=>0.0)
    if haskey(Params, "HouseShareFactor")
        HSparams["HouseShareFactor"] = Params["HouseShareFactor"]
    end
    if haskey(Params, "CarShareFactor")
        HSparams["CarShareFactor"] = Params["CarShareFactor"]
    end
    if HSparams["HouseShareFactor"] > 0 || HSparams["CarShareFactor"] > 0
        generate_car_share_and_house_share_graphs!(sim, HSparams)
    end
    
    apply_contact_mixing_params!(sim, Params)
    
    return sim
end

function init_pairs!(sim::Dict, PairParams::Dict)
    sim["fixed_job_pairings"] = Array{Int64,2}(undef,2,0)
    if PairParams["is_driver_pairs"]
        sim["driver_pairs"] = Array{Int64,2}(undef,2,0)
        if PairParams["fixed_driver_pairs"]
            sim["fixed_job_pairings"] = hcat(sim["fixed_job_pairings"],
                                        [transpose(collect(1:2:sim["N"][1]));
                                         transpose(collect(2:2:sim["N"][1]))])
        end
    end
    if PairParams["is_loader_pairs"]
        sim["loader_pairs"] = Array{Int64,2}(undef,2,0)
        if PairParams["fixed_loader_pairs"]
            sim["fixed_job_pairings"] = hcat(sim["fixed_job_pairings"],
                    [transpose(collect((sim["N"][1]+1):2:(sim["N"][1]+sim["N"][2])));
                     transpose(collect((sim["N"][1]+2):2:(sim["N"][1]+sim["N"][2])))])
        end
    end
end

"""

"""
function initialise(Params::Dict, PairParams::Dict, Incidence::Array{Float64,1},
                    Prevalence::Array{Float64,1})
    sim = init(Params, Incidence, Prevalence)
    if PairParams["is_driver_pairs"] || PairParams["is_loader_pairs"]
        init_pairs!(sim,  PairParams)
        sim["contact_times"] = BulkTimesPerDelivery
    else
        sim["contact_times"] = ParcelTimesPerDelivery
    end
    if Params["is_cohorts"]
        Nteams = [Params["NDteams"],Params["NLteams"],Params["NOteams"]]
        generate_cohort_graph!(sim, Nteams, Params["TeamTimes"],
                               Params["TeamsOutside"], Params["TeamDistances"])
    end
    
    sim["NTypes"] = Ninftypes
    
    return sim
end

"""

"""
function get_assignments!(sim::Dict, occ::Float64, Ncons::Int64, PairParams::Dict)
    
    at_work = zeros(Bool,sim["Ntot"])
    NAs = zeros(Int64,sim["Ntot"])
    
    
    if PairParams["is_driver_pairs"]
        at_work_drivers, NAdrivers = get_pair_assignments!(sim, occ, Ncons, PairParams, 1)
    else
        at_work_drivers, NAdrivers = get_individual_assignments!(sim, occ, Ncons, 1)
    end
    at_work[at_work_drivers] .= true
    NAs[at_work_drivers] .= NAdrivers
    
    
    if PairParams["is_loader_pairs"]
        at_work_loaders, NAloaders = get_pair_assignments!(sim, occ, Ncons, PairParams, 2)
    else
        at_work_loaders, NAloaders = get_individual_assignments!(sim, occ, Ncons, 2)
    end
    at_work[at_work_loaders] .= true
    NAs[at_work_loaders] .= NAloaders
    
    
    at_work_office, NAoffice = get_individual_assignments!(sim, occ, Ncons, 3)
    at_work[at_work_office] .= true
    NAs[at_work_office] .= NAoffice

    
    return at_work, NAs
end

"""

"""
function get_individual_assignments!(sim::Dict, occ::Float64, Ncons::Int64, job::Int)
    #without pairs
    in_job = sim["job_sorted_nodes"][job]
    Nj = length(in_job)
    
    available = in_job[sim["isolation_status"][in_job] .== false]
    if length(available)/Nj > occ
        w = sample(in_job, round(Int,occ*Nj), replace=false)
    else
        w = available
    end
    NW = length(w)
    NAs = zeros(Int64, NW)
    if job < 3
        nw = 1:NW
        NAs[nw] .= rand(Multinomial(Ncons,length(w)))
    end

    return w, NAs
end

"""

"""
function select_pairs!(sim::Dict, occ::Float64, fixed_pairs::Bool, job::Int64, Ncons::Int64)
    TotPairs = Int64(sim["N"][job]/2)
    jobgroup = sim["job_sorted_nodes"][job]
    is_available = .!(sim["isolation_status"][jobgroup])
    available = jobgroup[is_available]
    NJ = min(Int64(floor(0.5*length(available))*2), Int64(round(0.5*occ*sim["N"][job])*2))
    NP = Int64(NJ/2)
    pairs = zeros(Int64,(2,NP))
    if fixed_pairs
        p1 = 1:2:sim["N"][job]
        p2 = 2:2:sim["N"][job]
        is_pair_available = is_available[p1] .* is_available[p2]
        available_unpaired = jobgroup[vcat(p1[is_available[p1] .* (.!is_available[p2])],
                            p2[is_available[p2] .* (.!is_available[p1])])]
        available_pairs = (1:TotPairs)[is_pair_available]
        pairs_nos = vcat(transpose(jobgroup[p1]), transpose(jobgroup[p2]))
        apl = length(available_pairs)
#         if length(available_unpaired) > 0
#             print("Available unpaired: ", available_unpaired, '\n')
#         end
        if apl < NP
            pairs[1,1:apl] = pairs_nos[1,available_pairs]
            pairs[2,1:apl] = pairs_nos[2,available_pairs]
            if length(available_unpaired) > 2*(NP - apl)
                nos = sample(available_unpaired, 2*(NP - apl), replace = false)
                pairs[1,(apl+1):NP] = nos[1:(NP-apl)]
                pairs[2,(apl+1):NP] = nos[(NP-apl+1):(2*(NP - apl))]
            elseif length(available_unpaired) > 1
                #draw as many random pairs as are left, discard rest
                NPleft = Int64(floor(length(available_unpaired)/2))
                nos = sample(available_unpaired, 2*NPleft, replace = false)
                pairs[1,(apl+1):(apl+NPleft)] = nos[1:NPleft]
                pairs[2,(apl+1):(apl+NPleft)] = nos[(NPleft+1):(2*NPleft)]
                trimpairs = pairs[:, 1:(apl + NPleft)]
                pairs = pairs
            end
        else
            pair_nos = sample(available_pairs, NP, replace=false)
            pairs[1,:] = pairs_nos[1,pair_nos]
            pairs[2,:] = pairs_nos[2,pair_nos]
        end
    else
        nos = sample(available, NJ, replace=false)
        pairs[1,:] = nos[1:NP]
        pairs[2,:] = nos[(NP+1):NJ]
    end
    if size(pairs,2) > 0
        pair_cons = rand(Multinomial(Ncons,size(pairs,2)))
    else
        pair_cons = Array{Int64,1}(undef,0)
    end
    
    simple_pair_graph = MetaGraphs.MetaGraph(LightGraphs.SimpleGraph(sim["Ntot"]))
    NPfinal = size(pairs,2)
    pairs_out = Array{Array{Int64,1},1}(undef,NPfinal)
    for j in 1:size(pairs,2)
       add_edge!(simple_pair_graph, pairs[1,j], pairs[2,j]) 
       pairs_out[j] = pairs[:,j]
    end
    if job == 1
        sim["driver_pair_network"] = simple_pair_graph
    elseif job == 2
        sim["loader_pair_network"] = simple_pair_graph
    end
    
    return pairs_out, pair_cons
end

"""


"""
function get_pair_assignments!(sim::Dict, occ::Float64, Ncons::Int64, PairParams::Dict, job::Int)
    pairs = Array{Int64,2}(undef,2,0)
    NPassignments = Array{Int64,1}(undef,0)
    
    if job == 1
        pairs, NPassignments = select_pairs!(sim, occ, PairParams["fixed_driver_pairs"], job, Ncons)
        
        #component for delivery time
        exp_per_assign = return_infection_weight(1.0, 
            sim["contact_times"]["t_handling"] + sim["contact_times"]["t_doorstep"],
            true)
        
        #add component for cabin time
        exp_per_assign += return_infection_weight(1.0, sim["contact_times"]["t_cabin"], PairParams["is_window_open"])
        add_cohort_to_graph!.(Ref(sim["driver_pair_network"]), pairs, NPassignments .* exp_per_assign)
    end
    
    if job == 2
        pairs, NPassignments = select_pairs!(sim, occ, PairParams["fixed_loader_pairs"], job, Ncons)
        exp_per_assign = return_infection_weight(1.0, 
            sim["contact_times"]["t_picking"], true)
        add_cohort_to_graph!.(Ref(sim["loader_pair_network"]), pairs, NPassignments .* exp_per_assign)
    end
    
    
    return vcat(pairs...), vec(vcat(transpose(NPassignments),
            transpose(NPassignments)))
end



"""

"""
function get_in_work(sim::Dict, occ::Float64, Ncons::Int64,
                driver_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0),
                driver_cons::Array{Int64,1} = Array{Int64,1}(undef,0),
                loader_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0),
                loader_cons::Array{Int64,1} = Array{Int64,1}(undef,0))
    #wih driver pairs
    at_work = zeros(Bool,sim["Ntot"])
    NAs = zeros(Int64,sim["Ntot"])
    to_select = []
    if length(driver_pairs) > 0
        at_work[vec(driver_pairs)] .= true
        NAs[driver_pairs[1,:]] .= driver_cons
        NAs[driver_pairs[2,:]] .= driver_cons
    else
        push!(to_select,1)
    end
    if length(loader_pairs) > 0
        at_work[vec(loader_pairs)] .= true
        NAs[loader_pairs[1,:]] .= loader_cons
        NAs[loader_pairs[2,:]] .= loader_cons
    else
        push!(to_select,2)
    end
    push!(to_select,3)
    for j in to_select
        job = sim["job_sorted_nodes"][j]
        available = job[sim["isolation_status"][job] .== false]
        if length(available)/length(job) > occ
            w = sample(job, round(Int,occ*length(job)), replace=false)
        else
            w = available
        end
        if j < 3
            NAs[w] .= rand(Multinomial(Ncons,length(w)))
        end
        at_work[w] .= true
    end
    
    return at_work, NAs
end

# function get_infectivities(sim::Dict, i_day::Int)
#     nr = 1:sim["Ntot"]
#     infected = (sim["infection_status"] .!= Susc) .* (sim["infection_status"] .!= Recd)
#     inf = nr[infected .* sim["at_work"]]
#     inf_scales = zeros(length(inf))
#     j = 1
#     for i in inf
#         inf_scales[j] = sim["infection_profiles"][i][i_day - sim["inf_time"][i] + 1]
#         j += 1
#     end
#     return inf, inf_scales
# end

function reset_daily_contact_networks!(sim::Dict)
    sim["driver_pair_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    sim["loader_pair_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    sim["package_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    sim["rand_contact_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
    sim["room_contact_network"] = MetaGraphs.MetaGraph(SimpleGraph(sim["Ntot"]))
end

function add_to_pair_contacts_network!(sim::Dict, source::Int64, dest::Int64, weight::Float64)
    if sim["job"][source] == 1
        set_prop!(sim["driver_pair_network"], source, dest, :weight, weight)
    else
        set_prop!(sim["loader_pair_network"], source, dest, :weight, weight)
    end
end

#send list of packages that could be infected to this function
#with a list of source and destination nodes for each
#may contain repeats
function add_package_contacts_to_network!(sim::Dict, source::Int64, 
                                  dest::Int64, weight::Float64)
    if has_edge(sim["package_network"],source,dest)
        w = get_prop(sim["package_network"], source, dest, :weight)
        set_prop!(sim["package_network"], source, dest, :weight, w + weight)
    else
        add_edge!(sim["package_network"], source, dest)
        set_prop!(sim["package_network"], source, dest, :weight, weight)
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

function add_to_room_contact_network!(sim::Dict, source::Int64, 
                               dest::Int64, weight::Float64)
    if has_edge(sim["room_contact_network"],source,dest) == false
        add_edge!(sim["room_contact_network"], source, dest)
        set_prop!(sim["room_contact_network"], source, dest, :weight, weight)
    end
end

function generate_infected_packages_and_customers!(sim::Dict, NP::Int64,
        AllDrivers::Array{Array{Int64,1},1}, InfDrivers::Array{Int64,1},
        InfDriverPos::Array{Int64,1}, DriverInfs::Array{Float64,1},
        AllLoaders::Array{Array{Int64,1},1}, InfLoaders::Array{Int64,1}, 
        InfLoaderPos::Array{Int64,1}, LoaderInfs::Array{Float64,1},
        PkgParams::Dict, NADrivers::Array{Int64,1}, NALoaders::Array{Int64,1},
        modifiers::Dict)
    
    #assume package order is clustered, not random
    LAorder = cumsum(NALoaders)
    DAorder = cumsum(NADrivers)
    InfectorAtPickup = Array{Int64,1}(undef,0)
    LoadTimes = zeros(Float64,NP)
    DropTimes = zeros(Float64,NP)

    InfectorAtDropoff = Array{Int64,1}(undef,0)
    IDDriverPickup = Array{Int64,1}(undef,0)
    IDDriverDropoff = Array{Int64,1}(undef,0)

    InfProbAtDropoff = zeros(Float64,NP)
#     LoaderTeamAtPickup = vcat(fill.(collect(1:length(NALoaders)), NALoaders)...)
#     DriverTeamAtDropoff = vcat(fill.(collect(1:length(NADrivers)), NADrivers)...)
    AllLoadersAtPickup = vcat(fill.(AllLoaders, NALoaders)...)
    AllDriversAtDropoff = vcat(fill.(AllDrivers, NADrivers)...)
    
    #generate packages infected by loaders
    if PkgParams["p_fomite_trans"] > 0
        if length(InfLoaders) > 0
            #InfLoaders contains infectious loaders 
            #InfLoaderPos contains their position in AllLoaders which is an array of arrays 
            #NALoaders is the number of items assigned to each working group in AllLoaders
            for i in 1:length(InfLoaders)
                iLpkgs = NALoaders[InfLoaderPos[i]]
                PkgNos = collect(1:iLpkgs)
                if InfLoaderPos[i] > 1
                    PkgNos .= LAorder[InfLoaderPos[i]-1] .+ PkgNos 
                end
                LoadTimes[PkgNos] = PkgParams["Ltime"] .* rand(iLpkgs)  #generate time infected
                Lth = LoadTimes[PkgNos]
                #sterilisation time generated separately, if this loader deposits material, this is the time at which it will sterilise
                sth = Lth .+ rand(Exponential(PkgParams["PkgHlife"]/log(2)), iLpkgs)
                Iapcond = (sth .> PkgParams["Ltime"])  #check if could be infectious at pickup
                PDA = PkgNos[Iapcond]
                if sum(Iapcond) > 0
                    Nd = 1:length(AllDrivers)
                    dth = PkgParams["Ltime"] .+ PkgParams["Dtime"] * rand(sum(Iapcond))
                    Iadcond = (sth[Iapcond] .> dth)
                    for k in 1:length(PDA)  #for each package number generated, assign to driver based on cumsum
                        bool1 = (PDA[k] .> vcat([0],DAorder[1:(end-1)]))
                        bool2 = (PDA[k] .<= DAorder)
                        Dbool = bool1 .* bool2
                        for d in AllDrivers[Dbool][1]
                            InfectorAtPickup = push!(InfectorAtPickup,InfLoaders[i])  #push back number infectious at pickup
                            push!(IDDriverPickup, d)
                        end
                        if Iadcond[k]
                            InfProbAtDropoff[PDA[k]] += LoaderInfs[i]*PkgParams["p_fomite_trans"]
                            for d in AllDrivers[Dbool][1]
                                push!(InfectorAtDropoff,InfLoaders[i])  #push back number infectious at pickup
                                push!(IDDriverDropoff, d)
                            end
                        end
                    end
                    DropTimes[PDA] .= dth
                end
            end
            #print("Packages infected: ", PkgsInfected,"\n")
        end
    end

    #potential package infections by drivers here
    if length(InfDrivers) > 0
        #InfDrivers contains infectious drivers 
        #InfDriverPos contains their position in AllDrivers which is an array of arrays 
        #NADrivers is the number of items assigned to each working group in AllDrivers
        out_weight = return_infection_weight(1.0, 
                    sim["contact_times"]["t_doorstep"],true)
        in_weight = return_infection_weight(1.0, 
                    sim["contact_times"]["t_doorstep"],false)
        doorstep_weight = modifiers["outdoor_contact_frac"] * out_weight
         + (1-modifiers["outdoor_contact_frac"]) * in_weight
        
        for i in 1:length(InfDrivers)       #loop over all drivers to find infectees
            iDpkgs = NADrivers[InfDriverPos[i]]     #number of packages handled by driver
            PkgNos = collect(1:iDpkgs)
            if InfDriverPos[i] > 1
                PkgNos .= DAorder[InfDriverPos[i]-1] .+ PkgNos 
            end
            if PkgParams["p_fomite_trans"] > 0
                sth = PkgParams["Ltime"] .+ rand(Exponential(PkgParams["PkgHlife"]/log(2)), iDpkgs)
                dt_required = (DropTimes[PkgNos] .== 0) #check if dt needed
                DropTimes[PkgNos[dt_required]] = PkgParams["Ltime"] .+ rand(sum(dt_required)) .* PkgParams["Dtime"]
                dth = DropTimes[PkgNos]
                #could be infected at time of dropoff (only infects customers, not other drivers)
                InfProbAtDropoff[PkgNos] .= InfProbAtDropoff[PkgNos] .+ 
                                            (DriverInfs[i] * PkgParams["p_fomite_trans"])
                #could be infected at pickup and remain infected
                Iadcond = (sth .> dth)
                InfProbAtDropoff[PkgNos[Iadcond]] .= InfProbAtDropoff[PkgNos[Iadcond]] .+ 
                                            (DriverInfs[i] * PkgParams["p_fomite_trans"])

                for d in AllDrivers[InfDriverPos[i]]
                    if d != InfDrivers[i]
                        push!(InfectorAtDropoff, InfDrivers[i])  #push back number infectious at pickup
                        push!(IDDriverDropoff, d)
                    end
                    #add infections due to driver interaction
                    InfProbAtDropoff[PkgNos] .= InfProbAtDropoff[PkgNos] .+ 
                                               (DriverInfs[i] * doorstep_weight)
                end
            end
        end
    end
    
    add_package_contacts_to_network!.(Ref(sim), InfectorAtPickup,
               IDDriverPickup, Ref(PkgParams["p_fomite_trans"]))
    add_package_contacts_to_network!.(Ref(sim), InfectorAtDropoff,
               IDDriverDropoff, Ref(PkgParams["p_fomite_trans"]))
    
    #returns vectors of possibly infected customers, and prob of transmission for each
    #convert to probability
    
    InfProbAtDropoff = 1.0 .- exp.(-InfProbAtDropoff)
    
    return InfProbAtDropoff[InfProbAtDropoff .> 0], AllDriversAtDropoff[InfProbAtDropoff .> 0]
end

function get_package_and_customer_infections!(sim::Dict, NP::Int64,
           i_day::Int, PkgParams::Dict, PairParams::Dict,
           NAs::Array{Int64,1}, CustModifiers::Dict)
    
    CustInfProb = Array{Float64,1}(undef,0)
    DriverDelivering = Array{Array{Int64,1},1}(undef,0)
    if NP > 0
        inf, inf_scales = get_infectivities(sim, i_day)
        d_in_work = sim["job_sorted_nodes"][1][sim["at_work"][sim["job_sorted_nodes"][1]]]
        l_in_work = sim["job_sorted_nodes"][2][sim["at_work"][sim["job_sorted_nodes"][2]]]
        
        inf_drivers = inf[(sim["job"][inf] .== 1)]
        dinfs = inf_scales[(sim["job"][inf] .== 1)]
        inf_loaders = inf[(sim["job"][inf] .== 2)]
        linfs = inf_scales[(sim["job"][inf] .== 2)]
        
        infd_pos = zeros(Int64,length(inf_drivers))
        infl_pos = zeros(Int64,length(inf_loaders))
        
        NDPs = ne(sim["driver_pair_network"])
        if NDPs > 0
            alld = Array{Array{Int64,1},1}(undef,NDPs)
            NAD = zeros(Int64,NDPs)
            for (k,e) in enumerate(edges(sim["driver_pair_network"]))
                d1 = src(e)
                d2 = dst(e)
                alld[k] = [d1,d2]
                if d1 in inf_drivers
                    infd_pos[inf_drivers .== d1] .= k
                end
                if d2 in inf_drivers
                    infd_pos[inf_drivers .== d2] .= k
                end
                NAD[k] = NAs[d1]
            end
        else
            NDs = length(d_in_work)
            alld = Array{Array{Int64,1},1}(undef,NDs)
            ndw = 1:NDs
            for (k,d) in enumerate(d_in_work)
                alld[k] = [d]
            end
            for (k,d) in enumerate(inf_drivers)
                infd_pos[k] = ndw[d_in_work .== d][1]
            end
            NAD = NAs[d_in_work]
        end

        NLPs = ne(sim["loader_pair_network"])
        if NLPs > 0
            alll = Array{Array{Int64,1},1}(undef,NLPs)
            NAL = zeros(Int64,NLPs)
            for (k,e) in enumerate(edges(sim["loader_pair_network"]))
                l1 = src(e)
                l2 = dst(e)
                alll[k] = [l1,l2]
                if l1 in inf_loaders
                    infl_pos[inf_loaders .== l1] .= k
                end
                if l2 in inf_loaders
                    infl_pos[inf_loaders .== l2] .= k
                end
                NAL[k] = NAs[l1]
            end
        else
            NLs = length(l_in_work)
            alll = Array{Array{Int64,1},1}(undef,NLs)
            nlw = 1:NLs
            for (k,l) in enumerate(l_in_work)
                alll[k] = [l]
            end
            for (k,l) in enumerate(inf_loaders)
                infl_pos[k] = nlw[l_in_work .== l][1]
            end
            NAL = NAs[l_in_work]
        end
        
        
        CustInfProb, DriverDelivering = 
        generate_infected_packages_and_customers!(sim, NP, alld, inf_drivers,
            infd_pos, dinfs, alll, inf_loaders, infl_pos, linfs, PkgParams, NAD, NAL, CustModifiers)
    end
    
    return CustInfProb, DriverDelivering
end

function generate_random_contact_networks!(sim::Dict, Params::Dict, i_day::Int)
    inf, inf_scales = get_infectivities(sim, i_day)
    #for those who are infectious and in work, generate contacts randomly with others
    #in work
    ipairs = Array{Int64,2}(undef,3,0)
    if length(inf) > 0
        nr = 1:sim["Ntot"]
        #nw = nr[sim["at_work"]]
        #nw in work
        #j0 = sim["job"][nw]          #j0 is jobs of in work
        w_rand = return_infection_weight(1.0, t_F2F_random, false)
        wl_room = return_infection_weight(4.0, t_lunch, false)
        wo_room = return_infection_weight(4.0, t_office + t_lunch, false)
        for (k, i) in enumerate(inf)
            j = sim["job"][i]
            contacts = Array{Int64,1}(undef,0)
            for j0 in 1:3
                p1 = sim["contact_prob_mat"][j,j0]
                nwh = sim["job_sorted_nodes"][j0]
                nwj0 = nwh[sim["at_work"][nwh]]
                new_cs = randsubseq(nwj0,p1) 
                contacts = push!(contacts,new_cs...)
                if (j==2 && j0 == 2) || (j == 2 && j0 == 3) || (j == 3 && j0 == 2)
                    add_to_room_contact_network!.(Ref(sim), Ref(i), nwj0, Ref(wl_room))
                elseif (j == 3 && j0 == 3)
                    add_to_room_contact_network!.(Ref(sim), Ref(i), nwj0, Ref(wo_room))
                end
            end
            #contacts = nw[rand(length(nw)) .<  p1]
            
            add_to_random_contact_network!.(Ref(sim), Ref(i), contacts, Ref(w_rand))
        end
    end
end

function get_customer_introductions(sim::Dict, i_day::Int,
                  NAssignments::Array{Int64,1}, modifiers::Dict)
    pinf = zeros(Float64, sim["N"][1])
    
    out_weight = return_infection_weight(1.0, 
                    sim["contact_times"]["t_doorstep"],true)
    in_weight = return_infection_weight(1.0, 
                    sim["contact_times"]["t_doorstep"],true)
    
    ppd = modifiers["outdoor_contact_frac"] * out_weight
         + (1-modifiers["outdoor_contact_frac"]) * in_weight
    
    nds = sim["job_sorted_nodes"][1]
    pinf = convert_weight_to_prob.(sim["Prevalence"][i_day] * ppd * NAssignments[nds],
                 Ref(1.0), sim["susceptibility"][nds])
    
    intros = nds[rand(sim["N"][1]) .< pinf]
    
    return intros
end

function find_partner(pairs::Array{Int64,2}, p::Int64)
    c = 0
    if p in pairs
        if p in pairs[1,:]
            c = pairs[2,pairs[1,:].==p][1]
        else
            c = pairs[1,pairs[2,:].==p][1]
        end
    end
    return c
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
    nets = ["cohort_network", "driver_pair_network", "loader_pair_network", "package_network", 
            "room_contact_network", "rand_contact_network", "house_share_network", "car_share_network"]
    work_only = [true, true, true, true, true, true, false, true]
    indices = [network_contact, pair_contact, pair_contact, package_contact, room_transmission,
               non_network_contact, house_share, car_share]
    
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

function create_isolation_network!(sim::Dict, IsolParams::Dict)
    sim["isolation_network"] = SimpleGraph(sim["Ntot"])
    nets = ["cohort_network", "house_share_network", "car_share_network"]
    keys = ["cohort_isolation", "house_share_isolation", "car_share_isolation"]
    for (i,net_name) in enumerate(nets)
        if IsolParams[keys[i]]
            for e in edges(sim[net_name])
                add_edge!(sim["isolation_network"], src(e),dst(e))
            end
        end  
    end
    
    if IsolParams["fixed_pair_isolation"]
        for i in 1:size(sim["fixed_job_pairings"],2)
            add_edge!(sim["isolation_network"], sim["fixed_job_pairings"][1,i], 
                      sim["fixed_job_pairings"][2,i])
        end
    end
end

function setup_delivery_wp_model!(sim::Dict, Params::Dict, TestParams::Dict, OccPerDay::Array{Float64,1})
    summary, i_day, Anyinf = setup_transmission_model!(sim, Params, TestParams, OccPerDay)
    summary["CustomersInfected"] = zeros(Int64,length(OccPerDay))
    
    return summary, i_day, Anyinf
end
    
function update_sim_summary_delivery_wp!(summary::Dict, sim::Dict, i_day::Int,
        inf_pairs::Array{Int64,2}; CinfDrivers::Array{Array{Int64,1},1} =
        Array{Array{Int64,1},1}(undef,0))
    update_sim_summary!(summary, sim, i_day)
    inf, inf_scales = get_infectivities(sim, i_day)
    summary["CustomersInfected"][i_day] = length(CinfDrivers)
end

function sim_loop_delivery_wp!(sim::Dict, sim_summary::Dict, i_day::Int, Occ::Float64, Ncons::Int64, Params::Dict, PkgParams::Dict, PairParams::Dict, TestParams::Dict, cust_modifiers::Dict)
    
    reset_daily_contact_networks!(sim)
    
    #do_testing
    if TestParams["is_testing"]
        new_isolators = do_testing!(sim, TestParams, i_day, sim["isolation_network"])
    end
    
    #update infectivity and isolation status
    update_all_statuses!(sim, i_day)
    
    #update_in_work
    at_work, NAssignments = get_assignments!(sim, Occ, Ncons, PairParams)
    update_in_work!(sim, at_work)
    
    infpairs = Array{Int64,2}(undef,3,0)
    
    #introductions
    intro_pairs = get_introductions(sim, i_day)
    intros = get_customer_introductions(sim, i_day, NAssignments, cust_modifiers)
    intro_pairs = hcat(intro_pairs,
        [-transpose(ones(Int64,length(intros))); transpose(intros); 
        transpose(customer_contact*ones(Int64,length(intros)))])
    
    #get all contacts
    generate_random_contact_networks!(sim, Params, i_day)
    
    #customer infections
    CustInfProb, DriverDelivering = get_package_and_customer_infections!(sim, Ncons, i_day, PkgParams, PairParams, NAssignments, cust_modifiers)
    bool_infd = (rand(length(CustInfProb)) .< CustInfProb)
    CinfDs = DriverDelivering[bool_infd]
    
    g = collate_networks(sim)
    update_contact_network!(sim, g)
    infpairs = get_network_infections(sim, i_day)
    
    all_infpairs = hcat(intro_pairs, infpairs)
    
    infpairs_final = do_infections_randomly!(all_infpairs, sim, i_day)
    update_sim_summary_delivery_wp!(sim_summary, sim, i_day,
        infpairs_final; CinfDrivers=CinfDs)
    
    return infpairs_final
end

#run a either single outbreak with 1 initial case, or a sim with no cases but introductions
function run_sim_delivery_wp(Params::Dict, OccPerDay::Array{Float64,1}, NPPerDay::Array{Int64,1}; 
        PkgParams::Dict=DefaultPkgParams, PairParams::Dict=DefaultPairParams,
        TestParams::Dict=DefaultTestParams, Incidence::Array{Float64,1} = zeros(length(OccPerDay)),
        Prevalence::Array{Float64,1} = zeros(length(OccPerDay)))
    
    #if empty args are given, revert to default
    if length(PkgParams) == 0
        PkgParams = DefaultPkgParams
    end
    if length(PairParams) == 0
        PairParams = DefaultPairParams
    end
    if length(TestParams) == 0
        TestParams = DefaultTestParams
    end
    
    sim = initialise(Params, PairParams, Incidence, Prevalence)
    sim_summary, i_day, Go = setup_delivery_wp_model!(sim, Params, TestParams, OccPerDay)
    #defaults
    
    IsolParams = Dict("fixed_pair_isolation"=>false, "cohort_isolation"=>false, 
                      "car_share_isolation"=>false, "house_share_isolation"=>false)
    for key in IsolParams
        if haskey(Params,key)
            IsolParams[key] = Params[key]
        end
        if haskey(PairParams,key)
            IsolParams[key] = PairParams[key]
        end
    end
    create_isolation_network!(sim, IsolParams)
    
    CustModifiers = Dict("outdoor_contact_frac"=>1.0, "sanitise_frequency"=>0.0, "mask_prob"=>0.0)
    #edit these based on params
    
    
    while Go && (i_day <= length(OccPerDay))
        infpairs = sim_loop_delivery_wp!(sim, sim_summary, i_day,
            OccPerDay[i_day], NPPerDay[i_day], Params, PkgParams,
            PairParams, TestParams, CustModifiers)
        
        
        
        
#         fname = Printf.@sprintf("infection_network_%03d.png",i_day)
#         if visualise
#             print_infection_network(sim, fname, infpairs, node_x, node_y, pairs)
#         end
    
#         if testing && (i_day == test_days[test_day_counter])
#             test_day_counter += 1
#         end
        
        if Params["SimType"] == Outbreak_sim
            Go = any(get_infected_bool(sim))
        end
        
        i_day += 1
    end
    trim_sim_summary!(sim_summary, i_day-1, length(OccPerDay))

    return sim_summary
end



# "PackagesInfectiousOnDelivery"=>zeros(Int64,Ndays),
#                 "FomiteInfs"=>zeros(Int64,(3,Ndays)),
#                 "CustomersInfectedByPkgs"=>zeros(Int64,Ndays),
#                 "CustomersInfectedByDrivers"=>zeros(Int64,Ndays),


#List of modifiers
# modifiers["indoor_contact"] (for package delivery)