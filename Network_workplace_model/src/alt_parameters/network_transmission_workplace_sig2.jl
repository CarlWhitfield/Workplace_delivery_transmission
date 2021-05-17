include("viral_load_infectivity_testpos_sig2.jl")

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
using LightGraphs
using MetaGraphs
using GraphPlot

using Printf
using Base.Threads


#infection type definitions
const network_contact = 1
const non_network_contact = 2
const room_transmission = 3
const package_contact = 4
const customer_contact = 5
const introduction = 6
const pair_contact = 7
const edge_colours = [colorant"red", colorant"red", colorant"green", colorant"orange",
                      colorant"white",colorant"white",colorant"blue",colorant"lightgrey"]


#sim type definitions
const Outbreak_sim = 1
const Scenario_sim = 2


#const per_person_office_volume = 20   #m^3
#const office_turnover = 4*24    #ventilation rate per day
#const cabin_volume = 5           #m^3
#const cabin_turnover_closed = 5*24    #ventilation rate per day
#const cabin_turnover_open = 50*24    #ventilation rate per day
#const breathing_rate = 0.72*24  #lung turnovers per day m^3
const isol_time = 10
const Susc = 0
const Expd = 1
const Symp = 2
const Recd = 3
#const shift_pattern = 7     #days to draw shifts -- to be included
const infection_rate_F2F = 0.06   #infectivity of F2F interactions per hour = microCOVID approx
const infection_rate_F2F_outside = 0.012   #infectivity of F2F interactions per hour = microCOVID approx
const infection_rate_room = 0.0003    #shared space risk per hour per infected person (assuming average 4m distance, no talking) = microCOVID approx
const t_F2F = 1.0/6.0          #Average face to face contact duration (hours)
const t_office = 7.0
const t_lunch = 0.75

#const t_breaks = 1/24       #Time in break/lunch areas (1 hr)
#around 10 hrs per day: 85 deliveries per day,
const ParcelTimesPerDelivery = Dict("t_cabin"=>1/(12),"t_doorstep"=>1/(120), "t_handling"=>1/(60),
                                "t_picking"=>1/(60))
const BulkTimesPerDelivery = Dict("t_cabin"=>1/(12),"t_doorstep"=>1/(60), "t_handling"=>1/(12),
                              "t_picking"=>1/(12))
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


# function gamma_ab(mean::Float64, var::Float64)
#     #get gamma alpha beta from mean & var
#     b =  mean/var
#     return mean*b, b
# end

# function gamma_ab_mode(mode::Float64, beta::Float64)
#     #get gamme alpha from mode and beta
#     return mode*beta + 1
# end

# function generate_incubation_periods(Ntot::Int)
#     #generate incubation periods from gamma distribution
#     ai,bi = gamma_ab(mean_inc,std_inc^2)
#     p = rand(Gamma(ai,1/bi),Ntot)
#     return p
# end

# function generate_shedding_rate(Ntot::Int)
#     #generate incubation periods from gamma distribution
#     p = rand(Uniform(min_shed_rate,max_shed_rate),Ntot)
#     return p
# end

function generate_graph!(sim::Dict, phi::Float64, degree_logmean::Float64,
                         degree_logstd::Float64)
    #generate graph with lognormal node distribution and weighted probability for node connections
    #generate degree distribution
    node_degrees = Int64.(round.(rand(LogNormal(degree_logmean,
                                     degree_logstd),sim["Ntot"])))
    while any(node_degrees .> sim["Ntot"] - 1)
        Nleft = sum(node_degrees .> sim["Ntot"] - 1)
        node_degrees[node_degrees .> sim["Ntot"] - 1] .= Int64.(
            round.(rand(LogNormal(degree_logmean,degree_logstd),Nleft)))
    end
    #lognormal truncated at N-1
    #connect all loose edges
    sim["social_graph"] = simple_graph(sim["Ntot"],is_directed=false)
    nr = 1:sim["Ntot"]
    while sum(node_degrees .> 0) > 1
        i = sample(nr, Weights(node_degrees))
        ijob = sim["job"][i]
        #phi gives relative chance of edges with people in different job
        NDweights = node_degrees .* ((sim["job"] .== ijob) + phi .* (sim["job"] .!= ijob))
        j = sample(nr[nr .!= i], Weights(NDweights[nr .!= i]))
        Graphs.add_edge!(sim["social_graph"],i,j)
        node_degrees[i] -= 1
        node_degrees[j] -= 1
    end
end

# function generate_graph!(sim::Dict, av_degree::Float64, phi::Float64)
#     #generating incidence matrix: rather than full incidence,
#     #then also store list of lists for edges connected to nodes
#     #need to think about edge numbering
#     sim["social_graph"] = simple_graph(sim["Ntot"],is_directed=false)
#     Nh = zeros(Int64,3)
#     Nh[1] = 1
#     for i in 2:3
#         Nh[i] = Nh[i-1] + sim["N"][i-1]
#     end
#     #norm = (sim["N"][1]*((sim["N"][1]-1)/2 + phi*(sim["N"][2] + sim["N"][3])) +
#     #        sim["N"][2]*((sim["N"][2]-1)/2 + phi*sim["N"][3]) +
#     #        sim["N"][3]*(sim["N"][3]-1)/2)
#     coeff = av_degree/(sim["Ntot"]-1)  #av degree if phi = 1
#     for i in 1:3
#         for j in 1:i
#             if i == j
#                 Ngen = Int64((sim["N"][i]*(sim["N"][i]-1))/2)
#                 Pe = min(1, coeff)
#                 e = randsubseq(1:Ngen, Pe)
#                 Nij = cumsum((sim["N"][i]-1):-1:1)
#                 for edge in e
#                     ni = sum(edge .> Nij)
#                     if ni > 0
#                         nj = sim["N"][i] - edge + Nij[ni]
#                     else
#                         nj = sim["N"][i] - edge
#                     end
#                     Graphs.add_edge!(sim["social_graph"], ni+Nh[i], nj+Nh[i])
#                 end
#             else
#                 Ngen = sim["N"][i]*sim["N"][j]
#                 nr = 1:Ngen
#                 Pe = min(1, phi*coeff)
#                 e = randsubseq(1:Ngen, Pe)
#                 ni = Int64.(floor.((e.-1)./sim["N"][j])) .+ Nh[i]
#                 nj = mod.((e.-1),sim["N"][j]) .+ Nh[j]
#                 for k in 1:length(ni)
#                     Graphs.add_edge!(sim["social_graph"],ni[k],nj[k])
#                 end
#             end
#         end
#     end
# end

function init(Params)
    Ntot = Params["ND"]+Params["NL"]+Params["NO"]
    sim = Dict("Ntot"=>Ntot, "N"=>[Params["ND"],Params["NL"],Params["NO"]],
               "job"=>zeros(Int8, Ntot),  "inf_time"=>(zeros(Int64, Ntot) .- 1),
               "symp_time"=>-ones(Int64, Ntot), "asymptomatic"=>zeros(Bool, Ntot),
               "isolation_time"=>zeros(Int64, Ntot),
               "isolation_status"=>zeros(Bool, Ntot),
               "isolation_status_true_test"=>zeros(Bool, Ntot),
               "isolation_status_false_test"=>zeros(Bool, Ntot),
               "isolation_status_partner_true_test"=>zeros(Bool, Ntot),
               "isolation_status_partner_false_test"=>zeros(Bool, Ntot),
               "isolating_due_to_true_test"=>zeros(Bool, Ntot),
               "isolating_due_to_false_test"=>zeros(Bool, Ntot),
               "isolating_due_to_partner_true_test"=>zeros(Bool, Ntot),
               "isolating_due_to_partner_false_test"=>zeros(Bool, Ntot),
               "will_isolate"=>zeros(Bool, Ntot),
               "inf_mag"=>zeros(Float64, Ntot),
               "infection_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "VL_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "testing_paused"=>zeros(Bool,Ntot),
               "resume_testing"=>zeros(Int64,Ntot),
               "at_work"=>zeros(Bool, Ntot),
               "infection_status"=>zeros(Int8,Ntot),
               "days_infectious" => zeros(Int64,Ntot),
               "contact_prob_mat" => Params["p_contact"]*ones(3,3),
               "fixed_job_pairings" => Array{Int64,2}(undef,2,0))
    sim["job"][1:Params["ND"]] .= 1
    sim["job"][(Params["ND"]+1):(Params["ND"]+Params["NL"])] .= 2
    sim["job"][(Params["ND"]+Params["NL"]+1):(Params["ND"]+Params["NL"]+Params["NO"])] .= 3
    build_viral_load_distributions!(sim)
    sim["symp_day"] = Int64.(round.(rand(0:1,Ntot) .+ sim["symp_time"]))
    sim["asymptomatic"][generate_asymptomatics(sim["VL_mag"])] .= true
    sim["inf_mag"] .= generate_peak_inf.(sim["VL_mag"], sim["asymptomatic"])
    sim["will_isolate"][generate_isolations(Ntot, Params["Pisol"])] .= true
    nr = 1:sim["Ntot"]
    sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]
    sim["will_isolate"][sim["asymptomatic"]] .= false #asymptomatics don't self isolate, but are not "non isolators"
    sim["infection_profiles"] .= infectivity.(sim["VL_profiles"], sim["inf_mag"], sim["VL_mag"])
    sim["days_infectious"] = length.(sim["infection_profiles"])

    return sim
end

#maybe just change log10 viral load to be linear with beta, for low SAR this is effectively the same as being linear with p

function apply_contact_mixing_params!(sim::Dict, Params::Dict)
    for i = 1:3, j=1:3
        if i == 1 || j == 1
            sim["contact_prob_mat"][i,j] = Params["tD"]*sim["contact_prob_mat"][i,j]
        end
        if i != j
            sim["contact_prob_mat"][i,j] = Params["phi"]*sim["contact_prob_mat"][i,j]
        end
    end
    #     norm = Params["tD"]*sim["N"][1]*((sim["N"][1]-1)/2 + Params["phi"]*
    #              (sim["N"][2] + sim["N"][3])) + sim["N"][2]*(sim["N"][2]-1)/2 +
    #               sim["N"][3]*(sim["N"][3]-1)/2 + Params["phi"]*sim["N"][2]*sim["N"][3]
    #     coeff = 0.5*sim["Ntot"]*(sim["Ntot"]-1)/norm
    #     sim["inf_prob_mat"] .*= coeff
end



function initialise(Params::Dict)
    sim = init(Params)
    apply_contact_mixing_params!(sim, Params)

    return sim
end

function initialise(Params::Dict, PairParams::Dict)
    sim = init(Params)
    apply_contact_mixing_params!(sim, Params)
    init_pairs!(sim, PairParams)
    return sim
end

function initialise(Params::Dict, degree_logmean::Float64, degree_logstd::Float64)
    sim = init(Params)
    generate_graph!(sim, Params["phi"], degree_logmean, degree_logstd)
    #print("Average degree: ", 2*Graphs.num_edges(sim["social_graph"])
    #                           /Graphs.num_vertices(sim["social_graph"]),"\n",)
    return sim
end

function initialise(Params::Dict, PairParams::Dict, degree_logmean::Float64,
                         degree_logstd::Float64)
    sim = init(Params)
    generate_graph!(sim, Params["phi"], degree_logmean, degree_logstd::Float64)
    #print("Average degree: ", 2*Graphs.num_edges(sim["social_graph"])
    #                           /Graphs.num_vertices(sim["social_graph"]),"\n",)
    init_pairs!(sim,  PairParams)
    return sim
end

function initialise_withvis(Params::Dict, PairParams::Dict, degree_logmean::Float64,
                    degree_logstd::Float64, is_network::Bool, is_pairs::Bool,
                    Incidence::Array{Float64,1}, Prevalence::Array{Float64,1})
    if is_network
        if is_pairs
            sim = initialise(Params, PairParams, degree_logmean, degree_logstd)
        else
            sim = initialise(Params, degree_logmean, degree_logstd)
        end
        lg = copy_graph_into_lightgraph(sim["social_graph"])
        node_x, node_y = spring_layout(lg)
    else
        if is_pairs
            sim = initialise(Params, PairParams)
        else
            sim = initialise(Params)
        end
        lg = LightGraphs.SimpleGraph(sim["Ntot"])
        node_x, node_y = circular_layout(lg)
    end
    sim["Incidence"] = Incidence
    sim["Prevalence"] = Prevalence

    return sim, node_x, node_y
end

function initialise_novis(Params::Dict, PairParams::Dict, degree_logmean::Float64,
                    degree_logstd::Float64, is_network::Bool, is_pairs::Bool,
                    Incidence::Array{Float64,1}, Prevalence::Array{Float64,1})
    if is_network
        if is_pairs
            sim = initialise(Params, PairParams, degree_logmean, degree_logstd)
        else
            sim = initialise(Params, degree_logmean, degree_logstd)
        end
    else
        if is_pairs
            sim = initialise(Params, PairParams)
        else
            sim = initialise(Params)
        end
    end
    sim["Incidence"] = Incidence
    sim["Prevalence"] = Prevalence

    return sim
end

function infect_node!(sim::Dict, i::Int, time::Int, partner::Int64 = 0)
    sim["infection_status"][i] = Expd
    sim["inf_time"][i] = time
    #Assume partner isolation is guaranteed if infected node isolates (i.e. enforced)
    if sim["will_isolate"][i]
        if sim["isolation_time"][i] > time + sim["symp_day"][i] ||
           sim["isolation_time"][i] == 0 #might already be isolating
               sim["isolation_time"][i] = time + sim["symp_day"][i]
        end
#         print("Symp isol: ", i, ' ', time, ' ', sim["isolation_time"][i], '\n')
        if partner > 0   #pair isolation
            if sim["isolation_time"][partner] > time + sim["symp_day"][i] ||
               sim["isolation_time"][partner] == 0   #if not isolating already or isolating later than this
                    sim["isolation_time"][partner] = time + sim["symp_day"][i]
            end
        end
    end
end

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

function update_sim_summary!(summary::Dict, sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    isolator_bool = (sim["isolation_status"] .== true)
    new_isolator_bool = isolator_bool .* (sim["isolation_time"] .== i_day)
    isol_true_test = (sim["isolation_status_true_test"] .== true)
    isol_false_test = (sim["isolation_status_false_test"] .== true)
    isol_partner_true_test = (sim["isolation_status_partner_true_test"] .== true)
    isol_partner_false_test = (sim["isolation_status_partner_false_test"] .== true)
    isol_no_test = (.!(isol_true_test) .* .!(isol_false_test) .*
                    .!(isol_partner_true_test) .* .!(isol_partner_false_test))
    asymp_bool = (sim["asymptomatic"] .== true)
    inf_bool = (sim["infection_status"] .== Symp)
    at_work_bool = (sim["at_work"] .== true)
    for j in 1:3
        jt = (sim["job"] .== j)

        summary["Susceptible"][j,i_day] = sum(jt .* (sim["infection_status"] .== Susc))
        summary["Recovered"][j,i_day] = sum(jt .* (sim["infection_status"] .== Recd))
        summary["Infectious"][j,i_day] = (sum(jt) - summary["Susceptible"][j,i_day]
                                         - summary["Recovered"][j,i_day])
        summary["Exposed"][j,i_day] = sum(jt .* (sim["infection_status"] .== Expd))
        summary["Isolated"][j,i_day] = sum(jt .* isolator_bool)
        summary["NewIsolators"][j,i_day] = sum(jt .* new_isolator_bool)


        summary["IsolatedDueToSymptoms"][j,i_day] = sum(jt .* isol_no_test)
        summary["IsolatedDueToTestAsymp"][j,i_day] = sum(jt .* isol_true_test .* asymp_bool)
        summary["IsolatedDueToTestSymp"][j,i_day] = sum(jt .* isol_true_test .* .!(asymp_bool))
        summary["IsolatedDueToFalsePos"][j,i_day] = sum(jt .* isol_false_test)
        summary["IsolatedDueToPartnerTruePos"][j,i_day] = sum(jt .* isol_partner_true_test)
        summary["IsolatedDueToPartnerFalsePos"][j,i_day] = sum(jt .* isol_partner_false_test)

        summary["NewSympIsolators"][j,i_day] = sum(jt .* isol_no_test .* new_isolator_bool)
        summary["NewTestAsympIsolators"][j,i_day] = sum(jt .* isol_true_test .* new_isolator_bool .* asymp_bool)
        summary["NewTestSympIsolators"][j,i_day] = sum(jt .* isol_true_test .* new_isolator_bool .* .!(asymp_bool))
        summary["NewFalseIsolators"][j,i_day] = sum(jt .* isol_false_test .* new_isolator_bool)
        summary["NewPartnerTruePosIsolators"][j,i_day] = sum(jt .* isol_partner_true_test .* new_isolator_bool)
        summary["NewPartnerFalsePosIsolators"][j,i_day] = sum(jt .* isol_partner_false_test .* new_isolator_bool)

        summary["Presenting"][j,i_day] = sum(jt .* inf_bool .* at_work_bool)
        summary["Asymptomatic"][j,i_day] = sum(jt .* inf_bool .* asymp_bool)
    end
end

function update_sim_summary!(summary::Dict, sim::Dict, i_day::Int,
                             inf_packages::Int, inf_pairs::Array{Int64,2};
                             CinfsD::Int = 0, CinfsP::Int = 0)
    update_sim_summary!(summary, sim, i_day)
    summary["PackagesInfectiousOnDelivery"][i_day] = inf_packages
    summary["CustomersInfectedByPkgs"][i_day] = CinfsP
    summary["CustomersInfectedByDrivers"][i_day] = CinfsD
    nr = 1:sim["Ntot"]
    if haskey(summary,"IndexCaseInfections")
        summary["IndexCaseInfections"] += sum(inf_pairs[1,:] .== summary["IndexCase"])
    end
    for j in 1:3
        nj = nr[sim["job"] .== j]
        ipj = inf_pairs[:, (inf_pairs[2,:] .>= min(nj...)) .* (inf_pairs[2,:] .<= max(nj...))]
        summary["NetworkInfs"][j,i_day] = sum(ipj[3,:] .== network_contact)
        summary["ContactInfs"][j,i_day] = sum(ipj[3,:] .== non_network_contact)
        summary["RoomInfs"][j,i_day] = sum(ipj[3,:] .== room_transmission)
        summary["FomiteInfs"][j,i_day] = sum(ipj[3,:] .== package_contact)
        summary["PairInfs"][j,i_day] = sum(ipj[3,:] .== pair_contact)
        summary["CustomerIntroductions"][j,i_day] = sum(ipj[3,:] .== customer_contact)
        summary["ExternalIntroductions"][j,i_day] = sum(ipj[3,:] .== introduction)
    end
end

function update_infectivity!(sim::Dict, i_day::Int)
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

function update_isolation!(sim::Dict, i_day::Int)
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

function select_pairs(sim::Dict, occ::Float64, fixed_pairs::Bool, job::Int64, Ncons::Int64)
    TotPairs = Int64(sim["N"][job]/2)
    nr = 1:sim["Ntot"]
    jobgroup = nr[sim["job"] .== job]
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

    return pairs, pair_cons
end

function update_in_work!(sim::Dict, occ::Float64, Ncons::Int64,
                driver_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0),
                driver_cons::Array{Int64,1} = Array{Int64,1}(undef,0),
                loader_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0),
                loader_cons::Array{Int64,1} = Array{Int64,1}(undef,0))
    #wih driver pairs
    nr = 1:sim["Ntot"]
    sim["at_work"] .= false
    NAs = zeros(Int64,sim["Ntot"])
    to_select = []
    if length(driver_pairs) > 0
        sim["at_work"][vec(driver_pairs)] .= true
        NAs[driver_pairs[1,:]] .= driver_cons
        NAs[driver_pairs[2,:]] .= driver_cons
    else
        push!(to_select,1)
    end
    if length(loader_pairs) > 0
        sim["at_work"][vec(loader_pairs)] .= true
        NAs[loader_pairs[1,:]] .= loader_cons
        NAs[loader_pairs[2,:]] .= loader_cons
    else
        push!(to_select,2)
    end
    push!(to_select,3)
    for j in to_select
        job = nr[sim["job"] .== j]

        available = job[sim["isolation_status"][job] .== false]
        if length(available)/length(job) > occ
            w = sample(job, round(Int,occ*length(job)), replace=false)
        else
            w = available
        end
        if j < 3
            NAs[w] .= rand(Multinomial(Ncons,length(w)))
        end
        sim["at_work"][w] .= true
    end
    return NAs
end

function update_in_work!(sim::Dict, occ::Float64, Ncons::Int64)
    #without driver pairs
    nr = 1:sim["Ntot"]
    sim["at_work"] = zeros(Bool,sim["Ntot"])
    NAs = zeros(Int64,sim["Ntot"])
    for j in 1:3
        job = nr[sim["job"] .== j]
        available = job[sim["isolation_status"][job] .== false]
        if length(available)/length(job) > occ
            w = sample(job, round(Int,occ*length(job)), replace=false)
        else
            w = available
        end
        if j < 3
            NAs[w] .= rand(Multinomial(Ncons,length(w)))
        end

        sim["at_work"][w] .= true
    end



    return NAs
end

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

function generate_pair_infections(sim::Dict, i_day::Int, t_pair::Array{Float64,1},
                                  pairs::Array{Int64,2}, job::Int,
                                  t_cabin::Array{Float64,1}=Array{Float64,1}(undef,0),
                                  cabin_window::Bool=false)

    all_inf, all_inf_scales = get_infectivities(sim, i_day)
    ipairs = Array{Int64,2}(undef,3,0)
    inf_job_cat = (sim["job"][all_inf] .== job)
    nr = collect(1:sim["Ntot"])
    nrh = nr[sim["job"] .== job]
    if sum(inf_job_cat) > 0
        inf = all_inf[inf_job_cat]
        inf_scales = all_inf_scales[inf_job_cat]

        p1 = zeros(length(inf))
        p2 = zeros(length(inf))

        for (i_inf, infp) in enumerate(inf)
            pair_col = Bool.((pairs[1,:] .== infp) + (pairs[2,:] .== infp))
            p1[i_inf] = infection_rate_F2F_outside * t_pair[pair_col][1]
            ir = infection_rate_F2F
            if cabin_window
                ir = infection_rate_F2F_outside
            end
            if length(t_cabin) > 0
                 p2[i_inf] = ir * t_cabin[pair_col][1]
             end
        end
        iprob = 1.0 .- exp.(-inf_scales .* (p1 + p2))
        #randomly select those that infect partner
        pair_infs = inf[rand(length(inf)) .< (iprob)]
        for p in pair_infs
            I = findall(x->x==p, pairs)
            j = mod(I[1][1],2) + 1
            k = I[1][2]
            if sim["infection_status"][pairs[j,k]] == Susc
                ipairs = hcat(ipairs,[pairs[I[1][1],k]; pairs[j,k]; pair_contact])
            end
        end
    end
    return ipairs
end

function get_pair_infections(infpairs::Array{Int64,2}, sim::Dict, occ::Float64,
                            i_day::Int, PairParams::Dict, Ncons::Int64)
    dpairs = Array{Int64,2}(undef,2,0)
    lpairs = Array{Int64,2}(undef,2,0)
    NDPassignments = Array{Int64,1}(undef,0)
    NLPassignments = Array{Int64,1}(undef,0)
    if PairParams["is_driver_pairs"]
        dpairs, NDPassignments = select_pairs(sim, occ, PairParams["fixed_driver_pairs"], 1, Ncons)
    end
    if PairParams["is_loader_pairs"]
        lpairs, NLPassignments = select_pairs(sim, occ, PairParams["fixed_loader_pairs"], 2, Ncons)
    end
    NAs = update_in_work!(sim, occ, Ncons, dpairs, NDPassignments, lpairs, NLPassignments)
    ips = generate_pair_infections(sim, i_day, NDPassignments .* (BulkTimesPerDelivery["t_handling"] +
          BulkTimesPerDelivery["t_doorstep"]), dpairs, 1,
          NDPassignments .* BulkTimesPerDelivery["t_cabin"],  PairParams["is_window_open"])
    InfPairsNew = hcat(infpairs,ips)
    ips = generate_pair_infections(sim, i_day, NLPassignments .* BulkTimesPerDelivery["t_picking"],
                                     lpairs, 2)


    InfPairsNew = hcat(InfPairsNew,ips)

    return InfPairsNew, dpairs, lpairs, NAs
end

function generate_infectious_packages(sim::Dict, NP::Int64, AllDrivers::Array{Int64,1},
    InfDrivers::Array{Int64,1}, InfDriverScales::Array{Float64,1},
    AllLoaders::Array{Int64,1}, InfLoaders::Array{Int64,1}, InfLoaderScales::Array{Float64,1},
    PkgParams::Dict, NAs::Array{Int64,1})

    DAorder = cumsum(NAs[AllDrivers])
    InfectorAtPickup = Array{Int64,1}(undef,0)
    StimesAtPickup = Array{Float64,1}(undef,0)

    InfectorAtDropoff = Array{Int64,1}(undef,0)
    IDDriverPickup = Array{Int64,1}(undef,0)
    IDDriverDropoff = Array{Int64,1}(undef,0)

    # l_assignment = zeros(Int64,NP) #assigns potentially infectious packages to loaders or loader pairs
    # d_assignment = zeros(Int64,NP)
    NPavailable = collect(1:NP)
    #generate packages infected by loaders
    #PkgsInfected = []
    if length(InfLoaders) > 0
        #now InfLoaders either contains infectious loaders or loader pairs
        #the latter are treated as a unit with infectivity equal to their combined infectivitiy

        LInfPkgs = sum(NAs[InfLoaders])
        if LInfPkgs > 0
            #assign these packages to loaders (or pairs) and then infected ones to drivers (or pairs)
            NLIPinfected = zeros(Int64,length(InfLoaders))
            nlip = 1:length(LInfPkgs)
            for i in 1:length(InfLoaders)
                pkg_inf_prob = 1 - exp(-PkgParams["p_fomite_trans"] * InfLoaderScales[i])
                iLpkgs = NAs[InfLoaders[i]]
                Infh = rand(Binomial(iLpkgs, pkg_inf_prob))    #get no, of infected packages
                NLIPinfected[i] = Infh
                if Infh > 0
                    Lth = PkgParams["Ltime"] .* rand(Infh)  #generate time infected
                    #generate strilisation time
                    sth = Lth .+ rand(Exponential(PkgParams["PkgHlife"]/log(2)), Infh)
                    Iapcond = (sth .> PkgParams["Ltime"])  #check if infectious at pickup
                    if sum(Iapcond) > 0
                        InfectorAtPickup = vcat(InfectorAtPickup,InfLoaders[i]*ones(Int64,sum(Iapcond)))  #push back number infectious at pickup
                        StimesAtPickup = vcat(StimesAtPickup, sth[Iapcond])
                        PDAgen = sample(NPavailable, sum(Iapcond), replace=false)  #generate package number
                        filter!(eh -> !(eh in PDAgen), NPavailable)  #remove PDAgen from NP available
                        #push!.(Ref(PkgsInfected), PDAgen)
                        Nd = 1:length(AllDrivers)
                        for p in PDAgen  #for each package number generated, assign to driver based on cumsum
                            bool1 = p .> vcat([0],DAorder[1:(end-1)])
                            bool2 = p .<= DAorder
                            Dbool = bool1 .* bool2
                            push!(IDDriverPickup, AllDrivers[Dbool][1])
                        end
                    end
                end
            end
        end
        #print("Packages infected: ", PkgsInfected,"\n")
    end

    #package infections by drivers here
    if length(InfDrivers) > 0
        #now InfDrivers either contains infectious drivers or driver pairs
        #the latter are treated as a unit with infectivity equal to their combined infectivitiy
        DInfPkgs = sum(NAs[InfDrivers])  #packages handled by infectious drivers or driver pairs
        if DInfPkgs > 0
            NDIPinfected = zeros(Int64,length(InfDrivers))   #number of packages each driver infects
            ndip = 1:length(DInfPkgs)        #index of all packages contacted by inf drivers
            for i in 1:length(AllDrivers)       #loop over all drivers to find infectees
                iDpkgs = NAs[AllDrivers[i]]     #number of packages handled by driver
                Dtimes = zeros(iDpkgs)
                Stimes = zeros(iDpkgs)  #store sterilisation time of all packages delivered (0 if uninfected)
                Ilast = zeros(Int64,iDpkgs)
                #add those infected by loaders (that are infectious at pickup)
                IAPbool = (IDDriverPickup .== AllDrivers[i])
                NLAPs = sum(IAPbool)
                if NLAPs > 0   #just assign these as first NLAPs parcels
                    Stimes[1:NLAPs] .= StimesAtPickup[IAPbool]
                    Dtimes[1:NLAPs] .= PkgParams["Ltime"] .+ rand(NLAPs) .* PkgParams["Dtime"]
                    Ilast[1:NLAPs] .= InfectorAtPickup[IAPbool]
                end
                #add infected by drivers
                if AllDrivers[i] in InfDrivers
                    j = collect(1:length(InfDrivers))[AllDrivers[i] .== InfDrivers][1]
                    pkg_inf_prob = 1 - exp(-PkgParams["p_fomite_trans"] * InfDriverScales[j])
                    Infh = rand(Binomial(iDpkgs, pkg_inf_prob))    #get no, of infected packages
                    NDIPinfected[j] = Infh
                    if Infh > 0
                        #generate package numbers of those infected by drivers
                        NPinfh = sample(1:iDpkgs, Infh, replace=false)
                        #generate sterilisation times
                        Dth = PkgParams["Dtime"] .* rand(Infh)  #generate time delivered
                        binrand = rand(0:1, Infh)    #random whether infected on pickup or dropoff
                        InfTime = PkgParams["Ltime"] .+ binrand .* Dth   #time infected by driver
                        sth = InfTime .+ rand(Exponential(PkgParams["PkgHlife"]/log(2)), Infh)
                        update_sts_bool = (sth .> Stimes[NPinfh])
                        if sum(update_sts_bool) > 0
                            Stimes[NPinfh[update_sts_bool]] .= sth[update_sts_bool]
                            Dtimes[NPinfh[update_sts_bool]] .= PkgParams["Ltime"] .+ rand(sum(update_sts_bool)) .* PkgParams["Dtime"]
                            Ilast[NPinfh[update_sts_bool]] .= AllDrivers[i]
                        end
                    end
                end
                IADbool = (Stimes .> 0) .* (Stimes .> Dtimes)
                if sum(IADbool) > 0
                    InfectorAtDropoff = vcat(InfectorAtDropoff,Ilast[IADbool])
                    IDDriverDropoff = vcat(IDDriverDropoff,AllDrivers[i]*ones(Int64,sum(IADbool)))
                end
            end
        end
    end

    #returns vectors for packages that are infectious at pickup and dropoff
    #first is identity of infectors
    #second is identity of drivers (or pairs) delivering the packages
    return InfectorAtPickup, IDDriverPickup, InfectorAtDropoff, IDDriverDropoff
end

function get_package_infections(infpairs::Array{Int64,2}, sim::Dict, NP::Int64,
                                i_day::Int, PkgParams::Dict, NAs::Array{Int64,1})
    #requires testing and commenting
    #packages: find infectious people in group 2 (loaders)
    if NP > 0 && PkgParams["p_fomite_trans"] > 0
        nr = 1:sim["Ntot"]
        inf, inf_scales = get_infectivities(sim, i_day)
        d_in_work = nr[(sim["job"] .== 1) .* (sim["at_work"])]
        l_in_work = nr[(sim["job"] .== 2) .* (sim["at_work"])]
        inf_loaders = inf[(sim["job"][inf] .== 2)]
        inf_loader_scales = inf_scales[(sim["job"][inf] .== 2)]
        inf_drivers = inf[(sim["job"][inf] .== 1)]
        inf_driver_scales = inf_scales[(sim["job"][inf] .== 1)]
        InfectorsAtPickup, PickupAssignment, InfectorsAtDropoff, DropoffAssignment =
            generate_infectious_packages(sim, NP, d_in_work, inf_drivers, inf_driver_scales,
                            l_in_work, inf_loaders, inf_loader_scales, PkgParams, NAs)
        #InfectorsAtPickup contains infector number
        #same for InfectorsAtDropoff
        #Assignments show which drivers these are assigned to
        DInfectors = Array{Int64,1}(undef,0)
        DriversInfected = Array{Int64,1}(undef,0)
        if length(InfectorsAtPickup) > 0
            ni = collect(1:length(InfectorsAtPickup))
            Pinfs =  randsubseq(ni, 0.5*PkgParams["p_fomite_contr"])
            for Pf in Pinfs
                push!(DriversInfected, PickupAssignment[Pf])
                push!(DInfectors, InfectorsAtPickup[Pf])
            end
        end
        if length(InfectorsAtDropoff) > 0
            ni = collect(1:length(InfectorsAtDropoff))
            Dinfs =  randsubseq(ni, 0.5*PkgParams["p_fomite_contr"])
            for Df in Dinfs
                push!(DriversInfected, DropoffAssignment[Df])
                push!(DInfectors, InfectorsAtDropoff[Df])
            end
        end
        InfPairsNew = vcat(transpose(DInfectors), transpose(DriversInfected),
                           transpose(package_contact .* ones(Int64, length(DriversInfected))))

        return hcat(infpairs, InfPairsNew), InfectorsAtDropoff, DropoffAssignment
    else
        return infpairs, Array{Int64,1}(undef,0), Array{Int64,1}(undef,0)
    end
end

function get_package_infections(infpairs::Array{Int64,2}, sim::Dict, NP::Int64,
                                i_day::Int, PkgParams::Dict, PairParams::Dict,
                                NAs::Array{Int64,2},
                                driver_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0),
                                loader_pairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))

    if NP > 0 && PkgParams["p_fomite_trans"] > 0
        nr = 1:sim["Ntot"]
        inf, inf_scales = get_infectivities(sim, i_day)
        d_in_work = nr[(sim["job"] .== 1) .* (sim["at_work"])]
        l_in_work = nr[(sim["job"] .== 2) .* (sim["at_work"])]
        inf_loaders = inf[(sim["job"][inf] .== 2)]
        inf_loader_scales = inf_scales[(sim["job"][inf] .== 2)]
        inf_drivers = inf[(sim["job"][inf] .== 1)]
        inf_driver_scales = inf_scales[(sim["job"][inf] .== 1)]

        if length(driver_pairs > 0)
            NDPairs = size(driver_pairs,2)
            #pairs bit
            ndp = 1:NDPairs
            dp_infs = zeros(bool,NDPairs)
            dp_scales = zeros(Float64,NDPairs)
            for (i, dp) in enumerate(driver_pairs)
                dp1inf = false
                if (dp[1] in inf_drivers)
                    dp_infs[i] = true
                    dp1inf = true
                    dp_scales[i] += inf_scales[inf .== dp[1]][1]
                end
                if (dp[2] in inf_loaders)
                    dp_infs[i] = true
                    dp_scales[i] += inf_scales[inf .== dp[2]][1]
                    if dp1inf   #correction for joint infctivity as package can only be infected once
                        dp_scales[i] -= (inf_scales[inf .== dp[1]][1]
                                      * inf_scales[inf .== dp[2]][1])*PkgParams["p_fomite_trans"]
                    end
                end
            end
            alld = ndp
            infd = ndp[dp_infs]
            infd_scales = dp_scales[dp_infs]
            NA1 = NAs[driver_pairs[1,alld]]
        else
            alld = nr[sim["job"].==1]
            infd = inf_drivers
            infd_scales = inf_driver_scales
            NA1 = NAs[alld]
        end

        if length(loader_pairs > 0)
            NLPairs = size(loader_pairs,2)
            #pairs bit

            nlp = (1 + length(alld)):(1 + length(alld) + NLPairs)
            lp_infs = zeros(bool,NLPairs)
            lp_scales = zeros(Float64,NLPairs)
            for (i, lp) in enumerate(loader_pairs)
                lp1inf = false
                if (lp[1] in inf_loaders)
                    lp_infs[i] = true
                    lp1inf = true
                    lp_scales[i] += inf_scales[inf .== lp[1]][1]
                end
                if (lp[2] in inf_loaders)
                    lp_infs[i] = true
                    lp_scales[i] += inf_scales[inf .== lp[2]][1]
                    if lp1inf   #correction for joint infctivity as package can only be infected once
                        lp_scales[i] -= (inf_scales[inf .== lp[1]][1]*inf_scales[inf .== lp[2]][1])*PkgParams["p_fomite_trans"]
                    end
                end
            end
            alll = nlp
            infl = nlp[lp_infs]
            infl_scales = lp_scales[lp_infs]
            NA2 = NAs[loader_pairs[1,alll]]
        else
            alll = nr[sim["job"].==2]
            infl = inf_loaders
            infl_scales = inf_loader_scales
            NA2 = NAs[alll]
        end
        NAsend = hcat(NA1,NA2)
        InfectorsAtPickup, PickupAssignment, InfectorsAtDropoff, DropoffAssignment =
            generate_infectious_packages(sim, NP, d_in_work, inf_drivers, inf_driver_scales,
                            l_in_work, inf_loaders, inf_loader_scales, PkgParams, NAsend)

        #do driver infections by packages here
        if length(InfectorsAtPickup) > 0
            nap = collect(1:length(InfectorsAtPickup))
            nad = collect(1:length(InfectorsAtDropoff))
            PkgIAP = zeros(Int64, length(InfectorsAtPickup))  #store actual infectors (not pair numbers)
            PkgIAD = zeros(Int64, length(InfectorsAtDropoff))  #store actual infectors
            #bool separation of
            nDIAP = nap[(InfAtPickup .<= length(alld))]
            nLIAP = nap[(InfAtPickup .> length(alld))]
            nDIAD = nad[(InfAtDropoff .<= length(alld))]
            nLIAD = nad[(InfAtDropoff .> length(alld))]

            if length(driver_pairs) > 0  #if drivers are in pairs, extract them first
                DriversExposedAtPickup = vec(driver_pairs[:,PickupAssignment])
                DriversExposedAtDropoff = vec(driver_pairs[:,DropoffAssignment])
                for n in nDIAP
                    pair = driver_pairs[:,InfectorsAtPickup[n]]
                    #decide which of pair caused infection
                    #check if one is not infectious
                    if sim["infection_status"][pair[1]] == Susc || sim["infection_status"][pair[1]] == Recd
                        PkgIAP[n] = pair[2]
                    elseif sim["infection_status"][pair[2]] == Susc || sim["infection_status"][pair[2]] == Recd
                        PkgIAP[n] = pair[1]
                    else
                        #otherwise choose based on relative probability they infected it
                        p = 1 .- exp.( -PkgParams["p_fomite_trans"] .* inf_driver_scales[pair])
                        p1rel = p[1]/(p[1] + p[2])
                        r = rand()
                        if r < p1rel
                            PkgIAP[n] = pair[1]
                        else
                            PkgIAP[n] = pair[2]
                        end
                    end
                end

                for n in nDIAD
                    pair = driver_pairs[:,InfectorsAtDropoff[n]]
                    #decide which of pair caused infection
                    #check if one is not infectious
                    if sim["infection_status"][pair[1]] == Susc || sim["infection_status"][pair[1]] == Recd
                        PkgIAD[n] = pair[2]
                    elseif sim["infection_status"][pair[2]] == Susc || sim["infection_status"][pair[2]] == Recd
                        PkgIAD[n] = pair[1]
                    else
                        #otherwise choose based on relative probability they infected it
                        p = 1 .- exp.( -PkgParams["p_fomite_trans"] .* inf_driver_scales[pair])
                        p1rel = p[1]/(p[1] + p[2])
                        r = rand()
                        if r < p1rel
                            PkgIAD[n] = pair[1]
                        else
                            PkgIAD[n] = pair[2]
                        end
                    end
                end
            else
                DriversExposedAtPickup = PickupAssignment
                DriversExposedAtDropoff= DropoffAssignment
                PkgIAP[nDIAP] .= InfectorsAtPickup[nDIAP]
                PkgIAD[nDIAD] .= InfectorsAtDropoff[nDIAD]
            end
            if length(loader_pairs) > 0
                for n in nLIAP
                    pair = loader_pairs[:,InfectorsAtPickup[n] - length(alld)]
                    #decide which of pair caused infection
                    #check if one is not infectious
                    if sim["infection_status"][pair[1]] == Susc || sim["infection_status"][pair[1]] == Recd
                        PkgIAP[n] = pair[2]
                    elseif sim["infection_status"][pair[2]] == Susc || sim["infection_status"][pair[2]] == Recd
                        PkgIAP[n] = pair[1]
                    else
                        #otherwise choose based on relative probability they infected it
                        p = 1 .- exp.( -PkgParams["p_fomite_trans"] .* inf_loader_scales[pair])
                        p1rel = p[1]/(p[1] + p[2])
                        r = rand()
                        if r < p1rel
                            PkgIAP[n] = pair[1]
                        else
                            PkgIAP[n] = pair[2]
                        end
                    end
                end
                for n in nLIAD
                    pair = loader_pairs[:,InfectorsAtDropoff[n] - length(alld)]
                    #decide which of pair caused infection
                    #check if one is not infectious
                    if sim["infection_status"][pair[1]] == Susc || sim["infection_status"][pair[1]] == Recd
                        PkgIAD[n] = pair[2]
                    elseif sim["infection_status"][pair[2]] == Susc || sim["infection_status"][pair[2]] == Recd
                        PkgIAD[n] = pair[1]
                    else
                        #otherwise choose based on relative probability they infected it
                        p = 1 .- exp.( -PkgParams["p_fomite_trans"] .* inf_loader_scales[pair])
                        p1rel = p[1]/(p[1] + p[2])
                        r = rand()
                        if r < p1rel
                            PkgIAD[n] = pair[1]
                        else
                            PkgIAD[n] = pair[2]
                        end
                    end
                end
            else
                PkgIAP[nLIAD] .= InfectorsAtPickup[nLIAP]
                PkgIAD[nLIAD] .= InfectorsAtDropoff[nLIAD]
            end

            DriversInfected = Array{Int64,1}(undef,0)
            DriversInfectors = Array{Int64,1}(undef,0)
            PDrange = 1:length(DriversExposedAtPickup)
            PDinfs = randsubseq(PDrange, 0.5*PkgParams["p_fomite_contr"])
            DDrange = 1:length(DriversExposedAtDropoff)
            DDinfs = randsubseq(DDrange, 0.5*PkgParams["p_fomite_contr"])
            if length(PDinfs) > 0
                DriversInfected = vcat(DriversInfected, DriversExposedAtPickup[PDinfs])
                DriversInfectors = vcat(DriversInfectors, PkgIAP[PDinfs])
            end
            if length(DDinfs) > 0
                DriversInfected = vcat(DriversInfected, DDinfs)
                DriversInfectors = vcat(DriversInfectors, PkgIAD[DDinfs])
            end
        end
        InfPairsNew = vcat(DriversInfectors, DriversInfected,
                           package_contact .* ones(Int64,length(DriversInfected))) #packages labelled with zeros

        return hcat(infpairs, InfPairsNew), PkgIAD, DropoffAssignment
    else
        return infpairs, Array{Int64,1}(undef,0), Array{Int64,1}(undef,0)
    end
end

function get_network_infections(infpairs::Array{Int64,2}, sim::Dict, Params::Dict,
                                i_day::Int)
    inf, inf_scales = get_infectivities(sim, i_day)
    ipairs = Array{Int64,2}(undef,3,0)
    #collect all unique edges
    #all have different infection rates
    #get indices of all infectious edges
    #print(inf,'\n')
    if length(inf) > 0
        eind = edge_index.(vcat(out_edges.(inf,Ref(sim["social_graph"]))...))
        #infectious nodes in
        nin = vcat(fill.(inf,length.(out_edges.(inf,Ref(sim["social_graph"]))))...)
        #nodes out
        nout1 = source.(Graphs.edges(sim["social_graph"])[eind])
        nout2 = target.(Graphs.edges(sim["social_graph"])[eind])
        nout = ((nout1 .!= nin) .* nout1) .+ ((nout2 .!= nin) .* nout2)
        #ecprob is probability of edge contact
        #drivers have contact rate reduced by tD
        ecprob = Params["tD"] .* ((sim["job"][nin] .== 1) .| (sim["job"][nout] .== 1)) .+
                                 ((sim["job"][nin] .!= 1) .* (sim["job"][nout] .!= 1))
        #p_friend_contact < 1 reduced probability of contacting friends in work
        ecprob = Params["p_friend_contact"] .* ecprob
        #need to check if they are in work (inf nodes are by definition)
        not_at_work_bool = sim["at_work"][nout] .== false
        ecprob[not_at_work_bool] .= 0.0
        #get infectivity associated with each edge
        beta_vals = vcat(fill.(inf_scales,length.(out_edges.(inf,Ref(sim["social_graph"]))))...)
        eprob = 1 .- exp.(-infection_rate_F2F .* beta_vals .* t_F2F)

        #draw which will be infectious (duplicates don't matter)
        #print(ecprob,"\n\n")
        einf = eind[rand(length(eind)) .< (eprob .* ecprob)]
        if length(einf) > 0
            inf_edges = Graphs.edges(sim["social_graph"])[einf]
            s_nodes = source.(inf_edges)
            t_nodes = target.(inf_edges)
            #list all infectious contacts - even duplicates, then any nodes with
            #multiple infectors, one will be picked at random
            sn1 = sim["infection_status"][s_nodes] .== Susc
            if sum(sn1) > 0
                ipairs = hcat(ipairs, [transpose(t_nodes[sn1]); transpose(s_nodes[sn1]);
                      network_contact * transpose(ones(Int64,sum(sn1)))])
            end
            tn1 = sim["infection_status"][t_nodes] .== Susc
            if sum(tn1) > 0
                ipairs = hcat(ipairs, [transpose(s_nodes[tn1]); transpose(t_nodes[tn1]);
                      network_contact * transpose(ones(Int64,sum(tn1)))])
            end
        end
    end
     #susceptibles
     return hcat(infpairs, ipairs)
end

function get_contact_infections(infpairs::Array{Int64,2}, sim::Dict, Params::Dict,
                                i_day::Int)
    inf, inf_scales = get_infectivities(sim, i_day)
    #for those who are infectious and in work, generate contacts randomly with susceptibles
    #with given probability, ignore contacts between infectious cases (although these
    #are implied, just not explicitly simulated)
    ipairs = Array{Int64,2}(undef,3,0)
    if length(inf) > 0
        nr = 1:sim["Ntot"]
        nS = nr[(sim["infection_status"] .== Susc) .* (sim["at_work"])]
        #nS contains not currently infectious and in work
        j0 = sim["job"][nS] #j0 is jobs of not currently infectious
        for (k, i) in enumerate(inf)
            j = sim["job"][i]
            p1 = sim["contact_prob_mat"][j,j0] * (1 - exp.(-inf_scales[k] *
                                         infection_rate_F2F * t_F2F))
            if j == 3
                p2a = zeros(length(j0))
                p2a[j0 .== 3] .= 1.0
                p2a = p2a .* (1 - exp(-inf_scales[k] * (t_office + t_lunch)*infection_rate_room))
                p2b = zeros(length(j0))
                p2b[j0 .!= 3] .= 1.0
                p2b = p2b .* (1 - exp(-inf_scales[k] * t_lunch*infection_rate_room))
                p2 = p2a .+ p2b
            else
                p2 = ones(length(j0)) .* (1 - exp(-inf_scales[k] * t_lunch*infection_rate_room))
            end
            infectees1 = nS[rand(length(nS)) .<  p1]
            infectees2 = nS[rand(length(nS)) .<  p2]
            if length(infectees1) > 0
                ipairs = hcat(ipairs, [i*transpose(ones(Int64,length(infectees1)));
                                      transpose(infectees1);
                                      non_network_contact*transpose(ones(Int64,length(infectees1)))])
            end
            if length(infectees2) > 0
                ipairs = hcat(ipairs, [i*transpose(ones(Int64,length(infectees2)));
                                      transpose(infectees2);
                                      room_transmission*transpose(ones(Int64,length(infectees2)))])
            end
        end
    end
    return hcat(infpairs, ipairs)
end

function get_customer_infections(infpairs::Array{Int64,2}, sim::Dict, Params::Dict,
                                 PkgParams::Dict, Nassignments::Array{Int64,1},
                                 t_del::Float64, InfPkgs::Array{Int64,2}, i_day::Int)
    ipairs = Array{Int64,2}(undef,3,0)
    contact_infections = Array{Int64,2}(undef,2,0)
    pkg_infections = Array{Int64,2}(undef,2,0)
    nr = 1:sim["Ntot"]
    d_in_work = nr[(sim["at_work"]) .* (sim["job"] .== 1)]
    p_pos = sim["Prevalence"][i_day] * (1 - exp(-infection_rate_F2F*t_del))
    for id in d_in_work
        #assmue every contact has the same probability
        dinfs = randsubseq(1:Nassignments[id], p_pos)
        if length(dinfs) > 0
            ipairs = hcat(ipairs, [0; id; customer_contact])
        end
        if (sim["infection_status"][id] .!= Susc) .* (sim["infection_status"][id] .!= Recd)
            inf_scale = sim["infection_profiles"][id][i_day - sim["inf_time"][id] + 1]
            p_inf = 1 - exp(-infection_rate_F2F*inf_scale*t_del)
            cinfs = randsubseq(1:Nassignments[id], p_inf)
            if length(cinfs) > 0
                contact_infections = hcat(contact_infections, [id; length(cinfs)])
            end
        end
    end
    # print(contact_infections,'\n')
    #infectious packages
    pkg_infs = randsubseq(1:length(InfPkgs), PkgParams["p_fomite_contr"])
    infectors = unique(InfPkgs[pkg_infs])
    if length(infectors) > 0
        pkg_infections = zeros(Int64,(2,length(infectors)))
        pkg_infections[1,:] .= infectors
        for j in 1:length(infectors)
            pkg_infections[2,j] = sum(InfPkgs[pkg_infs] .== infectors[j])
        end
    end
    infpairs = hcat(infpairs, ipairs)
    #infpairs contains incoming infections in usual format
    #contact infections: first row infectors, second row number of customers infected
    #pkg infections: first row infectors, second row number of customers infected
    return infpairs, contact_infections, pkg_infections
end

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

function do_infections_randomly!(infpairs::Array{Int,2}, sim::Dict, i_day::Int,
                                 workpairs::Array{Int64,2} = Array{Int64,2}(undef,2,0))
    N = size(infpairs,2)
    nkeep = Array{Int64,1}()
    is_pairs = (length(workpairs) > 0)
    if N > 0
        ind = randperm(N)
        for i in ind
            k = infpairs[2,i]   #node to be infected
            if sim["infection_status"][k] .== Susc  #if susceptible
                if is_pairs
                    partner = find_partner(workpairs, k)
                    infect_node!(sim, k, i_day, partner)
                else
                    infect_node!(sim, k, i_day)
                end
                push!(nkeep, i)
            end
        end
    end
    infpairs_kept = infpairs[:,nkeep]
    return infpairs_kept
end


function update_testing_state!(sim::Dict, i_day::Int)
    nr = 1:sim["Ntot"]
    testing_paused = nr[sim["testing_paused"] .== true]
    resume_testing = testing_paused[i_day .> sim["resume_testing"][testing_paused]]
    sim["testing_paused"][resume_testing] .= false
end

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

function copy_graph_into_lightgraph(g::Graphs.GenericGraph)
    ges = collect(Graphs.edges(g))
    lg = LightGraphs.SimpleGraph(Graphs.num_vertices(g))

    LightGraphs.add_edge!.(Ref(lg), Graphs.source.(ges), Graphs.target.(ges))
    return lg
end

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

function plot_infection_profiles(N::Int, beta::Float64)
    profiles = Array{Array{Float64,1},1}(undef,N)
    VLprofiles = Array{Array{Float64,1},1}(undef,N)
    ip = generate_incubation_periods(N)
    alphas = gamma_ab_mode.(ip,Ref(beta))
    mags = generate_log10_viral_load(N)
    x = zeros(N)
    tl = 0
    for i in 1:N
        profiles[i] = calculate_infectivity(mags[i], alphas[i], beta)
        tl += length(profiles[i])
    end
    tl /= N
    pp = Plots.plot(profiles[1], xlabel = "Days", ylabel = "Infectiousness")
    for i in 2:N
        pp = Plots.plot!(profiles[i])
    end

    Plots.display(pp)
    return profiles
end

function init_pairs!(sim::Dict, PairParams::Dict)
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

function sim_loop!(sim::Dict, sim_summary::Dict, i_day::Int, occ::Float64,  np::Int,
        Params::Dict, PkgParams::Dict, PairParams::Dict, TestParams::Dict,
        is_pairs::Bool, is_network::Bool, is_testing::Bool, next_test::Int)

    update_infectivity!(sim, i_day)
    if i_day == next_test
        if is_pairs && PairParams["pair_isolation"]
            do_testing!(sim, TestParams, i_day, sim["fixed_job_pairings"])
        else
            do_testing!(sim, TestParams, i_day)
        end
    end
    update_isolation!(sim, i_day)
    infpairs = Array{Int64,2}(undef,3,0)
    driver_pairs = Array{Int64,2}(undef,2,0)
    loader_pairs = Array{Int64,2}(undef,2,0)
    NAs = Array{Int64,1}(undef,0)
    if sim["Incidence"][i_day] > 0
       infpairs = get_introductions(infpairs, sim, i_day)
    end
    if is_pairs
        infpairs, driver_pairs, loader_pairs, NAs =
                get_pair_infections(infpairs, sim, occ, i_day, PairParams, np)
    else
        NAs = update_in_work!(sim, occ, np)
    end
    #insert parcels here
    infpairs, pkg_infectors,  drivers_delivering_ipkgs =
        get_package_infections(infpairs, sim, np, i_day, PkgParams, NAs)
    if is_network
        infpairs = get_network_infections(infpairs, sim, Params, i_day)
    end
    infpairs = get_contact_infections(infpairs, sim, Params, i_day)

    if length(driver_pairs) > 0
        t_del = BulkTimesPerDelivery["t_doorstep"]
    else
        t_del = ParcelTimesPerDelivery["t_doorstep"]
    end
    infpairs, cust_infections_drivers, cust_infections_packages =
         get_customer_infections(infpairs, sim, Params, PkgParams,
                                 NAs, t_del, vcat(transpose(drivers_delivering_ipkgs),
                                 transpose(pkg_infectors)), i_day)


    # print(sum(cust_infections_drivers[2,:]),'\n')
    if is_pairs && PairParams["pair_isolation"]
        infpairs = do_infections_randomly!(infpairs, sim, i_day, sim["fixed_job_pairings"])
    else
        infpairs = do_infections_randomly!(infpairs, sim, i_day)
    end

    update_sim_summary!(sim_summary, sim, i_day, length(pkg_infectors), infpairs;
                        CinfsD = sum(cust_infections_drivers[2,:]),
                        CinfsP = sum(cust_infections_packages[2,:]))

    return infpairs, hcat(driver_pairs, loader_pairs)
end

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

#run a either single outbreak with 1 initial case, or a sim with no cases but introductions
function run_sim(Params::Dict, OccPerDay::Array{Float64,1},
     PkgParams::Dict, NPPerDay::Array{Int64,1}; is_pairs::Bool=false,
     PairParams::Dict=Dict(), is_network::Bool=false, degree_logmean::Float64=1.67,
     degree_logstd::Float64=0.97, visualise::Bool = false, testing::Bool=false,
     TestParams::Dict=Dict(), Incidence::Array{Float64,1} = zeros(length(OccPerDay)),
     Prevalence::Array{Float64,1} = zeros(length(OccPerDay)))

    if visualise
        sim, node_x, node_y = initialise_withvis(Params,
          PairParams, degree_logmean, degree_logstd, is_network, is_pairs, Incidence, Prevalence)
    else
        sim = initialise_novis(Params,
          PairParams, degree_logmean, degree_logstd, is_network, is_pairs, Incidence, Prevalence)
    end
    i_day = rand(1:7)
    if testing
        if (TestParams["protocol"] == PCR_mass_protocol || TestParams["protocol"] == LFD_mass_protocol)
            test_days, test_day_counter = init_testing!(sim,TestParams,i_day,length(OccPerDay))
        end  #add options for other protocols here
    else
        test_day_counter = 1
        test_days = [0]
    end
    Anyinf = true
    if Params["SimType"] == Outbreak_sim
        sim_summary = sim_setup!(sim, Params["InfInit"], i_day, length(OccPerDay))
        Anyinf = any(((sim["infection_status"] .== Susc)
              .+ (sim["infection_status"] .== Recd)) .== 0)
    elseif Params["SimType"] == Scenario_sim
        sim_summary = scenario_sim_setup!(sim, Incidence, Prevalence, i_day, length(OccPerDay))
    end
    while Anyinf && (i_day <= length(OccPerDay))
        infpairs, pairs = sim_loop!(sim, sim_summary, i_day, OccPerDay[i_day],
                                    NPPerDay[i_day], Params, PkgParams, PairParams, TestParams,
                                    is_pairs, is_network, testing, test_days[test_day_counter])
        if testing && (i_day == test_days[test_day_counter])
            test_day_counter += 1
        end
        fname = Printf.@sprintf("infection_network_%03d.png",i_day)
        if visualise
            print_infection_network(sim, fname, infpairs, node_x, node_y, pairs)
        end
        if Params["SimType"] == Outbreak_sim
            Anyinf = any(((sim["infection_status"] .== Susc)
                  .+ (sim["infection_status"] .== Recd)) .== 0)
        end
        i_day += 1
    end
    trim_sim_summary!(sim_summary, i_day-1, length(OccPerDay))

    return sim_summary
end
