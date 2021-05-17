using DataFrames
using CSV

using Distributions
using Random
using StatsBase

#testing protocol definitions
const PCR_mass_protocol = "PCR_mass_testing"
const LFD_mass_protocol = "LFD_mass_testing"

#VL trajectory params (Kissler)
const alpha_onset = 1.1
const beta_onset = 0.7
const alpha_peak = 2.3
const beta_peak = 0.7
const alpha_decay = 2.4
const beta_decay = 0.3
const peakVL_mean = 7.533
const peakVL_sd = 1.164
#PCR swab params (invented)
const swab_peak_accuracy = 0.97
const swab_midpoint = 2.0
const swab_disp = 1.0
const VL_cutoff = 0.0
#LFD params (Porton Down)
const LFD_PDcurve_b0 = -3.60
const LFD_PDcurve_b1 = 1.22
#VL ~ infectivity relation (Marks)
const inf_dep = 1.3
const VL_ref = mean(Normal(peakVL_mean,peakVL_sd))
const PIsigma = 1.0
const PImu = -0.5*PIsigma^2
#Viral load where people stop being infectious
const inf_VL_cutoff = 0.0 #3.0 -- should only have minor effect
#correlation gradient between peak_inf and pasymp
const pasymp_vl7 = 0.62
const pasymp_vl89 = 0.5
const pasymp_vl10 = 0.34
const asymp_frac = cdf(Normal(peakVL_mean,peakVL_sd),7.0)*pasymp_vl7 +
    (cdf(Normal(peakVL_mean,peakVL_sd),9.0) - cdf(Normal(peakVL_mean,peakVL_sd),7.0))*pasymp_vl89 +
    (1 - cdf(Normal(peakVL_mean,peakVL_sd),9.0))*pasymp_vl10
#relative infectivity of asymptomatics
const mean_asymp = 0.5

function generate_peak_viral_loads(Ntot::Int)
    return rand(Normal(peakVL_mean,peakVL_sd),Ntot)
end

function generate_onset_times(Ntot::Int)
    return rand(Gamma(alpha_onset,1/beta_onset),Ntot)
end

function generate_peak_times(Ntot::Int)
    return rand(Gamma(alpha_peak,1/beta_peak),Ntot)
end

function generate_decay_times(Ntot::Int)
    return rand(Gamma(alpha_decay,1/beta_decay),Ntot)
end

function build_viral_load_distributions!(sim::Dict)
    Ntot = sim["Ntot"]
    PVLs = generate_peak_viral_loads(Ntot)
    OTs = generate_onset_times(Ntot)
    PTs = generate_peak_times(Ntot)
    DTs = generate_decay_times(Ntot)
    v = zeros.(Int64.(ceil.(OTs .+ PTs .+ DTs .+ 1)))
    for n in 1:Ntot
        m_up = PVLs[n]/PTs[n]
        m_down = PVLs[n]/DTs[n]
        T = length(v[n])
        i = 1:T
        t = i .- 1
        i1 = Int64(ceil(OTs[n]))
        i2 = Int64(ceil(OTs[n] + PTs[n]))
        i3 = Int64(ceil(OTs[n] + PTs[n] + DTs[n]))
        cond1 = ((i .> i1) .* (i .<= i2))
        v[n][cond1] = m_up .* (t[cond1] .- OTs[n])
        cond2 = ((i .> i2) .* (i .<= i3))
        v[n][cond2] = PVLs[n] .- m_down .* (t[cond2] .- OTs[n] .- PTs[n])
    end
    sim["VL_mag"] = PVLs
    sim["symp_time"] = OTs .+ PTs
    sim["VL_profiles"] = v
    sim["decay_time"] = DTs
end

function generate_peak_inf(peak_VL::Float64, asymp::Bool)
    #assume reference person has average infectivity of 1 (so peak of 2)
    value = 2.0 * inf_dep^(peak_VL - VL_ref)
    ascale = mean_asymp*asymp_frac + (1 - asymp_frac)
    r = rand(LogNormal(PImu,PIsigma))
    if asymp
        r = mean_asymp * r
    end
    return r*value/ascale
end

function infectivity(VL::Array{Float64,1}, peak_inf::Float64, peak_VL::Float64)
    #assume infectivity scales linearly over time with log10 copies/ml above cutoff
    j = zeros(length(VL))
    j[VL .> inf_VL_cutoff] =  peak_inf .* (VL[VL .> inf_VL_cutoff] .- inf_VL_cutoff) ./ (peak_VL - inf_VL_cutoff)
    return j
end

function logistic_function(x::Float64, x0::Float64, k::Float64, ymax::Float64)
    return (ymax / (1  + exp(-k*(x-x0))))
end

function swab_positive_prob(VL::Float64)
    #given viral load VL, assume this is the probability that some virus makes it onto the swab -- this could be improved by a more mechanistic testing model
    return logistic_function(VL, swab_midpoint, swab_disp, swab_peak_accuracy)
end

function PCRtest_positive_prob(VL::Float64)
    #assume PCR is 100% accurate at detecting viral material on swab below Ct = 40
    if VL > VL_cutoff
        return swab_positive_prob(VL)
    else
        return 0
    end
end

function LFDtest_positive_prob(VL::Float64, sens_rel::Float64)
    #sens rel is 1 in Peto, could scale this for self-administering effects (sens rel < 1)
    if VL > VL_cutoff
        return (sens_rel * swab_positive_prob(VL) /
               (1 + exp(-LFD_PDcurve_b1*VL - LFD_PDcurve_b0)))
    else
        return 0
    end
end

function generate_isolations(Ntot::Int, Pisol::Float64)
    #randomly generate if people obey isolation or not
    return randsubseq(1:Ntot, Pisol)
end

function generate_asymptomatics(VL::Array{Float64,1})
    #randomly generate if people are asymptomatic
    nr = 1:length(VL)
    nVL7 = nr[VL .< 7]
    nVL89 = nr[(VL .>= 7) .* (VL .<= 10)]
    nVL10 = nr[VL .> 10]
    return vcat(randsubseq(nVL7,pasymp_vl7),randsubseq(nVL89,pasymp_vl89),
                randsubseq(nVL10,pasymp_vl10))
end

function init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)
    #simulate Ntot people with baseline isolation probability Pisol
    sim = Dict("Ntot"=>Ntot, "symp_time"=>-ones(Int64, Ntot),
               "asymptomatic"=>zeros(Bool, Ntot),
               "isolation_time"=>zeros(Int64, Ntot),
               "will_isolate"=>zeros(Bool, Ntot),
               "inf_mag"=>zeros(Float64, Ntot),
               "infection_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "VL_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "days_infectious" => zeros(Int64,Ntot))
    build_viral_load_distributions!(sim)
    sim["symp_day"] = Int64.(round.(rand(0:1,Ntot) .+ sim["symp_time"]))
    sim["asymptomatic"][generate_asymptomatics(sim["VL_mag"])] .= true
    sim["inf_mag"] .= generate_peak_inf.(sim["VL_mag"], sim["asymptomatic"])
    sim["will_isolate"][generate_isolations(Ntot, Pisol)] .= true
    nr = 1:sim["Ntot"]
    sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]
    sim["will_isolate"][sim["asymptomatic"]] .= false #asymptomatics don't self isolate, but are not "non isolators"
    sim["infection_profiles"] .= infectivity.(sim["VL_profiles"], sim["inf_mag"], sim["VL_mag"])
    sim["days_infectious"] = length.(sim["infection_profiles"])

    return sim
end

function init_testing!(sim::Dict, testing_params::Dict, i_day::Int, Ndays::Int)
    #add test positivity profiles to simulation Dict, i_day is start day, Ndays is total length of sim
    #testing params contains "tperiod" (days between test)
    #testing params contains "protocol" (LFD_mass_protocol or PCR_mass_protocol)
    #optional: sens_rel: relative sensitivity of LFD as VL -> infinity
    sim["will_isolate_with_test"] = ones(Bool,sim["Ntot"])
    sim["will_isolate_with_test"][sim["non_isolators"]] .= false
    new_compliers = randsubseq(sim["non_isolators"], testing_params["new_comply_prob"])
    if length(new_compliers) > 0
        sim["will_isolate_with_test"][new_compliers] .= true
    end
    test_day0 = rand(1:testing_params["tperiod"])
    test_days = collect(test_day0:testing_params["tperiod"]:Int64(ceil(Ndays)))
    test_days = push!(test_days, test_days[end] + testing_params["tperiod"])
    test_day_counter = 1 + sum(test_days .< i_day)
    sim["test_pos_profiles"] = Array{Array{Float64,1},1}(undef,sim["Ntot"])
    if testing_params["protocol"] == PCR_mass_protocol
        for i in 1:sim["Ntot"]
            sim["test_pos_profiles"][i] = PCRtest_positive_prob.(sim["VL_profiles"][i])
        end
    elseif testing_params["protocol"] == LFD_mass_protocol
        if haskey(testing_params,"sens_rel")
            sr = testing_params["sens_rel"]
        else
            sr = 1.0
        end
        for i in 1:sim["Ntot"]
            sim["test_pos_profiles"][i] = LFDtest_positive_prob.(sim["VL_profiles"][i], Ref(sr))
        end
    end

    return test_days, test_day_counter
end
