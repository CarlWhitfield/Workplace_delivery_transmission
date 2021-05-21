using DataFrames
using CSV

using Distributions
using Random
using StatsBase
using SpecialFunctions

#testing protocol definitions
const PCR_mass_protocol = "PCR_mass_testing"
const LFD_mass_protocol = "LFD_mass_testing"

#VL trajectory params (Kissler)
const beta_onset = 0.7
const alpha_onset = 1.5*beta_onset
const beta_peak = 0.7
const alpha_peak = 3.5*beta_peak
const beta_decay = 0.3
const alpha_decay_symp = 10.5*beta_decay
const alpha_decay_asymp = 6.7*beta_decay
const peakVL_mean = 7.533
const peakVL_sd = 1.164
const VL_LOD_Kissler = 2.66  #starting VL in log10 copies/ml
#PCR swab params (Smith et al. JCMB)
const PCR_VL0 = 8.522/4.408
const PCR_VLdisp = 4.408
const PCR_sens_max = 0.85 #assumed
#LFD params (Porton Down)
const LFD_VLmean = 10.836/2.680
const LFD_VLstd = 1.0/2.680
const LFD_sens_max = 0.75 #assumed
#VL ~ infectivity relation (Marks)
const inf_dep = 1.3
const VL_ref = peakVL_mean
const PIsigma = 0.5
const PImu = -0.5*PIsigma^2
#Viral load where people stop being infectious
const inf_VL_cutoff = 6.0 #3.0 -- should only have minor effect
#this scaling means that, on average, p_inf is approx 0.06 for hour long face-to-face contacts within 14 days of infection.
const j0scale = 4.1 * 1.1 *(4.0/(exp(log(inf_dep)^2*peakVL_sd^2/2)*(1 +
                erf((peakVL_mean - inf_VL_cutoff + log(inf_dep)*peakVL_sd^2)/(sqrt(2)*peakVL_sd)))))
#correlation gradient between peak_inf and pasymp
# const pasymp_vl7 = 0.62
# const pasymp_vl89 = 0.5
# const pasymp_vl10 = 0.34
p_asymp = 0.5
#const asymp_frac = cdf(Normal(peakVL_mean,peakVL_sd),7.0)*pasymp_vl7 +
#    (cdf(Normal(peakVL_mean,peakVL_sd),9.0) - cdf(Normal(peakVL_mean,peakVL_sd),7.0))*pasymp_vl89 +
#    (1 - cdf(Normal(peakVL_mean,peakVL_sd),9.0))*pasymp_vl10
#relative infectivity of asymptomatics
const mean_asymp = 0.5
#symp time distribution
const symp_beta = 4.84 / 2.6^2
const symp_alpha = 4.84 * symp_beta

function generate_peak_viral_loads(Ntot::Int)
    return rand(Normal(peakVL_mean,peakVL_sd),Ntot)
end

function generate_onset_times(Ntot::Int)
    return rand(Gamma(alpha_onset,1/beta_onset),Ntot)
end

function generate_peak_times(Ntot::Int)
    return rand(Gamma(alpha_peak,1/beta_peak),Ntot)
end

function generate_asymptomatics(Ntot::Int)
    #randomly generate if people are asymptomatic
    return randsubseq(1:Ntot,p_asymp)
end

function generate_decay_times(Ntot::Int, Asymp::Array{Int64,1})
    DTs = rand(Gamma(alpha_decay_symp,1/beta_decay),Ntot)
    DTs[Asymp] = rand(Gamma(alpha_decay_asymp,1.0/0.3),length(Asymp))
    return DTs
end


function generate_peak_inf(peak_VL::Float64)
    #assume reference person has average infectivity of 1 (so peak of 2)
    value = j0scale * inf_dep^(peak_VL - VL_ref)
    r = rand(LogNormal(PImu,PIsigma))
    return r*value
end

function build_viral_load_distributions!(sim::Dict)
    Ntot = sim["Ntot"]
    PVLs = generate_peak_viral_loads(Ntot)
    OTs = generate_onset_times(Ntot)
    PTs = generate_peak_times(Ntot)
    Asymp = generate_asymptomatics(Ntot)
    DTs = generate_decay_times(Ntot, Asymp)
    PInfs = generate_peak_inf.(PVLs)
    STs = zeros(Ntot)
    v = zeros.(Int64.(ceil.(OTs .+ PTs .+ DTs .+ 10)))
    inf = zeros.(Int64.(ceil.(OTs .+ PTs .+ DTs .+ 1)))
    for n in 1:Ntot
        m_up = (PVLs[n] - VL_LOD_Kissler)/PTs[n]
        m_down = (PVLs[n] - VL_LOD_Kissler)/DTs[n]
        T = length(v[n])
        i = 1:T
        t = i .- 1
        i1 = Int64(ceil(OTs[n]))
        i2 = Int64(ceil(OTs[n] + PTs[n]))
        i3 = Int64(ceil(OTs[n] + PTs[n] + DTs[n]))
        cond1 = (i .<= i2)
        v[n][cond1] = VL_LOD_Kissler .+ m_up .* (t[cond1] .- OTs[n])
        cond2 = (i .> i2)
        v[n][cond2] = PVLs[n] .- m_down .* (t[cond2] .- OTs[n] .- PTs[n])
        Strunc = truncated(Gamma(symp_alpha,1.0/symp_beta), OTs[n], OTs[n] + PTs[n] + DTs[n])
        STs[n] = rand(Strunc)
        if PVLs[n] > inf_VL_cutoff
            tinf_start = OTs[n] + (inf_VL_cutoff - VL_LOD_Kissler)/m_up
            tinf_end = OTs[n] + PTs[n] + (PVLs[n] - inf_VL_cutoff)/m_down
            cum_inf = zeros(T+1)
            t_inf = collect(-1:(T-1)) .+ 0.5
            
            cond1 = (t_inf .>= tinf_start) .* (t_inf .<  OTs[n] + PTs[n])
            cum_inf[cond1] = 0.5 * PInfs[n] .* (t_inf[cond1] .- tinf_start).^2 / 
                               (OTs[n] + PTs[n] - tinf_start)
            
            cond2 = (t_inf .>= OTs[n] + PTs[n]) .* (t_inf .<=  tinf_end)
            c_inf_peak = 0.5 * PInfs[n] * (OTs[n] + PTs[n] - tinf_start)
            cum_inf[cond2] = c_inf_peak .+ PInfs[n] .* (t_inf[cond2] .- OTs[n] .- PTs[n]) .*
                    (1.0 .- 0.5 .* (t_inf[cond2] .- OTs[n] .- PTs[n]) ./ (tinf_end  - OTs[n] - PTs[n]))
            cum_inf[t_inf .>  tinf_end] .= 0.5*PInfs[n]*(tinf_end - tinf_start)
            inf[n] = cum_inf[2:(T+1)] .- cum_inf[1:T]
            nel = 1:length(t)
            deleteat!(inf[n], nel[t .> Int64(ceil(tinf_end))])
        end
    end
    sim["VL_mag"] = PVLs
    sim["symp_time"] = STs
    sim["VL_profiles"] = v
    sim["infection_profiles"] = inf
    sim["asymptomatic"][Asymp] .= true
    sim["days_infectious"] = length.(sim["infection_profiles"])
end


# function infectivity(VL::Array{Float64,1}, peak_inf::Float64, peak_VL::Float64)
#     j = zeros(length(VL))
#     t = 0:(length(VL)-1)
#     cond1 = (VL .> inf_VL_cutoff)
#     j[cond1] = peak_inf .* (VL[cond1] .-  inf_VL_cutoff) ./ (peak_VL - inf_VL_cutoff)
#     return j
# end

# # function infectivity_alt(VL::Array{Float64,1}, peak_inf::Float64, peak_VL::Float64, peak_VL_time::Float64, decay_time::Float64)
# #     #assume infectivity scales linearly over time with log10 copies/ml above cutoff
# #     j = zeros(length(VL))
# #     t = 0:(length(VL)-1)
# #     cond1 = (VL .> inf_VL_cutoff)
# #     j[cond1] = peak_inf .* (VL[cond1] .- inf_VL_cutoff) ./ (peak_VL - inf_VL_cutoff)
# #     cond2 = (VL .> inf_VL_cutoff) .* (t .>= peak_VL_time)
# #     j[cond2] = peak_inf .* (1 .+ (peak_VL_time .- t[cond2])) ./ (0.38 * decay_time)

# #     return j
# # end

function logistic_function(x::Float64, x0::Float64, k::Float64, ymax::Float64)
    return (ymax / (1  + exp(-k*(x-x0))))
end

function probit_function(x::Float64, x0::Float64, xstd::Float64, ymax::Float64)
    return ymax * cdf(Normal(x0,xstd),x)
end

function PCRtest_positive_prob(VL::Float64)
    #assume PCR is 100% accurate at detecting viral material on swab below Ct = 40
    if VL > VL_LOD_Kissler
        return logistic_function(VL, PCR_VL0, PCR_VLdisp, PCR_sens_max)
    else
        return 0
    end
end

function LFDtest_positive_prob(VL::Float64, sens_rel::Float64)
    #sens rel is 1 in Peto, could scale this for self-administering effects (sens rel < 1)
    return probit_function(VL, LFD_VLmean, LFD_VLstd, LFD_sens_max)
end

function generate_isolations(Ntot::Int, Pisol::Float64)
    #randomly generate if people obey isolation or not
    return randsubseq(1:Ntot, Pisol)
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
    sim["will_isolate"][generate_isolations(Ntot, Pisol)] .= true
    nr = 1:sim["Ntot"]
    sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]
    sim["will_isolate"][sim["asymptomatic"]] .= false #asymptomatics don't self isolate, but are not "non isolators"

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
