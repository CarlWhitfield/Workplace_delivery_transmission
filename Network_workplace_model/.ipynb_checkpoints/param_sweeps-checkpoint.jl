include("network_transmission_workplace.jl")

function init_results_dataframe(Nrows::Int, AllParams::Dict)
    results = DataFrame([String, Int64, Int64, Float64,
                         Int64, Int64, Float64, Int64],
                         [:Group, :NStaff, :Iteration, :FracRecovered,
                         :TotInfPackagesDelivered, :IndexCaseInfections,
                         :FomiteInfectionFrac, :OverallOutbreakLength], Nrows)
    for p in keys(AllParams)
        if p != "ND" && p != "NO" && p != "NL"
            results[!,p] = fill(AllParams[p], Nrows)
        end
    end

    return results
end

function add_to_results_dataframe!(results::DataFrame, Params::Dict, SimOutput::Dict,
                                   irow_start::Int, Niteration::Int)

    #params and results common to all staff
    for p in keys(Params)
        if p in names(results)
            results[(irow_start):(irow_start+3),p] .= Params[p]
        end
    end
    results[(irow_start):(irow_start+3),"Iteration"] .= Niteration
    results[(irow_start):(irow_start+3), "TotInfPackagesDelivered"] .=
                                 sum(SimOutput["PackagesInfectiousOnDelivery"])
    results[(irow_start):(irow_start+3),"IndexCaseInfections"] .= SimOutput["IndexCaseInfections"]
    results[(irow_start):(irow_start+3),"OverallOutbreakLength"] .= length(SimOutput["time"])
    #CONTINUE HERE
    #params and results for specific staff groups
    results[(irow_start),"Group"] = "Drivers"
    results[(irow_start), "NStaff"] = Params["ND"]
    results[(irow_start),"FracRecovered"] = SimOutput["Recovered"][1,end] / Params["ND"]
    results[(irow_start),"FomiteInfectionFrac"] = SimOutput["NFomiteInfs"][1]/Params["ND"]
    results[(irow_start+1),"Group"] = "Pickers"
    results[(irow_start+1), "NStaff"] = Params["NL"]
    results[(irow_start+1),"FracRecovered"] = SimOutput["Recovered"][2,end] / Params["NL"]
    results[(irow_start+1),"FomiteInfectionFrac"] = SimOutput["NFomiteInfs"][2]/Params["NL"]
    results[(irow_start+2),"Group"] = "Office"
    results[(irow_start+2), "NStaff"] = Params["NO"]
    results[(irow_start+2),"FracRecovered"] = SimOutput["Recovered"][3,end] / Params["NO"]
    results[(irow_start+2),"FomiteInfectionFrac"] = SimOutput["NFomiteInfs"][3] / Params["NO"]
    results[(irow_start+3),"Group"] = "All"
    NStot = Params["ND"] + Params["NL"] + Params["NO"]
    results[(irow_start+3), "NStaff"] = NStot
    results[(irow_start+3),"FracRecovered"] = sum(SimOutput["Recovered"][:,end]) / NStot
    results[(irow_start+3),"FomiteInfectionFrac"] = sum(SimOutput["NFomiteInfs"]) / NStot
end

#needs updating to new format
function run_param_sweeps_simplest(ND::Int, NL::Int, NO::Int, p_asymp::Float64,
        OccPerDay::Array{Float64,1}, PIsol::Array{Float64,1}, PInf::Array{Float64,1},
        PContact::Array{Float64,1}, TD::Array{Float64,1}, Phi::Array{Float64,1},
        Nrepeats::Int, InfInit::Int; filename="output.csv", pf_trans::Float64=0.0,
        pf_contr::Float64=0.0, Dtime::Float64=1/6, Ltime::Float64=1/6, Hlife::Float64=0.5,
        NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay)))
     NPIsol = length(PIsol)
     NPInf = length(PInf)
     NPC = length(PContact)
     NPhi = length(Phi)
     NTD = length(TD)
     Nrows = NPIsol*NPInf*NPC*NPhi*NTD*Nrepeats*3
     results = DataFrame([Float64,Float64,Float64,Float64,Float64,Int64,String,Float64,
                          Float64],
                         [:Isolation_Prob, :Infection_Prob, :Contact_Prob,
                          :Driver_Mixing_Time, :Mixing_Param, :Iteration,
                          :Group, :Frac_recovered, :Inf_packages_delivered], Nrows)
     Params = Dict("ND"=>ND, "NL"=>NL, "NO"=>NO, "InfInit"=>InfInit)
     PkgParams = Dict("p_fomite_contr"=>pf_contr, "p_fomite_trans"=>pf_trans,
                       "PkgHlife"=>Hlife, "Dtime"=>Dtime, "Ltime"=>Ltime)
     m_step = 3*Nrepeats
     l_step = NPhi*m_step
     k_step = NTD*l_step
     j_step = NPC*k_step
     i_step = NPInf*j_step
     for (i, p_isol) in enumerate(PIsol)
         Params["Pisol"] = p_isol
         i_ind_start = (i-1)*i_step
         results[(i_ind_start + 1):(i_ind_start + i_step),"Isolation_Prob"] .= p_isol
         for (j, p_inf) in enumerate(PInf)
             Params["p_inf"] = p_inf
             j_ind_start = i_ind_start + (j-1)*j_step
             results[(j_ind_start+1):(j_ind_start+j_step), "Infection_Prob"] .= p_inf
             for (k, p_contact) in enumerate(PContact)
                 Params["p_contact"] = p_contact
                 k_ind_start = j_ind_start + (k-1)*k_step
                 results[(k_ind_start+1):(k_ind_start + k_step), "Contact_Prob"] .= p_contact
                 for (l, tD) in enumerate(TD)
                     Params["tD"] = tD
                     l_ind_start = k_ind_start + (l-1)*l_step
                     results[(l_ind_start+1):(l_ind_start + l_step), "Driver_Mixing_Time"] .= tD
                     for (m, phi) in enumerate(Phi)
                         Params["phi"] = phi
                         m_ind_start = l_ind_start + (m-1)*m_step
                         print(Int64(m_ind_start/m_step)+1,'/', Int64(Nrows/m_step),'\n')
                         results[(m_ind_start+1):(m_ind_start+m_step), "Mixing_Param"] .= phi
                         @threads for n in 1:Nrepeats
                             index_start = m_ind_start + (n-1)*3
                             out = run_outbreak_sim(Params, OccPerDay, PkgParams,
                                                NPPerDay)
                             results[(index_start+1):(index_start+3), "Iteration"] .= n
                             results[index_start+1, "Frac_recovered"] = out["Recovered"][1,end]/ND
                             results[index_start+2, "Frac_recovered"] = out["Recovered"][2,end]/NL
                             results[index_start+3, "Frac_recovered"] = out["Recovered"][3,end]/NO
                             if pf_trans > 0
                                results[(index_start+1):(index_start+3),
                                      "Inf_packages_delivered"] .=
                                          sum(out["PackagesInfectiousOnDelivery"])
                             else
                                results[(index_start+1):(index_start+3),
                                      "Inf_packages_delivered"] .= 0
                             end
                         end
                    end
                end
            end
        end
     end
     results[:, "Group"] .= repeat(job_names, Int64(Nrows/3))

     CSV.write(filename, results)
     return results
end




function run_many_sims(ParamsVec::Array{Dict,1}, Nrepeats::Int,
                OccPerDay::Array{Float64,1}, PkgParams::Array{Dict,1};
                NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay))
                IsNetwork::Array{Bool,1} = ones(Bool,length(ParamsVec)),
                IsPairs::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                PairParams::Array{Dict,1} = fill(Dict(),length(ParamsVec)),
                IsTesting::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                TestingParams::Array{Dict,1} = fill(Dict(),length(ParamsVec)),
                filename="output.csv")

     NParamSets = length(ParamsVec)
     Nrows = 4*NParamSets*Nrepeats
     TestParams = merge(ParamsVec[1],PkgParams[1],PairParams[1],TestingParams[1])
     results = init_results_dataframe(Nrows, TestParams)

     i_step = 4*Nrepeats
     for (i, Params) in enumerate(ParamsVec)
         i_ind_start = (i-1)*i_step + 1
         print(i,'\n')
         @threads for n in 1:Nrepeats
             index_start = i_ind_start + (n-1)*4
             out = run_sim(Params, OccPerDay, PkgParams[i], NPPerDay;
                   is_network=IsNetwork[i], is_pairs=IsPairs[i], PairParams=PairParams[i],
                   testing=IsTesting[i], TestParams=TestingParams[i])
             AllParams = merge(ParamsVec[i],PkgParams[i],PairParams[i],TestingParams[i])
             add_to_results_dataframe!(results, Params, out, index_start, n)
         end
     end

     CSV.write(filename, results)
     return results
end



function run_param_sweep_outbreak()
    NWeeks= 52
    NPh = 3000
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(round(NPh/80))
    NLh = Int64(round(NPh/150))
    NOh = Int64(round(NPh/300))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #other params
    PIsol = [0.25,0.5,0.75,1.0]
    Phi = [0.05,0.25,1.0]
    PFC = [0.5,1.0]
    II = [1,2,3]
    tD = [0.1, 0.5, 1.0]
    Nrepeats = 1000


    ParamVec = Array{Dict,1}(undef,0)
    for pi in PIsol
        for pf in PFC
            for ii in II
                for i in 1:3
                    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                         "p_contact"=>pc, "Pisol"=>pi,
                                         "InfInit"=>ii, "tD"=>tD[i], "phi"=>Phi[i], 
                                         "p_friend_contact"=>pf, "type"=>Outbreak_sim))
                    PkgParams = Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5)
                end
            end
        end
    end

    run_param_sweeps_network(ParamVec, Nrepeats, OccPattern; filename="param_sweep.csv")
end

function run_testing_script()
    NWeeks= 52
    NPh = 3000
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(round(NPh/80))
    NLh = Int64(round(NPh/150))
    NOh = Int64(round(NPh/300))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #other params
    pasymp = 0.3
    PIsol = [0.25,0.5,1.0]
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14
    Nrepeats = 1000


    Delay = [0,0,0,1,2]
    TestType = ["NoTest","LFD","PCR","PCR","PCR"]
    VL_ref = [4.5, 4.5, 1.5, 1.5, 1.5]
    VL_disp = [0.5, 0.5, 3.0, 3.0, 3.0]
    Sens_max = [0, 0.89, 0.907, 0.907, 0.907]
    Test_pause = [14, 14,90,90,90]

    ParamVec = Array{Dict,1}(undef,0)
    TestParamVec = Array{Dict,1}(undef,0)
    for pi in PIsol
        for ncp in NCP
            for tp in Tperiod
                for i in 1:length(TestType)
                    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                            "p_contact"=>pc, "Pasymp"=>pasymp, "Pisol"=>pi,
                            "InfInit"=>1, "tD"=>1.0, "phi"=>1.0, "p_friend_contact"=>1.0))
                    push!(TestParamVec, Dict("new_comply_prob"=>ncp, "tperiod"=>tp,
                          "type"=>TestType[i], "VL_ref"=>VL_ref[i], "VL_disp"=>VL_disp[i],
                          "sens_max"=>Sens_max[i], "specificity"=>0.995,
                          "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                end
            end
        end
    end

    run_param_sweeps_testing_network(ParamVec, Nrepeats, TestParamVec, OccPattern; filename="testing_loops.csv")
end
