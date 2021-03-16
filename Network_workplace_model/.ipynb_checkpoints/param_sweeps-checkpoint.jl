include("network_transmission_workplace.jl")

function init_results_dataframe(Nrows::Int, AllParams::Dict)
    results = DataFrame([String, Int64, Int64, Float64,
                         Int64, Float64],
                         [:Group, :NStaff, :Iteration, :FracRecovered,
                         :TotInfPackagesDelivered, :FomiteInfectionFrac, 
                         :NetworkInfectionFrac, :ContactInfectionFrac,
                         :PairInfectionFrac, :RoomInfectionFrac], Nrows)
    
    if AllParams["sim_type"] == Outbreak_sim      
        insertcols!(results,[:IndexCaseInfections, :OverallOutbreakLength] => 
                             zeros(Int64,Nrows))         
    end
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
    if Params["sim_type"] == Outbreak_sim
        results[(irow_start):(irow_start+3),"IndexCaseInfections"] .= SimOutput["IndexCaseInfections"]
        results[(irow_start):(irow_start+3),"OverallOutbreakLength"] .= length(SimOutput["time"])
    end
    
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
    results[(irow_start+3),"FomiteInfectionFrac"] = sum(SimOutput["FomiteInfs"]) / NStot
    results[(irow_start+3),"NetworkInfectionFrac"] = sum(SimOutput["NetworkInfs"]) / NStot
    results[(irow_start+3),"ContactInfectionFrac"] = sum(SimOutput["ContactInfs"]) / NStot
    results[(irow_start+3),"PairInfectionFrac"] = sum(SimOutput["PairInfs"]) / NStot
    results[(irow_start+3),"RoomInfectionFrac"] = sum(SimOutput["RoomInfs"]) / NStot
    results[(irow_start+3),"CustomerIntroFrac"] = sum(SimOutput["CustomerIntroductions"]) / NStot
    results[(irow_start+3),"ExtIntroFrac"] = sum(SimOutput["ExternalIntroductions"]) / NStot
    results[(irow_start+3),"TotInfPackagesDelivered"] = sum(SimOutput["PackagesInfectiousOnDelivery"])
    results[(irow_start+3), "CustomerInfectionsPkgs"] .= sum(SimOutput["CustomersInfectedByPkgs"])
    results[(irow_start+3), "CustomerInfectionsDrivers"] .= sum(SimOutput["CustomersInfectedByDrivers"])
    results[(irow_start+3), "IsolFrac"] .= sum(SimOutput["Isolated"])/(NStot * isol_time)
    results[(irow_start+3), "IsolTestFrac"] .= sum(SimOutput["IsolatedDueToTest"])/(NStot * isol_time)
end


#Check this works
function run_many_sims(ParamsVec::Array{Dict,1}, Nrepeats::Int,
                OccPerDay::Array{Float64,1}, PkgParams::Array{Dict,1};
                NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay)),
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

function run_param_sweep_outbreak_parcel()
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
    PkgParams = Array{Dict,1}(undef,0)
    for pi in PIsol
        for pf in PFC
            for ii in II
                for i in 1:3
                    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                         "p_contact"=>pc, "Pisol"=>pi,
                                         "InfInit"=>ii, "tD"=>tD[i], "phi"=>Phi[i], 
                                         "p_friend_contact"=>pf, "type"=>Outbreak_sim))
                    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
                end
            end
        end
    end

    run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; filename="param_sweep.csv")
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
    PkgParams = Array{Dict,1}(undef,0)
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
                    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
                end
            end
        end
    end

    run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, filename="testing_loops.csv")
end
