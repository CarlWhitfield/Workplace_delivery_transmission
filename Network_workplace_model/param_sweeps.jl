include("network_transmission_workplace.jl")

function init_results_dataframe(Nrows::Int, AllParams::Dict)
    results = DataFrame([String, Int64, Int64, Float64,
                         Int64, Float64, Float64, Float64, Float64,
                         Float64, Float64, Float64, Int64, Int64,
                         Float64, Float64, Float64, Float64, Float64],
                         [:Group, :NStaff, :Iteration, :FracRecovered,
                         :TotInfPackagesDelivered, :FomiteInfectionFrac,
                         :NetworkInfectionFrac, :ContactInfectionFrac,
                         :PairInfectionFrac, :RoomInfectionFrac,
                         :ExtIntroFrac, :CustIntroFrac, :CustomersInfectedByPkg,
                         :CustomersInfectedByDrivers, :IsolatorsFrac, :SympIsolatorsFrac,
                         :FPIsolatorsFrac, :TPSympIsolatorsFrac, :TPAsympIsolatorsFrac], Nrows)

    if AllParams["SimType"] == Outbreak_sim
        results.IndexCaseInfections = zeros(Int64,Nrows)
        results.IndexCaseViralLoad = zeros(Float64,Nrows)
        results.IndexCaseInfectivity = zeros(Float64,Nrows)
        results.OverallOutbreakLength = zeros(Int64,Nrows)
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
    if Params["SimType"] == Outbreak_sim
        results[(irow_start):(irow_start+3),"IndexCaseInfections"] .= SimOutput["IndexCaseInfections"]
        results[(irow_start):(irow_start+3),"IndexCaseViralLoad"] .= SimOutput["IndexCaseViralLoad"]
        results[(irow_start):(irow_start+3),"IndexCaseInfectivity"] .= SimOutput["IndexCaseInfectivity"]
        results[(irow_start):(irow_start+3),"OverallOutbreakLength"] .= length(SimOutput["time"])
    end

    #CONTINUE HERE
    #params and results for specific staff groups

    Names = ["Drivers","Pickers","Office"]
    Nh = [Params["ND"], Params["NL"], Params["NO"]]
    for i in 1:3
        results[(irow_start + i - 1),"Group"] = Names[i]
        results[(irow_start + i - 1), "NStaff"] = Nh[i]
        results[(irow_start + i - 1),"FracRecovered"] = SimOutput["Recovered"][i,end] / Nh[i]
        results[(irow_start + i - 1),"FomiteInfectionFrac"] = sum(SimOutput["FomiteInfs"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"NetworkInfectionFrac"] = sum(SimOutput["NetworkInfs"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"ContactInfectionFrac"] = sum(SimOutput["ContactInfs"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"PairInfectionFrac"] = sum(SimOutput["PairInfs"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"RoomInfectionFrac"] = sum(SimOutput["RoomInfs"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"CustIntroFrac"] = sum(SimOutput["CustomerIntroductions"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"ExtIntroFrac"] = sum(SimOutput["ExternalIntroductions"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"IsolatorsFrac"] = sum(SimOutput["NewIsolators"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"SympIsolatorsFrac"] = sum(SimOutput["NewSympIsolators"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"FPIsolatorsFrac"] = sum(SimOutput["NewFalseIsolators"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"TPSympIsolatorsFrac"] = sum(SimOutput["NewTestSympIsolators"][i,:])/ Nh[i]
        results[(irow_start + i - 1),"TPAsympIsolatorsFrac"] = sum(SimOutput["NewTestAsympIsolators"][i,:])/ Nh[i]
    end

    results[(irow_start+3),"Group"] = "All"
    NStot = Params["ND"] + Params["NL"] + Params["NO"]
    results[(irow_start+3), "NStaff"] = NStot
    results[(irow_start+3),"FracRecovered"] = sum(SimOutput["Recovered"][:,end]) / NStot
    results[(irow_start+3),"FomiteInfectionFrac"] = sum(SimOutput["FomiteInfs"]) / NStot
    results[(irow_start+3),"NetworkInfectionFrac"] = sum(SimOutput["NetworkInfs"]) / NStot
    results[(irow_start+3),"ContactInfectionFrac"] = sum(SimOutput["ContactInfs"]) / NStot
    results[(irow_start+3),"PairInfectionFrac"] = sum(SimOutput["PairInfs"]) / NStot
    results[(irow_start+3),"RoomInfectionFrac"] = sum(SimOutput["RoomInfs"]) / NStot
    results[(irow_start+3),"CustIntroFrac"] = sum(SimOutput["CustomerIntroductions"]) / NStot
    results[(irow_start+3),"ExtIntroFrac"] = sum(SimOutput["ExternalIntroductions"]) / NStot
    results[(irow_start+3),"TotInfPackagesDelivered"] = sum(SimOutput["PackagesInfectiousOnDelivery"])
    results[(irow_start+3), "CustomersInfectedByPkg"] = sum(SimOutput["CustomersInfectedByPkgs"])
    results[(irow_start+3), "CustomersInfectedByDrivers"] = sum(SimOutput["CustomersInfectedByDrivers"])
    results[(irow_start+3), "IsolatorsFrac"] = sum(SimOutput["NewIsolators"])/(NStot)
    results[(irow_start+3), "SympIsolatorsFrac"] = sum(SimOutput["NewSympIsolators"])/(NStot)
    results[(irow_start+3), "FPIsolatorsFrac"] = sum(SimOutput["NewFalseIsolators"])/(NStot)
    results[(irow_start+3), "TPSympIsolatorsFrac"] = sum(SimOutput["NewTestSympIsolators"])/(NStot)
    results[(irow_start+3), "TPAsympIsolatorsFrac"] = sum(SimOutput["NewTestAsympIsolators"])/(NStot)
end


#Check this works
function run_many_sims(ParamsVec::Array{Dict{Any,Any},1}, Nrepeats::Int,
                OccPerDay::Array{Float64,1}, PkgParams::Array{Dict{Any,Any},1};
                NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay)),
                IsNetwork::Array{Bool,1} = ones(Bool,length(ParamsVec)),
                IsPairs::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                PairParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
                IsTesting::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                TestingParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
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
    PIsol = 0.2
    PFC = 1.0
    II = [1,2,3]
    tD = 0.05:0.05:1.0
    Phi = 0.05:0.05:1.0
    Nrepeats = 1000


    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for phi in Phi
        for ii in II
            for td in tD
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                    "p_contact"=>pc, "Pisol"=>PIsol,
                                    "InfInit"=>ii, "tD"=>td, "phi"=>phi,
                                    "p_friend_contact"=>PFC, "SimType"=>Outbreak_sim))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; filename="param_sweep.csv")
end

function run_presenteeism_param_sweep_outbreak_parcel()
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
    PIsol = 0.05:0.05:1.0
    PFC = [0.25,0.5,1.0]
    II = [1,2,3]
    tD = 0.5
    Phi = 0.1
    Nrepeats = 10000


    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for pfc in PFC
        for ii in II
            for pi in PISol
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                    "p_contact"=>pc, "Pisol"=>pi,
                                    "InfInit"=>ii, "tD"=>tD, "phi"=>Phi,
                                    "p_friend_contact"=>pfc, "SimType"=>Outbreak_sim))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; filename="param_sweep.csv")
end

function run_testing_sweep_outbreak_parcel()
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
                            "InfInit"=>1, "tD"=>1.0, "phi"=>1.0, "p_friend_contact"=>1.0,
                            "SimType"=>Outbreak_sim))
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
