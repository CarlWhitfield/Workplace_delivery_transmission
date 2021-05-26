include("network_transmission_workplace.jl")
include("dataframe_write.jl")

#Check this works
function run_many_sims(ParamsVec::Array{Dict{Any,Any},1}, Nrepeats::Int,
                OccPerDay::Array{Float64,1}, PkgParams::Array{Dict{Any,Any},1};
                NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay)),
                IsNetwork::Array{Bool,1} = ones(Bool,length(ParamsVec)),
                IsPairs::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                PairParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
                IsTesting::Array{Bool,1} = zeros(Bool,length(ParamsVec)),
                TestingParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
                filename="output.csv", output::Bool = true,
                Incidence::Array{Float64,1} = zeros(Float64, length(OccPerDay)),
                Prevalence::Array{Float64,1} = zeros(Float64, length(OccPerDay)))

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
                   testing=IsTesting[i], TestParams=TestingParams[i], Incidence = Incidence,
                   Prevalence = Prevalence)
             AllParams = merge(ParamsVec[i],PkgParams[i],PairParams[i],TestingParams[i])
             add_to_results_dataframe!(results, AllParams, out, index_start, n)
         end
     end

     CSV.write(filename, results)
     return results
end

function run_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
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
    PIsol = 0.5
    PFC = 1.0
    II = [1,2,3]
    tD = 0.1:0.1:1.0
    Phi = [0.05,0.1,0.25,0.5,1.0]


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

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; NPPerDay = NPvec, filename="param_sweep.csv")
    return df
end

function run_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #other params
    PIsol = 0.5
    PFC = 1.0
    II = [1,2,3]
    tD = 0.1:0.1:1.0
    Phi = 0.05
    iswo = [true, false]
    fp = [true, false]
    pair_isol = [true, false]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for ii in II
        for td in tD
            for wo in iswo
                for fix in fp
                    for pih in pair_isol
                        push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                            "p_contact"=>pc, "Pisol"=>PIsol,
                                            "InfInit"=>ii, "tD"=>td, "phi"=>Phi,
                                            "p_friend_contact"=>PFC, "SimType"=>Outbreak_sim))
                        push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                               "fixed_driver_pairs"=>fix, "fixed_loader_pairs"=>fix,
                                               "is_window_open"=>wo, "pair_isolation"=>pih))
                        push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                            "Ltime"=>1/6, "PkgHlife"=>0.5))
                    end
                end
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; NPPerDay = NPvec,
                  filename="param_sweep_pairs.csv", IsPairs = ones(Bool,length(PairParams)),
                  PairParams = PairParams)
    return df
end

function run_presenteeism_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
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
    tD = 0.25
    Phi = 0.05

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for pfc in PFC
        for ii in II
            for pi in PIsol
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                    "p_contact"=>pc, "Pisol"=>pi,
                                    "InfInit"=>ii, "tD"=>tD, "phi"=>Phi,
                                    "p_friend_contact"=>pfc, "SimType"=>Outbreak_sim))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, filename="presenteeism_param_sweep.csv")
    return df
end

function run_presenteeism_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #other params
    PIsol = 0.05:0.05:1.0
    PFC = [0.25,0.5,1.0]
    II = [1,2,3]
    tD = 0.25
    Phi = 0.05

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for pfc in PFC
        for ii in II
            for pi in PIsol
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                    "p_contact"=>pc, "Pisol"=>pi,
                                    "InfInit"=>ii, "tD"=>tD, "phi"=>Phi,
                                    "p_friend_contact"=>pfc, "SimType"=>Outbreak_sim))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                                 "Ltime"=>1/6, "PkgHlife"=>0.5))
                push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                       "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                                        "is_window_open"=>false, "pair_isolation"=>true))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; NPPerDay = NPvec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  filename="presenteeism_param_sweep_pairs.csv")
    return df
end

function run_testing_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 3000
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(round(NPh/80))
    NLh = Int64(round(NPh/150))
    NOh = Int64(round(NPh/300))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol,
                        "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Outbreak_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Outbreak_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  output = false)

    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    dfout = vcat(df,df2)

    CSV.write("testing_sweep.csv", dfout)

    return df
end

function run_testing_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)

                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Outbreak_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
                push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                        "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                                         "is_window_open"=>false, "pair_isolation"=>true))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsTesting=ones(Bool,length(ParamVec)), TestingParams=TestParamVec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Outbreak_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))
    push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                            "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                            "is_window_open"=>false, "pair_isolation"=>true))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false)
    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    df = vcat(df,df2)

    CSV.write("testing_sweep_pairs.csv", df)
    return df
end


function run_testing_sweep_fixedprev_scenario_parcel(Prev_val::Float64, Nrepeats::Int = 10000)
    NWeeks= 26
    NPh = 3000
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(round(NPh/80))
    NLh = Int64(round(NPh/150))
    NOh = Int64(round(NPh/300))
    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol,
                        "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, output = false, Incidence = Inc,
                  Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  output = false, Incidence = Inc, Prevalence = Prev)

    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    dfout = vcat(df,df2)

    fname = string("testing_scenario_parcel_prev",string(100*Prev_val),".csv")

    CSV.write(fname, dfout)

    return df
end

function run_testing_sweep_fixedprev_scenario_pairs(Prev_val::Float64, Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)

                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
                push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                        "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                                         "is_window_open"=>false, "pair_isolation"=>true))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsTesting=ones(Bool,length(ParamVec)), TestingParams=TestParamVec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false, Incidence = Inc, Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))
    push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                            "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                            "is_window_open"=>false, "pair_isolation"=>true))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false, Incidence = Inc, Prevalence = Prev)
    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    df = vcat(df,df2)

    fname = string("testing_scenario_pairs_prev",string(100*Prev_val),".csv")
    CSV.write(fname, df)
    return df
end

function run_param_sweep_outbreak_fomite_parcel(Nrepeats::Int = 10000)
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
    PIsol = 0.5
    FCR = 1.0
    II = [1,2,3]
    tD = 0.25
    Phi = [0.05, 1.0]
    PFT = [0.001,0.01,0.1]
    PFC = [0.001,0.01,0.1]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for phi in Phi
        for ii in II
            for pft in PFT
                for pfc in PFC
                    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                    "p_contact"=>pc, "Pisol"=>PIsol,
                                    "InfInit"=>ii, "tD"=>tD, "phi"=>phi,
                                    "p_friend_contact"=>FCR, "SimType"=>Outbreak_sim))
                    push!(PkgParams, Dict("p_fomite_contr"=>pfc, "p_fomite_trans"=>pft, "Dtime"=>4,
                                            "Ltime"=>4, "PkgHlife"=>3))
                end
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; NPPerDay = NPvec, filename="param_sweep_fomite.csv")
    return df
end

function run_param_sweep_outbreak_fomite_pairs(Nrepeats::Int = 10000)
    NWeeks= 52
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #other params
    PIsol = 0.5
    FCR = 1.0
    II = [1,2,3]
    Phi = 0.05
    iswo = [true, false]
    fp = [true, false]
    pair_isol = [true, false]
    PFT = [0.001,0.01,0.1]
    PFC = [0.001,0.01,0.1]
    tD = 0.25
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for ii in II
        for pft in PFT
            for pfc in PFC
                for wo in iswo
                    for fix in fp
                        for pih in pair_isol
                            push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                                                "p_contact"=>pc, "Pisol"=>PIsol,
                                                "InfInit"=>ii, "tD"=>tD, "phi"=>Phi,
                                                "p_friend_contact"=>FCR, "SimType"=>Outbreak_sim))
                            push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                                   "fixed_driver_pairs"=>fix, "fixed_loader_pairs"=>fix,
                                                   "is_window_open"=>wo, "pair_isolation"=>pih))
                            push!(PkgParams, Dict("p_fomite_contr"=>pfc, "p_fomite_trans"=>pft, "Dtime"=>4,
                                                "Ltime"=>4, "PkgHlife"=>3))
                        end
                    end
                end
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams; NPPerDay = NPvec,
                  filename="param_sweep_pairs_fomite.csv", IsPairs = ones(Bool,length(PairParams)),
                  PairParams = PairParams)
    return df
end

function run_testing_fixedprev_wpsize_sens_parcel(Prev_val::Float64, Nrepeats::Int = 10000)
    NWeeks= 26
    NPh = [250,500,1000,1500,6000,9000,18000]
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7
    RandContacts = 2.0

    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    Tperiod = 2:2:14
    NCP = 1.0
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    df = DataFrame()
    for (i, NP) in enumerate(NPh)
        NPvec = Int64.(round.(NP*OccPattern))
        NDh = Int64(round(NP/80))
        NLh = Int64(round(NP/150))
        NOh = Int64(round(NP/300))
        pc = RandContacts/(NDh + NLh + NOh)
        for tp in Tperiod
            for k in 1:length(TestType)
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol,
                        "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
                      "protocol"=>TestType[k], "specificity"=>0.999,
                      "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, output = false, Incidence = Inc,
                  Prevalence = Prev)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                      NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                      TestingParams=TestParamVec, output = false, Incidence = Inc,
                      Prevalence = Prev)
            df = vcat(df,dfh)
        end

        #run baseline case
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        PairParams = Array{Dict{Any,Any},1}(undef,0)
        push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                            "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                            "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                            "SimType"=>Scenario_sim))
        push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                            "Ltime"=>1/6, "PkgHlife"=>0.5))

        df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                      output = false, Incidence = Inc, Prevalence = Prev)

        df2["new_comply_prob"] = zeros(nrow(df2))
        df2["tperiod"] = zeros(nrow(df2))
        df2["protocol"] = fill("No testing", nrow(df2))
        df2["specificity"] = 0.999 * ones(nrow(df2))
        df2["delay"] = zeros(nrow(df2))
        df2["test_pause"] = zeros(nrow(df2))
        df = vcat(df,df2)
    end

    fname = string("testing_scenario_parcel_wpsize_prev",string(100*Prev_val),".csv")

    CSV.write(fname, df)

end

function run_testing_fixedprev_wpsize_sens_pairs(Prev_val::Float64, Nrepeats::Int = 10000)
    NWeeks= 26
    NPh = [25,50,100,150,600,900,1800]
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7
    RandContacts = 2.0

    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    Tperiod = 2:2:14
    NCP = 1.0
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    df = DataFrame()
    for (i, NP) in enumerate(NPh)
        NPvec = Int64.(round.(NP*OccPattern))
        NDh = Int64(2*round(NP/30))
        NLh = Int64(2*round(NP/40))
        NOh = Int64(round(NP/40))
        pc = RandContacts/(NDh + NLh + NOh)
        for tp in Tperiod
            for k in 1:length(TestType)
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol,
                        "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
                      "protocol"=>TestType[k], "specificity"=>0.999,
                      "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, output = false, Incidence = Inc,
                  Prevalence = Prev)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                      NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                      TestingParams=TestParamVec, output = false, Incidence = Inc,
                      Prevalence = Prev)
            df = vcat(df,dfh)
        end

        #run baseline case
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        PairParams = Array{Dict{Any,Any},1}(undef,0)
        push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                            "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                            "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                            "SimType"=>Scenario_sim))
        push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                            "Ltime"=>1/6, "PkgHlife"=>0.5))

        df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                      output = false, Incidence = Inc, Prevalence = Prev)

        df2["new_comply_prob"] = zeros(nrow(df2))
        df2["tperiod"] = zeros(nrow(df2))
        df2["protocol"] = fill("No testing", nrow(df2))
        df2["specificity"] = 0.999 * ones(nrow(df2))
        df2["delay"] = zeros(nrow(df2))
        df2["test_pause"] = zeros(nrow(df2))
        df = vcat(df,df2)
    end

    fname = string("testing_scenario_pairs_wpsize_prev",string(100*Prev_val),".csv")

    CSV.write(fname, df)

end

function run_testing_sweep_variableprev_scenario_parcel(Prev::Array{Float64,1},
        Inc::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    NPh = 3000
    OccPattern = repeat([0.87,1.0,1.0,0.98,0.91,0.55,0],NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(round(NPh/80))
    NLh = Int64(round(NPh/150))
    NOh = Int64(round(NPh/300))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol,
                        "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
                  NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
                  TestingParams=TestParamVec, output = false, Incidence = Inc,
                  Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  output = false, Incidence = Inc, Prevalence = Prev)

    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    dfout = vcat(df,df2)

    fname = string("testing_varscenario_parcel_prev",string(100*Prev_val),".csv")

    CSV.write(fname, dfout)

    return df
end

function run_testing_sweep_variableprev_scenario_pairs(Prev::Array{Float64,1},
        Inc::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    NPh = 300
    OccPattern = repeat([0.80,0.94,0.95,0.94,1.0,0.96,0.53],NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    NPvec = Int64.(round.(NPh*OccPattern))
    NDh = Int64(2*round(NPh/30))
    NLh = Int64(2*round(NPh/40))
    NOh = Int64(round(NPh/40))
    ContactsPerDay = 2
    pc = ContactsPerDay/(NDh+NLh+NOh)
    #sweep these two in unison
    phi = 0.05
    tD = 0.25
    PIsol = 0.5
    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14

    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)

                push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>0.999,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                             "Ltime"=>1/6, "PkgHlife"=>0.5))
                push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                        "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                                         "is_window_open"=>false, "pair_isolation"=>true))
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsTesting=ones(Bool,length(ParamVec)), TestingParams=TestParamVec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false, Incidence = Inc, Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
                        "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
                        "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
                        "SimType"=>Scenario_sim))
    push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
                        "Ltime"=>1/6, "PkgHlife"=>0.5))
    push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                            "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                            "is_window_open"=>false, "pair_isolation"=>true))

    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
                  IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
                  output = false, Incidence = Inc, Prevalence = Prev)
    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    df = vcat(df,df2)

    fname = string("testing_varscenario_pairs_prev",string(100*Prev_val),".csv")
    CSV.write(fname, df)
    return df
end