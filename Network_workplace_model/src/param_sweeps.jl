include("dataframe_write.jl")

NweeksDefault = 52

#check and update
BasicParcelParams = Dict("ND"=>38, "NL"=>20, "NO"=>10, "NDteams"=>3, "NLteams"=>2, "NOteams"=>2,
                 "is_cohorts"=>true, "Pisol"=>0.5, "Psusc"=>1.0, "p_contact"=>(2.0/(38 + 20 + 10)),
                 "tD"=>0.05,"phi"=>0.1, "InfInit"=>0, "SimType"=>Outbreak_sim,
                 "TeamTimes"=>[0.25,1.0,1.0], "TeamsOutside"=>[true,true,false],
                 "TeamDistances"=>[2.0,2.0,2.0], "HouseShareFactor"=>0.5, "CarShareFactor"=>0.5)
ParcelOccPattern = 0.95 .* [0.90,1.0,1.0,0.99,0.91,0.55,0.0]
ParcelPkgPattern = [0.74,1.0,0.95,0.92,0.84,0.31,0.0]
ParcelPkgPattern = (6/7)*ParcelPkgPattern/mean(ParcelPkgPattern)
NPparcel = 3000


DefaultPkgParams = Dict("p_fomite_trans"=>0.0, "Dtime"=>8, "Ltime"=>4, "PkgHlife"=>3)

BasicBulkParams = Dict("ND"=>20, "NL"=>16, "NO"=>8, "NDteams"=>2, "NLteams"=>2, "NOteams"=>2,
                 "is_cohorts"=>true, "Pisol"=>0.5, "Psusc"=>1.0, "p_contact"=>(2.0/(20+16+8)),
                 "tD"=>0.05,"phi"=>0.1, "TeamTimes"=>[0.25,1.0,1.0], "TeamsOutside"=>[true,true,false],
                 "TeamDistances"=>[2.0,2.0,2.0], "HouseShareFactor"=>0.5, "CarShareFactor"=>0.5)
DefaultPairParams = Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                  "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                  "is_window_open"=>false, "pair_isolation"=>true)
BulkOccPattern = 0.95 .* [0.82, 0.98, 0.97, 0.99, 1.0, 0.84, 0.47]
BulkPkgPattern = [0.80, 0.94, 0.95, 0.94,  1.0, 0.81, 0.44]
BulkPkgPattern = BulkPkgPattern/mean(BulkPkgPattern)
NPbulk = 300

SpecDefault = 0.99   #Specificity

#Check this works
function run_many_sims(ParamsVec::Array{Dict{Any,Any},1}, Nrepeats::Int,
                OccPerDay::Array{Float64,1};
                PkgParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
                NPPerDay::Array{Int64,1} = zeros(Int64, length(OccPerDay)),
                PairParams::Array{Dict{Any,Any},1} = fill(Dict(),length(ParamsVec)),
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
            out = run_sim_delivery_wp(Params, OccPerDay, NPPerDay;
                PkgParams=PkgParams[i], PairParams=PairParams[i],
                TestParams=TestingParams[i], Incidence=Incidence,
                Prevalence=Prevalence)
            AllParams = merge(ParamsVec[i],PkgParams[i],PairParams[i],TestingParams[i])
            add_to_results_dataframe!(results, AllParams, out, index_start, n)
         end
     end

     CSV.write(filename, results)
     return results
end

function run_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))

    II = [1,2,3]
    Phis = 0.1:0.1:1.0
    NDteams = [2,5,10,3,3,3,3,3,3]
    NLteams = [2,2,2,2,4,8,2,2,2]
    NOteams = [2,2,2,2,2,2,1,3,6]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)

    for it in 1:length(NDteams)
        for ii in II
            for phi in Phis
                PP = copy(BasicParcelParams)
                PP["NDteams"] = NDteams[it]
                PP["NLteams"] = NLteams[it]
                PP["NOteams"] = NOteams[it]
                PP["phi"] = phi
                PP["InfInit"] = ii
                push!(ParamVec, PP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec, filename="param_sweep.csv")
    return df
end

function run_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    II = [1,2,3]
    Phis = 0.1:0.1:1.0
    iswo = [true, false]
    fp = [true, false]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for ii in II
        for phi in Phis
            for wo in iswo
                for fix in fp
                    PP = copy(BasicBulkParams)
                    PP["phi"] = phi
                    PP["InfInit"] = ii
                    PairPs = Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                  "fixed_driver_pairs"=>fix,"fixed_loader_pairs"=>fix,
                                  "is_window_open"=>wo)
                    push!(ParamVec, PP)
                    push!(PairParams, PairPs)
                end
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec,
                  filename="param_sweep_pairs.csv", IsPairs = ones(Bool,length(PairParams)),
                  PairParams = PairParams)
    return df
end

function run_presenteeism_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))
    #other params
    PIsol = 0.05:0.05:1.0
    CohortDistance = [1.0,1.5,2.0,3.0]
    II = [1,2,3]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    for cd in CohortDistance
        for ii in II
            for pi in PIsol
                PP = copy(BasicParcelParams)
                PP["Pisol"] = pi
                PP["InfInit"] = ii
                PP["TeamDistances"] = [cd,cd,cd]
                push!(ParamVec, PP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, filename="presenteeism_param_sweep.csv")
    return df
end

function run_presenteeism_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    PIsol = 0.05:0.05:1.0
    CohortDistance = [1.0,1.5,2.0,3.0]
    II = [1,2,3]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PairPs = copy(DefaultPairParams)
    for cd in CohortDistance
        for ii in II
            for pi in PIsol
                PP = copy(BasicBulkParams)
                PP["Pisol"] = pi
                PP["InfInit"] = ii
                PP["TeamDistances"] = [cd,cd,cd]
                push!(ParamVec, PP)
                push!(PairParams, PairPs)
            end
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec,
                  PairParams = PairParams, filename="presenteeism_param_sweep_pairs.csv")
    return df
end

function run_testing_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))

    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicParcelParams)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, PP)
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>SpecDefault,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))

            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec, output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, PP)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
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
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))

    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicBulkParams)
    PairPs = copy(DefaultPairParams)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)

                push!(ParamVec, PP)
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>SpecDefault,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PairParams, PairPs)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;  NPPerDay = NPvec,
             TestingParams=TestParamVec,  PairParams = PairParams, output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, PP)
    push!(PairParams, PairPs)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
        PairParams = PairParams, output = false)
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
    OccPattern = repeat(ParcelOccPattern,Nweeks)
    PkgPattern = repeat(ParcelPkgPattern,Nweeks)
    NPvec = Int64.(round.(NPparcel*PkgPattern))

    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicParcelParams)
    PP["SimType"] = Scenario_sim
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, PP)
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>SpecDefault,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))

            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec, output = false,
                  Incidence = Inc, Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, PP)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
                  output = false, Incidence = Inc, Prevalence = Prev)

    df2["new_comply_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill("No testing", nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    dfout = vcat(df,df2)

    CSV.write("testing_scenario_parcel_prev.csv", dfout)

    return df
end

function run_testing_sweep_fixedprev_scenario_pairs(Prev_val::Float64, Nrepeats::Int = 10000)
    NWeeks= 26
    OccPattern = repeat(BulkOccPattern,Nweek)
    PkgPattern = repeat(BulkPkgPattern,Nweeks)
    NPvec = Int64.(round.(NPbulk*PkgPattern))

    NCP = [0.0,0.5,1.0]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    Prev = Prev_val*ones(NWeeks*7)
    Inc = Prev_val*ones(NWeeks*7)/7

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicBulkParams)
    PP["SimType"] = Scenario_sim
    PairPs = copy(DefaultPairParams)
    for j in 1:3
        for tp in Tperiod
            for i in 1:length(TestType)
                push!(ParamVec, PP)
                push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
                      "protocol"=>TestType[i], "specificity"=>SpecDefault,
                      "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
                push!(PairParams, PairPs)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;  NPPerDay = NPvec,
             TestingParams=TestParamVec,  PairParams = PairParams, output = false,
             Incidence = Inc, Prevalence = Prev)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    push!(ParamVec, PP)
    push!(PairParams, PairPs)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
        PairParams = PairParams, output = false, Incidence = Inc, Prevalence = Prev)
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
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPh = [2000,3000,4000,5000]
    II = [1,2,3]
    PFT = [0.0001,0.001,0.01]

    df = DataFrame()
    PP = copy(BasicParcelParams)
    PkgP = copy(DefaultPkgParams)
    for (i,np) in enumerate(NPh)
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        NPvec = Int64.(round.(np*PkgPattern))
        for ii in II
            PP["InfInit"] = ii
            for pft in PFT
                PkgP["p_fomite_trans"] = pft
                push!(ParamVec, PP)
                push!(PkgParams, PkgP)
            end
        end
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                  NPPerDay = NPvec, TestingParams=TestParamVec, output = false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                  NPPerDay = NPvec, TestingParams=TestParamVec, output = false)
            df = vcat(df,dfh)
        end

    end
    CSV.write("fomite_param_sweep_parcel.csv", df)

    return df
end

function run_param_sweep_outbreak_fomite_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPh = [2000,3000,4000,5000]
    II = [1,2,3]
    PFT = [0.0001,0.001,0.01]

    df = DataFrame()
    PP = copy(BasicBulkParams)
    PairPs = copy(DefaultPairParams)
    PkgP = copy(DefaultPkgParams)
    for (i,np) in enumerate(NPh)
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        PairParams = Array{Dict{Any,Any},1}(undef,0)
        NPvec = Int64.(round.(np*PkgPattern))
        for ii in II
            PP["InfInit"] = ii
            for pft in PFT
                PkgP["p_fomite_trans"] = pft
                push!(ParamVec, PP)
                push!(PkgParams, PkgP)
                push!(PairParams, PairPs)
            end
        end
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                  PairParams = PairParams, NPPerDay = NPvec, TestingParams=TestParamVec,
                  output = false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                  PairParams = PairParams, NPPerDay = NPvec, TestingParams=TestParamVec,
                  output = false)
            df = vcat(df,dfh)
        end

    end
    CSV.write("fomite_param_sweep_pairs.csv", df)

    return df
end


# function run_testing_fixedprev_wpsize_sens_parcel(Prev_val::Float64, Nrepeats::Int = 10000)
#     NWeeks= 26
#     NPh = [250,500,1000,1500,6000,9000,18000]
#     OccPattern = repeat(ParcelOccPattern,NWeeks)
#     PkgPattern = repeat(ParcelPkgPattern,NWeeks)

#     Prev = Prev_val*ones(NWeeks*7)
#     Inc = Prev_val*ones(NWeeks*7)/7
#     RandContacts = 2.0

#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     Tperiod = 2:2:14
#     NCP = 1.0
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)

#     for (i, NP) in enumerate(NPh)
#         NPvec = Int64.(round.(NP*PkgPattern))
#         NDh = Int64(round(NP/80))
#         NLh = Int64(round(NP/150))
#         NOh = Int64(round(NP/300))
#         pc = RandContacts/(NDh + NLh + NOh)
#         for tp in Tperiod
#             for k in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
#                       "protocol"=>TestType[k], "specificity"=>0.999,
#                       "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#         if i == 1
#             df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)
#         else
#             dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                       NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                       TestingParams=TestParamVec, output = false, Incidence = Inc,
#                       Prevalence = Prev)
#             df = vcat(df,dfh)
#         end

#         #run baseline case
#         ParamVec = Array{Dict{Any,Any},1}(undef,0)
#         PkgParams = Array{Dict{Any,Any},1}(undef,0)
#         PairParams = Array{Dict{Any,Any},1}(undef,0)
#         push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                             "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                             "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                             "SimType"=>Scenario_sim))
#         push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                             "Ltime"=>1/6, "PkgHlife"=>0.5))

#         df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                       output = false, Incidence = Inc, Prevalence = Prev)

#         df2["new_comply_prob"] = zeros(nrow(df2))
#         df2["tperiod"] = zeros(nrow(df2))
#         df2["protocol"] = fill("No testing", nrow(df2))
#         df2["specificity"] = 0.999 * ones(nrow(df2))
#         df2["delay"] = zeros(nrow(df2))
#         df2["test_pause"] = zeros(nrow(df2))
#         df = vcat(df,df2)
#     end

#     fname = string("testing_scenario_parcel_wpsize_prev",string(100*Prev_val),".csv")

#     CSV.write(fname, df)

# end

# function run_testing_fixedprev_wpsize_sens_pairs(Prev_val::Float64, Nrepeats::Int = 10000)
#     NWeeks= 26


#     Prev = Prev_val*ones(NWeeks*7)
#     Inc = Prev_val*ones(NWeeks*7)/7
#     RandContacts = 2.0

#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     Tperiod = 2:2:14
#     NCP = 1.0
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     df = DataFrame()
#     for (i, NP) in enumerate(NPh)
#         NPvec = Int64.(round.(NP*OccPattern))
#         NDh = Int64(2*round(NP/30))
#         NLh = Int64(2*round(NP/40))
#         NOh = Int64(round(NP/40))
#         pc = RandContacts/(NDh + NLh + NOh)
#         for tp in Tperiod
#             for k in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
#                       "protocol"=>TestType[k], "specificity"=>0.999,
#                       "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#         if i == 1
#             df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)
#         else
#             dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                       NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                       TestingParams=TestParamVec, output = false, Incidence = Inc,
#                       Prevalence = Prev)
#             df = vcat(df,dfh)
#         end

#         #run baseline case
#         ParamVec = Array{Dict{Any,Any},1}(undef,0)
#         PkgParams = Array{Dict{Any,Any},1}(undef,0)
#         PairParams = Array{Dict{Any,Any},1}(undef,0)
#         push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                             "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                             "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                             "SimType"=>Scenario_sim))
#         push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                             "Ltime"=>1/6, "PkgHlife"=>0.5))

#         df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                       output = false, Incidence = Inc, Prevalence = Prev)

#         df2["new_comply_prob"] = zeros(nrow(df2))
#         df2["tperiod"] = zeros(nrow(df2))
#         df2["protocol"] = fill("No testing", nrow(df2))
#         df2["specificity"] = 0.999 * ones(nrow(df2))
#         df2["delay"] = zeros(nrow(df2))
#         df2["test_pause"] = zeros(nrow(df2))
#         df = vcat(df,df2)
#     end

#     fname = string("testing_scenario_pairs_wpsize_prev",string(100*Prev_val),".csv")

#     CSV.write(fname, df)

# end

# function run_testing_sweep_variableprev_scenario_parcel(Prev::Array{Float64,1},
#         Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = 3000
#     OccPattern = repeat(ParcelOccPattern,NWeeks)
#     PkgPattern = repeat(ParcelPkgPattern,NWeeks)
#     OccPattern = OccPattern[1:length(Prev)]
#     PkgPattern = PkgPattern[1:length(Prev)]
#     NPvec = Int64.(round.(NPh*PkgPattern))
#     NDh = Int64(round(NPh/80))
#     NLh = Int64(round(NPh/150))
#     NOh = Int64(round(NPh/300))
#     ContactsPerDay = 2
#     pc = ContactsPerDay/(NDh+NLh+NOh)
#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     NCP = [0.0,0.5,1.0]
#     Tperiod = 2:2:14

#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>0.999,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#     push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                         "Ltime"=>1/6, "PkgHlife"=>0.5))

#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                   output = false, Incidence = Inc, Prevalence = Prev)

#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     dfout = vcat(df,df2)

#     fname = string("testing_varscenario_parcel.csv")

#     CSV.write(fname, dfout)

#     return df
# end

# function run_testing_sweep_variableprev_scenario_pairs(Prev::Array{Float64,1},
#         Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = 300


#     OccPattern = OccPattern[1:length(Prev)]


#     NDh = Int64(2*round(NPh/30))
#     NLh = Int64(2*round(NPh/40))
#     NOh = Int64(round(NPh/40))
#     ContactsPerDay = 2
#     pc = ContactsPerDay/(NDh+NLh+NOh)
#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     NCP = [0.0,0.5,1.0]
#     Tperiod = 2:2:14

#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)

#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>0.999,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#                 push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
#                                         "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
#                                          "is_window_open"=>false, "pair_isolation"=>true))
#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;  NPPerDay = NPvec,
#                   IsTesting=ones(Bool,length(ParamVec)), TestingParams=TestParamVec,
#                   IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
#                   output = false, Incidence = Inc, Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#     push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                         "Ltime"=>1/6, "PkgHlife"=>0.5))
#     push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
#                             "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
#                             "is_window_open"=>false, "pair_isolation"=>true))

#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                   IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
#                   output = false, Incidence = Inc, Prevalence = Prev)
#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     df = vcat(df,df2)

#     fname = string("testing_varscenario_pairs.csv")
#     CSV.write(fname, df)
#     return df
# end

# function run_testing_variableprev_wpsize_sens_parcel(Prev::Array{Float64,1},
#         Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = [250,500,1000,1500,6000,9000,18000]
#     OccPattern = repeat(ParcelOccPattern,NWeeks)
#     PkgPattern = repeat(ParcelPkgPattern,NWeeks)
#     OccPattern = OccPattern[1:length(Prev)]
#     PkgPattern = PkgPattern[1:length(Prev)]

#     RandContacts = 2.0

#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     Tperiod = 2:2:14
#     NCP = 1.0
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     df = DataFrame()
#     for (i, NP) in enumerate(NPh)
#         NPvec = Int64.(round.(NP*PkgPattern))
#         NDh = Int64(round(NP/80))
#         NLh = Int64(round(NP/150))
#         NOh = Int64(round(NP/300))
#         pc = RandContacts/(NDh + NLh + NOh)
#         for tp in Tperiod
#             for k in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
#                       "protocol"=>TestType[k], "specificity"=>0.999,
#                       "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#         if i == 1
#             df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)
#         else
#             dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                       NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                       TestingParams=TestParamVec, output = false, Incidence = Inc,
#                       Prevalence = Prev)
#             df = vcat(df,dfh)
#         end

#         #run baseline case
#         ParamVec = Array{Dict{Any,Any},1}(undef,0)
#         PkgParams = Array{Dict{Any,Any},1}(undef,0)
#         PairParams = Array{Dict{Any,Any},1}(undef,0)
#         push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                             "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                             "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                             "SimType"=>Scenario_sim))
#         push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                             "Ltime"=>1/6, "PkgHlife"=>0.5))

#         df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                       output = false, Incidence = Inc, Prevalence = Prev)

#         df2["new_comply_prob"] = zeros(nrow(df2))
#         df2["tperiod"] = zeros(nrow(df2))
#         df2["protocol"] = fill("No testing", nrow(df2))
#         df2["specificity"] = 0.999 * ones(nrow(df2))
#         df2["delay"] = zeros(nrow(df2))
#         df2["test_pause"] = zeros(nrow(df2))
#         df = vcat(df,df2)
#     end

#     fname = string("testing_varscenario_parcel_wpsize.csv")

#     CSV.write(fname, df)
# end

# function run_testing_variableprev_wpsize_sens_pairs(Prev::Array{Float64,1},
#         Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = [25,50,100,150,600,900,1800]

#     OccPattern = repeat(PairsOccPattern,NWeeks)
#     PkgPattern = repeat(PairsPkgPattern,NWeeks)
#     OccPattern = OccPattern[1:length(Prev)]
#     PkgPattern = PkgPattern[1:length(Prev)]

#     RandContacts = 2.0

#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     Tperiod = 2:2:14
#     NCP = 1.0
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     df = DataFrame()
#     for (i, NP) in enumerate(NPh)
#         NPvec = Int64.(round.(NP*PkgPattern))
#         NDh = Int64(2*round(NP/30))
#         NLh = Int64(2*round(NP/40))
#         NOh = Int64(round(NP/40))
#         pc = RandContacts/(NDh + NLh + NOh)
#         for tp in Tperiod
#             for k in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP, "tperiod"=>tp,
#                       "protocol"=>TestType[k], "specificity"=>0.999,
#                       "delay"=>Delay[k], "test_pause"=>Test_pause[k]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#         if i == 1
#             df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)
#         else
#             dfh = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                       NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                       TestingParams=TestParamVec, output = false, Incidence = Inc,
#                       Prevalence = Prev)
#             df = vcat(df,dfh)
#         end

#         #run baseline case
#         ParamVec = Array{Dict{Any,Any},1}(undef,0)
#         PkgParams = Array{Dict{Any,Any},1}(undef,0)
#         PairParams = Array{Dict{Any,Any},1}(undef,0)
#         push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                             "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                             "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                             "SimType"=>Scenario_sim))
#         push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                             "Ltime"=>1/6, "PkgHlife"=>0.5))

#         df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                       output = false, Incidence = Inc, Prevalence = Prev)

#         df2["new_comply_prob"] = zeros(nrow(df2))
#         df2["tperiod"] = zeros(nrow(df2))
#         df2["protocol"] = fill("No testing", nrow(df2))
#         df2["specificity"] = 0.999 * ones(nrow(df2))
#         df2["delay"] = zeros(nrow(df2))
#         df2["test_pause"] = zeros(nrow(df2))
#         df = vcat(df,df2)
#     end

#     fname = string("testing_varscenario_pairs_wpsize.csv")

#     CSV.write(fname, df)
# end

# function run_all_interventions_variableprev_scenario_parcel(Prev::Array{Float64,1},
#                                    Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = 3000
#     OccPattern = repeat(ParcelOccPattern,NWeeks)
#     PkgPattern = repeat(ParcelPkgPattern,NWeeks)
#     OccPattern = OccPattern[1:length(Prev)]
#     PkgPattern = PkgPattern[1:length(Prev)]
#     NPvec = Int64.(round.(NPh*PkgPattern))
#     NDh = Int64(round(NPh/80))
#     NLh = Int64(round(NPh/150))
#     NOh = Int64(round(NPh/300))
#     ContactsPerDay = 2
#     pc = ContactsPerDay/(NDh+NLh+NOh)
#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = [0.0,1.0]
#     PFC = [0.5,1.0]
#     Mask_wearing = [true,false]
#     Fixed_pairs = [true,false]
#     Distancing = [true,false]
#     NCP = 0.0
#     Tperiod = 3.5

#     Delay = 0
#     TestType = LFD_mass_protocol
#     Test_pause = 21

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)
#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol,
#                         "InfInit"=>0, "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>0.999,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;
#                   NPPerDay = NPvec, IsTesting=ones(Bool,length(ParamVec)),
#                   TestingParams=TestParamVec, output = false, Incidence = Inc,
#                   Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#     push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                         "Ltime"=>1/6, "PkgHlife"=>0.5))

#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                   output = false, Incidence = Inc, Prevalence = Prev)

#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     dfout = vcat(df,df2)

#     fname = string("testing_varscenario_parcel.csv")

#     CSV.write(fname, dfout)

#     return df
# end

# function run_testing_sweep_variableprev_scenario_pairs(Prev::Array{Float64,1},
#         Inc::Array{Float64,1}, Nrepeats::Int = 10000)
#     NWeeks = Int64(ceil(length(Prev)/7))
#     NPh = 300
#     OccPattern = repeat(PairsOccPattern,NWeeks)
#     PkgPattern = repeat(PairsPkgPattern,NWeeks)
#     OccPattern = OccPattern[1:length(Prev)]
#     PkgPattern = PkgPattern[1:length(Prev)]
#     NPvec = Int64.(round.(NPh*PkgPattern))


#     NDh = Int64(2*round(NPh/30))
#     NLh = Int64(2*round(NPh/40))
#     NOh = Int64(round(NPh/40))
#     ContactsPerDay = 2
#     pc = ContactsPerDay/(NDh+NLh+NOh)
#     #sweep these two in unison
#     phi = 0.05
#     tD = 0.25
#     PIsol = 0.5
#     NCP = [0.0,0.5,1.0]
#     Tperiod = 2:2:14

#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)

#                 push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>0.999,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
#                 push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                              "Ltime"=>1/6, "PkgHlife"=>0.5))
#                 push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
#                                         "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
#                                          "is_window_open"=>false, "pair_isolation"=>true))
#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern, PkgParams;  NPPerDay = NPvec,
#                   IsTesting=ones(Bool,length(ParamVec)), TestingParams=TestParamVec,
#                   IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
#                   output = false, Incidence = Inc, Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgParams = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     push!(ParamVec, Dict("ND"=>NDh, "NL"=>NLh, "NO"=>NOh,
#                         "p_contact"=>pc, "Pisol"=>PIsol, "InfInit"=>0,
#                         "tD"=>tD, "phi"=>phi, "p_friend_contact"=>1.0,
#                         "SimType"=>Scenario_sim))
#     push!(PkgParams, Dict("p_fomite_contr"=>0.0, "p_fomite_trans"=>0.0, "Dtime"=>1/6,
#                         "Ltime"=>1/6, "PkgHlife"=>0.5))
#     push!(PairParams, Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
#                             "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
#                             "is_window_open"=>false, "pair_isolation"=>true))

#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern, PkgParams;  NPPerDay = NPvec,
#                   IsPairs = ones(Bool,length(PairParams)), PairParams = PairParams,
#                   output = false, Incidence = Inc, Prevalence = Prev)
#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     df = vcat(df,df2)

#     fname = string("testing_varscenario_pairs.csv")
#     CSV.write(fname, df)
#     return df
# end
