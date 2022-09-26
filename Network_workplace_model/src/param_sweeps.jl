using Distributed

include("dataframe_write.jl")

NweeksDefault = 52

#check and update
RandomAbsenceRate = 0.03
ParcelOccPattern = [0.90,1.0,1.0,0.99,0.91,0.55,0.0]
ParcelPkgPattern = [0.89,1.20,1.14,1.10,1.0,0.37,0.0]
NPparcel = 3000
NDparcel_def = Int64(ceil(mean(NPparcel*ParcelPkgPattern[1:6] ./
                     (85.0*(1.0-RandomAbsenceRate)*ParcelOccPattern[1:6]))))
NLparcel_def = Int64(ceil(0.5*NDparcel_def))
NOparcel_def = Int64(ceil(0.3*NDparcel_def))
NStaffparcel_def = NDparcel_def + NLparcel_def + NOparcel_def
BasicParcelParams = Dict("ND"=>NDparcel_def, "NL"=>NLparcel_def, "NO"=>NOparcel_def, "NDteams"=>3, "NLteams"=>2,
        "NOteams"=>1, "is_cohorts"=>true, "Pisol"=>0.9, "Psusc"=>1.0,
        "p_contact"=>(2.0/(NStaffparcel_def)), "tD"=>0.05,"phi"=>1.0, "InfInit"=>0,
        "SimType"=>Outbreak_sim, "TeamTimes"=>[0.25,1.0,1.0],
        "TeamsOutside"=>[false,false,false], "TeamDistances"=>[1.0,1.0,1.0],
         "HouseShareFactor"=>0.05, "CarShareFactor"=>0.05, "BreakContactProb"=>0.25,
         "CohortChangeRate"=>(1.0/(NStaffparcel_def)), "AbsenceRate"=>RandomAbsenceRate,
         "CustomerOutdoorFrac"=>1.0)

SpecDefault = 0.999 #Specificity
BasicTestingParams = Dict("is_testing"=>true, "test_miss_prob"=>0.4,
             "testing_enforced"=>false, "tperiod"=>3, "protocol"=>LFD_mass_protocol,
             "specificity"=>SpecDefault, "delay"=>0.0,
             "test_pause"=>21.0)
BasicPkgParams = Dict("p_fomite_trans"=>0.0001, "Dtime"=>8, "Ltime"=>4, "PkgHlife"=>3)

BulkOccPattern = [0.82, 0.98, 0.97, 0.99, 1.0, 0.84, 0.47]
BulkPkgPattern = [0.80, 0.94, 0.95, 0.94,  1.0, 0.81, 0.44]
NPbulk = 210
NDbulk_def = Int64(ceil(mean(NPbulk*BulkPkgPattern ./
                             (15.0*(1.0-RandomAbsenceRate)*BulkOccPattern))))
NLbulk_def = Int64(ceil(0.8*NDbulk_def))
NObulk_def = Int64(ceil(0.4*NDbulk_def))
NStaffbulk_def = NDbulk_def + NLbulk_def + NObulk_def
BasicBulkParams = Dict("ND"=>NDbulk_def, "NL"=>NLbulk_def, "NO"=>NObulk_def,
                       "NDteams"=>2, "NLteams"=>2, "NOteams"=>1,
                       "is_cohorts"=>true, "Pisol"=>0.9, "Psusc"=>1.0,
                       "p_contact"=>(2.0/(NStaffbulk_def)),
                       "tD"=>0.1, "phi"=>1.0, "InfInit"=>0, "SimType"=>Outbreak_sim,
                       "TeamTimes"=>[0.25,1.0,1.0], "TeamsOutside"=>[false,false,false],
                       "TeamDistances"=>[1.0,1.0,1.0], "HouseShareFactor"=>0.05,
                       "CarShareFactor"=>0.05, "BreakContactProb"=>0.25,
                       "CohortChangeRate"=>(1.0/(NStaffbulk_def)),
                       "AbsenceRate"=>RandomAbsenceRate, "CustomerOutdoorFrac"=>0.5)
BasicPairParams = Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                  "fixed_driver_pairs"=>true, "fixed_loader_pairs"=>true,
                  "is_window_open"=>false, "PairIsolation"=>true)


##RNA viral load model
VL_model = ke_model_no     #default option
#VL_model = kissler_model_no
#VL_model = HCS_model_no

##Infectiousness model
Inf_model = ke_inf_model_no
#Inf_model = flat_inf_model_no
#Inf_model = linear_inf_model_no

##LFD Test positivity model
LFD_model = porton_down_p3b
# LFD_model = porton_down
# LFD_model = SC_data_2021
# LFD_model = SC_data_Oct

##Model parameters and extra options
PCR_sens_max = 0.95            #maximum PCR sensitivity (0.95 - usually assumed, 0.83 - Ferretti et al 2021)

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
     IntArray, FloatArray, BoolArray, IntColMap, FloatColMap, BoolColMap =
                                    init_results_dataframe(Nrows, TestParams)
     print(Nrows,'\n')
     i_step = 4*Nrepeats
     for (i, Params) in enumerate(ParamsVec)
         i_ind_start = (i-1)*i_step + 1
         print(i,'/',length(ParamsVec),'\n')
         @distributed for n in 1:Nrepeats
            index_start = i_ind_start + (n-1)*4
            Ph = deepcopy(Params)
            PkgPh = deepcopy(PkgParams[i])
            PairPh = deepcopy(PairParams[i])
            TPh = deepcopy(TestingParams[i])
            out = run_sim_delivery_wp(Ph, deepcopy(OccPerDay), deepcopy(NPPerDay);
                PkgParams=PkgPh, PairParams=PairPh, TestParams=TPh,
                Incidence=deepcopy(Incidence), Prevalence=deepcopy(Prevalence))
            AllParams = merge(Ph,PkgPh,PairPh,TPh)
            add_to_results_dataframe!(IntArray, FloatArray, BoolArray, IntColMap,
                                      FloatColMap, BoolColMap, AllParams, out, index_start, n)
         end
         print(i_ind_start,'\n')
     end
    
     results = create_dataframe_from_arrays(IntArray, FloatArray, BoolArray, IntColMap,
                                            FloatColMap, BoolColMap)

     CSV.write(filename, results)
     return results
end

function run_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))

    II = [1,2,3]
    CR = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0] * BasicParcelParams["CohortChangeRate"]
    NDteams = [3,3,3,3,3,3,4,6,8]
    NLteams = [2,2,2,3,4,5,2,2,2]
    NOteams = [1,2,4,1,1,1,1,1,1]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for it in 1:length(NDteams)
        for ii in II
            for cr in CR
                PP = copy(BasicParcelParams)
                PP["NDteams"] = NDteams[it]
                PP["NLteams"] = NLteams[it]
                PP["NOteams"] = NOteams[it]
                PP["CohortChangeRate"] = cr
                PP["InfInit"] = ii
                push!(ParamVec, PP)
                push!(PkgVec, PkgP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec,
                       PkgParams = PkgVec, filename="param_sweep.csv")
    return df
end

function run_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    II = [1,2,3]
    NDteams = [3,3,3,3,3,3,4,6,8]
    NLteams = [2,2,2,3,4,5,2,2,2]
    NOteams = [1,2,4,1,1,1,1,1,1]
    iswo = [true, false]
    fp = [true, false]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for ii in II
        for it in 1:length(NDteams)
            for wo in iswo
                for fix in fp
                    PP = copy(BasicBulkParams)
                    PP["NDteams"] = NDteams[it]
                    PP["NLteams"] = NLteams[it]
                    PP["NOteams"] = NOteams[it]
                    PP["InfInit"] = ii
                    PairPs = Dict("is_driver_pairs"=>true, "is_loader_pairs"=>true,
                                  "fixed_driver_pairs"=>fix,"fixed_loader_pairs"=>fix,
                                  "is_window_open"=>wo, "PairIsolation"=>fix)
                    push!(ParamVec, PP)
                    push!(PairParams, PairPs)
                    push!(PkgVec, PkgP)
                end
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec,
                       PkgParams = PkgVec, filename="param_sweep_pairs.csv",
                       PairParams = PairParams)
    return df
end

function run_contact_sweeps_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))
    PC = [0, 0.5, 1.0, 2.0, 5.0] * BasicParcelParams["p_contact"]
    CR = [0, 0.1, 0.5, 1.0, 2.0, 5.0] * BasicParcelParams["CohortChangeRate"]
    II = [1,2,3]
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for pc in PC
        for ii in II
            for cr in CR
                PP = copy(BasicParcelParams)
                PP["p_contact"] = pc
                PP["CohortChangeRate"] = cr
                PP["InfInit"] = ii
                push!(ParamVec, PP)
                push!(PkgVec, PkgP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgVec,
                             NPPerDay = NPvec, filename="contact_sweep.csv")
    return df
end

function run_contacts_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    II = [1,2,3]
    PC = [0, 0.5, 1.0, 2.0, 5.0] * BasicBulkParams["p_contact"]
    CR = [0, 0.1, 0.5, 1.0, 2.0, 5.0] * BasicBulkParams["CohortChangeRate"]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
     for pc in PC
        for ii in II
            for cr in CR
                PP = copy(BasicBulkParams)
                PP["p_contact"] = pc
                PP["CohortChangeRate"] = cr
                PP["InfInit"] = ii
                PairPs = copy(BasicPairParams)
                push!(ParamVec, PP)
                push!(PairParams, PairPs)
                push!(PkgVec, PkgP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec, PkgParams = PkgVec,
                  filename="contact_sweep_pairs.csv", PairParams = PairParams)
    return df
end


function run_presenteeism_param_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))
    #other params
    PIsol = 0.1:0.1:1.0
    HHsharing = [0.05, 0.2, 0.5, 1.0]
    II = [1,2,3]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for hh in HHsharing
        for ii in II
            for pi in PIsol
                PP = copy(BasicParcelParams)
                PP["Pisol"] = pi
                PP["InfInit"] = ii
                PP["HouseShareFactor"] = hh
                push!(ParamVec, PP)
                push!(PkgVec, PkgP)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgVec,
                  NPPerDay = NPvec, filename="presenteeism_param_sweep.csv")
    return df
end

function run_presenteeism_param_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    PIsol = 0.1:0.1:1.0
    FPs = [false, true]
    II = [1,2,3]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for fp in FPs
        for ii in II
            for pi in PIsol
                PP = copy(BasicBulkParams)
                PairPs = copy(BasicPairParams)
                PP["Pisol"] = pi
                PP["InfInit"] = ii
                PairPs["fixed_driver_pairs"] = fp
                PairPs["fixed_loader_pairs"] = fp
                PairPs["PairIsolation"] = fp
                push!(ParamVec, PP)
                push!(PairParams, PairPs)
                push!(PkgVec, PkgP)
            end
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec, PkgParams = PkgVec,
                  PairParams = PairParams, filename="presenteeism_param_sweep_pairs.csv")
    return df
end

function run_car_house_share_sweep_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))
    #other params
    CarShareF = 0.0:0.2:1.0
    HouseShareF = 0.0:0.2:1.0
    CarShareIsol = [true, false]
    HouseShareIsol = [true, false]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for csf in CarShareF
        for hsf in HouseShareF
            for csisol in CarShareIsol
                for hsisol in HouseShareIsol
                    PP = copy(BasicParcelParams)
                    PP["CarShareFactor"] = csf
                    PP["HouseShareFactor"] = hsf
                    PP["CarShareIsolation"] = csisol
                    PP["HouseShareIsolation"] = hsisol
                    push!(ParamVec,PP)
                    push!(PkgVec, PkgP)
                end
            end
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec,
               PkgParams = PkgVec, filename="car_house_share_param_sweep_parcel.csv")
    return df

end

function run_car_house_share_sweep_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    #other params
    CarShareF = 0.0:0.2:1.0
    HouseShareF = 0.0:0.2:1.0
    CarShareIsol = [true, false]
    HouseShareIsol = [true, false]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PairPs = copy(BasicPairParams)
    PkgP = copy(BasicPkgParams)
    for csf in CarShareF
        for hsf in HouseShareF
            for csisol in CarShareIsol
                for hsisol in HouseShareIsol
                    PP = copy(BasicBulkParams)
                    PP["CarShareFactor"] = csf
                    PP["HouseShareFactor"] = hsf
                    PP["CarShareIsolation"] = csisol
                    PP["HouseShareIsolation"] = hsisol
                    push!(ParamVec,PP)
                    push!(PairParams,PairPs)
                    push!(PkgVec,PkgP)
                end
            end
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; NPPerDay = NPvec, PkgParams = PkgVec,
                  PairParams = PairParams, filename="car_house_share_param_sweep_pairs.csv")
    return df

end

function run_testing_sweep_outbreak_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))

    Enforced = [false,true]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicParcelParams)
    PkgP = copy(BasicPkgParams)
    for j in 1:2
        for tp in Tperiod
            for i in 1:length(TestType)
                TPh = copy(BasicTestingParams)
                push!(ParamVec, PP)
                push!(PkgVec, PkgP)
                TPh["testing_enforced"] = Enforced[j]
                TPh["tperiod"] = tp
                TPh["protocol"] = TestType[i]
                TPh["delay"] = Delay[i]
                TPh["test_pause"] = Test_pause[i]
                push!(TestParamVec, TPh)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgVec,
                  NPPerDay = NPvec, TestingParams=TestParamVec, output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    push!(PkgVec, PkgP)
    push!(ParamVec, PP)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
                  PkgParams = PkgVec, output = false)

    df2["is_testing"] = zeros(Bool,nrow(df2))
    df2["testing_enforced"] = zeros(Bool,nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill(0, nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    df2["test_miss_prob"] = zeros(nrow(df2))
    dfout = vcat(df,df2)

    CSV.write("testing_sweep.csv", dfout)

    return df
end

function run_testing_sweep_outbreak_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))

    Enforced = [false,true]
    Tperiod = 2:2:14
    Delay = [0,0,1,2]
    TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
    Test_pause = [21,90,90,90]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PP = copy(BasicBulkParams)
    PairPs = copy(BasicPairParams)
    PkgP = copy(BasicPkgParams)
    for j in 1:2
        for tp in Tperiod
            for i in 1:length(TestType)
                TPh = copy(BasicTestingParams)
                push!(ParamVec, PP)
                push!(PkgVec, PkgP)
                TPh["testing_enforced"] = Enforced[j]
                TPh["tperiod"] = tp
                TPh["protocol"] = TestType[i]
                TPh["delay"] = Delay[i]
                TPh["test_pause"] = Test_pause[i]
                push!(TestParamVec, TPh)
                push!(PairParams, PairPs)
            end
        end
    end

    df = run_many_sims(ParamVec, Nrepeats, OccPattern;  NPPerDay = NPvec, PkgParams = PkgVec,
             TestingParams=TestParamVec,  PairParams = PairParams, output = false)

    #run baseline case
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    push!(PkgVec, PkgP)
    push!(ParamVec, PP)
    push!(PairParams, PairPs)
    df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
                        PkgParams = PkgVec, PairParams = PairParams, output = false)
    df2["is_testing"] = zeros(Bool,nrow(df2))
    df2["testing_enforced"] = zeros(Bool,nrow(df2))
    df2["test_miss_prob"] = zeros(nrow(df2))
    df2["tperiod"] = zeros(nrow(df2))
    df2["protocol"] = fill(0, nrow(df2))
    df2["specificity"] = 0.999 * ones(nrow(df2))
    df2["delay"] = zeros(nrow(df2))
    df2["test_pause"] = zeros(nrow(df2))
    df = vcat(df,df2)

    CSV.write("testing_sweep_pairs.csv", df)
    return df
end

# function run_testing_sweep_fixedprev_scenario_parcel(Prev_val::Float64, Nrepeats::Int = 10000)
#     NWeeks= 26
#     OccPattern = repeat(ParcelOccPattern,Nweeks)
#     PkgPattern = repeat(ParcelPkgPattern,Nweeks)
#     NPvec = Int64.(round.(NPparcel*PkgPattern))

#     NCP = [0.0,0.5,1.0]
#     Tperiod = 2:2:14
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     Prev = Prev_val*ones(NWeeks*7)
#     Inc = Prev_val*ones(NWeeks*7)/7

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgVec = Array{Dict{Any,Any},1}(undef,0)
#     PP = copy(BasicParcelParams)
#     PkgP = copy(BasicPkgParams)
#     PP["SimType"] = Scenario_sim
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)
#                 push!(ParamVec, PP)
#                 push!(PkgVec, PkgP)
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>SpecDefault,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))

#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern;
#                   NPPerDay = NPvec, TestingParams=TestParamVec, output = false,
#                   PkgParams = PkgVec, Incidence = Inc, Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PkgVec = Array{Dict{Any,Any},1}(undef,0)
#     push!(PkgVec, PkgP)
#     push!(ParamVec, PP)
#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
#                   PkgParams = PkgVec, output = false, Incidence = Inc, Prevalence = Prev)

#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     dfout = vcat(df,df2)

#     CSV.write("testing_scenario_parcel_prev.csv", dfout)

#     return df
# end

# function run_testing_sweep_fixedprev_scenario_pairs(Prev_val::Float64, Nrepeats::Int = 10000)
#     NWeeks= 26
#     OccPattern = repeat(BulkOccPattern,NWeeks)
#     PkgPattern = repeat(BulkPkgPattern,NWeeks)
#     NPvec = Int64.(round.(NPbulk*PkgPattern))

#     NCP = [0.0,0.5,1.0]
#     Tperiod = 2:2:14
#     Delay = [0,0,1,2]
#     TestType = [LFD_mass_protocol,PCR_mass_protocol,PCR_mass_protocol,PCR_mass_protocol]
#     Test_pause = [21,90,90,90]

#     Prev = Prev_val*ones(NWeeks*7)
#     Inc = Prev_val*ones(NWeeks*7)/7

#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     TestParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     PkgVec = Array{Dict{Any,Any},1}(undef,0)
#     PP = copy(BasicBulkParams)
#     PP["SimType"] = Scenario_sim
#     PairPs = copy(BasicPairParams)
#     PkgP = copy(BasicPkgParams)
#     for j in 1:3
#         for tp in Tperiod
#             for i in 1:length(TestType)
#                 push!(ParamVec, PP)
#                 push!(PkgVec, PkgP)
#                 push!(TestParamVec, Dict("new_comply_prob"=>NCP[j], "tperiod"=>tp,
#                       "protocol"=>TestType[i], "specificity"=>SpecDefault,
#                       "delay"=>Delay[i], "test_pause"=>Test_pause[i]))
#                 push!(PairParams, PairPs)
#             end
#         end
#     end

#     df = run_many_sims(ParamVec, Nrepeats, OccPattern;  NPPerDay = NPvec,
#              TestingParams=TestParamVec,  PairParams = PairParams, output = false,
#              PkgParams = PkgVec, Incidence = Inc, Prevalence = Prev)

#     #run baseline case
#     ParamVec = Array{Dict{Any,Any},1}(undef,0)
#     PairParams = Array{Dict{Any,Any},1}(undef,0)
#     PkgVec = Array{Dict{Any,Any},1}(undef,0)
#     push!(PkgVec, PkgP)
#     push!(ParamVec, PP)
#     push!(PairParams, PairPs)
#     df2 = run_many_sims(ParamVec, Nrepeats*length(Tperiod), OccPattern;  NPPerDay = NPvec,
#              PkgParams = PkgVec, PairParams = PairParams, output = false,
#              Incidence = Inc, Prevalence = Prev)
#     df2["new_comply_prob"] = zeros(nrow(df2))
#     df2["tperiod"] = zeros(nrow(df2))
#     df2["protocol"] = fill("No testing", nrow(df2))
#     df2["specificity"] = 0.999 * ones(nrow(df2))
#     df2["delay"] = zeros(nrow(df2))
#     df2["test_pause"] = zeros(nrow(df2))
#     df = vcat(df,df2)

#     fname = string("testing_scenario_pairs_prev",string(100*Prev_val),".csv")
#     CSV.write(fname, df)
#     return df
# end

function run_param_sweep_outbreak_fomite_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPh = [2000,3000,4000,5000]
    II = [1,2,3]
    PFT = [0.0001,0.001,0.01]

    df = DataFrame()
    PP = copy(BasicParcelParams)
    PkgP = copy(BasicPkgParams)
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
                   NPPerDay = NPvec, output = false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                   NPPerDay = NPvec, output = false)
            df = vcat(df,dfh)
        end

    end
    CSV.write("fomite_param_sweep_parcel.csv", df)

    return df
end

function run_param_sweep_outbreak_fomite_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPh = [150,210,270,330]
    II = [1,2,3]
    PFT = [0.0001,0.001,0.01]

    df = DataFrame()
    PP = copy(BasicBulkParams)
    PairPs = copy(BasicPairParams)
    PkgP = copy(BasicPkgParams)
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
                   PairParams = PairParams, NPPerDay = NPvec, output = false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                   PairParams = PairParams, NPPerDay = NPvec, output = false)
            df = vcat(df,dfh)
        end

    end
    CSV.write("fomite_param_sweep_pairs.csv", df)

    return df
end

function run_all_interventions_separately_scenario_parcel(Prev::Array{Float64,1},
            Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(ParcelOccPattern,NWeeks)
    PkgPattern = repeat(ParcelPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    HHsharing = [0.05, 0.5]
    CarSharing = [0.05, 0.5]
    Adherence = [0.5, 0.9, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    TeamDistance = [1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    HouseShareIsolation = [false, false, false, true, false, false, false, false, false]
    Office_WFH = [false, false, false, false, true, false, false, false, false]
    Testing = [false, false, false, false, false, true, true, false, false]
    EnforcedTesting = [false, false, false, false, false, false, true, false, false]
    CarShareIsolation = [false, false, false, false, false, false, false, true, false]
    CohortIsolation = [false, false, false, false, false, false, false, false, true]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for j in 1:length(HHsharing)
        for i in 1:length(Office_WFH)
            PP = copy(BasicParcelParams)
            TP = copy(BasicTestingParams)
            PP["HouseShareFactor"] = HHsharing[j]
            PP["CarShareFactor"] = CarSharing[j]
            PP["Pisol"] = Adherence[i]
            PP["Office_WFH"] = Office_WFH[i]
            PP["TeamDistances"] = fill(TeamDistance[i],3)
            PP["HouseShareIsolation"] = HouseShareIsolation[i]
            PP["CarShareIsolation"] = CarShareIsolation[i]
            PP["CohortIsolation"] = CohortIsolation[i]
            PP["SimType"] = Scenario_sim
            TP["is_testing"] = Testing[i]
            TP["testing_enforced"] = EnforcedTesting[i]
            push!(ParamVec, PP)
            push!(TestParamVec, TP)
            push!(PkgVec, PkgP)
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev, PkgParams = PkgVec,
                  filename="all_interventions_sep_parcel.csv")
    return df
end

function run_all_interventions_separately_scenario_pairs(Prev::Array{Float64,1},
            Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(BulkOccPattern,NWeeks)
    PkgPattern = repeat(BulkPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    HHsharing = [0.05, 0.5]
    CarSharing = [0.05, 0.5]
    Adherence = [0.5, 0.9, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    TeamDistance = [1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    OpenWindows = [false, false, false, true, false, false, false, false, false, false, false]
    HouseShareIsolation = [false, false, false, false, true, false, false, false, false, false, false]
    FixedPairs = [false, false, false, false, false, true, false, false, false, false, false]
    PairIsolation = [false, false, false, false, false, true, false, false, false, false, false]
    Office_WFH = [false, false, false, false, false, false, true, false, false, false, false]
    Testing = [false, false, false, false, false, false, false, true, true, false, false]
    EnforcedTesting = [false, false, false, false, false, false, false, false, true, false, false]
    CarShareIsolation = [false, false, false, false, false, false, false, false, false, true, false]
    CohortIsolation = [false, false, false, false, false, false, false, false, false, false, true]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for j in 1:length(HHsharing)
        for i in 1:length(Office_WFH)
            PP = copy(BasicBulkParams)
            TP = copy(BasicTestingParams)
            PairPs = copy(BasicPairParams)
            PP["HouseShareFactor"] = HHsharing[j]
            PP["CarShareFactor"] = CarSharing[j]
            PP["Pisol"] = Adherence[i]
            PP["Office_WFH"] = Office_WFH[i]
            PP["TeamDistances"] = fill(TeamDistance[i],3)
            PP["HouseShareIsolation"] = HouseShareIsolation[i]
            PP["CarShareIsolation"] = CarShareIsolation[i]
            PP["CohortIsolation"] = CohortIsolation[i]
            PP["SimType"] = Scenario_sim
            TP["is_testing"] = Testing[i]
            TP["testing_enforced"] = EnforcedTesting[i]
            PairPs["fixed_driver_pairs"] = FixedPairs[i]
            PairPs["fixed_loader_pairs"] = FixedPairs[i]
            PairPs["PairIsolation"] = PairIsolation[i]
            push!(ParamVec, PP)
            push!(TestParamVec, TP)
            push!(PairParamVec, PairPs)
            push!(PkgVec, PkgP)
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev,
                  PairParams = PairParamVec, PkgParams = PkgVec,
                  filename="all_interventions_sep_pairs.csv")
    return df
end


function run_all_interventions_variableprev_scenario_parcel(Prev::Array{Float64,1},
            Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(ParcelOccPattern,NWeeks)
    PkgPattern = repeat(ParcelPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    #Non-isolation measures
    TeamDistance = [1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    NDteams = [3, 3, 6, 6, 6, 6, 6, 6, 6, 6]
    NLteams = [2, 2, 4, 4, 4, 4, 4, 4, 4, 4]
    NOteams = [1, 1, 3, 3, 3, 3, 3, 3, 3, 3]
    Office_WFH = [false, false, false, true, true, true, true, true, true, true]
    #isolation measures
    Adherence = [0.5, 0.5, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    HouseShareIsolation = [false, false, false, false, false, true, true, true, true, true]
    Testing = [false, false, false, false, false, false, true, true, true, true]
    EnforcedTesting = [false, false, false, false, false, false, false, true, true, true]
    CarShareIsolation = [false, false, false, false, false, false, false, false, true, true]
    CohortIsolation = [false, false, false, false, false, false, false, false, false, true]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for i in 1:length(Office_WFH)
        PP = copy(BasicParcelParams)
        TP = copy(BasicTestingParams)
        PP["Pisol"] = Adherence[i]
        PP["Office_WFH"] = Office_WFH[i]
        PP["TeamDistances"] = fill(TeamDistance[i],3)
        PP["HouseShareIsolation"] = HouseShareIsolation[i]
        PP["CarShareIsolation"] = CarShareIsolation[i]
        PP["CohortIsolation"] = CohortIsolation[i]
        PP["SimType"] = Scenario_sim
        TP["is_testing"] = Testing[i]
        TP["testing_enforced"] = EnforcedTesting[i]
        PP["NDteams"] = NDteams[i]
        PP["NLteams"] = NLteams[i]
        PP["NOteams"] = NOteams[i]
        push!(ParamVec, PP)
        push!(TestParamVec, TP)
        push!(PkgVec, PkgP)
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev, PkgParams = PkgVec,
                  filename="all_interventions_cumul_parcel.csv")
    return df
end

function run_all_interventions_variableprev_scenario_parcel_isolfirst(Prev::Array{Float64,1},
            Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(ParcelOccPattern,NWeeks)
    PkgPattern = repeat(ParcelPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    #isolation measures
    Adherence = [0.5, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    HouseShareIsolation = [false, false, true, true, true, true, true, true, true, true]
    Testing = [false, false, false, true, true, true, true, true, true, true]
    EnforcedTesting = [false, false, false, false, true, true, true, true, true, true]
    CarShareIsolation = [false, false, false, false, false, true, true, true, true, true]
    CohortIsolation = [false, false, false, false, false, false, true, true, true, true]

    #Non-isolation measures
    TeamDistance = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
    NDteams = [3, 3, 3, 3, 3, 3, 3, 3, 6, 6]
    NLteams = [2, 2, 2, 2, 2, 2, 2, 2, 4, 4]
    NOteams = [1, 1, 1, 1, 1, 1, 1, 1, 3, 3]
    Office_WFH = [false, false, false, false, false, false, false, false, false, true]


    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for i in 1:length(Office_WFH)
        PP = copy(BasicParcelParams)
        TP = copy(BasicTestingParams)
        PP["Pisol"] = Adherence[i]
        PP["Office_WFH"] = Office_WFH[i]
        PP["TeamDistances"] = fill(TeamDistance[i],3)
        PP["HouseShareIsolation"] = HouseShareIsolation[i]
        PP["CarShareIsolation"] = CarShareIsolation[i]
        PP["CohortIsolation"] = CohortIsolation[i]
        PP["SimType"] = Scenario_sim
        TP["is_testing"] = Testing[i]
        TP["testing_enforced"] = EnforcedTesting[i]
        PP["NDteams"] = NDteams[i]
        PP["NLteams"] = NLteams[i]
        PP["NOteams"] = NOteams[i]
        push!(ParamVec, PP)
        push!(TestParamVec, TP)
        push!(PkgVec, PkgP)
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev, PkgParams = PkgVec,
                  filename="all_interventions_cumul_parcel_isolfirst.csv")
    return df

end

function run_all_interventions_variableprev_scenario_pairs(Prev::Array{Float64,1},
             Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(BulkOccPattern,NWeeks)
    PkgPattern = repeat(BulkPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    #Non-isolation measures
    TeamDistance = [1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
    Office_WFH = [false, false, true, true, true, true, true, true, true, true]
    #isolation measures
    Adherence = [0.5, 0.5, 0.5, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    HouseShareIsolation = [false, false, false, false, true, true, true, true, true, true]
    FixedPairs = [false, false, false, false, false, true, true, true, true, true]
    Testing = [false, false, false, false, false, false, true, true, true, true]
    EnforcedTesting = [false, false, false, false, false, false, false, true, true, true]
    CarShareIsolation = [false, false, false, false, false, false, false, false, true, true]
    CohortIsolation = [false, false, false, false, false, false, false, false, false, true]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for i in 1:length(Office_WFH)
        PP = copy(BasicBulkParams)
        TP = copy(BasicTestingParams)
        PairPs = copy(BasicPairParams)
        PP["Pisol"] = Adherence[i]
        PP["Office_WFH"] = Office_WFH[i]
        PP["TeamDistances"] = fill(TeamDistance[i],3)
        PP["HouseShareIsolation"] = HouseShareIsolation[i]
        PP["CarShareIsolation"] = CarShareIsolation[i]
        PP["CohortIsolation"] = CohortIsolation[i]
        PP["SimType"] = Scenario_sim
        TP["is_testing"] = Testing[i]
        TP["testing_enforced"] = EnforcedTesting[i]
        PairPs["fixed_driver_pairs"] = FixedPairs[i]
        PairPs["fixed_loader_pairs"] = FixedPairs[i]
        PairPs["PairIsolation"] = FixedPairs[i]
        push!(ParamVec, PP)
        push!(TestParamVec, TP)
        push!(PairParamVec, PairPs)
        push!(PkgVec, PkgP)
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev,
                  PairParams = PairParamVec, PkgParams = PkgVec,
                  filename="all_interventions_cumul_pairs.csv")
    return df
end

function run_all_interventions_variableprev_scenario_pairs_isolfirst(Prev::Array{Float64,1},
             Inc::Array{Float64,1}, Demand::Array{Float64,1}, Nrepeats::Int = 10000)
    NWeeks = Int64(ceil(length(Prev)/7))
    OccPattern = repeat(BulkOccPattern,NWeeks)
    PkgPattern = repeat(BulkPkgPattern,NWeeks)
    OccPattern = OccPattern[1:length(Prev)]
    PkgPattern = PkgPattern[1:length(Prev)]
    NPvec = Int64.(round.(Demand.*PkgPattern))

    #isolation measures
    Adherence = [0.5, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9]
    HouseShareIsolation = [false, false, true, true, true, true, true, true, true, true]
    FixedPairs = [false, false, false, true, true, true, true, true, true, true]
    Testing = [false, false, false, false, true, true, true, true, true, true]
    EnforcedTesting = [false, false, false, false, false, true, true, true, true, true]
    CarShareIsolation = [false, false, false, false, false, false, true, true, true, true]
    CohortIsolation = [false, false, false, false, false, false, false, true, true, true]
    #Non-isolation measures
    TeamDistance = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0]
    Office_WFH = [false, false, false, false, false, false, false, false, false, true]

    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    TestParamVec = Array{Dict{Any,Any},1}(undef,0)
    PairParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgVec = Array{Dict{Any,Any},1}(undef,0)
    PkgP = copy(BasicPkgParams)
    for i in 1:length(Office_WFH)
        PP = copy(BasicBulkParams)
        TP = copy(BasicTestingParams)
        PairPs = copy(BasicPairParams)
        PP["Pisol"] = Adherence[i]
        PP["Office_WFH"] = Office_WFH[i]
        PP["TeamDistances"] = fill(TeamDistance[i],3)
        PP["HouseShareIsolation"] = HouseShareIsolation[i]
        PP["CarShareIsolation"] = CarShareIsolation[i]
        PP["CohortIsolation"] = CohortIsolation[i]
        PP["SimType"] = Scenario_sim
        TP["is_testing"] = Testing[i]
        TP["testing_enforced"] = EnforcedTesting[i]
        PairPs["fixed_driver_pairs"] = FixedPairs[i]
        PairPs["fixed_loader_pairs"] = FixedPairs[i]
        PairPs["PairIsolation"] = FixedPairs[i]
        push!(ParamVec, PP)
        push!(TestParamVec, TP)
        push!(PairParamVec, PairPs)
        push!(PkgVec, PkgP)
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern;
                  NPPerDay = NPvec, TestingParams=TestParamVec,
                  Incidence = Inc, Prevalence = Prev,
                  PairParams = PairParamVec, PkgParams = PkgVec,
                  filename="all_interventions_cumul_pairs_isolfirst.csv")
    return df
end

function run_param_sweep_outbreak_transmod_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPparcel*PkgPattern))
    F2F_mod = [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]
    SS_mod = [0.3, 0.5, 1.0, 1.5, 2.0, 2.5]

    PP = copy(BasicParcelParams)
    PkgP = copy(BasicPkgParams)
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    for (i,f2f) in enumerate(F2F_mod)
        PP["F2F_mod"] = f2f
        for (j,ss) in enumerate(SS_mod)
            PP["Aerosol_mod"] = ss
            push!(ParamVec, PP)
            push!(PkgParams, PkgP)
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                   NPPerDay = NPvec, filename = "transmod_sweep_parcel.csv")

    return df
end

function run_param_sweep_outbreak_transmod_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    NPvec = Int64.(round.(NPbulk*PkgPattern))
    F2F_mod = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    SS_mod = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    PP = copy(BasicBulkParams)
    PairPs = copy(BasicPairParams)
    PkgP = copy(BasicPkgParams)
    ParamVec = Array{Dict{Any,Any},1}(undef,0)
    PkgParams = Array{Dict{Any,Any},1}(undef,0)
    PairParams = Array{Dict{Any,Any},1}(undef,0)
    for (i,f2f) in enumerate(F2F_mod)
        PP["F2F_mod"] = f2f
        for (j,ss) in enumerate(SS_mod)
            PP["Aerosol_mod"] = ss
            push!(ParamVec, PP)
            push!(PkgParams, PkgP)
            push!(PairParams, PairPs)
        end
    end
    df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
            PairParams = PairParams, NPPerDay = NPvec, filename = "transmod_sweep_pairs.csv")

    return df
end

function run_param_sweep_outbreak_wpsize_parcel(Nrepeats::Int = 10000)
    OccPattern = repeat(ParcelOccPattern,NweeksDefault)
    PkgPattern = repeat(ParcelPkgPattern,NweeksDefault)

    rel_size = 0.5 .+ 1.5.*(0:20)./20

    df = DataFrame()
    PP = copy(BasicParcelParams)
    PkgP = copy(BasicPkgParams)
    for (i,rs) in enumerate(rel_size)
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        NPvec = Int64.(round.(rs*NPparcel*PkgPattern))
        PP["ND"] = Int64(round(rs*BasicParcelParams["ND"]))
        PP["NL"] = Int64(round(rs*BasicParcelParams["NL"]))
        PP["NO"] = Int64(round(rs*BasicParcelParams["NO"]))
        PP["NDteams"] = max(1, Int64(round(rs*BasicParcelParams["NDteams"])))
        PP["NLteams"] = max(1, Int64(round(rs*BasicParcelParams["NLteams"])))
        PP["NOteams"] = max(1, Int64(round(rs*BasicParcelParams["NOteams"])))
        PP["p_contact"] = BasicParcelParams["p_contact"]*
                         (NStaffparcel_def)/(PP["ND"]+PP["NL"]+PP["NO"])
        PP["CohortChangeRate"] = BasicParcelParams["CohortChangeRate"]*
                                 (NStaffparcel_def)/(PP["ND"]+PP["NL"]+PP["NO"])
        push!(ParamVec, PP)
        push!(PkgParams, PkgP)
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                   NPPerDay = NPvec, output=false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                   NPPerDay = NPvec, output=false)
            df = vcat(df,dfh)
        end
    end

    CSV.write("wpsize_param_sweep_parcel.csv", df)

    return df
end

function run_param_sweep_outbreak_wpsize_pairs(Nrepeats::Int = 10000)
    OccPattern = repeat(BulkOccPattern,NweeksDefault)
    PkgPattern = repeat(BulkPkgPattern,NweeksDefault)
    rel_size = 0.5 .+ 2.5.*(0:20)./20

    df = DataFrame()
    PP = copy(BasicBulkParams)
    PairPs = copy(BasicPairParams)
    PkgP = copy(BasicPkgParams)
    for (i,rs) in enumerate(rel_size)
        ParamVec = Array{Dict{Any,Any},1}(undef,0)
        PkgParams = Array{Dict{Any,Any},1}(undef,0)
        PairParams = Array{Dict{Any,Any},1}(undef,0)
        NPvec = Int64.(round.(rs*NPbulk*PkgPattern))
        PP["ND"] = Int64(2*max(1,round(rs*BasicBulkParams["ND"]/2)))
        PP["NL"] = Int64(2*max(1,round(rs*BasicBulkParams["NL"]/2)))  #ensure even number
        PP["NO"] = Int64(round(rs*BasicBulkParams["NO"]))
        PP["NDteams"] = max(1, Int64(round(rs*BasicBulkParams["NDteams"])))
        PP["NLteams"] = max(1, Int64(round(rs*BasicBulkParams["NLteams"])))
        PP["NOteams"] = max(1, Int64(round(rs*BasicBulkParams["NOteams"])))
        PP["p_contact"] = BasicBulkParams["p_contact"]*
                         (NStaffbulk_def)/(PP["ND"]+PP["NL"]+PP["NO"])
        PP["CohortChangeRate"] = BasicBulkParams["CohortChangeRate"]*
                                 (NStaffbulk_def)/(PP["ND"]+PP["NL"]+PP["NO"])
        push!(ParamVec, PP)
        push!(PkgParams, PkgP)
        push!(PairParams, PairPs)
        if i == 1
            df = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                PairParams = PairParams, NPPerDay = NPvec, output=false)
        else
            dfh = run_many_sims(ParamVec, Nrepeats, OccPattern; PkgParams = PkgParams,
                PairParams = PairParams, NPPerDay = NPvec, output=false)
            df = vcat(df,dfh)
        end
    end

    CSV.write("wpsize_param_sweep_pairs.csv", df)

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
