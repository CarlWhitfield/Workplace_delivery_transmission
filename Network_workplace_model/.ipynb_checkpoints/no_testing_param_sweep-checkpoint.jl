include("network_transmission_workplace.jl")

function run_param_sweep()
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
