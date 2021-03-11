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
    PIsol = [0.25,0.5,0.75,1.0]
    Phi = [0.05,0.25,1.0]
    PFC = [0.25,0.5,0.75,1.0]
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
                                         "p_friend_contact"=>pf))
                end
            end
        end
    end

    run_param_sweeps_network(ParamVec, Nrepeats, OccPattern; filename="param_sweep.csv")
end
