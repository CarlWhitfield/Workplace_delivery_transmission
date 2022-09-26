include("../src/simple_workplace_model.jl")

#select code options
VL_model = ke_model_no         #choice of viral load model (Ke et al is default)
LFD_model = porton_down_p3b         #choice of LFD sensitivity model (social care binned 2021 data is default)
PCR_sens_max = 0.95            #max PCR sensitivity (Ferretti et al 2021)
Inf_model = ke_inf_model_no    #infectivity model

function run_one_component_workplace_example()
    Nreps = 1000  #simulation repeats
    NDays = 300   #sim days
    #set shift pattern
    shift_pattern = [true, true, true, true, false, false, false,
                     true, true, true, true, true, false, false];
                     #set parameters
    Params = Dict("N"=>[100], "Pisol"=>1.0, "Psusc"=>1.0, "InfInit"=>0,
               "contact_mat"=>[0.296], "SimType"=>Outbreak_sim,
               "random_start"=>true)

    #run simulations and extract final outbreak sizes
    FinalSizesNoTesting = zeros(Nreps)
    for i in 1:Nreps
        summary, infs = run_sim_wp(Params, shift_pattern, NDays)
        FinalSizesNoTesting[i] = Params["N"][1] - 1 - summary["Susceptible"][end]
        if mod(i,100)==0 || i == Nreps
            print("No testing: ",i,'/',Nreps,'\n')
        end
    end

    FinalSizes2LFDs = zeros(Nreps)
    #set testing params
    test_params = Dict("is_testing"=>[true], "test_miss_prob"=>0.3, "policy_adherence"=>1.0,
                        "testing_enforced"=>false, "protocol"=>LFD_pattern,
                        "Pattern"=>[true, false, true, false, false, false, false,
                                    true, false, false, true, false, false, false],
                        "test_pause"=>21.0, "specificity"=>0.999, "delay"=>0)
    #run same simulations with testing
    for i in 1:Nreps
        summary, infs = run_sim_wp(Params, shift_pattern, NDays; TestParams=test_params)
        FinalSizes2LFDs[i] = Params["N"][1] - 1 - summary["Susceptible"][end]
        if mod(i,100)==0 || i == Nreps
            print("Testing: ",i,'/',Nreps,'\n')
        end
    end

    return FinalSizesNoTesting, FinalSizes2LFDs
end

function run_two_component_workplace_example()
    Nrep = 1000
    NDays = 300
    #baseline contact rate
    pc = 0.296
    scale_mat = ones(2,2)
    alpha = 0.2:0.2:1.8
    Np = [50,30]
    fs = 9/14
    beta = 1 .- (alpha.-1).*(Np[1]*fs*(Np[1]*fs - 1))/((Np[2]*(Np[2] - 1)))
    #base parameters
    Params = Dict("N"=>Np, "Pisol"=>1.0, "Psusc"=>1.0, "InfInit"=>1,
              "phi"=>1.0, "SimType"=>Outbreak_sim, "random_start"=>true)

    #define testing regime
    test_params = Dict("is_testing"=>[true,false], "test_miss_prob"=>0.3, "policy_adherence"=>1.0,
             "testing_enforced"=>false, "protocol"=>LFD_pattern,
             "Pattern"=>[true, false, true, false, false, false, false,
                         true, false, false, true, false, false, false],
             "test_pause"=>21.0, "specificity"=>0.999, "delay"=>0)
    #define shift pattern, group 1 is workers, 2 is residents (always 'at work')
    shift_patterns = ones(Bool,14,2)
    shift_patterns[:,1] .= [true, true, true, true, false, false, false,
                        true, true, true, true, true, false, false]

    #interested in counting resident infections only
    resident_infections_no_test = zeros(Float64,Nrep,length(alpha))
    resident_infections_test = zeros(Float64,Nrep,length(alpha))
    for na in 1:length(alpha)
        #loop over different between-group contact rates
        scale_mat[1,1] = alpha[na]
        scale_mat[2,2] = beta[na]
        Params["contact_mat"] = pc .* scale_mat

        for n in 1:Nrep
            #simulate testing and no testing cases
            sim_summary_no_test, infs_no_test = run_sim_wp(Params, shift_patterns, NDays)
            sim_summary_test, infs_test = run_sim_wp(Params, shift_patterns,
                                                    NDays; TestParams=test_params)
            #collect resident and worker infections
            infs_by_role = zeros(Float64,2)
            for k in 1:3
                infs_by_role += sum(sim_summary_no_test["InfsByType"][k],dims=2)
            end
            #store resident infections
            resident_infections_no_test[n,na] = infs_by_role[2]

            #repeat for testing case
            infs_by_role = zeros(Float64,2)
            for k in 1:3
                infs_by_role += sum(sim_summary_test["InfsByType"][k],dims=2)
            end
            resident_infections_test[n,na] = infs_by_role[2]
            if mod(n,100)==0 || n == Nrep
                print("Alpha loop: ", na, '/', length(alpha),", repeat: ", n,
                                     '/',Nrep,'\n')
            end
        end
    end

    return alpha, beta, resident_infections_no_test, resident_infections_test
end
