 include("transmission_model_framework.jl")

function run_example_sim()

    #Example parameters for initialisation
    N_per_role = [10,10]   #number of people in each job role (number of elements gives number of roles)
    Pisol = 0.9           #probability that a person isolates if they develop symptoms of COVID-19 (Symptomatic fraction
                          #is set to 0.5 and is defined in definitions.jl)
    Psusc = 1.0          #Probability that an individual starts the simulation susceptible
    NDays = 300
    Incidence = 0.005*ones(NDays)   #daily probability of a person getting infected in the community
                                 #the length of this vector determines how many days to run the simulation for
    Prevalence = 14*Incidence    #community prevalence (not actually used here)

    #this returns the simulation
    sim = init_transmission_model(N_per_role, Pisol, Psusc, Incidence, Prevalence)

    """
    these parameters are handed in Dict form
    two_options for "SimType":
    1) Outbreak_sim, infect one case at the start, run until no infectious people remain or we hit NDays
    2) Scenarion_sim, don't infect anybody at the start, run until no infectious people remain

    options for "InfInit"
    Only used if "Outbreak_sim" is selected, "InfInit" => 0 means that initial infection is random
    "InfInit" => 1,...,Njobtypes initial infection occurs in that group of workers
    "NinfTypes" => 2 is the number of types of contact that you want to simulate, in this case contacts
     will need to be labelled as either type 1 or 2
     """
    Params = Dict("SimType" => Outbreak_sim, "InfInit" => 0, "NContactTypes" => 2)


    """parameters for simulating testing policies are also added here, but for this simple
    example I will assume testing is turned off, which can be acheived as follows"""
    TestParams = Dict("is_testing" => false)

    sim_summary, start_day, Anyinf = setup_transmission_model!(sim, Params, TestParams, NDays)
    """sim_summary is another Dict for containing summary info about the simulation
    start_day is a random number from 1 to 7, determining the day of week that the simulation
    starts (can be largely ignored unless you want to include DOW effects)
    Anyinf is True or False, confirming whether there are any infected people in the population"""


    #create a very basic fixedcontact network
    contact_type = 1    #for illustration purposes only, you might use this to track different types of contacts
    weight = return_infection_weight(1.0, 0.25, false, true)  #15 minute contact @ 1m distance, inside, while talking
    for n1 in 1:sim["Ntot"]  #node into edge
        #loop over all possible edges
        n2vec = collect(1:(n1-1))
        n2selected = randsubseq(n2vec, 0.1) #select 10% of them to actually have edges
        for n2 in n2selected
            add_contact_to_contact_network!(sim,n1,n2,weight,contact_type)
        end
    end
    #add a second type of contact
    contact_type = 2
    weight = return_infection_weight(2.0, 1.0, false, false)   #1 hour contact @ 2m distance, inside, no talking
    for n1 in 1:sim["Ntot"]
        for n2 in 1:sim["Ntot"]
            add_contact_to_contact_network!(sim,n1,n2,weight,contact_type)
        end
    end

    #do sim loop
    i_day = start_day   #day counter
    Go = true
    while i_day < NDays && Go
        update_all_statuses!(sim, i_day)   #update infectious status etc.
        if TestParams["is_testing"]
            new_isolators = do_testing!(sim, TestParams, i_day, sim["isolation_network"])
        end
        update_in_work!(sim, ones(Bool, sim["Ntot"]))   #assume everyone is at work every day
        infpairs = get_network_infections(sim, i_day)
        infpairs_final = do_infections_randomly!(infpairs, sim, i_day)
        update_sim_summary!(sim_summary, sim, infpairs_final, i_day)

        i_day += 1
        if Params["SimType"] == Outbreak_sim  #trigger for outbreak sim to stop
            if !any(get_infected_bool(sim))
                Go = false
            end
        end
    end
    trim_sim_summary!(sim_summary,i_day-1,NDays)   #trim summary object down so it ends on last day of sim

    return sim_summary
end
