include("network_transmission_workplace.jl")

function init_results_dataframe(Nrows::Int, AllParams::Dict)
    results = DataFrame([String, Int64, Int64, Float64,
                         Float64, Float64, Float64, Float64,
                         Float64, Float64, Float64, Float64, Float64, Int64,
                         Float64, Float64, Float64, Float64, Float64],
                         [:Group, :NStaff, :Iteration, :FracRecovered,
                         :FomiteInfectionFrac, :CohortInfectionFrac,
                         :RandContactInfectionFrac,
                         :PairInfectionFrac, :RoomInfectionFrac,
                         :ExtIntroFrac, :CustIntroFrac, :CarShareInfFrac,
                         :HouseShareInfFrac, :CustomersInfected, :IsolatorsFrac,
                         :SympIsolatorsFrac, :FPIsolatorsFrac,
                         :TPSympIsolatorsFrac, :TPAsympIsolatorsFrac], Nrows)

    if AllParams["SimType"] == Outbreak_sim
        results.IndexCaseInfections = zeros(Int64,Nrows)
        results.IndexCaseViralLoad = zeros(Float64,Nrows)
        results.IndexCaseInfectivity = zeros(Float64,Nrows)
        results.OverallOutbreakLength = zeros(Int64,Nrows)
    end
    for p in keys(AllParams)
        if p != "ND" && p != "NO" && p != "NL"
            if length(AllParams[p]) > 1
                results[:,p] = fill(AllParams[p][1], Nrows)
            else
                results[:,p] = fill(AllParams[p], Nrows)
            end
        end
    end

    return results
end

function add_to_results_dataframe!(results::DataFrame, Params::Dict, SimOutput::Dict,
                                   irow_start::Int, Niteration::Int)

    #params and results common to all staff
    for p in keys(Params)
        if p in names(results)
            if length(Params[p]) > 1
                results[(irow_start):(irow_start+2),p] .= Params[p]
                results[(irow_start+3),p] = Params[p][1]
            else
                results[(irow_start):(irow_start+3),p] .= Params[p]
            end
        end
    end
    results[(irow_start):(irow_start+3),"Iteration"] .= Niteration
    if Params["SimType"] == Outbreak_sim
        results[(irow_start):(irow_start+3),"IndexCaseInfections"] .= SimOutput["IndexCaseInfections"]
        results[(irow_start):(irow_start+3),"IndexCaseViralLoad"] .= SimOutput["IndexCasePeakVL"]
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
        results[(irow_start + i - 1),"FracRecovered"] =
                            SimOutput["Recovered"][i,end] / Nh[i]
        results[(irow_start + i - 1),"FomiteInfectionFrac"] =
                      sum(SimOutput["InfsByType"][package_contact][i,:])/ Nh[i]
        results[(irow_start + i - 1),"CohortInfectionFrac"] =
                     sum(SimOutput["InfsByType"][network_contact][i,:])/ Nh[i]
        results[(irow_start + i - 1),"RandContactInfectionFrac"] =
                 sum(SimOutput["InfsByType"][non_network_contact][i,:])/ Nh[i]
        results[(irow_start + i - 1),"PairInfectionFrac"] =
                 sum(SimOutput["InfsByType"][pair_contact][i,:])/ Nh[i]
        results[(irow_start + i - 1),"RoomInfectionFrac"] =
                 sum(SimOutput["InfsByType"][room_transmission][i,:])/ Nh[i]
        results[(irow_start + i - 1),"CustIntroFrac"] =
                 sum(SimOutput["InfsByType"][customer_contact][i,:])/ Nh[i]
        results[(irow_start + i - 1),"ExtIntroFrac"] =
                 sum(SimOutput["InfsByType"][introduction][i,:])/ Nh[i]
        results[(irow_start + i - 1),"CarShareInfFrac"] =
                 sum(SimOutput["InfsByType"][car_share][i,:])/ Nh[i]
        results[(irow_start + i - 1),"HouseShareInfFrac"] =
                 sum(SimOutput["InfsByType"][house_share][i,:])/ Nh[i]
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
    results[(irow_start+3),"FomiteInfectionFrac"] =
                        sum(SimOutput["InfsByType"][package_contact]) / NStot
    results[(irow_start+3),"CohortInfectionFrac"] =
                        sum(SimOutput["InfsByType"][network_contact]) / NStot
    results[(irow_start+3),"RandContactInfectionFrac"] =
                    sum(SimOutput["InfsByType"][non_network_contact]) / NStot
    results[(irow_start+3),"PairInfectionFrac"] =
                    sum(SimOutput["InfsByType"][pair_contact]) / NStot
    results[(irow_start+3),"RoomInfectionFrac"] =
                    sum(SimOutput["InfsByType"][room_transmission]) / NStot
    results[(irow_start+3),"CustIntroFrac"] =
                    sum(SimOutput["InfsByType"][customer_contact]) / NStot
    results[(irow_start+3),"ExtIntroFrac"] =
                    sum(SimOutput["InfsByType"][introduction]) / NStot
    results[(irow_start+3),"CarShareInfFrac"] =
                    sum(SimOutput["InfsByType"][car_share]) / NStot
    results[(irow_start+3),"HouseShareInfFrac"] =
                    sum(SimOutput["InfsByType"][house_share]) / NStot
    results[(irow_start+3), "IsolatorsFrac"] = sum(SimOutput["NewIsolators"])/(NStot)
    results[(irow_start+3), "SympIsolatorsFrac"] = sum(SimOutput["NewSympIsolators"])/(NStot)
    results[(irow_start+3), "FPIsolatorsFrac"] = sum(SimOutput["NewFalseIsolators"])/(NStot)
    results[(irow_start+3), "TPSympIsolatorsFrac"] = sum(SimOutput["NewTestSympIsolators"])/(NStot)
    results[(irow_start+3), "TPAsympIsolatorsFrac"] = sum(SimOutput["NewTestAsympIsolators"])/(NStot)

    #multi-liners
    results[(irow_start):(irow_start+3), "CustomersInfected"] .=
                sum(SimOutput["CustomersInfected"])
end
