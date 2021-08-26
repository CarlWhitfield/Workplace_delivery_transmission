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
            if length(Params[p]) > 1
                results[(irow_start):(irow_start+3),p] .= Params[p]
            else
                results[(irow_start):(irow_start+2),p] = Params[p]
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
    results[(irow_start+3), "IsolatorsFrac"] = sum(SimOutput["NewIsolators"])/(NStot)
    results[(irow_start+3), "SympIsolatorsFrac"] = sum(SimOutput["NewSympIsolators"])/(NStot)
    results[(irow_start+3), "FPIsolatorsFrac"] = sum(SimOutput["NewFalseIsolators"])/(NStot)
    results[(irow_start+3), "TPSympIsolatorsFrac"] = sum(SimOutput["NewTestSympIsolators"])/(NStot)
    results[(irow_start+3), "TPAsympIsolatorsFrac"] = sum(SimOutput["NewTestAsympIsolators"])/(NStot)

    #multi-liners
    results[(irow_start):(irow_start+3),"TotInfPackagesDelivered"] .=
                sum(SimOutput["PackagesInfectiousOnDelivery"])
    results[(irow_start):(irow_start+3), "CustomersInfectedByPkg"] .=
                sum(SimOutput["CustomersInfectedByPkgs"])
    results[(irow_start):(irow_start+3), "CustomersInfectedByDrivers"] .=
                sum(SimOutput["CustomersInfectedByDrivers"])
end
