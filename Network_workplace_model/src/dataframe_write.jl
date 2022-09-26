using SharedArrays

include("network_transmission_workplace.jl")


function init_results_dataframe(Nrows::Int, AllParams::Dict)

    IntColNames = [:NStaff, :Iteration, :CustomersInfected, :GroupNo]
    FloatColNames = [:FracRecovered, :FomiteInfectionFrac, :CohortInfectionFrac,
    :RandContactInfectionFrac, :PairInfectionFrac, :RoomInfectionFrac,
    :ExtIntroFrac, :CustIntroFrac, :CarShareInfFrac, :HouseShareInfFrac,
    :IsolatorsFrac, :SympIsolatorsFrac, :FPIsolatorsFrac, :TPSympIsolatorsFrac,
    :TPAsympIsolatorsFrac]

    if AllParams["SimType"] == Outbreak_sim
        push!(IntColNames,:IndexCaseInfections,:OverallOutbreakLength)
        push!(FloatColNames,:IndexCaseViralLoad ,:IndexCaseInfectivity)
    end
    for p in keys(AllParams)
        if p != "ND" && p != "NO" && p != "NL"
            p0 = 0
            if length(AllParams[p]) > 1
                p0 = AllParams[p][1]
                # results[:,p] = fill(AllParams[p][1], Nrows)
            else
                p0 = AllParams[p]
                # results[:,p] = fill(AllParams[p], Nrows)
            end
            if (typeof(p0) <: Int || typeof(p0) <: Bool)
                push!(IntColNames,Symbol(p))
            elseif typeof(p0) <: AbstractFloat
                push!(FloatColNames,Symbol(p))
            end
        end
    end
    IntArray = SharedArray{Int64}(Nrows,length(IntColNames))
    FloatArray = SharedArray{Float64}(Nrows,length(FloatColNames))

    IntColMap = Dict{Symbol,Int}()
    for (i,col) in enumerate(IntColNames)
        IntColMap[col] = i
    end
    FloatColMap = Dict{Symbol,Int}()
    for (i,col) in enumerate(FloatColNames)
        FloatColMap[col] = i
    end

    return IntArray, FloatArray, IntColMap, FloatColMap
end

function add_to_shared_array!(ResultsArray::SharedArray, Value, icol::Int,
                              irow_start::Int)
    if length(Value) > 1
        ResultsArray[(irow_start):(irow_start+2),icol] .= Value
        ResultsArray[(irow_start+3),icol] = Value[1]
    else
        ResultsArray[(irow_start):(irow_start+3),icol] .= Value
    end
end

function add_to_results_dataframe!(IntArray::SharedArray{Int64,2},
    FloatArray::SharedArray{Float64,2}, IntColMap::Dict, FloatColMap::Dict,
    Params::Dict, SimOutput::Dict, irow_start::Int, Niteration::Int)

    for p in keys(Params)
        colname = Symbol(p)
        if colname in keys(FloatColMap)
            add_to_shared_array!(FloatArray,Params[p],FloatColMap[colname],irow_start)
        elseif colname in keys(IntColMap)
            add_to_shared_array!(IntArray,Params[p],IntColMap[colname],irow_start)
        end
    end

    if Params["SimType"] == Outbreak_sim
        IntArray[(irow_start):(irow_start+3),IntColMap[:IndexCaseInfections]] .=
                        SimOutput["IndexCaseInfections"]
        IntArray[(irow_start):(irow_start+3),IntColMap[:OverallOutbreakLength]] .=
                        length(SimOutput["time"])
        FloatArray[(irow_start):(irow_start+3),FloatColMap[:IndexCaseViralLoad]] .=
                        SimOutput["IndexCasePeakVL"]
        FloatArray[(irow_start):(irow_start+3),FloatColMap[:IndexCaseInfectivity]] .=
                        SimOutput["IndexCaseInfectivity"]
    end

    #params and results for specific staff groups

    Nh = [Params["ND"], Params["NL"], Params["NO"]]
    for i in 1:3
        IntArray[(irow_start + i - 1),IntColMap[:GroupNo]] = i
        IntArray[(irow_start + i - 1), IntColMap[:NStaff]] = Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:FracRecovered]] =
                            SimOutput["Recovered"][i,end] / Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:FomiteInfectionFrac]] =
                      sum(SimOutput["InfsByType"][package_contact][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:CohortInfectionFrac]] =
                     sum(SimOutput["InfsByType"][network_contact][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:RandContactInfectionFrac]] =
                 sum(SimOutput["InfsByType"][non_network_contact][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:PairInfectionFrac]] =
                 sum(SimOutput["InfsByType"][pair_contact][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:RoomInfectionFrac]] =
                 sum(SimOutput["InfsByType"][room_transmission][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:CustIntroFrac]] =
                 sum(SimOutput["InfsByType"][customer_contact][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:ExtIntroFrac]] =
                 sum(SimOutput["InfsByType"][introduction][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:CarShareInfFrac]] =
                 sum(SimOutput["InfsByType"][car_share][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:HouseShareInfFrac]] =
                 sum(SimOutput["InfsByType"][house_share][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:IsolatorsFrac]] =
                 sum(SimOutput["NewIsolators"][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:SympIsolatorsFrac]] =
                 sum(SimOutput["NewSympIsolators"][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:FPIsolatorsFrac]] =
                 sum(SimOutput["NewFalseIsolators"][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:TPSympIsolatorsFrac]] =
                 sum(SimOutput["NewTestSympIsolators"][i,:])/ Nh[i]
        FloatArray[(irow_start + i - 1),FloatColMap[:TPAsympIsolatorsFrac]] =
                 sum(SimOutput["NewTestAsympIsolators"][i,:])/ Nh[i]
    end

    IntArray[(irow_start+3),IntColMap[:GroupNo]] = 4
    NStot = Params["ND"] + Params["NL"] + Params["NO"]
    IntArray[(irow_start+3), IntColMap[:NStaff]] = NStot
    FloatArray[(irow_start+3),FloatColMap[:FracRecovered]] =
                    sum(SimOutput["Recovered"][:,end]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:FomiteInfectionFrac]] =
                    sum(SimOutput["InfsByType"][package_contact]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:CohortInfectionFrac]] =
                    sum(SimOutput["InfsByType"][network_contact]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:RandContactInfectionFrac]] =
                    sum(SimOutput["InfsByType"][non_network_contact]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:PairInfectionFrac]] =
                    sum(SimOutput["InfsByType"][pair_contact]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:RoomInfectionFrac]] =
                    sum(SimOutput["InfsByType"][room_transmission]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:CustIntroFrac]] =
                    sum(SimOutput["InfsByType"][customer_contact]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:ExtIntroFrac]]=
                    sum(SimOutput["InfsByType"][introduction]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:CarShareInfFrac]] =
                    sum(SimOutput["InfsByType"][car_share]) / NStot
    FloatArray[(irow_start+3),FloatColMap[:HouseShareInfFrac]] =
                    sum(SimOutput["InfsByType"][house_share]) / NStot
    FloatArray[(irow_start+3), FloatColMap[:IsolatorsFrac]] =
                    sum(SimOutput["NewIsolators"])/(NStot)
    FloatArray[(irow_start+3), FloatColMap[:SympIsolatorsFrac]] =
                    sum(SimOutput["NewSympIsolators"])/(NStot)
    FloatArray[(irow_start+3), FloatColMap[:FPIsolatorsFrac]] =
                    sum(SimOutput["NewFalseIsolators"])/(NStot)
    FloatArray[(irow_start+3), FloatColMap[:TPSympIsolatorsFrac]] =
                    sum(SimOutput["NewTestSympIsolators"])/(NStot)
    FloatArray[(irow_start+3), FloatColMap[:TPAsympIsolatorsFrac]] =
                    sum(SimOutput["NewTestAsympIsolators"])/(NStot)

    #multi-liners
    IntArray[(irow_start):(irow_start+3), IntColMap[:CustomersInfected]] .=
                sum(SimOutput["CustomersInfected"])
end

function create_dataframe_from_arrays(IntArray::SharedArray{Int64,2},
    FloatArray::SharedArray{Float64,2}, IntColMap::Dict{Symbol,Int},
    FloatColMap::Dict{Symbol,Int})

    df = DataFrame()
    for colname in keys(FloatColMap)
        df[!,colname] = FloatArray[:,FloatColMap[colname]]
    end
    for colname in keys(IntColMap)
        df[!,colname] = IntArray[:,IntColMap[colname]]
    end
    #fill in group numbers
    Names = ["Drivers","Pickers","Office","All"]
    df[!,:Group] = Names[IntArray[:,IntColMap[:GroupNo]]]

    return df
end
