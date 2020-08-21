using Markdown      # 引入Markdown支持
using Plots
using LaTeXStrings
using Printf
using DelimitedFiles


include("Methods.jl")


function _Measure!(
    Status::Array{Int64, 2}, 
    l::Int64, 
    n::Int64, 
    timeNum::Int64,
    measureCount::Int64,
    Tobs::Array{Array{Float64, 1}},
    Data::Array{Array{Float64, 1}}
)
    e::Float64 = 0
    m::Float64 = 0

    for xIndex = 1:l, yIndex = 1:l
        XNei, YNei = _Nei(xIndex, yIndex, l)
        xNext = XNei[2]
        yNext = YNei[2]

        e -= Status[xIndex, yIndex]*(Status[xNext, yIndex] + Status[xIndex, yNext])
        m += Status[xIndex, yIndex]
    end
    e /= n
    m = abs(m)/n

    if measureCount < 1
        Tobs[timeNum + measureCount] = [e, m]
    else
        Tobs[2:end] .= Tobs[1:end - 1]
        Tobs[1] = [e, m]

        Data[1] .+= [e, m]
        for i = 2:timeNum + 1
            @. Data[i] += Tobs[1] * Tobs[i - 1]
        end
    end

    return nothing
end


function _Main()
    #=
    读取配置
    =#
    Config = readdlm(
        "Config.in",
        ',',
        comments = true
    )
    Fns = Config[1, :] .|> (Fn -> "_MC"*Fn*"!") .|> Symbol
    lPowMax, lPowMin = Config[2, :]
    TS = Config[3, 1:2]
    bins, binSteps = Config[4, :]
    timeNum, tFreq = Config[5, :]

    for (fnCount, fn) in enumerate(Fns)
        println("\n", "The method used now is:\t", Config[1, fnCount], "\n", "*"^30)

        @eval _MC! = $fn

        LPows = lPowMin:lPowMax |> collect
        TPlot = plot(
            size = (1600, 800),
            xlabel = L"L",
            ylabel = L"\Theta_{int}",
            yscale = :log10,
            minorticks = 10
        )

        for T in TS
            println("\n", "The temperature is:\t", T, "\n", "*"^20)

            Aints = similar(LPows, Array{Float64, 1})
            for index in eachindex(Aints)
                Aints[index] = zeros(2)
            end
            Eints = similar(Aints)
            for index in eachindex(Eints)
                Eints[index] = zeros(2)
            end
            lPlot = plot(
                title = "T = $(T)" |> LaTeXString,
                size = (1600, 800),
                xlabel = L"\tau \quad [MC \quad steps]",
                ylabel = L"A_{\left| M \right|}(\tau)",
                yscale = :log10,
                minorticks = 10
            )

            for (lCount, lPow) in enumerate(LPows)
                l = 2^lPow
                n = l^2
                RawStatus = rand((-1, 1), (l, l))
                ADatas = Array{Array{Float64, 1}}(undef, bins, timeNum + 1)
                println("\n", "The length of eage is:\t", l, "\n", "*"^10)

                Threads.@threads for i = 1:bins
                    Status = deepcopy(RawStatus)
                    measureCount = -timeNum
                    Tobs = Array{Array{Float64, 1}}(undef, timeNum)
                    for index in eachindex(Tobs)
                        Tobs[index] = zeros(2)
                    end
                    Data = Array{Array{Float64, 1}}(undef, timeNum + 1)
                    for index in eachindex(Data)
                        Data[index] = zeros(2)
                    end

                    for j = 1:binSteps
                        _MC!(Status, l, n, T)

                        if j % tFreq == 0
                            measureCount += 1
                            _Measure!(Status, l, n, timeNum, measureCount, Tobs, Data)
                        end
                    end

                    println(Threads.threadid(), "\t", i, "\t", Data[1])
                    @. ADatas[i, :] = Data/measureCount
                end

                Avt = Array{Array{Float64, 1}}(undef, timeNum)
                for index in eachindex(Avt)
                    Avt[index] = zeros(2)
                end
                @inbounds for i = 1:bins, j = 1:timeNum
                    @. @fastmath Avt[j] += ADatas[i, j + 1] - ADatas[i, 1]^2
                end
                @fastmath @inbounds Avt = Avt .|> (X -> X./Avt[1])
                rightLimit = findfirst((X -> (X |> last) < 1e-2), Avt)

                Ert = Array{Array{Float64, 1}}(undef, timeNum)
                for index in eachindex(Ert)
                    Ert[index] = zeros(2)
                end
                @fastmath @inbounds for i = 1:bins
                    Abin = ADatas[i, 2:end] .|> (X -> X - ADatas[i, 1].^2)
                    Abin = Abin .|> (X -> X./Abin[1])

                    for j = 1:timeNum
                        Ert[j] += (Avt[j] - Abin[j]).^2
                    end
                end
                @fastmath @inbounds Ert = Ert .|> (X -> X./(bins*(bins - 1))) .|> (X -> sqrt.(X))

                timeFlag1 = 0
                for i = 1:timeNum
                    if Ert[i][1] < Avt[i][1]
                        Aints[lCount][1] += Avt[i][1]
                    else
                        timeFlag1 = i
                        break
                    end
                end
                timeFlag2 = 0
                for i = 1:timeNum
                    if Ert[i][2] < Avt[i][2]
                        Aints[lCount][2] += Avt[i][2]
                    else
                        timeFlag2 = i
                        break
                    end
                end
                Aints[lCount] .-= 0.5

                for i = 1:bins
                    Temp = zeros(2)

                    Abin = ADatas[i, 2:end] .|> (X -> X - ADatas[i, 1].^2)
                    Abin = Abin .|> (X -> X./Abin[1])

                    Temp[1] = Abin[1:timeFlag1] .|> first |> sum
                    Temp[2] = Abin[1:timeFlag2] .|> last |> sum
                    Temp .-= 0.5

                    @. Eints[lCount] += (Aints[lCount] - Temp)^2
                end
                @. Eints[lCount] = Eints[lCount]/(bins*(bins - 1)) |> sqrt

                Avt = Avt .|> (X -> tFreq*X)
                Ert = Ert .|> (X -> tFreq*X)
                Aints *= tFreq
                Eints *= tFreq

                if !isnothing(rightLimit)
                    rightLimit -= 1
                    Avt = Avt[1:rightLimit]
                    Ert = Ert[1:rightLimit]
                else
                    rightLimit = timeNum
                end
                plot!(
                    lPlot,
                    Avt[1:rightLimit] .|> last,
                    label = "L = $(l)",
                    yerror = Ert[1:rightLimit] .|> last 
                )
            end

            savefig(lPlot, "Auto/auto-"*Config[1, fnCount]*"-$(T).png")

            plot!(
                TPlot,
                LPows .|> (lPow -> 2^lPow),
                Aints .|> last,
                label = "T = $(T)" |> LaTeXString,
                yerror = Eints .|> last
            )
        end

        savefig(TPlot, "Int/int-"*Config[1, fnCount]*".png")
    end
end