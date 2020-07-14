using Markdown      # 引入Markdown支持
using Plots         # 引入绘图包
using Printf


const J₁ = 1    # x方向相邻的相互作用
const J₂ = 1    # y方向相邻的相互作用


function _Nei(xIndex::Int64, yIndex::Int64, n::Int64)
    # 对于临近索引赋初值，保证正常遍历过程不会得到
    XNei = fill(n + 1, 2)
    YNei = fill(n + 1, 2)

    # 处理近邻格点位置
    xIndex == 1 && (XNei[1] = n; true) || (XNei[1] = xIndex - 1; true)
    xIndex == n && (XNei[2] = 1; true) || (XNei[2] = xIndex + 1; true)
    yIndex == 1 && (YNei[1] = n; true) || (YNei[1] = yIndex - 1; true)
    yIndex == n && (YNei[2] = 1; true) || (YNei[2] = yIndex + 1; true)

    return XNei, YNei
end


function _MC!(Status::Array{Int64}, n::Int64,  pCluster::Float64, steps::Int64, measure::Bool)
    if !measure
        σASTemp = zeros(steps + 1)
        Acorr = zeros(steps + 1)
    else
        σASTemp = zeros(steps)
        Σ2ASTemp = zeros(steps)
    end
    σA::Float64 = 0
    Σ2A::Float64 = 0
    σMse::Float64 = 0
    Σ2Mse::Float64 = 0

    # 时间序列循环
    for time = 1:steps
        index = CartesianIndex(rand(1:n), rand(1:n))
        state = Status[index]

        StatusCheck = falses(n, n)

        CenterIndices = Channel{CartesianIndex{2}}(n^2)
        put!(CenterIndices, index)
        Cluster = Channel{CartesianIndex{2}}(n^2) # 用于存储本次生成的cluster的格点位置的数组
        put!(Cluster, index)
        StatusCheck[index] = true

        while isready(CenterIndices)
            index = take!(CenterIndices)

            xIndex = index[1]
            yIndex = index[2]
            XNei, YNei = _Nei(xIndex, yIndex, n)

            for i in XNei
                indexTemp = CartesianIndex(i, yIndex)
                if !StatusCheck[indexTemp] && Status[indexTemp] == state
                    if rand() < pCluster
                        put!(CenterIndices, indexTemp)
                        put!(Cluster, indexTemp)
                        StatusCheck[indexTemp] = true
                    end
                end
            end
            for j in YNei
                indexTemp = CartesianIndex(xIndex, j)
                if !StatusCheck[indexTemp] && Status[indexTemp] == state
                    if rand() < pCluster
                        put!(CenterIndices, indexTemp)
                        put!(Cluster, indexTemp)
                        StatusCheck[indexTemp] = true
                    end
                end
            end
        end
        close(CenterIndices)
        close(Cluster)

        for index in Cluster
            Status[index] *= -1
        end

        if !measure
            σA = sum(Status)/n^2 |> abs # 求平均后绝对值
            σASTemp[2:time + 1] .= σASTemp[1:time]
            σASTemp[1] = σA
            for j = 1:time + 1
                Acorr[j] += σA * σASTemp[j]
            end
            
            if Acorr[time + 1] |> iszero
                break
            end
        else
            temp = sum(Status) |> abs   # 计算总自旋
            Σ2A += temp^2
            Σ2ASTemp[time] = Σ2A
            temp /= n^2
            σA += temp
            σASTemp[time] = temp
        end
    end

    if measure
        σA /= steps
        Σ2A /= steps
        
        σMse = sum(@. (σASTemp - σA)^2)/steps
        Σ2Mse = sum(@. (Σ2ASTemp - Σ2A)^2)/steps
    end

    return σA, σMse, Σ2A, Σ2Mse
end


@doc md"""
    _Main(
        nPowMax::Int64 = 5, 
        timeScale::Int64 = 1000, 
        TMax::Int64 = 5
    )

    主函数，
    其中`nPowMax`用于迭代产生单方向节点个数，
    `timeScale`用于确定最大的迭代步数，
    `TMax`是最大的温度
""" ->
function _Main(
    nPowMax::Int64 = 8, 
    timeScale::Int64 = 10000, 
    measureScale::Int64 = 50000, 
    TMin::Float64 = 1.5,
    TMax::Float64 = 3.5,
    nTicks::Int64 = 50, 
    eps::Float64 = 1e-3
)
    TS = LinRange(TMin, TMax, nTicks) |> collect    # 温度序列
    σAS = similar(TS)  # 存储温度序列对应平均自旋的数组
    σMSE = similar(TS)
    Σ2AS = similar(TS)
    Σ2MSE = similar(TS)
    χAS = zeros(measureScale)

    σPlot = plot(
        size = (1600, 800), 
        title = "σAverage-nodes",
        xlabel = "kT/J",
        ylabel = "σ Average"
    )
    χPlot = plot(
        size = (1600, 800),
        title = "χAverage-nodes",
        xlabel = "kT/J",
        ylabel = "χ Average"
    )
    # 不同格点
    for nPow = 3:nPowMax
        n = 2^nPow  # 单方向格点数
        
        # 生成每次实验在不同温度下通用的初始状态
        # RawStatus::Array{Int64} = rand((-1, 1), (n, n))
        RawStatus = ones(Int64, (n, n))
        println("\nThe number of nodes is $(n)")
        println("*"^20)
        
        # 对于给定的单方向格点数遍历各个不同的温度
        Threads.@threads for (i, T) in enumerate(TS) |> collect
            Status = deepcopy(RawStatus)  # 赋初值
            
            pCluster = 1 - exp(-(J₁ + J₂)/T)

            # 时间序列循环
            _MC!(Status, n, pCluster, timeScale, false)
            σA, σMse, Σ2A, Σ2Mse = _MC!(Status, n, pCluster, measureScale, true)

            @printf(
                "Temperature\t %.5f \tσAverage\t %10.f \tSigma^2Average\t %.3f \n",
                T, σA, Σ2A
            )
            σAS[i] = σA    # 存储遍历数据
            σMSE[i] = σMse
            Σ2AS[i] = Σ2A
            Σ2MSE[i] = Σ2Mse
        end

        χAS = @. (Σ2AS - (σAS * n^2)^2)/(n^2 * TS)

        # 绘制不同实验的图像
        plot!(
            σPlot,
            TS, 
            σAS, 
            label = "$(n)-nodes", 
            yerror = σMSE,
            markershape = :x
        )
        plot!(
            χPlot, 
            TS,
            χAS,
            label = "$(n)-nodes",
            markershape = :x
        )
    end
    savefig(σPlot, "sigmaAverage-nodes.png")
    savefig(χPlot, "chiAverage-nodes.png")


    tempT = TS[findmax(χAS)[2]]
    println("\n", "*"^20)
    println("Solution is:\t", tempT)    # 将解打印
    println()

    "Everything goes fine."
end