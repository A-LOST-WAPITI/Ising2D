using Markdown      # 引入Markdown支持
using Plots         # 引入绘图包
using LaTeXStrings  # 引入LaTeX支持


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
    end
    σA::Float64 = 0
    mse::Float64 = 0

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
            temp = sum(Status)/n^2 |> abs # 求平均后绝对值
            σASTemp[time] = temp
            σA += temp
        end
    end

    if measure
        σA /= steps
        mse = sum(@. (σASTemp - σA)^2)/steps
    end

    return σA, mse
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
    nPowMax::Int64 = 6, 
    timeScale::Int64 = 20000, 
    measureScale::Int64 = 2000, 
    TMin::Int64 = 1,
    TMax::Int64 = 4,
    TTick::Float64 = 0.05
)
    TS = collect(TMin:TTick:TMax)    # 温度序列
    σAS = similar(TS)  # 存储温度序列对应平均自旋的数组
    MSE = similar(TS)

    plot(size = (1600, 800))
    # 不同格点
    for nPow = 3:nPowMax
        n = 2^nPow  # 单方向格点数
        
        # 生成每次实验在不同温度下通用的初始状态
        # RawStatus::Array{Int64} = rand((-1, 1), (n, n))
        RawStatus = ones(Int64, (n, n))
        println("\nThe number of nodes is $(n)")
        println("*"^20)
        
        # 对于给定的单方向格点数遍历各个不同的温度
        for (i, T) in enumerate(TS)
            Status = RawStatus  # 赋初值
            
            pCluster = 1 - exp(-(J₁ + J₂)/T)

            # 时间序列循环
            _MC!(Status, n, pCluster, timeScale, false)
            σA, mse = _MC!(Status, n, pCluster, measureScale, true)

            println("Temperature\t", T, "\tσAverage\t", σA)
            σAS[i] = σA    # 存储遍历数据
            MSE[i] = mse
        end

        # 绘制不同实验的图像
        plot!(
            TS, 
            σAS, 
            label = "$(n)-nodes", 
            yerror = MSE
        )
    end
    xlabel!("kT/J")     # x轴名称
    ylabel!("⟨σ⟩")       # y轴名称
    title!("σAverage-nodes")  # 图名
    savefig("sigmaAverage-nodes.png")

    #=
    solutionIndex = findfirst(iszero, σAS) # 找到第一个平均自旋降至0的温度
    println("\n", "*"^20)
    println("Solution is:\t", TS[solutionIndex])    # 将解打印
    println()
    =#

    "Everything goes fine."
end