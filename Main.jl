using Markdown      # 引入Markdown支持
using Plots         # 引入绘图包
using Printf


const J₁ = 1    # x方向相邻的相互作用
const J₂ = 1    # y方向相邻的相互作用


@doc md"""
    _Nei(xIndex::Int64, yIndex::Int64, n::Int64)

    获取满足周期性条件下的相邻位置,
    `xIndex`为当前横坐标,
    `yIndex`为当前纵坐标,
    `n`为单方向格点个数
""" ->
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


@doc md"""
    _MC!(
        Status::Array{Int64}, 
        n::Int64, 
        pCluster::Float64, 
        steps::Int64, 
        measure::Bool
    )

    使用Wolff Cluster进行蒙特卡洛模拟,
    `Status`存储当前状态的矩阵,
    `pCluster`当前温度下连入Cluster的概率,
    `steps`最大进行模拟的步数,
    `measure`是否进行统计测量
""" ->
function _MC!(
    Status::Array{Int64}, 
    n::Int64, 
    pCluster::Float64, 
    steps::Int64, 
    measure::Bool
)
    if !measure
        σASTemp = zeros(steps + 1)
        Σ2ASTemp = zeros(steps + 1)
        Acorr = zeros(steps + 1)
    else
        σASTemp = zeros(steps)
        Σ2ASTemp = zeros(steps)
    end
    σA::Float64 = 0 # 平均磁矩

    # 时间序列循环
    for time = 1:steps
        index = CartesianIndex(rand(1:n), rand(1:n))
        state = Status[index]

        StatusCheck = falses(n, n)

        # 需要检测位置的队列
        CenterIndices = Channel{CartesianIndex{2}}(n^2)
        put!(CenterIndices, index)
        # 用于存储本次生成的cluster的格点位置的数组
        Cluster = Channel{CartesianIndex{2}}(n^2)
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
            # 计算Acorr用于提前结束已稳定的状态
            for j = 1:time + 1
                Acorr[j] += σA * σASTemp[j]
            end
            
            if Acorr[time] |> iszero
                break
            end
        else
            temp = sum(Status) |> abs   # 计算总自旋
            Σ2ASTemp[time] = temp^2
            σASTemp[time] = temp/n^2
        end
    end

    return σASTemp, Σ2ASTemp
end


@doc md"""
    _Main(
        nPowMax::Int64 = 8, 
        timeScale::Int64 = 10000, 
        measureScale::Int64 = 50000, 
        TMin::Float64 = 1.5,
        TMax::Float64 = 3.5,
        nTicks::Int64 = 50
    )

    主函数，
    其中`nPowMax`用于迭代产生单方向节点个数,
    `timeScale`为弛豫步长,
    `measureScale`为用于统计的步长,
    `TMax`是温度下限,
    `TMax`是温度上限,
    `nTicks`为等间距温度序列的长度
""" ->
function _Main(
    nPowMax::Int64 = 8, 
    timeScale::Int64 = 10000, 
    measureScale::Int64 = 50000, 
    TMin::Float64 = 1.5,
    TMax::Float64 = 3.5,
    nTicks::Int64 = 50
)
    TS = LinRange(TMin, TMax, nTicks) |> collect    # 温度序列
    σAS = similar(TS)   # 存储温度序列对应平均磁矩的数组
    σMSE = similar(TS)  # 存储温度序列对应平均磁矩的均方误差的数组
    χAS = similar(TS)   # 存储温度序列对应磁化率的数组
    χMSE = similar(TS)  # 存储温度序列对应磁化率的均方误差的数组

    # 生成平均磁矩绘图的画布
    σPlot = plot(
        size = (1600, 800), 
        title = "σAverage-nodes",
        xlabel = "kT/J",
        ylabel = "σ Average"
    )
    # 生成磁化率绘图的画布
    χPlot = plot(
        size = (1600, 800),
        title = "χAverage-nodes",
        xlabel = "kT/J",
        ylabel = "χ Average"
    )
    # 不同格点
    for nPow = 3:nPowMax
        n = 2^nPow  # 单方向格点数
        
        println("\nThe number of nodes is $(n)")
        println("*"^20)
        
        # 对于给定的单方向格点数并行遍历各个不同的温度
        Threads.@threads for (i, T) in enumerate(TS) |> collect
            Status = ones(Int64, (n, n))  # 赋初值
            
            pCluster = 1 - exp(-(J₁ + J₂)/T)

            # 时间序列循环
            _MC!(Status, n, pCluster, timeScale, false)
            σASTemp, Σ2ASTemp = _MC!(Status, n, pCluster, measureScale, true)

            σA = sum(σASTemp)/measureScale
            Σ2A = sum(Σ2ASTemp)/measureScale
            χA = (Σ2A - (σA * n^2)^2)/(n^2 * T)
        
            σMse = sum(@. (σASTemp - σA)^2)/measureScale
            Σ2Mse = sum(@. (Σ2ASTemp - Σ2A)^2)/measureScale
            χMse = sum(@. (Σ2ASTemp - (σASTemp * n^2)^2)/(n^2 * T))/measureScale
            
            @printf(
                "Temperature\t %.5f \tσAverage\t %10.f \tchiAverage\t %.3f \n",
                T, σA, χA
            )
            # 存储遍历数据
            σAS[i] = σA
            σMSE[i] = σMse
            χAS[i] = χA
            χMSE[i] = χMse
        end

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
            yerror = χMSE,
            markershape = :x
        )
    end
    savefig(σPlot, "sigmaAverage-nodes.png")
    savefig(χPlot, "chiAverage-nodes.png")


    tempT = TS[findmax(χAS)[2]] # 寻找磁化率最高温度作为相变临界温度打印
    println("\n", "*"^20)
    println("Solution is:\t", tempT)    # 将解打印
    println()

    "Everything goes fine."
end