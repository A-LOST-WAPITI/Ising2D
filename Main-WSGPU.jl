using Markdown      # 引入Markdown支持
using Plots         # 引入绘图包
using Printf
using CUDA


const J₁ = 1f0    # x方向相邻的相互作用
const J₂ = 1f0    # y方向相邻的相互作用


function _Init!(n, Flags)
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    Flags[xIndex, yIndex] = (xIndex - 1) * n + yIndex
    return nothing
end


function _Cluster!(n, Flags, Right, Below, NewFlags)
    # 获取行索引
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    # 获取列索引
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    xIndex == n && (xNext = 1; true) || (xNext = xIndex + 1; true)
    xIndex == 1 && (xPrev = n; true) || (xPrev = xIndex - 1; true)
    yIndex == n && (yNext = 1; true) || (yNext = yIndex + 1; true)
    yIndex == 1 && (yPrev = n; true) || (yPrev = yIndex - 1; true)

    tempFlag = Flags[xIndex, yIndex]
    if Right[xIndex, yIndex]
        tempFlag = CUDA.min(tempFlag, Flags[xIndex, yNext])
    end
    if Right[xIndex, yPrev]
        tempFlag = CUDA.min(tempFlag, Flags[xIndex, yPrev])
    end
    if Below[xIndex, yIndex]
        tempFlag = CUDA.min(tempFlag, Flags[xNext, yIndex])
    end
    if Below[xPrev, yIndex]
        tempFlag = CUDA.min(tempFlag, Flags[xPrev, yIndex])
    end

    NewFlags[xIndex, yIndex] = tempFlag

    return nothing
end


function _Flip!(Flags, Probability, Status)
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if Probability[Flags[xIndex, yIndex]] < 0.5f0
        Status[xIndex, yIndex] = -Status[xIndex, yIndex]
    end

    return nothing
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
    TMin::Float32 = 3.0f0,
    TMax::Float32 = 3.5f0,
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
        Forward = [2:n |> collect; 1]
        thread = min(32, n)
        threads = (thread, thread)
        block = ceil(Int, n/thread)
        blocks = (block, block)
        
        println("\nThe number of nodes is $(n)")
        println("*"^20)
        
        Flags = CUDA.zeros(Int32, (n, n))
        NewFlags = similar(Flags)
        # 对于给定的单方向格点数并行遍历各个不同的温度
        for (i, T) in enumerate(TS)
            Status = CUDA.ones(Int32, (n, n)) # 赋初值

            σA::Float32 = 0f0 # 平均磁矩
            pCluster = 1 - exp(-(J₁ + J₂)/T)

            # 时间序列循环
            σASTemp = zeros(timeScale + 1)
            Acorr = zeros(timeScale + 1)
            for time = 1:timeScale
                @cuda blocks = blocks threads = threads _Init!(n, Flags)
                Probability = CUDA.rand(Float32, n^2)

                # 确定当前格点是否与右侧最近临格点相连
                Right = ((Status .== Status[:, Forward]) .& (CUDA.rand(Float32, (n, n)) .< pCluster))
                # 确定当前格点是否与下方最近临格点相连
                Below = ((Status .== Status[Forward, :]) .& (CUDA.rand(Float32, (n, n)) .< pCluster))

                for i = 1:n^2
                    @cuda blocks = blocks threads = threads _Cluster!(n, Flags, Right, Below, NewFlags)

                    if Flags == NewFlags
                        break
                    end
                    Flags .= NewFlags
                end
                @cuda blocks = blocks threads = threads _Flip!(Flags, Probability, Status)

                σA = sum(Status)/n^2 |> abs # 求平均后绝对值
                σASTemp[2:time + 1] .= σASTemp[1:time]
                σASTemp[1] = σA
                # 计算Acorr用于提前结束已稳定的状态
                for j = 1:time
                    Acorr[j] += σA * σASTemp[j]
                end
                
                if Acorr[time] |> iszero
                    break
                end
            end

            σASTemp = zeros(Float32, measureScale)
            Σ2ASTemp = zeros(Float32, measureScale)
            for time = 1:measureScale
                @cuda blocks = blocks threads = threads _Init!(n, Flags)
                Probability = CUDA.rand(Float32, n^2)

                # 确定当前格点是否与右侧最近临格点相连
                Right = ((Status .== Status[:, Forward]) .& (CUDA.rand(Float32, (n, n)) .< pCluster))
                # 确定当前格点是否与下方最近临格点相连
                Below = ((Status .== Status[Forward, :]) .& (CUDA.rand(Float32, (n, n)) .< pCluster))

                for i = 1:n
                    @cuda blocks = blocks threads = threads _Cluster!(n, Flags, Right, Below, NewFlags)
                
                    Flags .= NewFlags
                end
                @cuda blocks = blocks threads = threads _Flip!(Flags, Probability, Status)
                # display(Flags)
                # display(Status)
                # readline(stdin)

                Σ = Status |> sum |> abs
                # println(Σ)
                σ = Σ/n^2 
                Σ2 = Σ^2

                σASTemp[time] = σ
                Σ2ASTemp[time] = Σ2
            end

            # display(σASTemp)
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
    savefig(σPlot, "sigmaAverage-nodes-gpu.png")
    savefig(χPlot, "chiAverage-nodes-gpu.png")


    tempT = TS[findmax(χAS)[2]] # 寻找磁化率最高温度作为相变临界温度打印
    println("\n", "*"^20)
    println("Solution is:\t", tempT)    # 将解打印
    println()

    "Everything goes fine."
end