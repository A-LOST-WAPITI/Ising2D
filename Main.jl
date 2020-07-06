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


@doc md"""
    _ΔE(Status::Array{Int64}, xIndex::Int64, yIndex::Int64)

    用于计算给定所有点状态`Status`时，
    位于(`xIndex`,`yIndex`)的点发生反转时所引起的能量变化
""" ->
function _ΔE(Status::Array{Int64}, xIndex::Int64, yIndex::Int64)
    n = size(Status)[1]
    state = Status[xIndex, yIndex]  # 所选点的状态

    XNei, YNei = _Nei(xIndex, yIndex, n)

    # 计算原能量
    rawE = 0.0
    for i in XNei
        rawE += J₁ * state * Status[i, yIndex]
    end
    for j in YNei
        rawE += J₂ * state * Status[xIndex, j]
    end

    # 反转带来的能量差
    ΔE = 2rawE
end


@doc md"""
    _TheoryF(T::Float64, σ::Float64)

    用于确定温度为`T`时平均自旋`σ`的平均场近似理论解
""" ->
function _TheoryF(T::Float64, σ::Float64)
    tanh((2J₁ + 2J₂)σ/T) - σ
end


@doc md"""
    _Findσ(
        T::Float64, 
        TMax::Int64, 
        iterMax::Int64 = 1000, 
        eps::Float64 = 1e-8
    )

    使用二分法确定$_TheoryF = 0$的解来确定
    平均场近似下$\sigma$的值，
    其中`T`为给定温度，
    `TMax`为求解的最大温度，
    `iterMax`为迭代的最大步数，
    `eps`为允许误差最大值
""" ->
function _Findσ(
    T::Float64, 
    TMax::Int64, 
    iterMax::Int64 = 1000, 
    eps::Float64 = 1e-8
)
    upper = Float64(TMax)
    lower = 0.0
    mid = (lower + upper)/2

    for iter = 1:iterMax
        temp = _TheoryF(T, mid)
        if abs(temp) < eps
            return mid
        else
            temp > 0 && (lower = mid; true)
            temp < 0 && (upper = mid; true)
            mid = (lower + upper)/2
        end
    end

    return mid
end


@doc md"""
    _Main(
        nPowMax::Int64 = 6, 
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
    timeScale::Int64 = 1000, 
    TMax::Int64 = 5
)
    TS = collect(0.05:0.05:TMax)    # 温度序列
    σAS = similar(TS)  # 存储温度序列对应平均自旋的数组

    plot(size = (1200, 800))
    # 不同格点
    for nPow = 3:nPowMax
        n = 2^nPow  # 单方向格点数
        
        # 生成每次实验在不同温度下通用的初始状态
        RawStatus::Array{Int64} = rand((-1, 1), (n, n))
        println("\nThe number of nodes is $(n)")
        println("*"^5)
        
        for (i, T) in enumerate(TS)
            Status = RawStatus  # 赋初值
            Ratios = Dict()
            σASTemp = Array{Float64}([])
            Acorr = Array{Float64}([])
            σA::Float64 = 0
            coin = 1 - exp(-(J₁ + J₂)/T)

            # 时间序列循环
            for time = 1:timeScale
                StatusCheck = falses(n, n)
                Indices = Array{CartesianIndex}(
                    [CartesianIndex(rand(1:n), rand(1:n))]
                )
                

                while sum(StatusCheck) != n^2
                    Cluster = Array{CartesianIndex}([])

                    while lastindex(Indices) > 0
                        index = popfirst!(Indices)
                        state = Status[index]
                        StatusCheck[index] = true
                        push!(Cluster, index)

                        xIndex = index[1]
                        yIndex = index[2]

                        XNei, YNei = _Nei(xIndex, yIndex, n)

                        for i in XNei
                            indexTemp = CartesianIndex(i, yIndex)
                            if Status[indexTemp] == state && !StatusCheck[indexTemp]
                                if rand() < coin
                                    push!(Indices, indexTemp)
                                    StatusCheck[indexTemp] = true
                                    push!(Cluster, indexTemp)
                                end
                            end
                        end
                        for j in YNei
                            indexTemp = CartesianIndex(xIndex, j)
                            if Status[indexTemp] == state && !StatusCheck[indexTemp]
                                if rand() < coin
                                    push!(Indices, indexTemp)
                                    StatusCheck[indexTemp] = true
                                    push!(Cluster, indexTemp)
                                end
                            end
                        end
                    end

                    if rand() < 0.5
                        for index in Cluster
                            Status[index] *= -1
                        end
                    end

                    index = findfirst((x -> !x), StatusCheck)
                    if !isnothing(index)
                        push!(Indices, index)
                    end
                end

                σA = sum(Status)/n^2 |> abs # 求平均后绝对值
                pushfirst!(σASTemp, σA)
                push!(Acorr, 0)
                for j in eachindex(σASTemp)
                    Acorr[j] += σA * σASTemp[j]
                end
                
                if last(Acorr) == 0
                    break
                end
            end

            println("Temperature\t", T, "\tσAverage\t", σA)
            σAS[i] = σA    # 存储遍历数据
        end

        # 绘制不同实验的图像
        plot!(
            TS, 
            σAS, 
            label = "$(n)-nodes"
        )
    end
    xlabel!("kT/J") # x轴名称
    ylabel!("⟨σ⟩")  # y轴名称
    title!("σAverage-nodes")  # 图名
    savefig("sigmaAverage-nodes.png")

    solutionIndex = findfirst((x -> x < 0.1), σAS) # 找到第一个平均自旋降至最大值的0.1的温度
    #=
    TODO:
    这里通过图形来看确实在2-3之间便已经出现了相变，但0.1是否过大？
    在2-3之间仍有较大波动，
    如果将解的条件限制到0.01则会到达4左右与理论值差距减小，
    但又与整个图像的走势不同。
    =#
    println("\n", "*"^10)
    println("Solution is:\t", TS[solutionIndex])    # 将解打印
    println()

    σAST = TS .|> (x -> _Findσ(x, TMax))    # 使用二分法得到的平均场近似理论解存于σAST中

    plot(size = (1200, 800))
    # 最多节点数值解绘图
    plot!(
        TS, 
        σAS .|> abs, 
        yerror = σAS .- σAST,
        markerstrokecolor = :black,
        label = "computional"
    )
    # 理论解绘图
    plot!(
        TS, 
        σAST,
        label = "theoretical"
    )
    xlabel!("kT/J") # x轴名称
    ylabel!("⟨σ⟩")  # y轴名称
    title!("σAverage")  # 图名
    savefig("sigmaAverage.png")

    "Everything goes fine."
end