using CUDA
using Printf


CUDA.allowscalar(false)


const J₁ = 1f0    # x方向相邻的相互作用
const J₂ = 1f0    # y方向相邻的相互作用


function _Stablize!(n, Right, Below, Relation, NewRelation)
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    
    xIndex == n && (xNext = 1; true) || (xNext = xIndex + 1; true)
    xIndex == 1 && (xPre = n; true) || (xPre = xIndex - 1; true)
    yIndex == n && (yNext = 1; true) || (yNext = yIndex + 1; true)
    yIndex == 1 && (yPre = n; true) || (yPre = yIndex - 1; true)

    temp = 0f0
    numCount = 0f0
    if Right[xIndex, yIndex] 
        temp += Relation[xNext, yIndex]
        numCount += 1f0
    end
    if Right[xPre, yIndex]
        temp += Relation[xPre, yIndex]
        numCount += 1f0
    end
    if Below[xIndex, yIndex]
        temp += Relation[xIndex, yNext]
        numCount += 1f0
    end
    if Below[xIndex, yPre]
        temp += Relation[xIndex, yPre]
        numCount += 1f0
    end

    if numCount == 0f0
        NewRelation[xIndex, yIndex] = Relation[xIndex, yIndex]
    else
        NewRelation[xIndex, yIndex] = 0.9 * temp / numCount
        NewRelation[xIndex, yIndex] += 0.1 * Relation[xIndex, yIndex]
    end

    return nothing
end


function _Flip!(n, Right, Below, Status, Choices)
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    xIndex == 1 && (xPre = n; true) || (xPre = xIndex - 1; true)
    yIndex == 1 && (yPre = n; true) || (yPre = yIndex - 1; true)

    if Choices[xIndex, yIndex] 
        Status[xIndex, yIndex] = -Status[xIndex, yIndex]
    elseif Choices[xPre, yIndex] && Right[xPre, yIndex]
        Status[xIndex, yIndex] = -Status[xIndex, yIndex]
    elseif Choices[xIndex, yPre] && Below[xIndex, yPre]
        Status[xIndex, yIndex] = -Status[xIndex, yIndex]
    end

    return nothing
end


function _Main(
    TMin::Float32 = 1.5f0,
    TMax::Float32 = 3.5f0,
    TTicks::Int64 = 50,
    timeScale::Int64 = 10000
)
    TS = LinRange(TMin, TMax, TTicks) |> collect
    σAS = similar(TS)
    n = 2^4
    thread = min(32, n)
    threads = (thread, thread)
    block = ceil(Int, n/thread)
    blocks = (block, block)
    Forward = [2:n |> collect; 1]

    for (i, T) in enumerate(TS)
        σASTemp = zeros(Float32, timeScale + 1)
        Acorr = zeros(Float32, timeScale + 1)

        pCluster::Float32 = 1 - exp(-(J₁ + J₂)/T)
        σA::Float32 = 0

        Status = rand((-1, 1) .|> Int32, (n, n)) |> cu
        for time = 1:timeScale
            RightSign = (Status .== Status[:, Forward]) # 与右侧最近临格点是否自旋相同
            BelowSign = (Status .== Status[Forward, :]) # 与下方最近临格点是否自旋相同
            RightChances = (CUDA.rand(Float32, (n, n)) .< pCluster)
            BelowChances = (CUDA.rand(Float32, (n, n)) .< pCluster)
            Right = (RightSign .& RightChances)
            Below = (BelowSign .& BelowChances)

            Relation = CUDA.rand(Float32, (n, n))
            NewRelation = similar(Relation)
            for i = 1:1000
                @cuda blocks = blocks threads = threads _Stablize!(n, Right, Below, Relation, NewRelation)
                if NewRelation == Relation
                    # println(i)
                    break
                else
                    Relation .= NewRelation
                end
            end

            Choices = (Relation .< 0.5f0)

            @cuda blocks = blocks threads = threads _Flip!(n, Right, Below, Status, Choices)
            # display(Choices)
            # println(sum(Choices), "\t", sum(Choices)/n^2)
            # display(Status)
            # println(sum(Status))
            # readline(stdin)

            σA = sum(Status)/n^2 |> abs
            σASTemp[2:time + 1] .= σASTemp[1:time]
            σASTemp[1] = σA
            for j = 1:time + 1
                Acorr[j] += σA * σASTemp[j]
            end
            
            if Acorr[time] |> iszero
                # println(time)
                break
            end
        end

        @printf("Temperature\t %.5f \tσAverage\t %10.f\n", T, σA)

        σAS[i] = σA
    end
end