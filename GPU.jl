using CUDA


function _Stablize!(n, Right, Below, Relation)
    xIndex = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    yIndex = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    
    xIndex == n && (xNext = 1; true) || (xNext = xIndex + 1; true)
    yIndex == n && (yNext = 1; true) || (yNext = yIndex + 1; true)

    if Right[xIndex, yIndex] && Below[xIndex, yIndex]
        @inbounds Relation[xIndex, yIndex] = 
            0.1 * Relation[xIndex, yIndex] + 
            0.45 * Relation[xNext, yIndex] + 
            0.45 * Relation[xIndex, yNext]
    elseif Right[xIndex, yIndex]
        @inbounds Relation[xIndex, yIndex] = 
            0.1 * Relation[xIndex, yIndex] + 
            0.9 * Relation[xNext, yIndex]
    elseif Below[xIndex, yIndex]
        @inbounds Relation[xIndex, yIndex] = 
            0.1 * Relation[xIndex, yIndex] + 
            0.9 * Relation[xIndex, yNext]
    end

    return nothing
end


function _Main()
    n = 2^3
    pCluster = 0.5f0
    Forward = [2:n |> collect; 1]
    Backward = [n; 1:n - 1 |> collect]

    Status = ones(Int32, (n, n)) |> cu

    RightSign = (Status .== Status[:, Forward]) # 与右侧最近临格点是否自旋相同
    BelowSign = (Status .== Status[Forward, :]) # 与下方最近临格点是否自旋相同
    RightChances = (CUDA.rand(Float32, (n, n)) .< pCluster)
    BelowChances = (CUDA.rand(Float32, (n, n)) .< pCluster)
    Right = (RightSign .& RightChances)
    Below = (BelowSign .& BelowChances)

    Relation = CUDA.rand(Float32, (n, n))
    for i = 1:1000
        @cuda threads = (n, n) _Stablize!(n, Right, Below, Relation)
    end

    return Relation |> Array
end