using Markdown


const J₁ = 1    # x方向相邻的相互作用
const J₂ = 1    # y方向相邻的相互作用


@doc md"""
    _Nei(xIndex::Int64, yIndex::Int64, n::Int64)

    获取满足周期性条件下的相邻位置,
    `xIndex`为当前横坐标,
    `yIndex`为当前纵坐标,
    `n`为单方向格点个数
""" ->
function _Nei(xIndex::Int64, yIndex::Int64, l::Int64)
    # 对于临近索引赋初值，保证正常遍历过程不会得到
    XNei = fill(l + 1, 2)
    YNei = fill(l + 1, 2)

    # 处理近邻格点位置
    @inbounds xIndex == 1 && (XNei[1] = l; true) || (XNei[1] = xIndex - 1; true)
    @inbounds xIndex == l && (XNei[2] = 1; true) || (XNei[2] = xIndex + 1; true)
    @inbounds yIndex == 1 && (YNei[1] = l; true) || (YNei[1] = yIndex - 1; true)
    @inbounds yIndex == l && (YNei[2] = 1; true) || (YNei[2] = yIndex + 1; true)

    return XNei, YNei
end


function _MCMetropolis!(Status::Array{Int64, 2}, l::Int64, n::Int64, T::Float64)
    for i = 1:n
        xIndex = rand(1:l)
        yIndex = rand(1:l)
        
        XNei, YNei = _Nei(xIndex, yIndex, l)
        xPrev, xNext = XNei
        yPrev, yNext = YNei


        @inbounds ΔE::Float64 = 
            Status[xIndex, yPrev] +
            Status[xIndex, yNext] +
            Status[xPrev, yIndex] +
            Status[xNext, yIndex]
        ΔE *= (J₁ + J₂)/T * Status[xIndex, yIndex]


        if ΔE < 0 || rand() < exp(- ΔE)
            Status[xIndex, yIndex] = - Status[xIndex, yIndex]
        end
    end


    return nothing
end


function _MCWolff!(Status::Array{Int64, 2}, l::Int64, n::Int64, T::Float64)
    pCluster = 1 - exp(-(J₁ + J₂)/T)

    StatusCheck = falses(l, l)

    index = CartesianIndex(rand(1:l), rand(1:l))
    state = Status[index]
    centerCount = 1
    clusterCount = 1

    # 需要检测位置的队列
    CenterIndices = Array{CartesianIndex{2}}(undef, n)
    CenterIndices[clusterCount] = index
    # 用于存储本次生成的cluster的格点位置的数组
    Cluster = Array{CartesianIndex{2}}(undef, n)
    Cluster[clusterCount] = index
    StatusCheck[index] = true

    while centerCount > 0
        index = CenterIndices[1]
        centerCount -= 1
        CenterIndices[1:centerCount] .= @view CenterIndices[2:centerCount + 1]

        xIndex = index[1]
        yIndex = index[2]
        XNei, YNei = _Nei(xIndex, yIndex, l)

        for i in XNei
            indexTemp = CartesianIndex(i, yIndex)
            if !StatusCheck[indexTemp] && Status[indexTemp] == state
                if rand() < pCluster
                    centerCount += 1
                    clusterCount += 1

                    CenterIndices[centerCount] = indexTemp
                    Cluster[clusterCount] = indexTemp

                    StatusCheck[indexTemp] = true
                end
            end
        end
        for j in YNei
            indexTemp = CartesianIndex(xIndex, j)
            if !StatusCheck[indexTemp] && Status[indexTemp] == state
                if rand() < pCluster
                    centerCount += 1
                    clusterCount += 1

                    CenterIndices[centerCount] = indexTemp
                    Cluster[clusterCount] = indexTemp

                    StatusCheck[indexTemp] = true
                end
            end
        end
    end

    for i = 1:clusterCount
        index = Cluster[i]
        Status[index] *= -1
    end

    return nothing
end


function _MCWS!(Status::Array{Int64, 2}, l::Int64, n::Int64, T::Float64)
    pCluster = 1 - exp(-(J₁ + J₂)/T)

    StatusCheck = falses(l, l)

    while sum(StatusCheck) < n
        index = findfirst(!, StatusCheck)
        state = Status[index]
        centerCount = 1
        clusterCount = 1

        # 需要检测位置的队列
        CenterIndices = Array{CartesianIndex{2}}(undef, n)
        CenterIndices[clusterCount] = index
        # 用于存储本次生成的cluster的格点位置的数组
        Cluster = Array{CartesianIndex{2}}(undef, n)
        Cluster[clusterCount] = index
        StatusCheck[index] = true

        while centerCount > 0
            index = CenterIndices[1]
            centerCount -= 1
            CenterIndices[1:centerCount] .= @view CenterIndices[2:centerCount + 1]

            xIndex = index[1]
            yIndex = index[2]
            XNei, YNei = _Nei(xIndex, yIndex, l)

            for i in XNei
                indexTemp = CartesianIndex(i, yIndex)
                if !StatusCheck[indexTemp] && Status[indexTemp] == state
                    if rand() < pCluster
                        centerCount += 1
                        clusterCount += 1

                        CenterIndices[centerCount] = indexTemp
                        Cluster[clusterCount] = indexTemp

                        StatusCheck[indexTemp] = true
                    end
                end
            end
            for j in YNei
                indexTemp = CartesianIndex(xIndex, j)
                if !StatusCheck[indexTemp] && Status[indexTemp] == state
                    if rand() < pCluster
                        centerCount += 1
                        clusterCount += 1

                        CenterIndices[centerCount] = indexTemp
                        Cluster[clusterCount] = indexTemp

                        StatusCheck[indexTemp] = true
                    end
                end
            end
        end

        if rand() < 0.5
            for i = 1:clusterCount
                index = Cluster[i]
                Status[index] *= -1
            end
        end
    end

    return nothing
end