using MatrixNetworks
using Combinatorics

immutable SymTensor3
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    V::Vector{Float64}
    dimension::Int64
    order::Int64
end

immutable SymTensor4
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    I4::Vector{Int64}    
    V::Vector{Float64}
    dimension::Int64
    order::Int64
end

immutable SymTensor5
    I1::Vector{Int64}
    I2::Vector{Int64}
    I3::Vector{Int64}
    I4::Vector{Int64}
    I5::Vector{Int64}
    V::Vector{Float64}
    dimension::Int64
    order::Int64
end
const SymTensor = Union{SymTensor3, SymTensor4, SymTensor5}
const SpFltMat = SparseMatrixCSC{Float64,Int64}

function apply(T::SymTensor3, x::Vector{Float64})
    dim = T.dimension
    if length(x) != dim
        throw("Vector dimension does not match tensor dimension")
    end
    y = zeros(Float64, dim)
    for (i1, i2, i3, v) in zip(T.I1, T.I2, T.I3, T.V)
        x1, x2, x3 = x[i1], x[i2], x[i3]
        y[i1] += v * x2 * x3
        y[i2] += v * x1 * x3
        y[i3] += v * x1 * x2
    end
    return y * 2
end

function apply(T::SymTensor4, x::Vector{Float64})
    dim = T.dimension
    if length(x) != dim
        throw("Vector dimension does not match tensor dimension")
    end
    y = zeros(Float64, dim)
    for (i1, i2, i3, i4, v) in zip(T.I1, T.I2, T.I3, T.I4, T.V)
        x1, x2, x3, x4 = x[i1], x[i2], x[i3], x[i4]
        y[i1] += v * x2 * x3 * x4
        y[i2] += v * x1 * x3 * x4
        y[i3] += v * x1 * x2 * x4
        y[i4] += v * x1 * x2 * x3        
    end
    return y * 6
end

function apply(T::SymTensor5, x::Vector{Float64})
    dim = T.dimension
    if length(x) != dim
        throw("Vector dimension does not match tensor dimension")
    end
    y = zeros(Float64, dim)
    for (i1, i2, i3, i4, i5, v) in zip(T.I1, T.I2, T.I3, T.I4, T.I5, T.V)
        x1, x2, x3, x4, x5 = x[i1], x[i2], x[i3], x[i4], x[i5]
        y[i1] += v * x2 * x3 * x4 * x5
        y[i2] += v * x1 * x3 * x4 * x5
        y[i3] += v * x1 * x2 * x4 * x5
        y[i4] += v * x1 * x2 * x3 * x5
        y[i5] += v * x1 * x2 * x3 * x4        
    end
    return y * 24
end

function reduce(T::SymTensor3, x::Vector{Float64})
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, v) in zip(T.I1, T.I2, T.I3, T.V)
        vx1, vx2, vx3 = v * [x[i1], x[i2], x[i3]]
        push!(I, i1,  i1,  i2)
        push!(J, i2,  i3,  i3)
        push!(V, vx3, vx2, vx1)
    end
    n = T.dimension
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    return A + A'
end

function reduce(T::SymTensor4, x::Vector{Float64})
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, v) in zip(T.I1, T.I2, T.I3, T.I4, T.V)
        vx12, vx13, vx14 = (v * x[i1]) * [x[i2], x[i3], x[i4]]
        vx23, vx24       = (v * x[i2]) * [x[i3], x[i4]]
        vx34             = (v * x[i3]) * x[i4]
        push!(I, i1,   i1,   i1,   i2,   i2,   i3)
        push!(J, i2,   i3,   i4,   i3,   i4,   i4)
        push!(V, vx34, vx24, vx23, vx14, vx23, vx12)
    end
    n = T.dimension
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    return A + A'
end

function reduce(T::SymTensor5, x::Vector{Float64})
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, i5, v) in zip(T.I1, T.I2, T.I3, T.I4, T.I5, T.V)
        x1, x2, x3, x4, x5 = x[i1], x[i2], x[i3], x[i4], x[i5]
        vx123, vx124, vx125 = (v * x1 * x2) * [x3, x4, x5]
        vx134, vx135        = (v * x1 * x3) * [x4, x5]
        vx145               = (v * x1 * x4) * x5
        vx234, vx235        = (v * x2 * x3) * [x4, x5]
        vx245               = (v * x2 * x4) * x5
        vx345               = (v * x3 * x4) * x5
        push!(I, i1,    i1,    i1,    i1,    i2,    i2,    i2,    i3,    i3,    i4)
        push!(J, i2,    i3,    i4,    i5,    i3,    i4,    i5,    i4,    i5,    i5)
        push!(V, vx345, vx245, vx235, vx234, vx145, vx135, vx134, vx125, vx124, vx123)
    end
    n = T.dimension
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    return A + A'
end

function squeeze_tensor(T::SymTensor3)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, v) in zip(T.I1, T.I2, T.I3, T.V)
        push!(I, i1, i1, i2)
        push!(J, i2, i3, i3)
        push!(V,  v,  v,  v)
    end
    n = T.dimension    
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    A + A'
end

function squeeze_tensor(T::SymTensor4)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, v) in zip(T.I1, T.I2, T.I3, T.I4, T.V)
        push!(I, i1, i1, i1, i2, i2, i3)
        push!(J, i2, i3, i4, i3, i4, i4)
        push!(V,  v,  v,  v,  v,  v,  v)
    end
    n = T.dimension    
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    return A + A'
end

function squeeze_tensor(T::SymTensor5)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, i5, v) in zip(T.I1, T.I2, T.I3, T.I4, T.I5, T.V)
        push!(I, i1, i1, i1, i1, i2, i2, i2, i3, i3, i4)
        push!(J, i2, i3, i4, i5, i3, i4, i5, i4, i5, i5)
        push!(V,  v,  v,  v,  v,  v,  v,  v,  v,  v,  v)
    end
    n = T.dimension    
    A = convert(SpFltMat, sparse(I, J, V, n, n))
    return A + A'
end

function largest_component(T::SymTensor3)
    lcc = MatrixNetworks.largest_component(squeeze_tensor(T))[2]
    ind_map = zeros(Int64, length(lcc))
    new_dim = sum(lcc)
    ind_map[find(lcc)] = collect(1:new_dim)
    I1, I2, I3 = Int64[], Int64[], Int64[]
    V = Float64[]
    for (i1, i2, i3, v) in zip(T.I1, T.I2, T.I3, T.V)
        if all(lcc[[i1, i2, i3]])
            push!(I1, ind_map[i1])
            push!(I2, ind_map[i2])
            push!(I3, ind_map[i3])
            push!(V, v)
        end
    end
    return (SymTensor3(I1, I2, I3, V, new_dim, 3), lcc, ind_map)
end

function largest_component(T::SymTensor4)
    lcc = MatrixNetworks.largest_component(squeeze_tensor(T))[2]
    ind_map = zeros(Int64, length(lcc))
    new_dim = sum(lcc)
    ind_map[find(lcc)] = collect(1:new_dim)
    I1, I2, I3, I4 = Int64[], Int64[], Int64[], Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, v) in zip(T.I1, T.I2, T.I3, T.I4, T.V)    
        if all(lcc[[i1, i2, i3, i4]])
            push!(I1, ind_map[i1])
            push!(I2, ind_map[i2])
            push!(I3, ind_map[i3])
            push!(I4, ind_map[i4])
            push!(V, v)
        end
    end
    return (SymTensor4(I1, I2, I3, I4, V, new_dim, 4), lcc, ind_map)
end

function largest_component(T::SymTensor5)
    lcc = MatrixNetworks.largest_component(squeeze_tensor(T))[2]
    ind_map = zeros(Int64, length(lcc))
    new_dim = sum(lcc)
    ind_map[find(lcc)] = collect(1:new_dim)
    I1, I2, I3, I4, I5 = Int64[], Int64[], Int64[], Int64[], Int64[]
    V = Float64[]
    for (i1, i2, i3, i4, i5, v) in zip(T.I1, T.I2, T.I3, T.I4, T.I5, T.V)
        if all(lcc[[i1, i2, i3, i4, i5]])
            push!(I1, ind_map[i1])
            push!(I2, ind_map[i2])
            push!(I3, ind_map[i3])
            push!(I4, ind_map[i4])
            push!(I5, ind_map[i5])            
            push!(V, v)
        end
    end
    return (SymTensor5(I1, I2, I3, I4, I5, V, new_dim, 5), lcc, ind_map)
end