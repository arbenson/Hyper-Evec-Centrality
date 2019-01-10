include("centrality.jl")
using Combinatorics
using Random
using Test

function test3()
    dim = 11
    Random.seed!(1234)
    
    entries = collect(combinations(1:dim, 3))
    vals = rand(length(entries))
    T = SymTensor3([e[1] for e in entries],
                   [e[2] for e in entries],
                   [e[3] for e in entries],
                   vals)
    entry_vals = Dict{NTuple{3,Int64},Float64}()
    for (e, v) in zip(entries, vals)
        entry_vals[NTuple{3, Int64}(e)] = v
    end

    y = zeros(Float64, dim)
    x = rand(dim)
    for i in 1:dim, j in 1:dim, k in 1:dim
        key = NTuple{3,Int64}(sort([i, j, k]))
        if haskey(entry_vals, key)
            Tijk = entry_vals[key]
            y[i] += Tijk * x[j] * x[k]
        end
    end
    @test y ≈ apply(T, x)

    Y = zeros(Float64, dim, dim)
    for i in 1:dim, j in 1:dim, k in 1:dim
        key = NTuple{3,Int64}(sort([i, j, k]))
        if haskey(entry_vals, key)
            Tijk = entry_vals[key]
            Y[i, j] += Tijk * x[k]
        end
    end
    @test Y ≈ Matrix(reduce(T, x))
end

function test4()
    dim = 10
    Random.seed!(1234)
    
    entries = collect(combinations(1:dim, 4))
    vals = rand(length(entries))
    T = SymTensor4([e[1] for e in entries],
                   [e[2] for e in entries],
                   [e[3] for e in entries],
                   [e[4] for e in entries],                   
                   vals)
    entry_vals = Dict{NTuple{4,Int64},Float64}()
    for (e, v) in zip(entries, vals)
        entry_vals[NTuple{4,Int64}(e)] = v
    end

    y = zeros(Float64, dim)
    x = rand(dim)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
        key = NTuple{4,Int64}(sort([i, j, k, l]))
        if haskey(entry_vals, key)
            Tijkl = entry_vals[key]
            y[i] += Tijkl * x[j] * x[k] * x[l]
        end
    end
    @test y ≈ apply(T, x)

    Y = zeros(Float64, dim, dim)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
        key = NTuple{4,Int64}(sort([i, j, k, l]))
        if haskey(entry_vals, key)
            Tijkl = entry_vals[key]
            Y[i, j] += Tijkl * x[k] * x[l]
        end
    end
    @test Y ≈ Matrix(reduce(T, x))
end

function test5()
    dim = 9
    Random.seed!(1234)
    
    entries = collect(combinations(1:dim, 5))
    vals = rand(length(entries))
    T = SymTensor5([e[1] for e in entries],
                   [e[2] for e in entries],
                   [e[3] for e in entries],
                   [e[4] for e in entries],
                   [e[5] for e in entries],                   
                   vals)
    entry_vals = Dict{NTuple{5,Int64},Float64}()
    for (e, v) in zip(entries, vals)
        entry_vals[NTuple{5,Int64}(e)] = v
    end

    y = zeros(Float64, dim)
    x = rand(dim)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim, m in 1:dim
        key = NTuple{5,Int64}(sort([i, j, k, l, m]))
        if haskey(entry_vals, key)
            Tijklm = entry_vals[key]
            y[i] += Tijklm * x[j] * x[k] * x[l] * x[m]
        end
    end
    @test y ≈ apply(T, x)

    Y = zeros(Float64, dim, dim)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim, m in 1:dim
        key = NTuple{5,Int64}(sort([i, j, k, l, m]))
        if haskey(entry_vals, key)
            Tijklm = entry_vals[key]
            Y[i, j] += Tijklm * x[k] * x[l] * x[m]
        end
    end
    @test Y ≈ Matrix(reduce(T, x))
end

function main()
    test3()
    test4()
    test5()
end
;
