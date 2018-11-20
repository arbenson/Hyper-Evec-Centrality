using Combinatorics
using DelimitedFiles


function read_data_unweighted(dataset::String, order::Int64, exact_match::Bool=false)
    if order < 3 || order > 5; error("Only support for order-3,4,5 tensors"); end
    
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    simplices = read("data/$(dataset)/$(dataset)-simplices.txt")
    nverts = read("data/$(dataset)/$(dataset)-nverts.txt")

    curr_ind = 1
    S = Set{NTuple{order,Int64}}()
    for nvert in nverts
        simp = sort(simplices[curr_ind:(curr_ind + nvert - 1)])
        curr_ind += nvert
        if (exact_match && nvert != order); continue; end
        for hedge in combinations(simp, order)
            push!(S, NTuple{order,Int64}(hedge))
        end
    end

    I = Vector{Vector{Int64}}(undef, order)
    for j in 1:order; I[j] = Vector{Int64}(); end
    for hedge in S
        for (j, v) in enumerate(hedge)
            push!(I[j], v)
        end
    end
    dim = maximum([maximum(I[j]) for j in 1:order])
    V = ones(Float64, length(I[1])) / length(I[1])
    if order == 3; return SymTensor3(I[1], I[2], I[3], V, dim); end
    if order == 4; return SymTensor4(I[1], I[2], I[3], I[4], V, dim); end
    if order == 5; return SymTensor5(I[1], I[2], I[3], I[4], I[5], V, dim); end    
end

function read_node_labels(dataset::String)
    labels = String[]
    open("data/$(dataset)/$(dataset)-node-labels.txt") do f
        for line in eachline(f)
            data = split(strip(line))
            push!(labels, join(data[2:end], " "))
        end
    end
    return labels
end
;
