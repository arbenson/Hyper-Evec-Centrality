include("tensor.jl")
using MatrixNetworks

function read_data_unweighted(dataset::String, order::Int64, exact_match::Bool=false)
    if order < 3 || order > 5; throw("Only support for order-3,4,5 tensors"); end
    
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    simplices = read("data/$(dataset)/$(dataset)-simplices.txt")
    nverts = read("data/$(dataset)/$(dataset)-nverts.txt")

    curr_ind = 1
    S = Set{NTuple{order,Int64}}()
    for nvert in nverts
        simp = simplices[curr_ind:(curr_ind + nvert - 1)]
        curr_ind += nvert
        if (exact_match && nvert != order); continue ;end
        for hedge in combinations(simp, order)
            push!(S, NTuple{order,Int64}(hedge))
        end
    end

    I = Vector{Vector{Int64}}(order)
    for j in 1:order; I[j] = Vector{Int64}(); end
    for hedge in S
        for (j, v) in enumerate(hedge)
            push!(I[j], v)
        end
    end
    dim = maximum([maximum(I[j]) for j in 1:order])
    V = ones(Float64, length(I[1])) / length(I[1])
    if order == 3; return SymTensor3(I[1], I[2], I[3], V, dim, 3); end
    if order == 4; return SymTensor4(I[1], I[2], I[3], I[4], V, dim, 4); end
    if order == 5; return SymTensor5(I[1], I[2], I[3], I[4], I[5], V, dim, 5); end    
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

function read_data_weighted(dataset::String, order::Int64, exact_match::Bool=false)
    if order < 3 || order > 5; throw("Only support for order-3,4,5 tensors"); end
    
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    simplices = read("data/$(dataset)/$(dataset)-simplices.txt")
    nverts = read("data/$(dataset)/$(dataset)-nverts.txt")

    I = Vector{Vector{Int64}}(order)
    for j in 1:order; I[j] = Vector{Int64}(); end
    curr_ind = 1
    for nvert in nverts
        simp = simplices[curr_ind:(curr_ind + nvert - 1)]
        curr_ind += nvert
        if (exact_match && nvert != order); continue ;end
        for hedge in combinations(simp, order)
            for (j, v) in enumerate(hedge)
                push!(I[j], v)
            end
        end
    end
    
    dim = maximum([maximum(I[j]) for j in 1:order])
    V = ones(Float64, length(I[1])) / length(I[1])
    if order == 3; return SymTensor3(I[1], I[2], I[3], V, dim, 3); end
    if order == 4; return SymTensor4(I[1], I[2], I[3], I[4], V, dim, 4); end
    if order == 5; return SymTensor5(I[1], I[2], I[3], I[4], I[5], V, dim, 5); end
end

""" Compute largest real eigenvector of A with unit 1-norm """
function LR_evec(A::SpFltMat)
    evec = eigs(A, nev=1, which=:LR, tol=1e-5, maxiter=200)[2][:,1]
    if evec[1] < 0; evec = -evec; end
    return evec / norm(evec, 1)
end

# Benson and Gleich algorithm for computing the leading Z-eigenvector.
function Z_evec_dynsys(T::SymTensor, tol::Float64=1e-5, niter::Int64=200)
    f(u::Vector{Float64}) = LR_evec(reduce(T, u)) - u
    x_curr = ones(Float64, T.dimension) / T.dimension    
    h = 0.5
    converged = false
    for i = 1:niter
        print("$i of $niter \r")
        flush(STDOUT)
        x_next = x_curr + h * f(x_curr)

        s = x_next ./ x_curr
        converged = (maximum(s) - minimum(s)) / minimum(s) < tol
        if converged; break; end

        x_curr = x_next
    end

    evec = x_curr
    return (evec, converged)
end

# NQI algorithm for computing the leading H-eigenvector.
function H_evec_NQI(T::SymTensor, niter::Int64=2000, tol::Float64=1e-5)
    m = T.order
    converged = false
    x = ones(Float64, T.dimension) / T.dimension
    y = apply(T, x)
    for i in 1:niter
        print("$i of $niter \r")
        flush(STDOUT)
        y_scaled = y .^ (1.0 / (m - 1))
        x = y_scaled / norm(y_scaled, 1)
        y = apply(T, x)
        s = y ./ (x .^ (m - 1))
        converged = (maximum(s) - minimum(s)) / minimum(s) < tol
        if converged; break; end
    end
    return (x, converged)
end

function CMEC(T::SymTensor)
    W = squeeze_tensor(T)
    lcc = MatrixNetworks.largest_component(W)[2]
    inds = find(lcc)
    return (LR_evec(W[inds, inds]), true, lcc)
end

function ZEC(T::SymTensor)
    Tcc, lcc = largest_component(T)[1:2]
    (centrality, converged) = Z_evec_dynsys(Tcc)
    return (centrality, converged, lcc)
end

function HEC(T::SymTensor)
    Tcc, lcc, ind_map = largest_component(T)
    centrality, converged = H_evec_LZI(Tcc)
    return (centrality, converged, lcc)
end
;
