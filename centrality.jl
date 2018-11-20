include("tensor.jl")
include("data_io.jl")

using Arpack
using LinearAlgebra
using FileIO
using JLD2

# Compute largest real eigenvector of A with unit 1-norm
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
        flush(stdout)
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
function H_evec_NQI(T::SymTensor, m::Int64, niter::Int64=2000, tol::Float64=1e-5)
    converged = false
    x = ones(Float64, T.dimension) / T.dimension
    y = apply(T, x)
    for i in 1:niter
        print("$i of $niter \r")
        flush(stdout)
        y_scaled = y .^ (1.0 / (m - 1))
        x = y_scaled / norm(y_scaled, 1)
        y = apply(T, x)
        s = y ./ (x .^ (m - 1))
        converged = (maximum(s) - minimum(s)) / minimum(s) < tol
        if converged; break; end
    end
    return (x, converged)
end

function CEC(T::SymTensor)
    c = LR_evec(squeeze_tensor(T))
    return (c / norm(c, 1), true)
end

function ZEC(T::SymTensor)
    (c, converged) = Z_evec_dynsys(T)
    return (c / norm(c, 1), converged)
end

function HEC(T::SymTensor)
    m = 0
    if typeof(T) == SymTensor3; m = 3; end
    if typeof(T) == SymTensor4; m = 4; end
    if typeof(T) == SymTensor5; m = 5; end
    if m == 0; error("Only support for order-3,4,5 tensors"); end
    c, converged = H_evec_NQI(T, m)
    return (c / norm(c, 1), converged)
end

function collect_data(dataset::String, k::Int, exact::Bool)
    T = read_data_unweighted(dataset, k, exact)
    Tcc, lcc = largest_component(T)
    labels = read_node_labels(dataset)
    lcc_labels = labels[findall(lcc)]
    get_top_labels(c::Vector{Float64}) =
        lcc_labels[sortperm(c, rev=true)[1:20]]
    cec_c, cec_conv = CEC(Tcc)
    zec_c, zec_conv  = ZEC(Tcc)
    hec_c, hec_conv  = HEC(Tcc)
    
    save("results/$dataset-$k.jld2",
         Dict("cec_c"         => cec_c,
              "cec_converged" => cec_conv,
              "hec_c"         => hec_c,
              "hec_converged" => hec_conv,
              "zec_c"         => zec_c,
              "zec_converged" => zec_conv,                     
              "order"         => k,
              "lcc_labels"    => lcc_labels,
              "nnz"           => length(Tcc.I1),
              "nnodes"        => Tcc.dimension))
    
    cec_labels = get_top_labels(cec_c)
    zec_labels = get_top_labels(zec_c)
    hec_labels = get_top_labels(hec_c)
    println("k = $k")
    for (ind, (x, y, z)) in enumerate(zip(cec_labels, zec_labels, hec_labels))
        println("$ind: $x\t\t$y\t\t$z")
    end
    println("--------------\n")
end
;
