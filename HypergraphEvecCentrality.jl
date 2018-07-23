include("tensor.jl")

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

function CEC(T::SymTensor)
    c = LR_evec(squeeze_tensor(T))
    return (c / norm(c, 1), true)
end

function ZEC(T::SymTensor)
    (c, converged) = Z_evec_dynsys(T)
    return (c / norm(c, 1), converged)
end

function HEC(T::SymTensor)
    c, converged = H_evec_NQI(T)
    return (c / norm(c, 1), converged)
end
;
