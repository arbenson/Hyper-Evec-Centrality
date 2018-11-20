include("centrality.jl")

using Statistics

function main()
    m = 3
    dim = 7
    T = SymTensor3([1, 1, 3, 5],
                   [2, 2, 5, 6],
                   [3, 4, 6, 7],
                   ones(Float64, 4),
                   dim)
    
    z, conv = Z_evec_dynsys(T)
    println("dynamical systems algorithm converged: $conv")
    z2 = z / norm(z, 2)
    λ1 = mean(apply(T, z2) ./ z2)
    println("tensor z-eigenvalue: $(λ1)")
    println("tensor z-eigenvector (2-norm normalized): $(z2)")

    # Check for stability
    Ustar = Array(qr((I - z2 * z2') * randn(dim, dim - 1)).Q)
    maxval = maximum(z2' * Ustar )
    println("orthogonality test (should be close to 0): $(maxval)")
    diff = Ustar' * Ustar - I
    diffnorm = [norm(vec(diff[:, ind]), 2) for ind in 1:size(diff)[2]]
    println("orthonormality test (should be close to 0): $(diffnorm)")    

    C = Ustar' * ((m - 1) * reduce(T, z2) - λ1 * I) * Ustar
    U = triu(C, 1)  # for symmetrization
    C_sym = Matrix(U + U' + Diagonal(C))
    spectrum = eigen(C_sym).values
    println("Spectrum of projected Hessian of the Lagrangian: $spectrum")
end
;
