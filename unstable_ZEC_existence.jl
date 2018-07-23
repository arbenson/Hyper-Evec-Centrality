include("HypergraphEvecCentrality.jl")

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
    Ustar = qr((I - z2 * z2') * randn(dim, dim - 1), thin=true)[1]
    maxval = maximum(z2' * Ustar )
    println("orthogonality test (should be close to 0): $(maxval)")
    diff = vecnorm(Ustar' * Ustar - eye(dim - 1))
    println("orthonormality test (should be close to 0): $(diff)")    

    C = Ustar' * ((m - 1) * reduce(T, z2) - λ1 * eye(dim)) * Ustar
    U = triu(C, 1) # for symmetrization
    C_sym = full(U + U' + diagm(diag(C)))
    spectrum = eig(C_sym)[1]
    println("Spectrum of projected Hessian of the Lagrangian: $spectrum")
end
