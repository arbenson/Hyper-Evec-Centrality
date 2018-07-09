include("HypergraphEvecCentrality.jl")

using Combinatorics

function main()
    dim = 7
    m = 3

    hedges = collect(combinations(1:dim, m))
    total = 2^length(hedges)
    for test_hgraph in combinations(hedges, 4)
    #for ind in shuffle(collect(0:(total-1)))
        #keeps = [parse(Bool, b) for b in bin(ind, length(hedges))]
        #test_hgraph = hedges[keeps]
        #=
        T = SymTensor4([h[1] for h in test_hgraph],
                       [h[2] for h in test_hgraph],
                       [h[3] for h in test_hgraph],
                       [h[4] for h in test_hgraph],                                   
                       ones(Float64, length(test_hgraph)),
                       dim, m)
        =#
        T = SymTensor3([h[1] for h in test_hgraph],
                       [h[2] for h in test_hgraph],
                       [h[3] for h in test_hgraph],
                       ones(Float64, length(test_hgraph)),
                       dim, m)

        # Skip cases that aren't connected
        lcc = largest_component(T)[2]
        if length(find(lcc)) < dim; continue; end
        
        # Get evec / eval
        z, conv = Z_evec_dynsys(T)
        if !conv; continue; end
        z2 = z / norm(z, 2)
        eval = mean(apply(T, z2) ./ z2)

        # Check for stability
        Ustar = qr((I - z2 * z2') * randn(dim, dim - 1), thin=true)[1]
        C = Ustar' * ((m - 1) * reduce(T, z2) - eval * eye(dim)) * Ustar
        U = triu(C, 1)
        λ = real.(eig(U + U' + diagm(diag(C)))[1])

        if maximum(λ) > 0 && minimum(λ) < 0; @show λ, test_hgraph, z; end
    end
end

function test2()
    # (λ, test_hgraph, z) =
    # ([-2.82842, -2.82842, -2.82842, -2.06111, -1.41421, 0.646902],
    # Array{Int64,1}[[1, 2, 3], [1, 2, 4], [3, 5, 6], [5, 6, 7]],
    # [0.158494, 0.158494, 0.183012, 0.091507, 0.158494, 0.158494, 0.091507])
    dim = 7
    m = 3
    T = SymTensor3([1, 1, 3, 5],
                   [2, 2, 5, 6],
                   [3, 4, 6, 7],
                   ones(Float64, 4),
                   dim, m)
    M = reduce(T, ones(Float64, dim))
    return M, T
end
