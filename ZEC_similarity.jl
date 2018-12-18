include("centrality.jl")

using Random
using StatsBase

function collect_data(dataset::String, k::Int64, exact::Bool)
    T = read_data_unweighted(dataset, k, exact)
    Tcc, lcc = largest_component(T)    
    ntrials = 100
    z_evecs = zeros(Float64, Tcc.dimension, ntrials)
    convergences = zeros(Bool, ntrials)
    Random.seed!(1234)
    for trial = 1:ntrials
        x_init = rand(Tcc.dimension)
        x_init /= sum(x_init)
        evec, conv = Z_evec_dynsys(Tcc, 1e-5, 500, x_init)
        println("trial $trial converged = $conv")
        z_evecs[:, trial] = evec
        convergences[trial] = conv
    end

    save("results/$dataset-$k-random.jld2",    
         Dict("z_evecs" => z_evecs,
              "conv"    => convergences))
end

function compute_similarity(dataset::String, k::Int64)
    data = load("results/$dataset-$k-random.jld2")
    z_evecs = data["z_evecs"]
    ntrials = size(z_evecs, 2)
    max_diff = 0.0
    min_cor  = 1.0
    for i = 1:ntrials
        zi = z_evecs[:, i]
        for j = (i + 1):ntrials
            zj = z_evecs[:, j]
            @show i, j, norm(zi - zj, 1)
            max_diff = max(max_diff, norm(zi - zj, 1))
            min_cor  = min(min_cor, corspearman(zi, zj))
        end
    end
    return max_diff, min_cor
end
