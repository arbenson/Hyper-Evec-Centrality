using FileIO
using JLD2
using Distances
using PyPlot
using StatsBase
using Statistics

function read_evecs(dataset::String, k::Int64)
    data1 = load("results/$dataset-$k.jld2")
    data2 = load("results/$dataset-$k-random.jld2")

    # Get mode Z-eigenvector
    z_evecs = data2["z_evecs"]
    dists = pairwise(Cityblock(), z_evecs, dims=2)
    cutoff = 1e-5
    n = size(z_evecs, 2)
    # group vectors by distance
    components = zeros(Int64, n)
    for i = 1:n
        if components[i] == 0
            components[i] = maximum(components) + 1
        end
        for j = 1:n
            if dists[i, j] < cutoff && components[j] == 0
                components[j] = components[i]
            end
        end
    end
    # find mode vector index
    biggest = argmax(counts(components))
    z_mode_ind = 0
    for i = 1:n
        if components[i] == biggest
            z_mode_ind = i
            break
        end
    end
    
    return data1["cec_c"], data1["hec_c"], z_evecs[:, z_mode_ind]
end

function summary_statistics(dataset::String, k::Int64)
    data = load("results/$dataset-$k.jld2")
    nnodes, nnz = data["nnodes"], data["nnz"]
    println("# nodes: $nnodes")
    println("nnz (ignoring symmetries): $nnz")
end

# Top ranked items
function top_ranked(dataset::String, numtop::Int64, singlerow::Bool=false)
    label_mat = Array{String}(undef, numtop, 9)
    col_ind = 1
    for k in [3, 4, 5]
        data = load("results/$dataset-$(k).jld2")
        lcc_labels = data["lcc_labels"]
        cec_c, hec_c, zec_c = read_evecs(dataset, k)

        get_top_labels(c::Vector{Float64}) =
            [string(s) for s in lcc_labels[sortperm(c, rev=true)[1:numtop]]]
        
        label_mat[:, col_ind] = get_top_labels(cec_c); col_ind += 1
        label_mat[:, col_ind] = get_top_labels(zec_c); col_ind += 1
        label_mat[:, col_ind] = get_top_labels(hec_c); col_ind += 1

        local_label_mat = Array{String}(undef, numtop, 3)
        local_label_mat[:, 1] = get_top_labels(cec_c)
        local_label_mat[:, 2] = get_top_labels(zec_c)
        local_label_mat[:, 3] = get_top_labels(hec_c)

        if !singlerow
            for i in 1:numtop
                #str = [join(split(s)[2:end], " ") for s in collect(local_label_mat[i, :])]
                #str = join([lowercase(s) for s in str], " & ")
                str = join([lowercase(s) for s in collect(local_label_mat[i, :])], " & ")
                str = "& $i & $(str) \\\\"
                println(str)
            end
            println("--------------------")
        end
    end

    if singlerow
        for i in 1:numtop
            str = join([lowercase(s) for s in collect(label_mat[i, :])], " & ")            
            str = "$i & $(str) \\\\"
            println(str)
        end
    end
end

# Rank correlation plot
function rank_corr(dataset::String, k::Int64, index_start::Int64)
    close()
    cec_c, hec_c, zec_c = read_evecs(dataset, k)
    @show length(cec_c), length(hec_c), length(zec_c)
    
    start = log10(index_start - 1)
    finish = log10(length(cec_c) - 1)
    ran = [convert(Int64, round(v)) for v in 10 .^ range(start, stop=finish, length=500)]
    
    sp = sortperm(cec_c)
    cec_c, hec_c, zec_c = cec_c[sp], hec_c[sp], zec_c[sp]
    cec_hec = [corspearman(cec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
    cec_zec = [corspearman(cec_c[(end-s):end], zec_c[(end-s):end]) for s in ran]
    sp = sortperm(zec_c)
    cec_c, hec_c, zec_c = cec_c[sp], hec_c[sp], zec_c[sp]
    zec_cec = [corspearman(zec_c[(end-s):end], cec_c[(end-s):end]) for s in ran]
    zec_hec = [corspearman(zec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
    sp = sortperm(hec_c)
    cec_c, hec_c, zec_c = cec_c[sp], hec_c[sp], zec_c[sp]
    hec_cec = [corspearman(hec_c[(end-s):end], cec_c[(end-s):end]) for s in ran]
    hec_zec = [corspearman(hec_c[(end-s):end], zec_c[(end-s):end]) for s in ran]

    xs = collect(ran) .+ 1
    semilogx(xs, cec_hec, linestyle="-",  lw=1,    label="CEC-HEC")
    semilogx(xs, cec_zec, linestyle="-",  lw=3.5,  label="CEC-ZEC")
    semilogx(xs, zec_cec, linestyle="--", lw=2.25, label="ZEC-CEC")
    semilogx(xs, zec_hec, linestyle="--", lw=1,    label="ZEC-HEC")
    semilogx(xs, hec_cec, linestyle=":",  lw=2.25, label="HEC-CEC")
    semilogx(xs, hec_zec, linestyle=":",  lw=3.5,  label="HEC-ZEC")

    fsz = 20
    ax = gca()
    ax[:tick_params](axis="both", labelsize=fsz-2, length=7, width=1.5)
    ax[:tick_params](axis="x",    which="minor", length=4, width=1)
    if k == 3; legend(fontsize=fsz-4, frameon=false); end
    title("$(dataset) ($(k)-uniform)", fontsize=fsz)
    xlabel("Number of top-ranked elements", fontsize=fsz)
    ylabel("Spearman's rank corr. coeff.", fontsize=fsz)
    savefig("$(dataset)-$(k)-rankcorr.eps", bbox_inches="tight")
end

function get_rank(dataset::String, k::Int, target_label::String)
    cec_c, hec_c, zec_c = read_evecs(dataset, k)    
    n = length(cec_c)
    sp_cec = zeros(Int64, n); sp_cec[sortperm(cec_c, rev=true)] = collect(1:n)
    sp_zec = zeros(Int64, n); sp_zec[sortperm(zec_c, rev=true)] = collect(1:n)
    sp_hec = zeros(Int64, n); sp_hec[sortperm(hec_c, rev=true)] = collect(1:n)    
    data = load("results/$dataset-$k.jld2")
    labels = data["lcc_labels"]

    for (label, cec_rk, zec_rk, hec_rk) in zip(labels, sp_cec, sp_zec, sp_hec)
        if occursin(target_label, string(label))
            println("$label $(cec_rk) $(zec_rk) $(hec_rk)")
        end
    end
end
;
