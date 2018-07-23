using MAT
using PyPlot
using StatsBase

function summary_statistics(dataset::String, k::Int64)
    data = matread("results/$dataset-$k.mat")    
    nnodes, nnz = data["nnodes"], data["nnz"]
    println("# nodes = $nnodes")
    println("nnz = $nnz")
end

function top_ranked(dataset::String, numtop::Int64, combined::Bool=false)
    label_mat = Array{String}(numtop, 9)
    col_ind = 1
    for k in [3, 4, 5]
        data = matread("results/$dataset-$(k).mat")
        lcc_labels = data["lcc_labels"]
        cec_c, hec_c, zec_c = data["cec_c"], data["hec_c"], data["zec_c"]

        get_top_labels(c::Vector{Float64}) =
            [string(s) for s in lcc_labels[sortperm(c, rev=true)[1:numtop]]]
        
        label_mat[:, col_ind] = get_top_labels(cec_c); col_ind += 1
        label_mat[:, col_ind] = get_top_labels(zec_c); col_ind += 1
        label_mat[:, col_ind] = get_top_labels(hec_c); col_ind += 1

        local_label_mat = Array{String}(numtop, 3)
        local_label_mat[:, 1] = get_top_labels(cec_c)
        local_label_mat[:, 2] = get_top_labels(zec_c)
        local_label_mat[:, 3] = get_top_labels(hec_c)

        if !combined
            for i in 1:numtop
                str = join(collect(local_label_mat[i, :]), " & ")
                str = "& $i & $(str) \\\\"
                println(str)
            end
        end
    end

    if combined
        for i in 1:numtop
            strs = [split(s, "[")[1] for s in collect(label_mat[i, :])]
            str = join(strs, " & ")
            str = "$i & $(str) \\\\"
            println(str)
        end
    end
end

function rank_corr(dataset::String, index_start::Int64, k::Int64)
    close()
    data = matread("results/$dataset-$k.mat")
    cec_c, zec_c, hec_c = data["cec_c"], data["zec_c"], data["hec_c"]
    
    start = log10(index_start - 1)
    finish = log10(length(cec_c) - 1)
    ran = [convert(Int64, round(v)) for v in logspace(start, finish, 500)]
    
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

    xs = collect(ran) + 1
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
    data = matread("results/$dataset-$k.mat")
    cec_c, hec_c, zec_c = data["cec_c"], data["hec_c"], data["zec_c"]
    n = length(cec_c)
    sp_cec = zeros(Int64, n); sp_cec[sortperm(cec_c, rev=true)] = collect(1:n)
    sp_zec = zeros(Int64, n); sp_zec[sortperm(zec_c, rev=true)] = collect(1:n)
    sp_hec = zeros(Int64, n); sp_hec[sortperm(hec_c, rev=true)] = collect(1:n)    
    labels = data["lcc_labels"]

    for (label, cec_rk, zec_rk, hec_rk) in zip(labels, sp_cec, sp_zec, sp_hec)
        if contains(string(label), target_label)
            println("$label $(cec_rk) $(zec_rk) $(hec_rk)")
        end
    end
end
