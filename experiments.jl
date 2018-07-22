include("HypergraphEvecCentrality.jl")
using MAT
using PyPlot
using StatsBase

function collect_data(dataset::String, k::Int, exact::Bool)
    T = read_data_unweighted(dataset, k, exact)
    Tcc, lcc = largest_component(T)
    labels = read_node_labels(dataset)
    lcc_labels = labels[find(lcc)]
    get_top_labels(c::Vector{Float64}) =
        lcc_labels[sortperm(c, rev=true)[1:20]]
    cec_c, cec_conv = CEC(Tcc)
    zec_c, zec_conv  = ZEC(Tcc)
    hec_c, hec_conv  = HEC(Tcc)
    
    matwrite("$dataset-$k.mat",
             Dict("cec_c"      => cec_c,
                  "cec_conv"   => cec_conv,
                  "hec_c"      => hec_c,
                  "hec_conv"   => hec_conv,
                  "zec_c"      => zec_c,
                  "zec_conv"   => zec_conv,                     
                  "order"      => kunif,
                  "lcc_labels" => lcc_labels))
    
    cec_labels = get_top_labels(cec_c)
    zec_labels = get_top_labels(zec_c)
    hec_labels = get_top_labels(hec_c)
    println("k = $k")
    for (ind, (x, y, z)) in enumerate(zip(cec_labels, zec_labels, hec_labels))
        println("$ind: $x\t\t$y\t\t$z")
    end
    println("--------------\n")
end

function rank_corr(dataset::String, start::Int64)
    close()
    for k in [3, 4, 5]
        data = matread("results/$dataset-$k.mat")
        cec_c, zec_c, hec_c = data["cec_c"], data["zec_c"], data["hec_c"]
        sp = sortperm(cec_c)
        cec_c, hec_c, zec_c = cec_c[sp], hec_c[sp], zec_c[sp]

        start = log10(start - 1)
        finish = log10(length(cec_c) - 1)
        ran = [convert(Int64, round(v)) for v in logspace(start, finish, 500)]
        cec_hec = [corspearman(cec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
        cec_zec = [corspearman(cec_c[(end-s):end], zec_c[(end-s):end]) for s in ran]
        hec_zec = [corspearman(hec_c[(end-s):end], zec_c[(end-s):end]) for s in ran]
        figure()
        xs = collect(ran) + 1
        semilogx(xs, cec_hec, linestyle="-",  lw=3, label="CEC/HEC")
        semilogx(xs, cec_zec, linestyle="-",  lw=1, label="CEC/ZEC")
        semilogx(xs, hec_zec, linestyle="--", lw=2, label="HEC/ZEC")
        fsz = 20
        ax = gca()
        ax[:tick_params](axis="both", labelsize=fsz-2, length=7, width=1.5)
        ax[:tick_params](axis="x",    which="minor", length=4, width=1)
        if k == 3; legend(fontsize=fsz-2, frameon=false); end
        title("$(dataset) ($(k)-uniform)", fontsize=fsz)
        xlabel("Number of top-ranked elements", fontsize=fsz)
        ylabel("Spearman's rank corr. coeff.", fontsize=fsz)
        savefig("$(dataset)-$(k)-rankcorr.eps", bbox_inches="tight")
    end
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
