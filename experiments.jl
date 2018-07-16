include("HypergraphEvecCentrality.jl")
using MAT
using PyPlot
using StatsBase

function collect_data(dataset::String, exact::Bool, kunif::Int)
    exact_int = convert(Int64, exact)

    T = read_data_unweighted(dataset, kunif, exact)
    Tcc, lcc = largest_component(T)
    labels = read_node_labels(dataset)
    lcc_labels = labels[find(lcc)]
    get_top_labels(c::Vector{Float64}) =
        lcc_labels[sortperm(c, rev=true)[1:10]]
    cec_c, cec_conv = CEC(Tcc)
    zec_c, zec_conv  = ZEC(Tcc)
    hec_c, hec_conv  = HEC(Tcc)
    
    matwrite("output/$dataset-$(kunif)-$(exact_int).mat",
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
    println("k = $(kunif)")
    for (ind, (x, y, z)) in enumerate(zip(cec_labels, zec_labels, hec_labels))
        println("$ind: $x\t\t$y\t\t$z")
    end
    println("--------------\n")
end

function rank_corr(dataset::String, exact::Bool)
    close()
    exact_int = convert(Int64, exact)
    for k in [3, 4, 5]
        data = matread("output/$dataset-$k-$(exact_int).mat")
        cec_c, hec_c, zec_c = data["cec_c"], data["hec_c"], data["zec_c"]
        sp = sortperm(cec_c)
        cec_c, hec_c, zec_c = cec_c[sp], hec_c[sp], zec_c[sp]

        ran = 9:(length(cec_c) - 1)
        cec_hec = [corspearman(cec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
        cec_zec = [corspearman(cec_c[(end-s):end], zec_c[(end-s):end]) for s in ran]
        hec_zec  = [corspearman(hec_c[(end-s):end],  zec_c[(end-s):end]) for s in ran]
        figure()
        xs = collect(ran) + 1
        semilogx(xs, cec_hec, linestyle="-",  lw=3, label="CEC/HEC")
        semilogx(xs, cec_zec, linestyle="-",  lw=1, label="CEC/ZEC")
        semilogx(xs, hec_zec,  linestyle="--", lw=2, label="HEC/ZEC")
        fsz = 18
        ax = gca()
        ax[:tick_params](axis="both", labelsize=fsz-2, length=7, width=1.5)
        ax[:tick_params](axis="x",    which="minor", length=4, width=1)
        legend(fontsize=fsz-2, frameon=false)
        title("$(dataset) ($(k)-uniform)", fontsize=fsz)
        xlabel("Number of top-ranked elements", fontsize=fsz)
        ylabel("Spearman's rank corr. coeff.", fontsize=fsz)
        savefig("output/$(dataset)-$(k)-rankcorr.eps", bbox_inches="tight")
    end
end

function tops(dataset::String, exact::Bool, kunif::Int64, numtop::Int64)
    exact_int = convert(Int64, exact)
    data = matread("output/$dataset-$(kunif)-$(exact_int).mat")
    lcc_labels = data["lcc_labels"]
    cec_c, hec_c, zec_c = data["cec_c"], data["hec_c"], data["zec_c"]

    get_top_labels(c::Vector{Float64}) =
        lcc_labels[sortperm(c, rev=true)[1:numtop]]

    cec_labels = get_top_labels(cec_c)
    zec_labels = get_top_labels(zec_c)
    hec_labels = get_top_labels(hec_c)
    println("k = $(kunif)")
    for (ind, (x, y, z)) in enumerate(zip(cec_labels, zec_labels, hec_labels))
        println("$ind: $x\t\t$y\t\t$z")
    end
    println("--------------\n")
end
