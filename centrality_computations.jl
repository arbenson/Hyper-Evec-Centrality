include("HypergraphEvecCentrality.jl")
using MAT

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
                  "lcc_labels" => lcc_labels,
                  "nnz"        => length(Tcc.I1),
                  "nnodes"     => Tcc.dimension))
    
    cec_labels = get_top_labels(cec_c)
    zec_labels = get_top_labels(zec_c)
    hec_labels = get_top_labels(hec_c)
    println("k = $k")
    for (ind, (x, y, z)) in enumerate(zip(cec_labels, zec_labels, hec_labels))
        println("$ind: $x\t\t$y\t\t$z")
    end
    println("--------------\n")
end
