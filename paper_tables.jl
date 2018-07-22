include("HypergraphEvecCentrality.jl")
using MAT

function summary_statistics(dataset::String, exact::Bool, kunif::Int)
    exact_int = convert(Int64, exact)
    T = read_data_unweighted(dataset, kunif, exact)
    Tcc, lcc = largest_component(T)
    nnodes = Tcc.dimension
    nnz = length(Tcc.I1)
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
