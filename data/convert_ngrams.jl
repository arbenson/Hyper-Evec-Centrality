# Convert n-grams data to same format.
function main(n::Int64)
    gram_map = Dict{AbstractString,Int64}()
    function get_value(key::AbstractString)
        if !haskey(gram_map, key)
            val = length(gram_map) + 1
            gram_map[key] = val
            return val
        end
        return gram_map[key]
    end
    
    uniques = Set{NTuple{n,Int64}}()
    open("ngrams-$(n)/w$(n)_.txt") do f
        for line in eachline(f)
            vals = NTuple{n,Int64}([get_value(k) for k in split(line)[2:end]])
            if length(vals) != n
                println(line)
                @show split(line)[2:end]
                @show length(vals)
                assert(length(vals) == n)
            end
            push!(uniques, vals)
        end
    end

    nverts = Int64[]
    simplices = Int64[]
    for ngram in uniques
        push!(nverts, n)
        append!(simplices, collect(ngram))
    end
    
    open("ngrams-$(n)/ngrams-$(n)-nverts.txt", "w") do f
        for nv in nverts; write(f, "$(nv)\n"); end
    end
    open("ngrams-$(n)/ngrams-$(n)-simplices.txt", "w") do f
        for node in simplices; write(f, "$(node)\n"); end
    end
    open("ngrams-$(n)/ngrams-$(n)-node-labels.txt", "w") do f
        for (val, key) in sort([(val, key) for (key, val) in gram_map])
            write(f, "$(val) $(key)\n")
        end
    end
end
