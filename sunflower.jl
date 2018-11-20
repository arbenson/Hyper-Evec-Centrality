include("centrality.jl")

# Generate 3-uniform, r-petal sunflower
function sunflower3(r::Int64)
    I1, I2, I3 = Int64[], Int64[], Int64[]
    core = 1
    curr_ind = 2
    for j in 1:r  # each petal
        push!(I1, core)
        push!(I2, curr_ind); curr_ind += 1
        push!(I3, curr_ind); curr_ind += 1
    end
    return SymTensor3(I1, I2, I3, ones(Float64, length(I1)), curr_ind - 1)
end

# Generate 4-uniform, r-petal sunflower
function sunflower4(r::Int64)
    I1, I2, I3, I4 = Int64[], Int64[], Int64[], Int64[]
    core = 1
    curr_ind = 2
    for j in 1:r  # each petal
        push!(I1, core)
        push!(I2, curr_ind); curr_ind += 1
        push!(I3, curr_ind); curr_ind += 1
        push!(I4, curr_ind); curr_ind += 1        
    end
    return SymTensor4(I1, I2, I3, I4, ones(Float64, length(I1)), curr_ind - 1)
end

# Generate 5-uniform, r-petal sunflower
function sunflower5(r::Int64)
    I1, I2, I3, I4, I5 = Int64[], Int64[], Int64[], Int64[], Int64[]
    core = 1
    curr_ind = 2
    for j in 1:r  # each petal
        push!(I1, core)
        push!(I2, curr_ind); curr_ind += 1
        push!(I3, curr_ind); curr_ind += 1
        push!(I4, curr_ind); curr_ind += 1
        push!(I5, curr_ind); curr_ind += 1        
    end
    return SymTensor5(I1, I2, I3, I4, I5, ones(Float64, length(I1)), curr_ind - 1)
end

# Test of ratios
function sf_test1()
    function ratios(m::Int64, r::Int64)
        cec_ratio = 2 * r * (m - 1) / (sqrt(m^2 + 4 * (m - 1) * (r - 1)) + m - 2)
        zec_ratio = sqrt(r)
        hec_ratio = r^(1/m)
        return (cec_ratio, zec_ratio, hec_ratio)
    end

    function get_sunflower(m::Int64, r::Int64)
        if m == 3; return sunflower3(r); end
        if m == 4; return sunflower4(r); end
        if m == 5; return sunflower5(r); end
    end
    
    for m in 3:5, r in 2:8
        println("m = $m; r = $r")
        S = get_sunflower(m, r)
        cec_ratio, zec_ratio, hec_ratio = ratios(m, r)
        cec_v = CEC(S)[1]
        zec_v = ZEC(S)[1]
        hec_v = HEC(S)[1]
        cec_normalized = cec_v[1] ./ cec_v
        zec_normalized = zec_v[1] ./ zec_v
        hec_normalized = hec_v[1] ./ hec_v        
        println("cec ratio = $(cec_ratio)")
        println("cec vec = $(cec_normalized)")
        println("zec ratio = $(zec_ratio)")
        println("zec vec = $(zec_normalized)")
        println("hec ratio = $(hec_ratio)")
        println("hec vec = $(hec_normalized)")
        println("-----")
    end
end

function sf_test2()
    r = 4
    S = sunflower3(r)
    c = zeros(2 * r + 1)
    cP = [1.0, 2.0, 3.0, 4.0]
    c[2:3] .= cP[1]
    c[4:5] .= cP[2]
    c[6:7] .= cP[3]
    c[8:9] .= cP[4]
    c[1] = sqrt(sum(cP .^ 2))
    rat = apply(S, c) ./ c
    println("Elements should be same: $(rat)")

    S = sunflower4(r)
    c = zeros(3 * r + 1)    
    cP = [1.0, 2.0, 3.0, 4.0]
    c[2:4]   .= cP[1]
    c[5:7]   .= cP[2]
    c[8:10]  .= cP[3]
    c[11:13] .= cP[4]
    c[1] = sqrt(sum(cP .^ 2))
    rat = apply(S, c) ./ c    
    println("Elements should be different: $(rat)")    
end
;
