# Hypergraph eigenvector centralities

This code and data repository accompanies the paper

- Three hypergraph eigenvector centralities. Austin R. Benson. [arXiv:1807.09644](http://arxiv.org/abs/1807.09644), 2018.

All of the code is written in Julia 1.0.

For questions, please email Austin at arb@cs.cornell.edu.

### Data

The tags and DAWN datasets are in the `data/` directory. These files were downloaded directly from [here](http://www.cs.cornell.edu/~arb/data/tags-ask-ubuntu/index.html) and [here](http://www.cs.cornell.edu/~arb/data/DAWN/index.html). The ngrams dataset is available from https://www.ngrams.info/, but I cannot include the data directly due to the usage terms. The script `data/convert_ngrams.jl ` can be used to convert the raw ngrams data into a format that can be automatically read by the existing code.

### Setup

First, download the code:

```bash
git clone https://github.com/arbenson/Hyper-Evec-Centrality.git
```

The code is written in Julia 1.0. To reproduce everything, you will need the following packages:

```julia
using Pkg
Pkg.add("Arpack")
Pkg.add("Combinatorics")
Pkg.add("FileIO")
Pkg.add("JLD2")
Pkg.add("MatrixNetworks")
Pkg.add("PyPlot")
Pkg.add("StatsBase")
```

Then you can test the code

```julia
include("test.jl")
main()
```

### Examples for computing eigenvectors

```julia
include("centrality.jl")
# Create the hypergraph adjacency tensor for the hypergraph in Proposition 2.8.
# Hyperedges: {{1, 2, 3}, {1, 2, 4}, {3, 5, 6}, {5, 6, 7}}
T = SymTensor3([1, 1, 3, 5], [2, 2, 5, 6], [3, 4, 6, 7], ones(Float64, 4))

# We might first want to take the largest component.
Tlcc, lcc_inds = largest_component(T)
# However, the hypergraph is connected here.

# Clique motif Eigenvector Centrality (CEC)
(cec_c, cec_converged) = CEC(T)
# Z-eigenvector Centrality (ZEC)
(zec_c, zec_converged) = ZEC(T)
# H-eigenvector Centrality (HEC)
(hec_c, hec_converged) = HEC(T)
```

### Computational experiments for Proposition 2.8

Proposition 2.8 says that Z-eigenvector centrality (ZEC) vectors can be unstable in the sense of the definitions [the paper](https://epubs.siam.org/doi/abs/10.1137/100801482) from Kolda and Mayo analyzing their symmetric shifted higher-order power method for computing Z-eigenvectors of symmetric tensors. The consequence is that we cannot rely on power-method like algorithms for computing ZEC vectors. This script contains the computational experiments demonstrating the proof of the proposition:

```julia
include("unstable_ZEC_existence.jl")
main()
```

The output should look something like the following

```
dynamical systems algorithm converged: true
tensor z-eigenvalue: 1.4142111552433783
tensor z-eigenvector (2-norm normalized): [0.408248, 0.408248, 0.471403, 0.235704, 0.408248, 0.408248, 0.235704]
orthogonality test (should be close to 0): 9.71445146547012e-16
orthonormality test (should be close to 0): [3.70907e-16, 2.66455e-16, 2.35813e-16, 2.87556e-16, 1.0074e-15, 9.14796e-16]
Spectrum of projected Hessian of the Lagrangian: [-2.82842, -2.82842, -2.82842, -2.06111, -1.41421, 0.646902]
```

Notice that the projected Hessian is indefinite since it has both a negative and positive eigenvalue.

### Computational experiments for sunflower hypergraphs.

The paper provides analytic results for sunflower hypergraphs where the core is a singleton. The script `sunflower.jl ` provides numerical evidence for both the results listed in the right of Figure 1 and Proposition 2.12.

```julia
include("sunflower.jl")
sf_test1()  # computational experiments related to the table within Figure 1.
sf_test2()  # computational experiments related to Proposition 2.12.
```

### Reproduce the figures and tables in the paper

First, we need to compute all of the centralities. Here is an exmaple for the 4-uniform tags dataset. This code is designed to use the data format of the temporal higher-order networks available [here](http://www.cs.cornell.edu/~arb/data/).

```julia
include("centrality.jl")
collect_data("tags-ask-ubuntu", 4, false)  # produces tags-ask-ubuntu-4.jld2
collect_data("DAWN", 5, true)  # produces DAWN-5.jld2
```

In the second command, the parameter "4" means 4-uniform hypergraph, and the parameter "false" means to not look for exact size matching in the dataset. In this case, A size-5 hyperedge becomes five different size-4 hyperedges. In the third command, we look at a 5-uniform hypergraph, and we only take hyperedges in the original data since the third parameter is set to "true".

We have stored all of the computational results in the `results/` directory:

-  `results/ngrams-{3,4,5}.jld2`
-  `results/tags-ask-ubuntu-{3,4,5}.jld2`
-  `results/DAWN-{3,4,5}.jld2`

##### Figures 2, 3, and 4.

These are the rank correlation plots.

```julia
include("paper_tables_and_figures.jl")

# Create ngrams-{3,4,5}-rankcorr.eps (x axis starts from 20)
for k in [3, 4, 5]; rank_corr("ngrams", k, 20); end 

# Create tags-ask-ubuntu-{3,4,5}-rankcorr.eps (x axis starts from 10)
for k in [3, 4, 5]; rank_corr("tags-ask-ubuntu", k, 10); end

# Create DAWN-{3,4,5}-rankcorr.eps (x axis starts from 10)
for k in [3, 4, 5]; rank_corr("DAWN", k, 10); end
```

##### Table 1

This table has the summary statistics of the datasets.

```julia
include("paper_tables_and_figures.jl")
for k in [3, 4, 5]; summary_statistics("ngrams", k); end
for k in [3, 4, 5]; summary_statistics("DAWN", k); end
for k in [3, 4, 5]; summary_statistics("tags-ask-ubuntu", k); end
```

##### Tables 2, 3, and 4

These are the top-k ranked elements according to the centrality vectors.

```julia
include("paper_tables_and_figures.jl")

# Table 2 
# Print top-20 words in ngrams datasets in a single row
top_ranked("ngrams", 20, true)

# Table 3
# Print top-10 tags in ask ubuntu dataset; one block for each hypergraph uniformity.
top_ranked("tags-ask-ubuntu", 10, false)

# Table 4
# Print top-10 drugs in DAWN dataset; one block for each hypergraph uniformity.
top_ranked("DAWN", 10, false)
```

