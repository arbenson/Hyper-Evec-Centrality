# Hypergraph eigenvector centralities

This code and data repository accompanies the paper

- Three hypergraph eigenvector centralities. Austin R. Benson. Submitted, 2018.

All of the code is written in Julia.

For questions, please email Austin at arb@cs.cornell.edu.

### Data

The tags and DAWN datasets are in the `data/` directory. These files were downloaded directly from [here](http://www.cs.cornell.edu/~arb/data/tags-ask-ubuntu/index.html) and [here](http://www.cs.cornell.edu/~arb/data/DAWN/index.html).

### Examples computing eigenvectors

### Computational experiment backing up the proof of Proposition 2.8

Proposition 2.8 says that Z-eigenvector centrality (ZEC) vectors can be unstable in the sense of the definitions [the paper](https://epubs.siam.org/doi/abs/10.1137/100801482) from Kolda and Mayo analyzing their symmetric shifted higher-order power method for computing Z-eigenvectors of symmetric tensors. The consequence is that we cannot rely on power-method like algorithms for computing ZEC vectors. This script contains the computational experiments demonstrating the proof of the proposition:

```julia
include("unstable_ZEC_existence.jl")
main()
```

The output should look something like the following

```
dynamical systems algorithm converged: true
tensor z-eigenvalue: 1.4142111552433778
tensor z-eigenvector (2-norm normalized): [0.408248, 0.408248, 0.471403, 0.235704, 0.408248, 0.408248, 0.235704]
orthogonality test (should be close to 0): 1.1102230246251565e-16
orthonormality test (should be close to 0): 1.4383797347324747e-15
Spectrum of projected Hessian of the Lagrangian: [-2.82842, -2.82842, -2.82842, -2.06111, -1.41421, 0.646902]
```

Notice that the projected Hessian is indefinite since it has both a negative and positive eigenvalue.

## Reproduce the figures and tables in the paper





##### Figures 2, 3, and 4.

```julia
include("paper_tables_and_figures.jl")
```

##### Table 1

```julia
include("paper_tables_and_figures.jl")
for k in [3, 4, 5]; summary_statistics("ngrams", k); end
for k in [3, 4, 5]; summary_statistics("DAWN", k); end
for k in [3, 4, 5]; summary_statistics("tags-ask-ubuntu", k); end
```

##### Tables 2, 3, and 4
