# Hypergraph eigenvector centralities

This code and data repository accompanies the paper

- Three hypergraph eigenvector centralities. Austin R. Benson. Submitted, 2018.

All of the code is written in Julia.

For questions, please email Austin at arb@cs.cornell.edu.

### Data

The tags and DAWN datasets are in the `data/` directory. These were downloaded directly from
[here](http://www.cs.cornell.edu/~arb/data/tags-ask-ubuntu/index.html) and
[here](http://www.cs.cornell.edu/~arb/data/DAWN/index.html).

### Reproduce the figures and tables in the paper

##### Figures 2, 3, and 4.

```julia
include("paper_tables_and_figures.jl")
```

##### Table 1

```julia
include("paper_tables_and_figures.jl")
```

##### Tables 2, 3, and 4
