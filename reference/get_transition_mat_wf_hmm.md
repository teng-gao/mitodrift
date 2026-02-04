# Get the transition matrix for WF model with HMM (with caching) TODO: add log option for small probabilities

Get the transition matrix for WF model with HMM (with caching) TODO: add
log option for small probabilities

## Usage

``` r
get_transition_mat_wf_hmm(k, eps, N, ngen, safe = FALSE)
```

## Arguments

- k:

  number of VAF bins

- eps:

  Variant detection error rate

- N:

  population size

- ngen:

  number of generations

- safe:

  whether to add small probability to avoid 0s

## Value

transition matrix
