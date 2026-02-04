# MitoDrift R6 Class

MitoDrift R6 Class

MitoDrift R6 Class

## Details

A class for mitochondrial drift analysis with tree-based inference

## Public fields

- `tree_init`:

  Initial phylogenetic tree (ape::phylo)

- `A`:

  Transition matrix

- `leaf_likelihoods`:

  List of leaf likelihoods

- `logA`:

  Vector of log transition probabilities

- `logP`:

  List of log probabilities

- `model_params`:

  Named vector containing model parameters (eps, err, npop, ngen, k)

- `amat`:

  Alternative allele count matrix

- `dmat`:

  Total depth matrix

- `vmat`:

  VAF matrix

- `mcmc_trace`:

  MCMC trace for parameter fitting

- `tree_list`:

  List of trees from optimization

- `tree_ml`:

  Maximum likelihood tree

- `tree_annot`:

  Annotated tree with clade frequencies

- `ml_trace_file`:

  File path for ML trace

- `mcmc_trace_file`:

  File path for MCMC trace

- `mcmc_diag_file`:

  File path for MCMC diagnostics

- `param_trace_file`:

  File path for parameter fitting trace

## Methods

### Public methods

- [`MitoDrift$new()`](#method-MitoDrift-new)

- [`MitoDrift$print()`](#method-MitoDrift-print)

- [`MitoDrift$make_model()`](#method-MitoDrift-make_model)

- [`MitoDrift$fit_params_em()`](#method-MitoDrift-fit_params_em)

- [`MitoDrift$optimize_tree()`](#method-MitoDrift-optimize_tree)

- [`MitoDrift$run_mcmc()`](#method-MitoDrift-run_mcmc)

- [`MitoDrift$annotate_tree()`](#method-MitoDrift-annotate_tree)

- [`MitoDrift$copy()`](#method-MitoDrift-copy)

- [`MitoDrift$clone()`](#method-MitoDrift-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize MitoDrift object

#### Usage

    MitoDrift$new(
      mut_dat = NULL,
      amat = NULL,
      dmat = NULL,
      model_params = NULL,
      build_tree = TRUE,
      ncores = 1
    )

#### Arguments

- `mut_dat`:

  Mutation data in long format (data.frame with columns: variant, cell,
  d, a)

- `amat`:

  Alternative allele count matrix (optional, if mut_dat not provided)

- `dmat`:

  Total depth matrix (optional, if mut_dat not provided)

- `model_params`:

  Model parameters (optional)

- `build_tree`:

  Whether to build an initial NJ tree (default: TRUE)

- `ncores`:

  Number of cores used when building the initial NJ tree (default: 1)

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print object summary

Prints a summary of the MitoDrift object including data dimensions, tree
information, model parameters, and available components.

#### Usage

    MitoDrift$print()

#### Returns

None (prints to console)

------------------------------------------------------------------------

### Method `make_model()`

Make model components

#### Usage

    MitoDrift$make_model(model_params = NULL, ncores = 1)

#### Arguments

- `model_params`:

  Model parameters (optional)

- `ncores`:

  Number of cores to use

#### Returns

None

------------------------------------------------------------------------

### Method `fit_params_em()`

Fit tree parameters using Expectation-Maximization (EM) algorithm

#### Usage

    MitoDrift$fit_params_em(
      tree_fit,
      initial_params = c(ngen = 100, log_eps = log(0.001), log_err = log(0.001)),
      lower_bounds = c(ngen = 1, log_eps = log(1e-12), log_err = log(1e-12)),
      upper_bounds = c(ngen = 1000, log_eps = log(0.2), log_err = log(0.2)),
      max_iter = 10,
      epsilon = 0.001,
      ncores = 1
    )

#### Arguments

- `tree_fit`:

  Tree to fit parameters to

- `initial_params`:

  Initial parameter values (ngen, log_eps, log_err)

- `lower_bounds`:

  Lower bounds for parameters in transformed space

- `upper_bounds`:

  Upper bounds for parameters in transformed space

- `max_iter`:

  Maximum number of EM iterations (default: 10)

- `epsilon`:

  Convergence threshold (default: 1e-3)

- `ncores`:

  Number of cores to use (default: 3)

- `trace`:

  Whether to return trace of parameter values (default: TRUE)

#### Returns

Fitted parameters or list of parameters and trace

------------------------------------------------------------------------

### Method `optimize_tree()`

Optimize tree using C++ implementation

#### Usage

    MitoDrift$optimize_tree(
      max_iter = 100,
      ncores = 1,
      outfile = NULL,
      resume = FALSE,
      trace_interval = 5
    )

#### Arguments

- `max_iter`:

  Maximum number of iterations (default: 100)

- `ncores`:

  Number of cores to use (default: 1)

- `outfile`:

  Output file for saving results (optional)

- `resume`:

  Whether to resume from existing file (default: FALSE)

- `trace_interval`:

  Interval for saving trace (default: 5)

#### Returns

Optimization result object

------------------------------------------------------------------------

### Method `run_mcmc()`

Run phylogenetic MCMC sampling

#### Usage

    MitoDrift$run_mcmc(
      max_iter = 10000,
      conv_thres = NULL,
      nchains = 1000,
      ncores = 1,
      ncores_qs = 1,
      batch_size = 1000,
      outfile = NULL,
      resume = FALSE,
      diag = TRUE,
      use_nj = FALSE,
      diagfile = NULL
    )

#### Arguments

- `max_iter`:

  Maximum number of MCMC iterations (default: 10000)

- `conv_thres`:

  Convergence threshold for ASDSF (default: NULL)

- `nchains`:

  Number of MCMC chains (default: 1000)

- `ncores`:

  Number of cores to use (default: 1)

- `ncores_qs`:

  Number of cores to use for QS operations (default: 1)

- `batch_size`:

  Batch size for MCMC (default: 1000)

- `outfile`:

  Output file for saving results (optional)

- `resume`:

  Whether to resume from existing file (default: FALSE)

- `diag`:

  Whether to compute ASDSF diagnostics (default: TRUE)

- `use_nj`:

  Whether to use NJ tree instead of ML tree as initial tree (default:
  FALSE)

- `diagfile`:

  File path for saving ASDSF diagnostics (RDS)

#### Returns

MCMC result object

------------------------------------------------------------------------

### Method `annotate_tree()`

Trim tree using MCMC results

#### Usage

    MitoDrift$annotate_tree(
      burnin = 0,
      max_iter = 1e+08,
      use_nj = FALSE,
      ncores = 1,
      ncores_qs = 1,
      mcmc_trace_file = NULL
    )

#### Arguments

- `burnin`:

  Number of burnin iterations (default: 0)

- `max_iter`:

  Maximum iteration to include (default: 1e8)

- `use_nj`:

  Whether to use NJ tree instead of ML tree (default: FALSE)

- `ncores`:

  Number of cores to use (default: 1)

- `ncores_qs`:

  Number of cores to use for QS operations (default: 1)

- `mcmc_trace_file`:

  MCMC result file (optional, uses stored file if NULL)

#### Returns

Trimmed tree with clade frequencies

------------------------------------------------------------------------

### Method `copy()`

Create a copy of a MitoDrift object with updated structure

This method creates a new MitoDrift instance with the same data but
updated class structure. Useful for refreshing objects created with
older versions of the code.

#### Usage

    MitoDrift$copy(old_obj, rebuild_tree = FALSE)

#### Arguments

- `old_obj`:

  The old MitoDrift object to copy from

- `rebuild_tree`:

  Whether to rebuild the initial tree from data (default: FALSE)

#### Returns

A new MitoDrift object with the same data

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MitoDrift$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
