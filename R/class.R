#' MitoDrift R6 Class
#' 
#' A class for mitochondrial drift analysis with tree-based inference
#' 
#' @field tree_init Initial phylogenetic tree (ape::phylo)
#' @field A Transition matrix
#' @field leaf_likelihoods List of leaf likelihoods
#' @field logA Vector of log transition probabilities
#' @field logP List of log probabilities
#' @field model_params Named vector containing model parameters (eps, err, npop, ngen, k)
#' @field amat Alternative allele count matrix
#' @field dmat Total depth matrix  
#' @field vmat VAF matrix
#' @field mcmc_trace MCMC trace for parameter fitting
#' @field tree_list List of trees from optimization
#' @field tree_ml Maximum likelihood tree
#' @field tree_annot Annotated tree with clade frequencies
#' @field ml_trace_file File path for ML trace
#' @field mcmc_trace_file File path for MCMC trace
#' @field param_trace_file File path for parameter fitting trace
#' @export
MitoDrift <- R6::R6Class("MitoDrift",
    public = list(
        # Member variables
        tree_init = NULL,
        A = NULL,
        leaf_likelihoods = NULL,
        logA = NULL,
        logP = NULL,
        amat = NULL,
        dmat = NULL,
        vmat = NULL,
        model_params = c('k' = 20, 'npop' = 600),
        mcmc_trace = NULL,
        tree_list = NULL,
        tree_ml = NULL,
        tree_annot = NULL,
        ml_trace_file = NULL,
        mcmc_trace_file = NULL,
        param_trace_file = NULL,
        
        #' @description Initialize MitoDrift object
        #' 
        #' @param mut_dat Mutation data in long format (data.frame with columns: variant, cell, d, a)
        #' @param amat Alternative allele count matrix (optional, if mut_dat not provided)
        #' @param dmat Total depth matrix (optional, if mut_dat not provided)
        #' @param model_params Model parameters (optional)
        initialize = function(
                mut_dat = NULL, amat = NULL, dmat = NULL, 
                model_params = NULL
            ) {
            
            # Set model parameters if complete
            if (!is.null(model_params) && all(c("eps", "err", "npop", "ngen", "k") %in% names(model_params))) {
                self$model_params <- model_params
            }

            # Process input data
            if (!is.null(mut_dat)) {
                # Convert long format to matrices
                if (!all(c("variant", "cell", "d", "a") %in% colnames(mut_dat))) {
                    stop("mut_dat must contain columns: variant, cell, d, a")
                }
                
                # Add VAF column if not present
                if (!"vaf" %in% colnames(mut_dat)) {
                    mut_dat$vaf <- mut_dat$a / mut_dat$d
                }
                
                # Convert to matrices
                self$vmat <- long_to_mat(mut_dat, 'vaf')
                self$dmat <- long_to_mat(mut_dat, 'd')
                self$amat <- long_to_mat(mut_dat, 'a')
                
            } else if (!is.null(amat) && !is.null(dmat)) {
                # Use provided matrices
                if (!identical(dim(amat), dim(dmat))) {
                    stop("amat and dmat must have the same dimensions")
                }
                self$amat <- amat
                self$dmat <- dmat
                self$vmat <- self$amat / self$dmat
            } else {
                stop("Either mut_dat or both amat and dmat must be provided")
            }
            
            # Check data validity
            n_muts <- nrow(self$vmat)
            n_cells <- ncol(self$vmat)
            
            if (n_cells == 0 || n_muts == 0) {
                stop('Zero cells or no shared mutations')
            }
            
            message('Building initial tree...')
            self$tree_init <- make_rooted_nj(self$vmat) %>% reorder_phylo()

        },
        
        #' @description Print object summary
        #' 
        #' Prints a summary of the MitoDrift object including data dimensions,
        #' tree information, model parameters, and available components.
        #' 
        #' @return None (prints to console)
        print = function() {
            cat("MitoDrift object\n")
            cat("================\n")
            cat("Data matrix: ", nrow(self$vmat), " mutations x ", ncol(self$vmat), " cells\n")
            cat("Initial tree: ", length(self$tree_init$tip.label), " tips, ", 
                self$tree_init$Nnode, " internal nodes\n")
            if (!is.null(self$A)) {
                cat("Transition matrix: ", nrow(self$A), "x", ncol(self$A), "\n")
            }
            if (!is.null(self$model_params)) {
                cat("Model parameters:\n")
                for (param in names(self$model_params)) {
                    cat("  ", param, ": ", self$model_params[param], "\n")
                }
            }
            if (!is.null(self$tree_ml)) {
                cat("ML tree: Available\n")
            }
            if (!is.null(self$tree_annot)) {
                cat("Annotated tree with clade frequencies: Available\n")
            }
            if (!is.null(self$param_trace_file)) {
                cat("Parameter fitting trace: ", self$param_trace_file, "\n")
            }
            if (!is.null(self$ml_trace_file)) {
                cat("ML tree optimization trace: ", self$ml_trace_file, "\n")
            }
            if (!is.null(self$mcmc_trace_file)) {
                cat("Phylogenetic MCMC trace: ", self$mcmc_trace_file, "\n")
            }
        },
        
        #' @description Make model components
        #' @param model_params Model parameters (optional)
        #' @param ncores Number of cores to use
        #' @return None
        make_model = function(model_params = NULL, ncores = 1) {

            if (!is.null(model_params)) {
                self$model_params <- model_params
            }

            if (!all(c("eps", "err", "npop", "ngen", "k") %in% names(self$model_params))) {
                stop("Model parameters must contain eps, err, npop, ngen, k")
            }

            # Get transition matrix
            self$A <- get_transition_mat_wf_hmm(
                k = self$model_params["k"], 
                eps = self$model_params["eps"], 
                N = self$model_params["npop"], 
                ngen = self$model_params["ngen"]
            )
            
            # Get leaf likelihoods
            self$leaf_likelihoods <- get_leaf_liks_mat(
                self$amat, 
                self$dmat, 
                get_vaf_bins(k = self$model_params["k"]), 
                eps = self$model_params["err"], 
                log = TRUE,
                ncores = ncores
            )
            
            # Convert to log probabilities
            self$logA <- as.vector(t(log(self$A)))
            self$logP <- convert_logliks_to_logP_list(self$leaf_likelihoods, self$tree_init)
        },
        
        #' @description Fit tree parameters using MCMC
        #' @param nsteps Number of MCMC steps (default: 500)
        #' @param nchains Number of MCMC chains (default: 1)
        #' @param ncores Number of cores to use (default: 1)
        #' @param keep Number of samples to keep at the end of the chain (default: 100)
        #' @param outfile Output file for MCMC results (optional)
        #' @param check_conv whether to check convergence of parameter fitting
        #' @return MCMC result object
        fit_params_mcmc = function(
            nsteps = 500,
            nchains = 1,
            ncores = 1,
            keep = 100,
            check_conv = FALSE,
            outfile = NULL
        ) {
            
            if (!is.null(outfile)) {
                self$param_trace_file <- normalizePath(outfile, mustWork = FALSE)
            }
            
            message("Fitting tree parameters using MCMC...")
            
            # Run MCMC parameter fitting
            mcmc_trace <- fit_params_mcmc(
                tree_fit = self$tree_init,
                amat = self$amat,
                dmat = self$dmat,
                nsteps = nsteps,
                nchains = nchains,
                ncores = ncores,
                outfile = self$param_trace_file,
                npop = self$model_params["npop"],
                k = self$model_params["k"],
                check_conv = check_conv
            )

            params_est = mcmc_trace %>%
                group_by(chain, variable) %>%
                slice_tail(n = keep) %>%
                ungroup() %>%
                group_by(variable) %>%
                summarise(est = mean(value), .groups = 'drop') %>%
                {setNames(.$est, .$variable)}
                        
            # Update model parameters with fitted values
            for (param in names(params_est)) {
                self$model_params[param] <- params_est[param]
            }

            message("MCMC parameter fitting completed!")
            message("Fitted parameters: ", paste(names(params_est), params_est, sep = " = ", collapse = ", "))
            if (!is.null(self$param_trace_file)) {
                message("Results saved to ", self$param_trace_file)
            }
            
            return(params_est)
        },
        
        #' @description Fit tree parameters using Expectation-Maximization (EM) algorithm
        #' @param initial_params Initial parameter values (ngen, log_eps, log_err)
        #' @param lower_bounds Lower bounds for parameters in transformed space
        #' @param upper_bounds Upper bounds for parameters in transformed space
        #' @param max_iter Maximum number of EM iterations (default: 10)
        #' @param epsilon Convergence threshold (default: 1e-3)
        #' @param ncores Number of cores to use (default: 3)
        #' @param trace Whether to return trace of parameter values (default: TRUE)
        #' @return Fitted parameters or list of parameters and trace
        fit_params_em = function(
            initial_params = c('ngen' = 100, 'log_eps' = log(1e-3), 'log_err' = log(1e-3)),
            lower_bounds = c('ngen' = 1, 'log_eps' = log(1e-12), 'log_err' = log(1e-12)),
            upper_bounds = c('ngen' = 1000, 'log_eps' = log(0.2), 'log_err' = log(0.2)),
            max_iter = 10,
            epsilon = 1e-3,
            ncores = 1
        ) {
            
            message("Fitting tree parameters using EM with ", ncores, " cores")
            
            # Run EM parameter fitting using the existing function
            params_est <- fit_params_em(
                tree_fit = self$tree_init,
                amat = self$amat,
                dmat = self$dmat,
                initial_params = initial_params,
                lower_bounds = lower_bounds,
                upper_bounds = upper_bounds,
                max_iter = max_iter,
                k = self$model_params["k"],
                npop = self$model_params["npop"],
                ncores = ncores,
                epsilon = epsilon,
                trace = FALSE
            )

            # Update model parameters with fitted values
            for (param in c('ngen', 'eps', 'err')) {
                self$model_params[param] <- params_est[param]
            }
            
            message("EM parameter fitting completed!")
            message("Fitted parameters: ", paste(names(params_est), params_est, sep = " = ", collapse = ", "))
            
            return(params_est)
        },
        
        
        #' @description Optimize tree using C++ implementation
        #' @param max_iter Maximum number of iterations (default: 100)
        #' @param ncores Number of cores to use (default: 1)
        #' @param outfile Output file for saving results (optional)
        #' @param resume Whether to resume from existing file (default: FALSE)
        #' @param trace_interval Interval for saving trace (default: 5)
        #' @return Optimization result object
        optimize_tree = function(
            max_iter = 100,
            ncores = 1,
            outfile = NULL,
            resume = FALSE,
            trace_interval = 5
        ) {
            
            # Check if model components are initialized
            if (is.null(self$model_params) || is.null(self$logP) || is.null(self$logA)) {
                stop("Model components must be initialized. Use make_model() first.")
            }

            if (!is.null(outfile)) {
                self$ml_trace_file <- normalizePath(outfile, mustWork = FALSE)
            } else {
                if (resume & is.null(self$ml_trace_file)) {
                    stop("Cannot resume from non-existent file. Provide a valid ml_trace_file.")
                }
            }
                
            tree_list <- optimize_tree_cpp(
                tree_init = self$tree_init,
                logP = self$logP,
                logA = self$logA,
                max_iter = max_iter,
                ncores = ncores,
                resume = resume,
                outfile = self$ml_trace_file,
                trace_interval = trace_interval
            )

            self$tree_ml = tree_list %>% .[[length(.)]]

            message("ML tree optimization completed!")
            if (!is.null(self$ml_trace_file)) {
                message("Results saved to ", self$ml_trace_file)
            }
        },
        
        #' @description Run phylogenetic MCMC sampling
        #' @param max_iter Maximum number of MCMC iterations (default: 10000)
        #' @param nchains Number of MCMC chains (default: 1000)
        #' @param ncores Number of cores to use (default: 1)
        #' @param outfile Output file for saving results (optional)
        #' @param resume Whether to resume from existing file (default: FALSE)
        #' @param use_nj Whether to use NJ tree instead of ML tree as initial tree (default: FALSE)
        #' @return MCMC result object
        run_mcmc = function(
            max_iter = 10000,
            nchains = 1000,
            ncores = 1,
            batch_size = 1000,
            outfile = NULL,
            resume = FALSE,
            use_nj = FALSE
        ) {

            if (!is.null(outfile)) {
                self$mcmc_trace_file <- normalizePath(outfile, mustWork = FALSE)
            }
            
            # Check if model components are initialized
            if (is.null(self$model_params) || is.null(self$logP) || is.null(self$logA)) {
                stop("Model components must be initialized. Use make_model() first.")
            }
            
            # Select initial tree
            if (use_nj) {
                message('Using NJ tree as initial tree')
                tree_init_mcmc <- self$tree_init
            } else {
                # Check if tree optimization has been run
                if (is.null(self$tree_ml)) {
                    stop("Tree optimization must be run first. Use optimize_tree() first.")
                }
                message('Using ML tree as initial tree')
                tree_init_mcmc <- self$tree_ml
            }
            
            message('Running phylogenetic MCMC with batch size ', batch_size, '...')
            
            # Run MCMC
            run_tree_mcmc_batch(
                phy_init = tree_init_mcmc,
                logP_list = self$logP,
                logA_vec = self$logA,
                max_iter = max_iter,
                nchains = nchains,
                ncores = ncores,
                outfile = self$mcmc_trace_file,
                resume = resume,
                batch_size = batch_size
            )
            
            message('Phylogenetic MCMC completed!')
            if (!is.null(self$mcmc_trace_file)) {
                message('Results saved to ', self$mcmc_trace_file)
            }
        },
        
        #' @description Trim tree using MCMC results
        #' @param burnin Number of burnin iterations (default: 0)
        #' @param max_iter Maximum iteration to include (default: 1e8)
        #' @param use_nj Whether to use NJ tree instead of ML tree (default: FALSE)
        #' @param mcmc_trace_file MCMC result file (optional, uses stored file if NULL)
        #' @param ncores Number of cores to use (default: 1)
        #' @return Trimmed tree with clade frequencies
        annotate_tree = function(
            burnin = 0,
            max_iter = 1e8,
            use_nj = FALSE,
            ncores = 1,
            mcmc_trace_file = NULL
        ) {
            
            # Select base tree
            if (use_nj) {
                message('Using NJ tree as base tree')
                tree <- self$tree_init
            } else {
                message('Using ML tree as base tree')
                if (is.null(self$tree_ml)) {
                    stop("ML tree optimization must be run first. Use optimize_tree() first.")
                }
                tree <- self$tree_ml
            }
            
            # Determine MCMC file to use
            if (is.null(mcmc_trace_file)) {
                if (is.null(self$mcmc_trace_file)) {
                    stop("MCMC file must be provided or MCMC must be run first")
                }
                mcmc_trace_file <- self$mcmc_trace_file
            }
            
            message('Loading MCMC results from ', mcmc_trace_file)
            res_mcmc <- readRDS(mcmc_trace_file)
            
            message('Collecting MCMC chains...')
            trees_mcmc <- collect_chains(res_mcmc, burnin = burnin, max_iter = max_iter)
            
            message('Adding clade frequencies to tree...')
            self$tree_annot <- add_clade_freq(tree, trees_mcmc, ncores = ncores)
            
            message('Tree annotation completed!')
        },
        
        #' @description Create a copy of a MitoDrift object with updated structure
        #' 
        #' This method creates a new MitoDrift instance with the same data but updated
        #' class structure. Useful for refreshing objects created with older versions
        #' of the code.
        #' 
        #' @param old_obj The old MitoDrift object to copy from
        #' @return A new MitoDrift object with the same data
        copy = function(old_obj) {
            
            # Create new instance with same data
            new_obj <- MitoDrift$new(
                amat = old_obj$amat,
                dmat = old_obj$dmat,
                model_params = old_obj$model_params
            )
            
            # List of all fields to copy
            all_fields <- c("tree_init", "A", "leaf_likelihoods", "logA", "logP", 
                           "tree_ml", "tree_annot", "mcmc_trace", "tree_list",
                           "ml_trace_file", "mcmc_trace_file", "param_trace_file")
            
            # Copy over all fields
            for (field in all_fields) {
                if (!is.null(old_obj[[field]])) {
                    new_obj[[field]] <- old_obj[[field]]
                }
            }
            
            message("MitoDrift object copied with updated structure")
            return(new_obj)
        }
    ),
    
    private = list(
        # Private methods can be added here if needed
    )
) 