#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppParallel;

/////////////////////////////////////// NNI ////////////////////////////////////////

// Below functions are modified from R-package `ape' by Emmanuel Paradis and Klaus Schliep
// Functions are made thread-safe using RcppArmadillo. 

// https://github.com/KlausVigo/phangorn/blob/master/src/phangorn_utils.cpp
// [[Rcpp::export]]
std::vector<std::vector<int>> allChildrenCPP(const arma::Mat<int> E) {

    arma::Col<int> parent = E.col(0);
    arma::Col<int> children = E.col(1);
    int m = max(parent);

    std::vector<std::vector<int>> out(m);

    for(int i = 0; i<parent.size(); i++) {
        out[parent(i)-1L].push_back(children(i));
    }

    return out;
}

void bar_reorderRcpp(int node, int nTips, const arma::Col<int> & e1,
    const arma::Col<int> & e2, std::vector<int> & neworder, const arma::Col<int> & L,
    const arma::Col<int> & xi, const arma::Col<int> & xj, int & iii)
{
    int i = node - nTips - 1, j, k;

    for (j = xj[i] -1; j >= 0; j--)
        neworder[iii--] = L[xi[i] + j ] + 1;

    for (j = 0; j < xj[i]; j++) {
        k = e2[L[xi[i] + j ]];
        if (k > nTips)
            bar_reorderRcpp(k, nTips, e1, e2, neworder, L, xi, xj, iii);
    }
}

// [[Rcpp::export]]
arma::Mat<int> reorder_rows(arma::Mat<int> x, arma::Col<int> y) {

    // Create an output matrix
    arma::Mat<int> out = x;

    // Loop through each row and copy the data. 
    for (int i = 0; i < y.n_elem; ++i) {
        out.row(i) = x.row(y[i]-1);
    }

    return out;
}

// [[Rcpp::export]]
arma::Mat<int> reorderRcpp(arma::Mat<int> E) {

    int n = E.n_rows;
    int nTips = n/2 + 1;
    int root = nTips + 1;

    arma::Col<int> e1 = E.col(0);
    arma::Col<int> e2 = E.col(1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
    
    arma::Col<int> L(n);
    std::vector<int> neworder(n);
    arma::Col<int> pos(nnode);
    arma::Col<int> xi(nnode);
    arma::Col<int> xj(nnode);
    for (int i = 0; i < n; i++) {
        xj[e1[i] - nTips - 1]++;
    }
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i-1] + xj[i - 1];
    }
    for (int i = 0; i < n; i++) {
        k = e1[i] - nTips - 1;
        j = pos[k]; /* the current 'column' position corresponding to k */
        L[xi[k] + j] = i;
        pos[k]++;
    }

    int iii = n - 1;

    bar_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj, iii);

    E = reorder_rows(E, neworder);

    return E;
}

// Modified from R-package `phangorn' by Klaus Schliep
// n goes from 1 to total number of edges
// [[Rcpp::export]]
std::vector<arma::Mat<int>> nnin_cpp(const arma::Mat<int> E, const int n) {

    arma::Mat<int> E1 = E;
    arma::Mat<int> E2 = E;
    arma::Col<int> parent = E.col(0);
    arma::Col<int> child = E.col(1);
    int k = min(parent) - 1;
    arma::uvec indvec = find(child > k);
    int ind = indvec[n-1];
    int p1 = parent[ind];
    int p2 = child[ind];
    arma::uvec ind1_vec = find(parent == p1);
    ind1_vec = ind1_vec.elem(find(ind1_vec != ind));
    int ind1 = ind1_vec[0];
    arma::uvec ind2 = find(parent == p2);
    
    int e1 = child[ind1];
    int e2 = child[ind2[0]];
    int e3 = child[ind2[1]];

    E1(ind1, 1) = e2;
    E1(ind2[0], 1) = e1;
    E2(ind1, 1) = e3;
    E2(ind2[1], 1) = e1;

    std::vector<arma::Mat<int>> res(2);

    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);

    return res;
}

// Based on https://github.com/cran/ape/blob/390386e67f9ff6cd8e6e523b7c43379a1551c565/src/plot_phylo.c
// [[Rcpp::export]]
NumericVector node_depth(int ntip, NumericVector e1, NumericVector e2,
        int nedge, NumericVector xx, int method)
/* method == 1: the node depths are proportional to the number of tips
   method == 2: the node depths are evenly spaced */
{

    int i;

    /* First set the coordinates for all tips */
    for (i = 0; i < ntip; i++) xx[i] = 1;

    /* Then compute recursively for the nodes; we assume `xx' has */
    /* been initialized with 0's which is true if it has been */
    /* created in R (the tree must be in pruningwise order) */
    if (method == 1) {
        for (i = 0; i < nedge; i++)
            xx[e1[i] - 1] = xx[e1[i] - 1] + xx[e2[i] - 1];
    } else { /* *method == 2 */
        for (i = 0; i < nedge; i++) {
            /* if a value > 0 has already been assigned to the ancestor
               node of this edge, check that the descendant node is not
               at the same level or more */
            if (xx[e1[i] - 1])
            if (xx[e1[i] - 1] >= xx[e2[i] - 1] + 1) continue;
            xx[e1[i] - 1] = xx[e2[i] - 1] + 1;
        }
    }
    return xx;
}


/////////////////////////////////////// MitoDrfit ////////////////////////////////////////

//' definitions for logSumExp function
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

//' logSumExp function for a vector
//'
//' @param x NumericVector
//' @return double logSumExp of x
//' @export
// [[Rcpp::export]]
double logSumExp(const arma::vec& x) {
    // https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
    unsigned int maxi = x.index_max();
    LDOUBLE maxv = x(maxi);
    if (!(maxv > -arma::datum::inf)) {
        return -arma::datum::inf;
    }
    LDOUBLE cumsum = 0.0;
    for (unsigned int i = 0; i < x.n_elem; i++) {
        if ((i != maxi) & (x(i) > -arma::datum::inf)) {
            cumsum += EXPL(x(i) - maxv);
        }
    }
  
    return maxv + log1p(cumsum);
}

// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
// [[Rcpp::export]]
double score_tree_bp(arma::Mat<int> E, const arma::mat& logP, const arma::mat logA) {

    E = E - 1;
    arma::mat tlogA = logA.t();

    int n = logP.n_cols; // Number of nodes
    int C = logP.n_rows; // Number of character states

    // Initialize messages
    arma::mat log_messages(C, n, arma::fill::zeros);

    // Step 2: Find the root node (assuming postorder)
    int root = E(E.n_rows - 1, 0);

    // Step 3: Process nodes in postorder as given by E
    for (size_t i = 0; i < E.n_rows; i++) {

        int node = E(i, 1);
        int par = E(i, 0);

        // Compute message from node → parent in log-space
        arma::vec log_sums(C, arma::fill::zeros);
        for (int c = 0; c < C; c++) { // parent state
            // arma::vec child_log_values(C);
            // for (int c_child = 0; c_child < C; c_child++) { // child state
            //     child_log_values(c_child) = logA(c, c_child) + logP(c_child, node) + log_messages(c_child, node);
            // }
            arma::vec child_log_values = tlogA.col(c) + logP.col(node) + log_messages.col(node);
            log_sums(c) = logSumExp(child_log_values);
        }

        log_messages.col(par) += log_sums; // Update parent message
    }

    // // Step 4: Compute log-partition function at root
    double log_partition = logSumExp(logP.col(root) + log_messages.col(root));

    return(log_partition);
}

// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
// [[Rcpp::export]]
arma::vec score_tree_bp_wrapper(arma::Mat<int> E, const arma::cube& logP, const arma::cube& logA) {
    
    int L = logP.n_slices; // Number of loci
    arma::vec logZ(L, arma::fill::zeros); // Store log-likelihood for each locus
    
    // Loop over loci
    for (int l = 0; l < L; l++) {
        logZ(l) = score_tree_bp(E, logP.slice(l), logA.slice(l));
    }

    return(logZ);
}

// Compute logSumExp over a specified axis (0 for columns, 1 for rows)
// [[Rcpp::export]]
arma::vec logSumExpMat(const arma::mat& X, int axis = 0) {
    if (axis == 0) {
        // Sum over columns (returns a row vector)
        arma::rowvec max_vals = max(X, 0);
        arma::rowvec sum_exp = sum(exp(X.each_row() - max_vals), 0);
        return (max_vals + log(sum_exp)).t(); // Convert row vector to column vector
    } else {
        // Sum over rows (returns a column vector)
        arma::colvec max_vals = max(X, 1);
        arma::colvec sum_exp = sum(exp(X.each_col() - max_vals), 1);
        return max_vals + log(sum_exp);
    }
}


// [[Rcpp::export]]
arma::vec score_tree_bp_multi(arma::Mat<int> E, const arma::cube& logP, const arma::cube& logA) {

    E = E - 1; // Convert to 0-based indexing

    int n = logP.n_cols;  // Number of nodes
    int C = logP.n_rows;  // Number of character states
    int L = logP.n_slices; // Number of loci

    arma::vec logZ(L, arma::fill::zeros); // Vector to store log-likelihoods for each locus

    // Initialize messages as a 3D cube (C × n × L)
    arma::cube log_messages(C, n, L, arma::fill::zeros);

    // Step 2: Find the root node (assuming postorder)
    int root = E(E.n_rows - 1, 0);

    // Step 3: Process nodes in postorder as given by E
    for (size_t i = 0; i < E.n_rows; i++) {

        int node = E(i, 1);
        int par = E(i, 0);

        // Compute message from node → parent in log-space for all loci
        for (int l = 0; l < L; l++) {
            arma::mat child_log_values(C, C, arma::fill::zeros);
            
            for (int c = 0; c < C; c++) { // Parent state
                for (int c_child = 0; c_child < C; c_child++) { // Child state
                    child_log_values(c, c_child) = logA(c, c_child, l) + logP(c_child, node, l) + log_messages(c_child, node, l);
                }
                log_messages(c, par, l) += logSumExp(child_log_values.row(c).t()); // Sum over child states
            }
        }

        // // Process all loci simultaneously
        // for (int l = 0; l < L; l++) {
        //     arma::mat child_log_values = logA.slice(l) + logP.slice(l).col(node) + log_messages.slice(l).col(node);
        //     log_messages.slice(l).col(par) += logSumExpMat(child_log_values, 1); // Sum over child states
        // }
    }

    // Step 4: Compute log-partition function at root for each locus
    for (int l = 0; l < L; l++) {
        logZ(l) = logSumExp(logP.slice(l).col(root) + log_messages.slice(l).col(root));
    }

    return logZ; // Return vector of log-likelihoods for all loci
}


struct score_neighbours_max: public Worker {

    // original tree
    const arma::Mat<int> E;
    const arma::cube logP;
    const arma::cube logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_max(const arma::Mat<int> E, const arma::cube logP, const arma::cube logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Mat<int>> Ep = nnin_cpp(E, i+1);
            arma::vec res1 = score_tree_bp_wrapper(Ep[0], logP, logA);
            arma::vec res2 = score_tree_bp_wrapper(Ep[1], logP, logA);
            scores[2*i] = sum(res1);
            scores[2*i+1] = sum(res2);
        }
    }
};

// Note that the tree has to be already in post-order
// [[Rcpp::export]]
NumericVector nni_cpp_parallel(const List tree, const arma::cube logP, const arma::cube logA) {
    
    arma::Mat<int> E = tree["edge"];

    E = reorderRcpp(E);

    int n = E.n_rows/2 - 1;

    NumericVector scores(2*n);

    score_neighbours_max score_neighbours_max(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours_max);

    return scores;

}

// doesn't work yet
// //' logSumExp function for a matrix with axis argument
// //'
// //' @param X arma::mat (C × n)
// //' @param dim int (0 = column-wise, 1 = row-wise)
// //' @return arma::rowvec or arma::colvec logSumExp along the specified dimension
// // [[Rcpp::export]]
// arma::vec logSumExpMat(const arma::mat& X, int dim = 0) {
//     if (dim == 0) {  // Column-wise
//         arma::rowvec maxv = max(X, 0);
//         arma::rowvec cumsum = sum(exp(X.each_row() - maxv), 0);
//         return maxv + log1p(cumsum);
//     } else if (dim == 1) {  // Row-wise
//         arma::colvec maxv = max(X, 1);
//         arma::colvec cumsum = sum(exp(X.each_col() - maxv), 1);
//         return maxv + log1p(cumsum);
//     } else {
//         Rcpp::stop("Invalid dimension for logSumExpMat. Use dim=0 for columns, dim=1 for rows.");
//     }
// }

// //' logSumExp function for a cube with axis argument
// //'
// //' @param X arma::cube (C × n × L)
// //' @param dim int (0 = over rows, 1 = over columns, 2 = over slices)
// //' @return arma::mat logSumExp along the specified dimension
// // [[Rcpp::export]]
// arma::mat logSumExpCube(const arma::cube& X, int dim = 0) {
//     if (dim == 0) {  // Sum over rows (C → 1)
//         arma::mat maxv = max(X, 0);
//         arma::mat cumsum = sum(exp(X.each_slice() - maxv), 0);
//         return maxv + log1p(cumsum);
//     } else if (dim == 1) {  // Sum over columns (n → 1)
//         arma::mat maxv = max(X, 1);
//         arma::mat cumsum = sum(exp(X.each_slice().t() - maxv.t()), 1).t();
//         return maxv + log1p(cumsum);
//     } else if (dim == 2) {  // Sum over slices (L → 1)
//         arma::mat maxv = max(X, 2);
//         arma::mat cumsum = sum(exp(X - maxv.each_slice()), 2);
//         return maxv + log1p(cumsum);
//     } else {
//         Rcpp::stop("Invalid dimension for logSumExpCube. Use dim=0 for rows, dim=1 for columns, dim=2 for slices.");
//     }
// }

// // [[Rcpp::export]]
// Rcpp::List score_tree_bp_multi(arma::Mat<int> E, const arma::cube logP, const arma::cube logA, bool report_beliefs = false) {

//     // Ensure E is in postorder
//     E = reorderRcpp(E);
//     E = E - 1; // Convert to 0-based indexing

//     int n = logP.n_cols;   // Number of nodes
//     int C = logP.n_rows;   // Number of character states
//     int L = logP.n_slices; // Number of loci

//     // Initialize messages (C × n × L)
//     arma::cube log_messages(C, n, L, arma::fill::zeros);
//     arma::Col<int> parent(n, arma::fill::value(-1));

//     // Step 1: Assign parents based on E (assuming E is in postorder)
//     for (size_t i = 0; i < E.n_rows; i++) {
//         int par = E(i, 0);
//         int node = E(i, 1);
//         parent[node] = par;
//     }

//     // Step 2: Find the root node (node that never appears in the second column)
//     int root = E(E.n_rows - 1, 0);

//     // Step 3: Process nodes in postorder as given by E
//     for (size_t i = 0; i < E.n_rows; i++) {
//         int node = E(i, 1);
//         int par = parent[node];

//         // Compute message from node → parent in log-space
//         arma::cube child_log_values(C, C, L, arma::fill::zeros);

//         // Broadcast logP + log_messages for each state transition (vectorized)
//         for (int c = 0; c < C; c++) {  // Parent state
//             child_log_values.row(c) = logA.slice(c) + logP.slice(node) + log_messages.slice(node);
//         }

//         // Compute log-sums over child states (loop over loci)
//         arma::mat log_sums(C, L, arma::fill::zeros);
//         for (int l = 0; l < L; l++) {
//             for (int c = 0; c < C; c++) {
//                 log_sums(c, l) = logSumExp(arma::vectorise(child_log_values.tube(c, 0, c, C - 1)));
//             }
//         }

//         // Update parent message (vectorized)
//         log_messages.slice(par) += log_sums;
//     }

//     // Step 4: Compute total log-partition function over all loci
//     double total_log_partition = 0.0;
//     for (int l = 0; l < L; l++) {
//         arma::vec root_log_values(C, arma::fill::zeros);
//         for (int c = 0; c < C; c++) {
//             root_log_values(c) = logP(c, root, l) + log_messages(c, root, l);
//         }
//         total_log_partition += logSumExp(root_log_values);
//     }

//     // Compute node beliefs if requested
//     arma::cube node_beliefs(C, n, L, arma::fill::zeros);
//     if (report_beliefs) {
//         node_beliefs = logP + log_messages; // Compute raw node beliefs
        
//         for (int l = 0; l < L; l++) {
//             for (int node = 0; node < n; node++) {
//                 node_beliefs.slice(l).col(node) -= logSumExp(node_beliefs.slice(l).col(node));
//             }
//         }
//     }

//     // Return log-partition function and (optionally) beliefs
//     if (report_beliefs) {
//         return Rcpp::List::create(
//             Rcpp::Named("logZ") = total_log_partition,
//             Rcpp::Named("node_beliefs") = node_beliefs
//         );
//     } else {
//         return Rcpp::List::create(
//             Rcpp::Named("logZ") = total_log_partition
//         );
//     }
// }

