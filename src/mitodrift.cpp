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
double logSumExp(const arma::vec x) {
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


std::vector<std::vector<double>> armaMatToStdVec(const arma::mat mat) {
    std::vector<std::vector<double>> vec;
    for (size_t i = 0; i < mat.n_rows; ++i) {
        std::vector<double> row;
        for (size_t j = 0; j < mat.n_cols; ++j) {
            row.push_back(mat(i, j));
        }
        vec.push_back(row);
    }
    return vec;
}

std::vector<std::vector<double>> armaMatToVector(const arma::mat mat) {

    int nrows = mat.n_rows;
    int ncols = mat.n_cols;

    std::vector<std::vector<double>> result(nrows, std::vector<double>(ncols));

    // Copy data from arma::mat to std::vector
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            result[i][j] = mat(i, j);
        }
    }
    
    return result;
}

// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp(const arma::Mat<int> E, const arma::mat logP, const arma::mat logA) {

    int n = logP.n_cols; // Number of nodes
    int m = E.n_rows; // Number of edges
    int C = logP.n_rows; // Number of character states

    // Use a single contiguous memory block for log_messages
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);

    // Step 3: Process nodes in postorder as given by E
    for (size_t i = 0; i < m; i++) {
        int node = E(i, 1);
        int par = E(i, 0);
        // Compute message from node → parent in log-space
        for (int c = 0; c < C; c++) { // Parent state
            for (int c_child = 0; c_child < C; c_child++) { // Child state
                state_log_values[c_child] = logA(c, c_child) + logP(c_child, node) + log_messages[c_child * n + node];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }

    // Step 4: Compute log-partition function at root
    int root = E(m - 1, 0);
    for (int c = 0; c < C; c++) {
        state_log_values[c] = logP(c, root) + log_messages[c * n + root];
    }

    return logSumExp(state_log_values);
}

// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp_wrapper(arma::Mat<int> E, const arma::cube logP, const arma::mat logA) {
    
    int L = logP.n_slices; // Number of loci
    double logZ = 0;
    E = E - 1;
    
    // Loop over loci
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp(E, logP.slice(l), logA);
    }

    return(logZ);
}

// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp2(const arma::Mat<int> E, const arma::mat logP, const std::vector<std::vector<double>> logA) {

    int n = logP.n_cols; // Number of nodes
    int m = E.n_rows; // Number of edges
    int C = logP.n_rows; // Number of character states

    // Use a single contiguous memory block for log_messages
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);

    // Step 3: Process nodes in postorder as given by E
    for (size_t i = 0; i < m; i++) {
        int node = E(i, 1);
        int par = E(i, 0);
        // Compute message from node → parent in log-space
        for (int c = 0; c < C; c++) { // Parent state
            for (int c_child = 0; c_child < C; c_child++) { // Child state
                state_log_values[c_child] = logA[c][c_child] + logP(c_child, node) + log_messages[c_child * n + node];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }

    // Step 4: Compute log-partition function at root
    int root = E(m - 1, 0);
    for (int c = 0; c < C; c++) {
        state_log_values[c] = logP(c, root) + log_messages[c * n + root];
    }

    return logSumExp(state_log_values);
}


// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp_wrapper2(arma::Mat<int> E, const arma::cube logP, const arma::mat logA) {
    
    int L = logP.n_slices; // Number of loci
    double logZ = 0;
    E = E - 1;

    std::vector<std::vector<double>> logAA = armaMatToStdVec(logA);
    
    // Loop over loci
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp2(E, logP.slice(l), logAA);
    }

    return(logZ);
}


// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp3(const arma::Mat<int> E, const arma::mat logP, const std::vector<double> logA) {

    int n = logP.n_cols; // Number of nodes
    int m = E.n_rows; // Number of edges
    int C = logP.n_rows; // Number of character states

    // Use a single contiguous memory block for log_messages
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);

    // Step 3: Process nodes in postorder as given by E
    for (size_t i = 0; i < m; i++) {
        int node = E(i, 1);
        int par = E(i, 0);
        // Compute message from node → parent in log-space
        for (int c = 0; c < C; c++) { // Parent state
            for (int c_child = 0; c_child < C; c_child++) { // Child state
                state_log_values[c_child] = logA[c * C + c_child] + logP(c_child, node) + log_messages[c_child * n + node];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }

    // Step 4: Compute log-partition function at root
    int root = E(m - 1, 0);
    for (int c = 0; c < C; c++) {
        state_log_values[c] = logP(c, root) + log_messages[c * n + root];
    }

    return logSumExp(state_log_values);
}


// E is edge matrix of the tree; rows must be in postorder
// P is the likelihood matrix (character state x node); must be ordered according to node indices in E
// output: total log likelihood for the tree graph (partition function) via the sum product algorithm
//' @export
// [[Rcpp::export]]
double score_tree_bp_wrapper3(arma::Mat<int> E, const arma::cube logP, const std::vector<double> logA) {
    
    int L = logP.n_slices; // Number of loci
    double logZ = 0;
    E = E - 1;
    
    // Loop over loci
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp3(E, logP.slice(l), logA);
    }

    return(logZ);
}




// [[Rcpp::export]]
double score_tree_bp_multi(arma::Mat<int> E, const arma::cube logP, const arma::mat logA) {

    E = E - 1;

    int n = logP.n_cols;   // Number of nodes
    int m = E.n_rows;      // Number of edges
    int C = logP.n_rows;   // Number of character states
    int L = logP.n_slices; // Number of loci
    int root = E(m - 1, 0);

    // Output vector for logZ values across loci
    double logZ = 0;
    std::vector<double> state_log_values(C);
    std::vector<double> log_messages(C * n, 0.0);

    // Process each locus separately
    for (int l = 0; l < L; l++) {

        std::fill(log_messages.begin(), log_messages.end(), 0.0);

        for (size_t i = 0; i < m; i++) {
            int node = E(i, 1);
            int par = E(i, 0);
            // Compute message from node → parent in log-space
            for (int c = 0; c < C; c++) { // Parent state
                for (int c_child = 0; c_child < C; c_child++) { // Child state
                    state_log_values[c_child] = logA(c, c_child) + logP(c_child, node, l) + log_messages[c_child * n + node];
                }
                log_messages[c * n + par] += logSumExp(state_log_values);
            }
        }

        // Step 4: Compute log-partition function at root
        for (int c = 0; c < C; c++) {
            state_log_values[c] = logP(c, root, l) + log_messages[c * n + root];
        }
        logZ += logSumExp(state_log_values);
    }

    return logZ;
}



struct score_neighbours_max: public Worker {

    // original tree
    const arma::Mat<int> E;
    const arma::cube logP;
    const arma::mat logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_max(const arma::Mat<int> E, const arma::cube logP, const arma::mat logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Mat<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper2(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper2(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel(arma::Mat<int> E, const arma::cube logP, const arma::mat logA) {

    E = reorderRcpp(E);

    int n = E.n_rows/2 - 1;

    NumericVector scores(2*n);

    score_neighbours_max score_neighbours_max(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours_max);

    return scores;

}


