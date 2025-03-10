#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppParallel;

/////////////////////////////////////// NNI ////////////////////////////////////////

// An inline helper that performs the recursive postorder traversal.
// Instead of using Armadillo’s element access repeatedly, we use raw pointers
// to speed up the recursion.
inline void bar_reorderRcpp_inline(int node, int nTips,
    const int* e1_ptr, const int* e2_ptr,
    std::vector<int>& neworder,
    const int* L_ptr, const int* xi_ptr, const int* xj_ptr,
    int &iii) {

    int i = node - nTips - 1;
    // First, add the current block in reverse order.
    for (int j = xj_ptr[i] - 1; j >= 0; j--) {
        neworder[iii--] = L_ptr[xi_ptr[i] + j] + 1; // +1 for R's 1-indexing
    }
    // Then, recursively process each child.
    for (int j = 0; j < xj_ptr[i]; j++) {
        int k = e2_ptr[ L_ptr[xi_ptr[i] + j] ];
    if (k > nTips)
        bar_reorderRcpp_inline(k, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    }
}

// A vectorized version of reorder_rows using Armadillo’s built-in indexing.
arma::Mat<int> reorder_rows(const arma::Mat<int>& x, const std::vector<int>& order) {
    // Convert the std::vector to an arma::uvec.
    arma::uvec idx = arma::conv_to<arma::uvec>::from(order);
    // Adjust for 0-indexing (order was stored 1-indexed).
    idx -= 1;
    return x.rows(idx);
}

// [[Rcpp::export]]
arma::Mat<int> reorderRcpp(arma::Mat<int> E) {
    int n = E.n_rows;
    int nTips = n / 2 + 1;
    int root = nTips + 1;

    // Get copies of the parent and child columns.
    arma::Col<int> e1 = E.col(0);
    arma::Col<int> e2 = E.col(1);

    // The maximum parent label tells us the total number of nodes.
    int m = e1.max(); 
    int nnode = m - nTips;

    // Allocate working arrays.
    std::vector<int> L(n);            // Will store indices corresponding to each internal node.
    std::vector<int> neworder(n);       // The final reordering.
    std::vector<int> pos(nnode, 0);     // Current fill position for each node.
    std::vector<int> xi(nnode, 0);      // Starting index for each node in L.
    std::vector<int> xj(nnode, 0);      // Count of children per internal node.

    // First pass: count children per node (using e1, which stores parent IDs).
    for (int i = 0; i < n; i++) {
        int idx = e1[i] - nTips - 1;
        xj[idx]++;
    }

    // Compute starting positions xi (cumulative sums).
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i - 1] + xj[i - 1];
    }

    // Fill L: For each edge, place its row index in the appropriate block.
    for (int i = 0; i < n; i++) {
        int k = e1[i] - nTips - 1;
        int j = pos[k];
        L[xi[k] + j] = i;
        pos[k]++;
    }

    // Reset the new order index.
    int iii = n - 1;

    // Use raw pointers to speed up inner loops.
    const int* e1_ptr = e1.memptr();
    const int* e2_ptr = e2.memptr();
    const int* L_ptr   = L.data();
    const int* xi_ptr  = xi.data();
    const int* xj_ptr  = xj.data();

    // Run the recursive postorder traversal.
    bar_reorderRcpp_inline(root, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);

    // Reorder the rows of E using the computed order.
    E = reorder_rows(E, neworder);

    return E;
}

// [[Rcpp::export]]
std::vector<arma::Mat<int>> nnin_cpp(const arma::Mat<int> E, const int n) {
    // Make copies for the two alternative topologies.
    arma::Mat<int> E1 = E;
    arma::Mat<int> E2 = E;
    
    // Get the parent and child columns.
    arma::Col<int> parent = E.col(0);
    arma::Col<int> child  = E.col(1);
    
    // Get raw pointers to data for speed.
    const int* p_parent = parent.memptr();
    const int* p_child  = child.memptr();
    
    int numEdges = child.n_elem;
    int nTips = numEdges / 2 + 1;
    
    // Find the nth internal edge
    int count = 0, ind = -1;
    for (int i = 0; i < numEdges; i++) {
        if (p_child[i] > nTips) {
            count++;
            if (count == n) {
                ind = i;
                break;
            }
        }
    }

    if (ind < 0) {
        Rcpp::stop("n is larger than the number of valid edges.");
    }
    
    int p1 = p_parent[ind];
    int p2 = p_child[ind];
    
    // Find first index (other than ind) where parent equals p1.
    int ind1 = -1;
    for (int i = 0; i < numEdges; i++) {
        if (i != ind && p_parent[i] == p1) {
            ind1 = i;
            break;
        }
    }

    if (ind1 < 0) {
        Rcpp::stop("No valid index found for p1 interchange.");
    }
    
    // Find two indices where parent equals p2.
    int ind2_0 = -1, ind2_1 = -1;
    for (int i = 0; i < numEdges; i++) {
        if (p_parent[i] == p2) {
            if (ind2_0 < 0) {
                ind2_0 = i;
            } else if (ind2_1 < 0 && i != ind) {
                ind2_1 = i;
                break;
            }
        }
    }

    if (ind2_0 < 0 || ind2_1 < 0) {
        Rcpp::stop("No valid indices found for p2 interchange.");
    }
    
    // Retrieve children values for the swapping.
    int e1 = p_child[ind1];
    int e2 = p_child[ind2_0];
    int e3 = p_child[ind2_1];
    
    // Perform the nearest neighbor interchanges:
    // In the first topology, swap child at ind1 with that at ind2_0.
    E1(ind1, 1)   = e2;
    E1(ind2_0, 1) = e1;
    
    // In the second topology, swap child at ind1 with that at ind2_1.
    E2(ind1, 1)   = e3;
    E2(ind2_1, 1) = e1;
    
    std::vector<arma::Mat<int>> res(2);
    // Assume reorderRcpp is a function that reorders the tree appropriately.
    res[0] = reorderRcpp(E1);
    res[1] = reorderRcpp(E2);
    
    return res;
}


/////////////////////////////////////// MitoDrfit ////////////////////////////////////////

//' definitions for logSumExp function
// https://github.com/helske/seqHMM/blob/master/src/logSumExp.cpp
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
    unsigned int maxi = x.index_max();
    LDOUBLE maxv = x(maxi);
    if (!(maxv > -arma::datum::inf)) {
        return -arma::datum::inf;
    }
    LDOUBLE cumsum = 1.0; // Include the max element's contribution (exp(0)=1)
    for (unsigned int i = 0; i < maxi; i++) {
        cumsum += EXPL(x(i) - maxv);
    }
    for (unsigned int i = maxi + 1; i < x.n_elem; i++) {
        cumsum += EXPL(x(i) - maxv);
    }
    return maxv + std::log(cumsum);
}


// bp: Belief-propagation function.
// logP is a flattened likelihood matrix (row-major; dimensions: C x n)
// logA is a flattened transition matrix (dimensions: C x C)
// E is a flattened edge list in postorder (each edge: parent, child)
// n: number of nodes
// C: number of states
// m: number of edges (E.size() / 2)
// root: index of the root node
// [[Rcpp::export]]
double score_tree_bp(std::vector<int> E,
                        std::vector<double> logP,
                        std::vector<double> logA,
                        int n, int C, int m, int root) {
    // Allocate memory for messages and temporary state values.
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);

    // Process nodes in postorder; E is assumed to be provided in postorder.
    for (int i = 0; i < m; i++) {
        int par  = E[2 * i];     // Parent node
        int node = E[2 * i + 1]; // Child node

        // For each possible state at the parent.
        for (int c = 0; c < C; c++) {
            // For each possible state at the child.
            for (int c_child = 0; c_child < C; c_child++) {
            // For row-major ordering, element (c_child, node) is at index: c_child * n + node.
            state_log_values[c_child] = logA[c * C + c_child] +
                                        logP[c_child * n + node] +
                                        log_messages[c_child * n + node];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }

    // Compute the final log-partition function at the root.
    for (int c = 0; c < C; c++) {
        state_log_values[c] = logP[c * n + root] + log_messages[c * n + root];
    }
    return logSumExp(state_log_values);
}

// wrapper: Wrapper function that precomputes dimensions and passes them to score_tree_bp3.
// E: Edge matrix (each row is a (parent, child) pair, 1-indexed from R)
// logP_list: List of flattened likelihood matrices (each in row-major order)
// logA: Flattened transition matrix (row-major order)
// Computes:
//   L: number of loci,
//   C: number of states (inferred from logA),
//   n: number of nodes (inferred from first likelihood matrix),
//   m: number of edges,
//   root: parent of the last edge in E.
// [[Rcpp::export]]
double score_tree_bp_wrapper(arma::Mat<int> E,
                              std::vector< std::vector<double> > logP_list,
                              std::vector<double> logA) {
    // Compute number of loci.
    int L = logP_list.size();
    // Infer number of states from the transition matrix.
    int C = std::sqrt(logA.size());
    // Infer number of nodes from the first likelihood matrix.
    int n = logP_list[0].size() / C;

    // Adjust the edge matrix from 1-indexing (R) to 0-indexing (C++).
    E = E - 1;

    // Flatten E into a 1D vector where each row contributes (parent, child).
    std::vector<int> E_vec;
    E_vec.reserve(E.n_rows * 2);
    for (int i = 0; i < E.n_rows; i++) {
        E_vec.push_back(E(i, 0));  // Parent
        E_vec.push_back(E(i, 1));  // Child
    }

    // Compute m once (number of edges) and root (parent from the last edge).
    int m = E_vec.size() / 2;
    int root = E_vec[2 * (m - 1)];

    double logZ = 0.0;
    // Loop over loci; each element in logP_list is a flattened likelihood matrix for one locus.
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp(E_vec, logP_list[l], logA, n, C, m, root);
    }
    return logZ;
}

struct score_neighbours: public Worker {

    // original tree
    const arma::Mat<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<double> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours(const arma::Mat<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Mat<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel(arma::Mat<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA) {

    E = reorderRcpp(E);

    int n = E.n_rows/2 - 1;

    NumericVector scores(2*n);

    score_neighbours score_neighbours(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours);

    return scores;

}


