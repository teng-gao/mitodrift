#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppParallel;

/////////////////////////////////////// NNI ////////////////////////////////////////

// Inline helper performing the recursive postorder traversal.
// Parameters:
//   node, nTips: as before.
//   e1_ptr and e2_ptr: raw pointers to the parent's and child's arrays.
//   neworder: vector (of length nEdges) that will be filled with the new row ordering (stored 1-indexed).
//   L_ptr, xi_ptr, xj_ptr: auxiliary arrays computed for the internal nodes.
//   iii: current index to fill in neworder (initialized to nEdges-1).
inline void bar_reorderRcpp_inline(int node, int nTips,
    const int* e1_ptr, const int* e2_ptr,
    std::vector<int>& neworder,
    const int* L_ptr, const int* xi_ptr, const int* xj_ptr,
    int &iii) {
    
    int i = node - nTips - 1;
    // Output the current node's block (edges for this internal node) in reverse order.
    for (int j = xj_ptr[i] - 1; j >= 0; j--) {
        neworder[iii--] = L_ptr[xi_ptr[i] + j] + 1; // +1 converts to R's 1-indexing
    }
    // Recursively process each child.
    for (int j = 0; j < xj_ptr[i]; j++) {
        int k = e2_ptr[ L_ptr[xi_ptr[i] + j] ];
        if (k > nTips)
            bar_reorderRcpp_inline(k, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    }
}

// Helper function to reorder the rows of a column-major edge matrix.
// E is a one-dimensional arma::Col<int> where the first nEdges entries are the parent column
// and the next nEdges entries are the child column.
// 'order' is a vector of 1-indexed row numbers in the desired order.
// The function returns a new arma::Col<int> with rows rearranged (still column-major).
arma::Col<int> reorder_rows(const arma::Col<int>& E, const std::vector<int>& order) {
    int nEdges = E.n_elem / 2;
    arma::Col<int> newE(E.n_elem);
    for (int i = 0; i < nEdges; i++) {
        int orig = order[i] - 1; // convert from 1-indexed to 0-indexed row number
        newE[i] = E[orig];             // parent's value (first column)
        newE[nEdges + i] = E[nEdges + orig];  // child's value (second column)
    }
    return newE;
}

// E is a one-dimensional vector representing an edge matrix in column-major order.
// [[Rcpp::export]]
arma::Col<int> reorderRcpp(const arma::Col<int>& E) {
    int nEdges = E.n_elem / 2;
    int nTips = nEdges / 2 + 1;
    int root = nTips + 1;
    
    arma::Col<int> parent = E.rows(0, nEdges - 1);
    arma::Col<int> child = E.rows(nEdges, E.n_elem - 1);
    
    // The maximum parent label tells us the total number of nodes.
    int m_val = parent.max();
    int nnode = m_val - nTips; // number of internal nodes

    // Allocate working arrays.
    std::vector<int> L(nEdges);            // Will store edge indices for each internal node.
    std::vector<int> neworder(nEdges);       // The final reordering (stored as 1-indexed row numbers).
    std::vector<int> pos(nnode, 0);          // Current fill position for each internal node.
    std::vector<int> xi(nnode, 0);           // Starting index for each internal node in L.
    std::vector<int> xj(nnode, 0);           // Count of children per internal node.

    // First pass: count children per internal node using parent's values.
    for (int i = 0; i < nEdges; i++) {
        int idx = parent[i] - nTips - 1;
        xj[idx]++;
    }
    
    // Compute starting positions xi as cumulative sums.
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i - 1] + xj[i - 1];
    }
    
    // Fill L: For each edge, assign its row index to the appropriate block.
    for (int i = 0; i < nEdges; i++) {
        int k = parent[i] - nTips - 1;
        int j = pos[k];
        L[xi[k] + j] = i;
        pos[k]++;
    }
    
    // Reset the new order index.
    int iii = nEdges - 1;
    
    // Get raw pointers for fast access.
    const int* e1_ptr = parent.memptr();
    const int* e2_ptr = child.memptr();
    const int* L_ptr   = L.data();
    const int* xi_ptr  = xi.data();
    const int* xj_ptr  = xj.data();
    
    // Run the recursive postorder traversal.
    bar_reorderRcpp_inline(root, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    
    // Use the computed new order to reorder the rows of E.
    arma::Col<int> newE = reorder_rows(E, neworder);
    
    return newE;
}

// [[Rcpp::export]]
std::vector<arma::Col<int>> nnin_cpp(const arma::Col<int>& E, const int n) {

    const int numEdges = E.n_elem / 2;
    const int nTips = numEdges / 2 + 1;
    
    // Create subviews for parent's and child's columns without copying.
    arma::subview_col<int> parent = E.subvec(0, numEdges - 1);
    arma::subview_col<int> child  = E.subvec(numEdges, E.n_elem - 1);
    
    // Find internal edges (child > nTips) using vectorized operations.
    arma::uvec internalEdges = arma::find(child > nTips);
    if (internalEdges.n_elem < (unsigned int)n)
        stop("n is larger than the number of valid internal edges.");
    
    // Select the nth internal edge (0-indexed)
    const int ind = internalEdges(n - 1);
    
    // Retrieve parent's value (p1) and child's value (p2) for the chosen edge.
    const int p1 = parent(ind);
    const int p2 = child(ind);
    
    // Find indices where parent equals p1.
    arma::uvec indices_p1 = arma::find(parent == p1);
    // Choose the index that is not equal to 'ind'
    const int ind1 = (indices_p1(0) == (unsigned int) ind) ? indices_p1(1) : indices_p1(0);
    
    // Find indices where parent equals p2.
    arma::uvec indices_p2 = arma::find(parent == p2);
    const int ind2_0 = indices_p2(0);
    const int ind2_1 = indices_p2(1);
    
    // Retrieve the child values for these swap candidates.
    const int e1_val = child(ind1);
    const int e2_val = child(ind2_0);
    const int e3_val = child(ind2_1);
    
    // Create copies for the two alternative topologies.
    arma::Col<int> E1 = E;  // copy of E
    arma::Col<int> E2 = E;  // another copy
    
    // Topology 1: swap child at ind1 with that at ind2_0.
    E1(numEdges + ind1) = e2_val;
    E1(numEdges + ind2_0) = e1_val;
    
    // Topology 2: swap child at ind1 with that at ind2_1.
    E2(numEdges + ind1) = e3_val;
    E2(numEdges + ind2_1) = e1_val;
    
    // Reorder each topology (e.g., into postorder) using the helper reorderRcpp.    
    std::vector<arma::Col<int>> res(2);
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
// E is a flattened edge list in postorder (each edge: parent, child) stored in column-major order.
// n: number of nodes, C: number of states, m: number of edges, root: index of the root node.
// [[Rcpp::export]]
double score_tree_bp(const arma::Col<int>& E,
                     const std::vector<double>& logP, 
                     const std::vector<double>& logA,
                     const int n, const int C, const int m, const int root) {
    // Allocate memory for messages and temporary state values.
    std::vector<double> log_messages(C * n, 0.0);
    std::vector<double> state_log_values(C);
    std::vector<double> temp(C);
    int idx;

    for (int i = 0; i < m; i++) {
        int par  = E(i);         // parent's value for edge i
        int node = E(m + i);     // child's value for edge i

        // Precompute common expression for each c_child.
        for (int c_child = 0; c_child < C; c_child++) {
            idx = c_child * n + node;
            temp[c_child] = logP[idx] + log_messages[idx];
        }
        
        for (int c = 0; c < C; c++) {
            for (int c_child = 0; c_child < C; c_child++) {
                state_log_values[c_child] = logA[c * C + c_child] + temp[c_child];
            }
            log_messages[c * n + par] += logSumExp(state_log_values);
        }
    }
  
    for (int c = 0; c < C; c++) {
        idx = c * n + root;
        state_log_values[c] = logP[idx] + log_messages[idx];
    }
    return logSumExp(state_log_values);
}

// wrapper: Precomputes dimensions and passes them to score_tree_bp.
// E: Edge matrix (each row is a (parent, child) pair, 1-indexed from R) provided as an arma::Col<int> in column-major order.
// logP_list: List of flattened likelihood matrices (each in row-major order)
// logA: Flattened transition matrix (row-major order)
// Computes:
//   - L: number of loci (length of logP_list),
//   - C: number of states (inferred from logA),
//   - n: number of nodes (inferred from the first likelihood matrix),
//   - m: number of edges (number of rows in the original edge matrix),
//   - root: parent's value from the last edge.
// [[Rcpp::export]]
double score_tree_bp_wrapper(arma::Col<int> E,
                             const std::vector< std::vector<double> >& logP_list,
                             const std::vector<double>& logA) {
  // Compute number of loci.
  int L = logP_list.size();
  // Infer number of states from the transition matrix.
  int C = std::sqrt(logA.size());
  // Infer number of nodes from the first likelihood matrix.
  int n = logP_list[0].size() / C;
  // Compute m: number of edges.
  int m = E.n_elem / 2;
  
  // Adjust the edge matrix from 1-indexing (R) to 0-indexing (C++).
  E = E - 1;
  int root = E(m - 1);
  
  double logZ = 0.0;
  for (int l = 0; l < L; l++) {
    logZ += score_tree_bp(E, logP_list[l], logA, n, C, m, root);
  }
  return logZ;
}

struct score_neighbours: public Worker {

    // original tree
    const arma::Col<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<double> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours(const arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Col<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel(arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA) {

    E = reorderRcpp(E);

    int n = E.n_elem / 4 - 1;

    NumericVector scores(2*n);

    score_neighbours score_neighbours(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours);

    return scores;

}


