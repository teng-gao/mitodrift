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
arma::Mat<int> reorder_rows(const arma::Mat<int> x, const std::vector<int> order) {
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

// Inline helper performing the recursive postorder traversal.
// Parameters:
//   node, nTips: as before.
//   e1_ptr and e2_ptr: raw pointers to the parent's and child's arrays.
//   neworder: vector (of length nEdges) that will be filled with the new row ordering (stored 1-indexed).
//   L_ptr, xi_ptr, xj_ptr: auxiliary arrays computed for the internal nodes.
//   iii: current index to fill in neworder (initialized to nEdges-1).
inline void bar_reorderRcpp_inline2(int node, int nTips,
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
            bar_reorderRcpp_inline2(k, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    }
}

// Helper function to reorder the rows of a column-major edge matrix.
// E is a one-dimensional arma::Col<int> where the first nEdges entries are the parent column
// and the next nEdges entries are the child column.
// 'order' is a vector of 1-indexed row numbers in the desired order.
// The function returns a new arma::Col<int> with rows rearranged (still column-major).
arma::Col<int> reorder_rows2(const arma::Col<int>& E, const std::vector<int>& order) {
    int nEdges = E.n_elem / 2;
    arma::Col<int> newE(E.n_elem);
    for (int i = 0; i < nEdges; i++) {
        int orig = order[i] - 1; // convert from 1-indexed to 0-indexed row number
        newE[i] = E[orig];             // parent's value (first column)
        newE[nEdges + i] = E[nEdges + orig];  // child's value (second column)
    }
    return newE;
}

// [[Rcpp::export]]
arma::Col<int> reorderRcpp2(arma::Col<int> E) {
    // E is a one-dimensional vector representing an edge matrix in column-major order.
    // Let nEdges be the number of rows (edges) in the original edge matrix.
    int nEdges = E.n_elem / 2;
    // For a fully bifurcating (binary) tree, the standard relationship is:
    // nTips = nEdges/2 + 1.
    int nTips = nEdges / 2 + 1;
    // The root is assumed to be nTips + 1.
    int root = nTips + 1;
    
    // Extract the parent's and child's columns from E.
    // In column-major storage, the first nEdges elements are the parent column.
    arma::Col<int> parent = E.rows(0, nEdges - 1);
    // The next nEdges elements form the child column.
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
    bar_reorderRcpp_inline2(root, nTips, e1_ptr, e2_ptr, neworder, L_ptr, xi_ptr, xj_ptr, iii);
    
    // Use the computed new order to reorder the rows of E.
    arma::Col<int> newE = reorder_rows2(E, neworder);
    
    return newE;
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


// [[Rcpp::export]]
std::vector<arma::Col<int>> nnin_cpp_vec(arma::Col<int> E, const int n) {
    // E is assumed to be in column-major order.
    // Let numEdges be the number of rows (edges) in the original edge matrix.
    int numEdges = E.n_elem / 2;
    // For a fully bifurcating binary tree, the standard relationship is:
    // nTips = numEdges/2 + 1.
    int nTips = numEdges / 2 + 1;
    
    // Extract parent's and child's columns using subvec (using parentheses for element access).
    arma::Col<int> parent = E.subvec(0, numEdges - 1);
    arma::Col<int> child  = E.subvec(numEdges, E.n_elem - 1);
    
    // Find internal edges (those with child > nTips) using vectorized operations.
    arma::uvec internalEdges = arma::find(child > nTips);
    if (internalEdges.n_elem < (unsigned int)n)
        stop("n is larger than the number of valid internal edges.");
    // Select the nth internal edge (0-indexed)
    int ind = internalEdges(n - 1);
    
    // Retrieve parent's value (p1) and child's value (p2) for the chosen edge.
    int p1 = parent(ind);
    int p2 = child(ind);
    
    // Find indices where parent equals p1.
    arma::uvec indices_p1 = arma::find(parent == p1);
    // Since indices_p1 has at most 2 elements, choose the one that is not 'ind'.
    int ind1 = (indices_p1(0) == (unsigned int) ind) ? indices_p1(1) : indices_p1(0);
    
    // Find indices where parent equals p2.
    arma::uvec indices_p2 = arma::find(parent == p2);
    int ind2_0 = indices_p2(0);
    int ind2_1 = indices_p2(1);
    
    // Retrieve the child values for these edges.
    int e1_val = child(ind1);      // child corresponding to p1 (other than the target edge)
    int e2_val = child(ind2_0);      // one child of p2
    int e3_val = child(ind2_1);      // the other child of p2
    
    // Create copies of E for the two alternative topologies.
    arma::Col<int> E1 = E;
    arma::Col<int> E2 = E;
    
    // In the first topology, swap the child at ind1 with that at ind2_0.
    E1(numEdges + ind1)   = e2_val;
    E1(numEdges + ind2_0) = e1_val;
    
    // In the second topology, swap the child at ind1 with that at ind2_1.
    E2(numEdges + ind1)   = e3_val;
    E2(numEdges + ind2_1) = e1_val;

    // Return both alternative topologies as a vector of arma::Col<int>.
    std::vector<arma::Col<int>> res(2);
    res[0] = reorderRcpp2(E1);
    res[1] = reorderRcpp2(E2);
    return res;
}


// [[Rcpp::export]]
std::vector<arma::Col<int>> nnin_cpp_vec2(arma::Col<int> E, const int n) {

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
    
    // Reorder each topology (e.g., into postorder) using the helper reorderRcpp2.    
    std::vector<arma::Col<int>> res(2);
    res[0] = reorderRcpp2(E1);
    res[1] = reorderRcpp2(E2);
    return res;
}

// // E is assumed to be a one-dimensional vector representing an edge matrix in row-major order
// // [[Rcpp::export]]
// std::vector<std::vector<int>> nnin_cpp_vec(std::vector<int> E, const int n) {
    
//     int numEdges = E.size() / 2;
//     int nTips = numEdges / 2 + 1;
    
//     // Find the nth internal edge (an edge with child > nTips).
//     int count = 0, ind = -1;
//     for (int i = 0; i < numEdges; i++) {
//         int childVal = E[2 * i + 1];
//         if (childVal > nTips) {
//             count++;
//             if (count == n) {
//                 ind = i;
//                 break;
//             }
//         }
//     }
//     if (ind < 0) {
//         Rcpp::stop("n is larger than the number of valid internal edges.");
//     }
    
//     // Get the parent's value (p1) and child’s value (p2) for the chosen edge.
//     int p1 = E[2 * ind];
//     int p2 = E[2 * ind + 1];
    
//     // Find the first index (other than ind) where the parent equals p1.
//     int ind1 = -1;
//     for (int i = 0; i < numEdges; i++) {
//         if (i != ind && E[2 * i] == p1) {
//             ind1 = i;
//             break;
//         }
//     }
//     if (ind1 < 0) {
//         Rcpp::stop("No valid index found for p1 interchange.");
//     }
    
//     // Find two indices where the parent equals p2.
//     int ind2_0 = -1, ind2_1 = -1;
//     for (int i = 0; i < numEdges; i++) {
//         if (E[2 * i] == p2) {
//             if (ind2_0 < 0) {
//                 ind2_0 = i;
//             } else if (ind2_1 < 0 && i != ind) {
//                 ind2_1 = i;
//                 break;
//             }
//         }
//     }
//     if (ind2_0 < 0 || ind2_1 < 0) {
//         Rcpp::stop("No valid indices found for p2 interchange.");
//     }
    
//     // Retrieve the child values at these indices.
//     int e1_val = E[2 * ind1 + 1];
//     int e2_val = E[2 * ind2_0 + 1];
//     int e3_val = E[2 * ind2_1 + 1];
    
//     // Make copies of E for the two alternative topologies.
//     std::vector<int> E1 = E;
//     std::vector<int> E2 = E;
    
//     // In the first topology, swap the child at ind1 with that at ind2_0.
//     E1[2 * ind1 + 1]   = e2_val;
//     E1[2 * ind2_0 + 1] = e1_val;
    
//     // In the second topology, swap the child at ind1 with that at ind2_1.
//     E2[2 * ind1 + 1]   = e3_val;
//     E2[2 * ind2_1 + 1] = e1_val;
    
//     // Return both alternative topologies.
//     std::vector<std::vector<int>> res(2);
//     res[0] = E1;
//     res[1] = E2;
//     return res;
// }

// Recursive helper that, given a node label, returns all edge indices (from L)
// in the subtree rooted at that node.
// If the node is a tip (node <= nTips), return an empty set.
// Parameters:
//   node: the node label (tip or internal)
//   nTips: number of tip nodes
//   L: an array (vector<int>) built from the edge matrix that groups edges by parent.
//   xi: starting indices for each internal node in L (indexed by node - nTips - 1)
//   xj: number of children (edges) for each internal node (again indexed by node - nTips - 1)
//   child: the vector of child labels (from the edge matrix)
// Returns an unordered_set<int> containing edge indices (0-indexed) in the subtree.
std::unordered_set<int> getSubtreeEdges(int node, int nTips, 
                                          const std::vector<int>& L, 
                                          const std::vector<int>& xi, 
                                          const std::vector<int>& xj,
                                          const arma::Col<int>& child) {
    std::unordered_set<int> subtree;
    // For a tip, no subtree edges.
    if (node <= nTips) return subtree;
    
    // For an internal node, compute its internal index.
    int internalIndex = node - nTips - 1;
    // The immediate edges for this node are in L starting at xi[internalIndex] of count xj[internalIndex].
    for (int j = 0; j < xj[internalIndex]; j++) {
        int edgeIndex = L[xi[internalIndex] + j];
        // Add the current edge.
        subtree.insert(edgeIndex);
        // Get the child of this edge.
        int childNode = child(edgeIndex);
        // If the child is internal, recursively add its subtree edges.
        if (childNode > nTips) {
            std::unordered_set<int> subSubtree = getSubtreeEdges(childNode, nTips, L, xi, xj, child);
            subtree.insert(subSubtree.begin(), subSubtree.end());
        }
    }
    return subtree;
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

// bp: Belief-propagation function.
// logP is a flattened likelihood matrix (row-major; dimensions: C x n)
// logA is a flattened transition matrix (dimensions: C x C)
// E is a flattened edge list in postorder (each edge: parent, child) stored in column-major order.
// n: number of nodes, C: number of states, m: number of edges, root: index of the root node.
// [[Rcpp::export]]
double score_tree_bp2(arma::Col<int> E, 
                     std::vector<double> logP, 
                     std::vector<double> logA,
                     int n, int C, int m, int root) {
  // Allocate memory for messages and temporary state values.
  std::vector<double> log_messages(C * n, 0.0);
  std::vector<double> state_log_values(C);
  
  // Process nodes in postorder; E is assumed to be provided in postorder.
  // In column-major order: parent's column occupies indices 0..(m-1) and child's column indices m..(2*m-1)
  for (int i = 0; i < m; i++) {
    int par  = E(i);        // parent's value for edge i
    int node = E(m + i);     // child's value for edge i
    
    for (int c = 0; c < C; c++) {
      for (int c_child = 0; c_child < C; c_child++) {
        // logP is row-major, so element (c_child, node) is at index: c_child * n + node.
        state_log_values[c_child] = logA[c * C + c_child] +
                                    logP[c_child * n + node] +
                                    log_messages[c_child * n + node];
      }
      log_messages[c * n + par] += logSumExp(state_log_values);
    }
  }
  
  for (int c = 0; c < C; c++) {
    state_log_values[c] = logP[c * n + root] + log_messages[c * n + root];
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
double score_tree_bp_wrapper2(arma::Col<int> E,
                             std::vector< std::vector<double> > logP_list,
                             std::vector<double> logA) {
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
    logZ += score_tree_bp2(E, logP_list[l], logA, n, C, m, root);
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


struct score_neighbours2: public Worker {

    // original tree
    const arma::Col<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<double> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours2(const arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Col<int>> Ep = nnin_cpp_vec2(E, i+1);
            scores[2*i] = score_tree_bp_wrapper2(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper2(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel2(arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA) {

    E = reorderRcpp2(E);

    int n = E.n_elem / 4 - 1;

    NumericVector scores(2*n);

    score_neighbours2 score_neighbours2(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours2);

    return scores;

}


