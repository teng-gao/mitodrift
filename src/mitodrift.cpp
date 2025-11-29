#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <limits>
#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>
#include <array>
#include <string>
#include <vector>
#include <sstream>
#include <numeric>
#include <mutex>
#include <memory>
#include <thread>

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

/////////////////////////////////////// ML Tree Search ////////////////////////////////////////

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
double logSumExp(const arma::vec& x) {
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


static inline double logsumexp_array(const double* x, int len) {
	double maxv = -std::numeric_limits<double>::infinity();
	for (int i = 0; i < len; ++i) {
		if (x[i] > maxv) maxv = x[i];
	}
	if (!(maxv > -std::numeric_limits<double>::infinity())) {
		return -std::numeric_limits<double>::infinity();
	}
	double sum = 0.0;
	for (int i = 0; i < len; ++i) {
		sum += std::exp(x[i] - maxv);
	}
	return maxv + std::log(sum);
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
    // Find the root
    int root = E(m - 1);

    double logZ = 0.0;
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp(E, logP_list[l], logA, n, C, m, root);
    }
    return logZ;
}




// Core routine for belief propagation with precomputed exp-shifted A and row maxes.
// Messages are stored in node-major layout: msg[node * C + c]
static inline double score_tree_bp2_core(const int* e_ptr, int m, int n, int C, int root,
                                         const double* logP_ptr,
                                         const double* expA_shifted,   // size C*C, row-major
                                         const double* row_maxA,       // size C
                                         double* msg_nm,               // size n*C, node-major
                                         double* temp,                 // size C (scratch)
                                         double* u) {                  // size C (scratch)
	for (int i = 0; i < m; ++i) {
		const int par  = e_ptr[i];
		const int node = e_ptr[m + i];

		// Build temp[c] = logP[c*n + node] + msg_nm[node*C + c], and find its max for stability.
		double max_t = -std::numeric_limits<double>::infinity();
		int offP = node;
		const int base_child = node * C;
		for (int c = 0; c < C; ++c, offP += n) {
			const double v = logP_ptr[offP] + msg_nm[base_child + c];
			temp[c] = v;
			if (v > max_t) max_t = v;
		}
		// u[c] = exp(temp[c] - max_t)  (shared across all rows of A)
		for (int c = 0; c < C; ++c) u[c] = std::exp(temp[c] - max_t);

		// For each parent state c, accumulate log-sum-exp as:
		// row_maxA[c] + max_t + log( dot( expA_shifted_row_c , u ) )
		const int base_par = par * C;
		const double* rowA = expA_shifted;
		for (int c = 0; c < C; ++c, rowA += C) {
			double s = 0.0;
			// dot product between rowA (length C) and u (length C)
			for (int j = 0; j < C; ++j) s += rowA[j] * u[j];
			msg_nm[base_par + c] += row_maxA[c] + max_t + std::log(s);
		}
	}

	// Aggregate at root over states: logsumexp_c( logP[c*n + root] + msg_nm[root*C + c] )
	int offP = root;
	const int base_root = root * C;
	for (int c = 0; c < C; ++c, offP += n) {
		temp[c] = logP_ptr[offP] + msg_nm[base_root + c];
	}
	return logsumexp_array(temp, C);
}

// bp: Belief-propagation function.
// logP is a flattened likelihood matrix (row-major; dimensions: C x n)
// logA is a flattened transition matrix (dimensions: C x C)
// E is a flattened edge list in postorder (each edge: parent, child) stored in column-major order.
// n: number of nodes, C: number of states, m: number of edges, root: index of the root node.
// [[Rcpp::export]]
double score_tree_bp2(const arma::Col<int>& E,
                      const std::vector<double>& logP, 
                      const std::vector<double>& logA,
                      const int n, const int C, const int m, const int root) {
	// Precompute per-row max of A and exp(A - max_row) once.
	std::vector<double> row_maxA(C);
	std::vector<double> expA_shifted(static_cast<size_t>(C) * C);
	for (int r = 0; r < C; ++r) {
		const double* Arow = logA.data() + static_cast<size_t>(r) * C;
		double mr = Arow[0];
		for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
		row_maxA[r] = mr;
		double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
		for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
	}

	// Node-major message buffer and small scratch.
	std::vector<double> msg_nm(static_cast<size_t>(n) * C, 0.0);
	std::vector<double> temp(C);
	std::vector<double> u(C);

	return score_tree_bp2_core(E.memptr(), m, n, C, root,
	                           logP.data(),
	                           expA_shifted.data(), row_maxA.data(),
	                           msg_nm.data(), temp.data(), u.data());
}

// bp2: Belief-propagation function.
// logP_list: List of flattened likelihood matrices (each in row-major order)
// logA: Flattened transition matrix (row-major order)
// [[Rcpp::export]]
double score_tree_bp_wrapper2(arma::Col<int> E,
                              const std::vector< std::vector<double> >& logP_list,
                              const std::vector<double>& logA) {
	const int L = static_cast<int>(logP_list.size());
	const int C = static_cast<int>(std::sqrt(logA.size()));
	const int n = static_cast<int>(logP_list[0].size() / C);
	const int m = static_cast<int>(E.n_elem / 2);

	// 1-indexed (R) -> 0-indexed (C++)
	E -= 1;
	const int root = E(m - 1);
	const int* e_ptr = E.memptr();

	// Precompute exp-shifted rows of A and per-row max once (shared across loci).
	std::vector<double> row_maxA(C);
	std::vector<double> expA_shifted(static_cast<size_t>(C) * C);
	for (int r = 0; r < C; ++r) {
		const double* Arow = logA.data() + static_cast<size_t>(r) * C;
		double mr = Arow[0];
		for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
		row_maxA[r] = mr;
		double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
		for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
	}

	// Reusable buffers
	std::vector<double> msg_nm(static_cast<size_t>(n) * C);
	std::vector<double> temp(C);
	std::vector<double> u(C);

	double logZ = 0.0;
	for (int l = 0; l < L; ++l) {
		std::fill(msg_nm.begin(), msg_nm.end(), 0.0);
		logZ += score_tree_bp2_core(e_ptr, m, n, C, root,
		                            logP_list[l].data(),
		                            expA_shifted.data(), row_maxA.data(),
		                            msg_nm.data(), temp.data(), u.data());
	}
	return logZ;
}

// ----------------------------
// Node & Edge Belief Computation (BP2)
// ----------------------------

static inline void compute_down_msg_fast(
	const double* S_parent,
	const double* expA_rowmajor,
	int C,
	double* out,
	double* su);

struct BeliefScratch {
	std::vector<double> Uin;
	std::vector<double> up_to_parent;
	std::vector<double> down_in;
	std::vector<double> temp;
	std::vector<double> u;
	std::vector<double> S_parent;
	std::vector<double> down_msg;
	std::vector<double> a;
	std::vector<double> b;
	std::vector<int> stack;

	BeliefScratch(int n, int C)
		: Uin(static_cast<std::size_t>(n) * static_cast<std::size_t>(C)),
		  up_to_parent(static_cast<std::size_t>(n) * static_cast<std::size_t>(C)),
		  down_in(static_cast<std::size_t>(n) * static_cast<std::size_t>(C)),
		  temp(C), u(C), S_parent(C), down_msg(C), a(C), b(C) {
		stack.reserve(static_cast<std::size_t>(n));
	}
};

static inline double process_belief_locus(
	const double* P,
	double* node_out,
	double* edge_out,
	BeliefScratch& scratch,
	int m,
	int n,
	int C,
	int root,
	const int* parent_ptr,
	const int* child_ptr,
	const std::vector< std::array<int,2> >& children_of,
	const double* expA_shifted_ptr,
	const double* expA_ptr,
	const double* row_maxA_ptr)
{
	const double neg_inf = -std::numeric_limits<double>::infinity();
	const std::size_t nC = static_cast<std::size_t>(n) * static_cast<std::size_t>(C);

	double* const Uin_data = scratch.Uin.data();
	double* const up_data = scratch.up_to_parent.data();
	double* const down_data = scratch.down_in.data();
	double* const temp_data = scratch.temp.data();
	double* const u_data = scratch.u.data();
	double* const S_parent_data = scratch.S_parent.data();
	double* const down_msg_data = scratch.down_msg.data();
	double* const a_data = scratch.a.data();
	double* const b_data = scratch.b.data();
	std::vector<int>& stack = scratch.stack;

	std::fill(Uin_data, Uin_data + nC, 0.0);
	std::fill(up_data, up_data + nC, 0.0);
	std::fill(down_data, down_data + nC, 0.0);

	for (int i = 0; i < m; ++i) {
		const int par = parent_ptr[i];
		const int node = child_ptr[i];
		double max_t = neg_inf;
		int offP = node;
		double* const Uin_node = Uin_data + static_cast<std::size_t>(node) * C;
		double* const up_node = up_data + static_cast<std::size_t>(node) * C;
		double* const Uin_parent = Uin_data + static_cast<std::size_t>(par) * C;
		for (int c = 0; c < C; ++c, offP += n) {
			const double v = P[offP] + Uin_node[c];
			temp_data[c] = v;
			if (v > max_t) max_t = v;
		}
		for (int c = 0; c < C; ++c) u_data[c] = std::exp(temp_data[c] - max_t);

		const double* rowA = expA_shifted_ptr;
		for (int r = 0; r < C; ++r, rowA += C) {
			double s = 0.0;
			for (int j = 0; j < C; ++j) s += rowA[j] * u_data[j];
			const double f = row_maxA_ptr[r] + max_t + std::log(s);
			up_node[r] = f;
			Uin_parent[r] += f;
		}
	}

	stack.clear();
	stack.push_back(root);
	while (!stack.empty()) {
		const int u_node = stack.back();
		stack.pop_back();
		double* const down_parent = down_data + static_cast<std::size_t>(u_node) * C;
		double* const Uin_parent = Uin_data + static_cast<std::size_t>(u_node) * C;
		const auto& ch = children_of[u_node];
		for (int t = 0; t < 2; ++t) {
			const int v = ch[t];
			if (v < 0) continue;
			double* const down_child = down_data + static_cast<std::size_t>(v) * C;
			double* const up_child = up_data + static_cast<std::size_t>(v) * C;
			for (int r = 0; r < C; ++r) {
				const double P_ur = (u_node == root)
					? (r == 0 ? 0.0 : neg_inf)
					: P[static_cast<std::size_t>(r) * n + u_node];
				S_parent_data[r] = P_ur + down_parent[r] + Uin_parent[r] - up_child[r];
			}
			compute_down_msg_fast(S_parent_data, expA_ptr, C, down_msg_data, u_data);
			for (int k = 0; k < C; ++k) down_child[k] = down_msg_data[k];
			stack.push_back(v);
		}
	}

	double logZ_root = 0.0;
	for (int v = 0; v < n; ++v) {
		double maxW = neg_inf;
		const double* const down_v = down_data + static_cast<std::size_t>(v) * C;
		const double* const Uin_v = Uin_data + static_cast<std::size_t>(v) * C;
		for (int c = 0; c < C; ++c) {
			const double P_vc = (v == root)
				? (c == 0 ? 0.0 : neg_inf)
				: P[static_cast<std::size_t>(c) * n + v];
			temp_data[c] = P_vc + Uin_v[c] + down_v[c];
			if (temp_data[c] > maxW) maxW = temp_data[c];
		}
		double denom = 0.0;
		double* const node_col = node_out + static_cast<std::size_t>(v);
		for (int c = 0; c < C; ++c) {
			const double val = std::exp(temp_data[c] - maxW);
			node_col[static_cast<std::size_t>(c) * n] = val;
			denom += val;
		}
		const double inv = 1.0 / denom;
		for (int c = 0; c < C; ++c) node_col[static_cast<std::size_t>(c) * n] *= inv;
		if (v == root) logZ_root = maxW + std::log(denom);
	}

	const std::size_t stride_parent = static_cast<std::size_t>(m);
	const std::size_t stride_child = stride_parent * static_cast<std::size_t>(C);
	for (int i = 0; i < m; ++i) {
		const int u_node = parent_ptr[i];
		const int v = child_ptr[i];
		const double* const down_u = down_data + static_cast<std::size_t>(u_node) * C;
		const double* const Uin_u = Uin_data + static_cast<std::size_t>(u_node) * C;
		const double* const up_v = up_data + static_cast<std::size_t>(v) * C;
		const double* const down_v = down_data + static_cast<std::size_t>(v) * C;
		const double* const Uin_v = Uin_data + static_cast<std::size_t>(v) * C;
		for (int r = 0; r < C; ++r) {
			const double P_ur = (u_node == root)
				? (r == 0 ? 0.0 : neg_inf)
				: P[static_cast<std::size_t>(r) * n + u_node];
			S_parent_data[r] = P_ur + down_u[r] + Uin_u[r] - up_v[r];
		}
		for (int k = 0; k < C; ++k) down_msg_data[k] = P[static_cast<std::size_t>(k) * n + v] + Uin_v[k];
		double maxWr = neg_inf;
		for (int r = 0; r < C; ++r) if (S_parent_data[r] > maxWr) maxWr = S_parent_data[r];
		for (int r = 0; r < C; ++r) a_data[r] = std::exp(S_parent_data[r] - maxWr);
		double maxY = neg_inf;
		for (int k = 0; k < C; ++k) if (down_msg_data[k] > maxY) maxY = down_msg_data[k];
		for (int k = 0; k < C; ++k) b_data[k] = std::exp(down_msg_data[k] - maxY);
		double denom = 0.0;
		for (int r = 0; r < C; ++r) {
			double s = 0.0;
			const double* const Arow = expA_ptr + static_cast<std::size_t>(r) * C;
			for (int k = 0; k < C; ++k) s += Arow[k] * b_data[k];
			denom += a_data[r] * s;
		}
		const double inv = 1.0 / denom;
		for (int r = 0; r < C; ++r) {
			const double* const Arow = expA_ptr + static_cast<std::size_t>(r) * C;
			const std::size_t base_idx_r = static_cast<std::size_t>(i) + static_cast<std::size_t>(r) * stride_parent;
			for (int k = 0; k < C; ++k) {
				edge_out[base_idx_r + static_cast<std::size_t>(k) * stride_child] = a_data[r] * Arow[k] * b_data[k] * inv;
			}
		}
	}

	return logZ_root;
}

struct ComputeNodeEdgeBeliefsWorker : public Worker {
	const int m;
	const int n;
	const int C;
	const int root;
	const int* parent_ptr;
	const int* child_ptr;
	const std::vector< std::array<int,2> >& children_of;
	const std::vector<double>& row_maxA;
	const std::vector<double>& expA_shifted;
	const std::vector<double>& expA;
	const std::vector< std::vector<double> >& logP_list;
	std::vector< std::vector<double> >& node_storage;
	std::vector< std::vector<double> >& edge_storage;
	std::vector<double>& logZ_vals;

	ComputeNodeEdgeBeliefsWorker(
		int m_, int n_, int C_, int root_,
		const int* parent_ptr_, const int* child_ptr_,
		const std::vector< std::array<int,2> >& children_of_,
		const std::vector<double>& row_maxA_,
		const std::vector<double>& expA_shifted_,
		const std::vector<double>& expA_,
		const std::vector< std::vector<double> >& logP_list_,
		std::vector< std::vector<double> >& node_storage_,
		std::vector< std::vector<double> >& edge_storage_,
		std::vector<double>& logZ_vals_)
		: m(m_), n(n_), C(C_), root(root_),
		  parent_ptr(parent_ptr_), child_ptr(child_ptr_),
		  children_of(children_of_),
		  row_maxA(row_maxA_), expA_shifted(expA_shifted_), expA(expA_),
		  logP_list(logP_list_),
		  node_storage(node_storage_), edge_storage(edge_storage_),
		  logZ_vals(logZ_vals_) {}

	void operator()(std::size_t begin, std::size_t end) {
		const double* const expA_shifted_ptr = expA_shifted.data();
		const double* const expA_ptr = expA.data();
		const double* const row_maxA_ptr = row_maxA.data();
		BeliefScratch scratch(n, C);

		for (std::size_t idx_l = begin; idx_l < end; ++idx_l) {
			logZ_vals[idx_l] = process_belief_locus(
				logP_list[idx_l].data(),
				node_storage[idx_l].data(),
				edge_storage[idx_l].data(),
				scratch,
				m, n, C, root,
				parent_ptr,
				child_ptr,
				children_of,
				expA_shifted_ptr,
				expA_ptr,
				row_maxA_ptr);
		}
	}
};

// Compute downward message m_{u->v}(c_v) using fast column-wise expA and S_parent.
// expA_rowmajor: exp(logA) in row-major order, S_parent: log domain, out: result, su: scratch (length C)
static inline void compute_down_msg_fast(
	const double* S_parent,          // length C
	const double* expA_rowmajor,     // length C*C, exp(logA) row-major
	int C,
	double* out,                     // length C
	double* su)                      // scratch length C (holds exp(S_parent - maxS))
{
	// m_{u->v}(k) = logsum_u [ S_parent[u] + logA[u,k] ]
	// = maxS + log( dot( exp(S_parent - maxS), expA[,k] ) )
	double maxS = -std::numeric_limits<double>::infinity();
	for (int u = 0; u < C; ++u) if (S_parent[u] > maxS) maxS = S_parent[u];
	for (int u = 0; u < C; ++u) su[u] = std::exp(S_parent[u] - maxS);
	for (int k = 0; k < C; ++k) {
		double dot = 0.0;
		// column k of row-major matrix: index (u*C + k)
		for (int u = 0; u < C; ++u) dot += expA_rowmajor[u * C + k] * su[u];
		out[k] = std::log(dot) + maxS;
	}
}


// Returns:
//  - node_beliefs: List of L matrices (n x C)
//  - edge_beliefs: List of L 3D arrays (m x C x C) with dims (edge, parent, child)
//  - logZ: total log-likelihood summed over loci
//  - Ecounts_parent_child: per-locus CxC matrices summing edge_beliefs over edges (parent rows, child cols)
//  - Ecounts_child_parent: per-locus transpose of the above (child rows, parent cols)
//  - Ecounts_parent_child_internal: same as parent->child but only counting edges with internal children
//  - Ecounts_low_high: per‑locus CxC matrices summing edge_beliefs with rows=smaller node id, cols=larger node id
// [[Rcpp::export]]
Rcpp::List compute_node_edge_beliefs_bp2(
	arma::Col<int> E,
	const std::vector< std::vector<double> >& logP_list,
	const std::vector<double>& logA)
{
	// --- Dimensions ---
	const int L = static_cast<int>(logP_list.size());
	const int C = static_cast<int>(std::sqrt(logA.size()));
	if (L <= 0) Rcpp::stop("logP_list must be non-empty");
	const int n = static_cast<int>(logP_list[0].size() / C);
	const int m = static_cast<int>(E.n_elem / 2);

	// Ensure postorder and 0-indexing for BP core loops
	E = reorderRcpp(E);
	E -= 1;
	int root = E(m - 1);

	// --- Topology helpers ---
	std::vector<int> parent_of(n, -1);
	std::vector< std::array<int,2> > children_of(n, std::array<int,2>{-1, -1});
	for (int i = 0; i < m; ++i) {
		const int p = E[i];
		const int c = E[m + i];
		parent_of[c] = p;
		if (children_of[p][0] == -1) children_of[p][0] = c; else children_of[p][1] = c;
	}

	// --- Precompute A rows: row_maxA & exp(A - row_max) ---
	std::vector<double> row_maxA(C);
	std::vector<double> expA_shifted(static_cast<size_t>(C) * C);
	for (int r = 0; r < C; ++r) {
		const double* Arow = logA.data() + static_cast<size_t>(r) * C;
		double mr = Arow[0];
		for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
		row_maxA[r] = mr;
		double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
		for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
	}
	// Also precompute exp(logA) once for fast column-wise ops (downward + edge beliefs).
	std::vector<double> expA(static_cast<size_t>(C) * C);
	for (size_t i = 0; i < expA.size(); ++i) expA[i] = std::exp(logA[i]);

	std::vector< std::vector<double> > node_storage(
		static_cast<std::size_t>(L),
		std::vector<double>(static_cast<std::size_t>(n) * static_cast<std::size_t>(C)));
	std::vector< std::vector<double> > edge_storage(
		static_cast<std::size_t>(L),
		std::vector<double>(static_cast<std::size_t>(m) * static_cast<std::size_t>(C) * static_cast<std::size_t>(C)));
	std::vector<double> logZ_vals(static_cast<std::size_t>(L), 0.0);

	const int* parent_ptr = E.memptr();
	const int* child_ptr = E.memptr() + m;

	ComputeNodeEdgeBeliefsWorker worker(
		m, n, C, root,
		parent_ptr, child_ptr,
		children_of,
		row_maxA,
		expA_shifted,
		expA,
		logP_list,
		node_storage,
		edge_storage,
		logZ_vals);

	parallelFor(0, static_cast<std::size_t>(L), worker);

	Rcpp::List node_list(L);
	Rcpp::List edge_list(L);
	Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(m, C, C);

	double logZ_total = 0.0;
	for (int l = 0; l < L; ++l) {
		Rcpp::NumericMatrix node_beliefs(n, C);
		std::copy(node_storage[static_cast<std::size_t>(l)].begin(),
		          node_storage[static_cast<std::size_t>(l)].end(),
		          node_beliefs.begin());

		Rcpp::NumericVector edge_beliefs(edge_storage[static_cast<std::size_t>(l)].size());
		std::copy(edge_storage[static_cast<std::size_t>(l)].begin(),
		          edge_storage[static_cast<std::size_t>(l)].end(),
		          edge_beliefs.begin());
		edge_beliefs.attr("dim") = dims;

		node_list[l] = node_beliefs;
		edge_list[l] = edge_beliefs;
		logZ_total += logZ_vals[static_cast<std::size_t>(l)];
	}

	return Rcpp::List::create(
		Rcpp::Named("node_beliefs") = node_list,
		Rcpp::Named("edge_beliefs") = edge_list,
		Rcpp::Named("logZ") = logZ_total
	);
}

/* ===========================
 * Message-caching NNI engine
 * ===========================
 *
 * We cache, for every locus ℓ and node v, the length‑C vector F_ℓ[v,·] which is the
 * child→parent contribution added to the parent's message (i.e., the quantity we
 * accumulate into msg_nm[parent,·] inside score_tree_bp2_core). For a parent p with two
 * children a and b, msg_nm[p,·] = F[a,·] + F[b,·]. An NNI around edge (p1→p2) only changes
 * the sets of children of p1 and p2, so we can recompute F at p2, then p1, and then walk
 * upward to the root, updating F along that single path. This yields O(height × C × L)
 * updates per proposal instead of O(n × C × L).
 */

struct NNICache {
	// static tree info
	int n;					// #nodes
	int m;					// #edges
	int C;					// #states
	int L;					// #loci
	int nTips;				// #tips
	int root;				// root node id (0-indexed)
	arma::Col<int> E;		// edge list, 0-indexed, postorder (parent col [0..m-1], child col [m..2m-1])

	// topology helpers
	std::vector<int> parent_of;					// size n
	std::vector< std::array<int,2> > children_of;		// for internal nodes: exactly two children
	// transition precompute
	std::vector<double> row_maxA;				// size C
	std::vector<double> expA_shifted;			// size C*C, row-major

	// data likelihoods
	std::vector< std::vector<double> > logP_list;	// L × (C*n), row-major per locus
	std::vector<int> locus_ids;

	// per-locus cached child→parent contributions F (n*C) and root logZ
    // F layout: [node][c][l] -> node * (C*L) + c*L + l
	std::vector<double> F;		// n * C * L
	std::vector<double> logZ;	// L
    std::vector<double> logP_storage; // n * C * L
    double current_total_loglik_val;

    // Scratch space for calculations to avoid reallocations
    // Size C * L
    mutable std::vector<double> scratch_temp;
    mutable std::vector<double> scratch_u;
    mutable std::vector<double> scratch_F_p2;
    mutable std::vector<double> scratch_F_p1;
    mutable std::vector<double> scratch_childF;
    mutable std::vector<double> scratch_prev_childF;
    mutable std::vector<double> scratch_F_v;
    mutable std::vector<double> scratch_curF;
    
    // Size L
    mutable std::vector<double> scratch_max_t;
    mutable std::vector<double> scratch_s;
	mutable std::mutex scratch_mutex;

	NNICache(arma::Col<int> E_in,
		const std::vector< std::vector<double> >& logP_in,
		const std::vector<double>& logA_in,
		bool reorder = true) {
		initialize_cache(E_in, logP_in, logA_in, reorder, nullptr);
	}

	NNICache(arma::Col<int> E_in,
		const std::vector< std::vector<double> >& logP_in,
		const std::vector<double>& logA_in,
		const std::vector<int>& locus_subset,
		bool reorder)
	{
		initialize_cache(E_in, logP_in, logA_in, reorder, &locus_subset);
	}

	void initialize_cache(
		arma::Col<int> E_in,
		const std::vector< std::vector<double> >& logP_in,
		const std::vector<double>& logA_in,
		bool reorder,
		const std::vector<int>* locus_subset)
	{
		if (logP_in.empty()) {
			stop("logP_in must contain at least one locus");
		}
		const int total_loci = static_cast<int>(logP_in.size());
		if (locus_subset) {
			if (locus_subset->empty()) {
				stop("locus_subset must not be empty");
			}
			locus_ids = *locus_subset;
			for (int idx : locus_ids) {
				if (idx < 0 || idx >= total_loci) {
					stop("locus_subset index out of range");
				}
			}
		} else {
			locus_ids.resize(total_loci);
			std::iota(locus_ids.begin(), locus_ids.end(), 0);
		}

		L = static_cast<int>(locus_ids.size());
		C = static_cast<int>(std::sqrt(logA_in.size()));
		if (C <= 0) {
			stop("logA_in must define at least one state");
		}
		n = static_cast<int>(logP_in[locus_ids[0]].size() / C);

		// Initialize scratch space
		int CL = C * L;
		scratch_temp.resize(CL);
		scratch_u.resize(CL);
		scratch_F_p2.resize(CL);
		scratch_F_p1.resize(CL);
		scratch_childF.resize(CL);
		scratch_prev_childF.resize(CL);
		scratch_F_v.resize(CL);
		scratch_curF.resize(CL);
		
		scratch_max_t.resize(L);
		scratch_s.resize(L);

		// Bring E to postorder and 0-indexed
		if (reorder) {
			E_in = reorderRcpp(E_in);
		}
		E = E_in - 1;

		m = static_cast<int>(E.n_elem / 2);
		nTips = m / 2 + 1;
		root = E(m - 1);

		// parent/children
		parent_of.assign(n, -1);
		children_of.assign(n, std::array<int,2>{-1,-1});
		for (int i = 0; i < m; ++i) {
			const int p = E[i];
			const int c = E[m + i];
			parent_of[c] = p;
			// fill two slots
			if (children_of[p][0] == -1) children_of[p][0] = c;
			else children_of[p][1] = c;
		}

		// A precompute
		row_maxA.resize(C);
		expA_shifted.resize(static_cast<size_t>(C) * C);
		for (int r = 0; r < C; ++r) {
			const double* Arow = logA_in.data() + static_cast<size_t>(r) * C;
			double mr = Arow[0];
			for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
			row_maxA[r] = mr;
			double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
			for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
		}

		// Transpose logP
		logP_storage.resize(n * C * L);
		for (int l = 0; l < L; ++l) {
			const auto& P_l = logP_in[locus_ids[l]];
			for (int node = 0; node < n; ++node) {
				for (int c = 0; c < C; ++c) {
					// Old: P_l[c * n + node]
					// New: logP_storage[node * C * L + c * L + l]
					logP_storage[node * C * L + c * L + l] = P_l[c * n + node];
				}
			}
		}

		// allocate caches
		F.assign(n * C * L, 0.0);
		logZ.assign(L, 0.0);

		// initialize F and logZ per locus with one postorder pass
		std::vector<double> msg_nm(n * C * L, 0.0); // Accumulator for children messages
		
		// Scratch buffers for initialization
		std::vector<double> init_temp(C * L);
		std::vector<double> init_u(C * L);
		std::vector<double> init_max_t(L);
		std::vector<double> init_s(L);
		
		for (int i = 0; i < m; ++i) {
			const int par  = E[i];
			const int node = E[m + i];
			
			// Compute F[node]
			// temp = P[node] + msg_nm[node]
			
			const double* P_ptr = logP_storage.data() + node * C * L;
			const double* msg_ptr = msg_nm.data() + node * C * L;
			double* F_ptr = F.data() + node * C * L;
			double* msg_par_ptr = msg_nm.data() + par * C * L;
			
			std::fill(init_max_t.begin(), init_max_t.end(), -std::numeric_limits<double>::infinity());
			
			for (int c = 0; c < C; ++c) {
				const double* p = P_ptr + c * L;
				const double* m_ = msg_ptr + c * L;
				double* t = init_temp.data() + c * L;
				for (int l = 0; l < L; ++l) {
					double val = p[l] + m_[l];
					t[l] = val;
					if (val > init_max_t[l]) init_max_t[l] = val;
				}
			}
			
			for (int c = 0; c < C; ++c) {
				double* t = init_temp.data() + c * L;
				double* u = init_u.data() + c * L;
				for (int l = 0; l < L; ++l) {
					if (init_max_t[l] == -std::numeric_limits<double>::infinity()) {
						u[l] = 0.0;
					} else {
						u[l] = std::exp(t[l] - init_max_t[l]);
					}
				}
			}
			
			const double* rowA = expA_shifted.data();
			for (int r = 0; r < C; ++r, rowA += C) {
				std::fill(init_s.begin(), init_s.end(), 0.0);
				for (int j = 0; j < C; ++j) {
					double A_val = rowA[j];
					double* u = init_u.data() + j * L;
					for (int l = 0; l < L; ++l) {
						init_s[l] += A_val * u[l];
					}
				}
				
				double row_max = row_maxA[r];
				double* f = F_ptr + r * L;
				double* mp = msg_par_ptr + r * L;
				for (int l = 0; l < L; ++l) {
					if (init_max_t[l] == -std::numeric_limits<double>::infinity()) {
						f[l] = -std::numeric_limits<double>::infinity();
						mp[l] += -std::numeric_limits<double>::infinity();
					} else {
						double val = row_max + init_max_t[l] + std::log(init_s[l]);
						f[l] = val;
						mp[l] += val;
					}
				}
			}
		}
		
		// Root marginal
		const double* P_root = logP_storage.data() + root * C * L;
		const double* msg_root = msg_nm.data() + root * C * L;
		
		std::fill(init_max_t.begin(), init_max_t.end(), -std::numeric_limits<double>::infinity());
		
		for (int c = 0; c < C; ++c) {
			const double* p = P_root + c * L;
			const double* m_ = msg_root + c * L;
			double* t = init_temp.data() + c * L;
			for (int l = 0; l < L; ++l) {
				double val = p[l] + m_[l];
				t[l] = val;
				if (val > init_max_t[l]) init_max_t[l] = val;
			}
		}
		
		current_total_loglik_val = 0.0;
		for (int l = 0; l < L; ++l) {
			if (init_max_t[l] == -std::numeric_limits<double>::infinity()) {
				logZ[l] = -std::numeric_limits<double>::infinity();
			} else {
				double sum_exp = 0.0;
				for (int c = 0; c < C; ++c) {
					sum_exp += std::exp(init_temp[c * L + l] - init_max_t[l]);
				}
				logZ[l] = init_max_t[l] + std::log(sum_exp);
			}
			current_total_loglik_val += logZ[l];
		}
	}

	inline int other_child(int parent, int child) const {
		const auto &ch = children_of[parent];
		return (ch[0] == child) ? ch[1] : ch[0];
	}

	inline void compute_F_vectorized(
        const double* F_c1, const double* F_c2, // pointers to [c=0][l=0]
        const double* P_base, // pointer to P[node][0][0]
        double* outF,
        // scratch buffers of size C*L
        double* u, double* temp, double* s_vec, double* max_t
    ) const {
        // 1. Compute temp = P + F1 + F2 and max_t
        std::fill(max_t, max_t + L, -std::numeric_limits<double>::infinity());
        
        for (int c = 0; c < C; ++c) {
            const double* p_ptr = P_base + c * L;
            const double* f1_ptr = F_c1 + c * L;
            const double* f2_ptr = F_c2 + c * L;
            double* t_ptr = temp + c * L;
            
            for (int l = 0; l < L; ++l) {
                double val = p_ptr[l] + f1_ptr[l] + f2_ptr[l];
                t_ptr[l] = val;
                if (val > max_t[l]) max_t[l] = val;
            }
        }
        
        // 2. Compute u = exp(temp - max_t)
        for (int c = 0; c < C; ++c) {
            double* t_ptr = temp + c * L;
            double* u_ptr = u + c * L;
            for (int l = 0; l < L; ++l) {
                if (max_t[l] == -std::numeric_limits<double>::infinity()) {
                    u_ptr[l] = 0.0;
                } else {
                    u_ptr[l] = std::exp(t_ptr[l] - max_t[l]);
                }
            }
        }
        
        // 3. Compute s = sum(A * u) and outF
        const double* rowA = expA_shifted.data();
        for (int r = 0; r < C; ++r, rowA += C) {
            double* out_ptr = outF + r * L;
            
            std::fill(s_vec, s_vec + L, 0.0);
            
            for (int j = 0; j < C; ++j) {
                double A_val = rowA[j];
                const double* u_ptr = u + j * L;
                for (int l = 0; l < L; ++l) {
                    s_vec[l] += A_val * u_ptr[l];
                }
            }
            
            double row_max = row_maxA[r];
            for (int l = 0; l < L; ++l) {
                if (max_t[l] == -std::numeric_limits<double>::infinity()) {
                    out_ptr[l] = -std::numeric_limits<double>::infinity();
                } else {
                    out_ptr[l] = row_max + max_t[l] + std::log(s_vec[l]);
                }
            }
        }
    }

    inline double compute_root_logZ_vectorized(
        const double* F_c1, const double* F_c2,
        const double* P_base,
        double* temp, double* max_t
    ) const {
        std::fill(max_t, max_t + L, -std::numeric_limits<double>::infinity());
        
        for (int c = 0; c < C; ++c) {
            const double* p_ptr = P_base + c * L;
            const double* f1_ptr = F_c1 + c * L;
            const double* f2_ptr = F_c2 + c * L;
            double* t_ptr = temp + c * L;
            
            for (int l = 0; l < L; ++l) {
                double val = p_ptr[l] + f1_ptr[l] + f2_ptr[l];
                t_ptr[l] = val;
                if (val > max_t[l]) max_t[l] = val;
            }
        }
        
        double total_logZ = 0.0;
        for (int l = 0; l < L; ++l) {
            if (max_t[l] == -std::numeric_limits<double>::infinity()) {
                total_logZ += -std::numeric_limits<double>::infinity();
            } else {
                double sum_exp = 0.0;
                for (int c = 0; c < C; ++c) {
                    sum_exp += std::exp(temp[c * L + l] - max_t[l]);
                }
                total_logZ += max_t[l] + std::log(sum_exp);
            }
        }
        if (std::isnan(total_logZ)) return -std::numeric_limits<double>::infinity();
        return total_logZ;
    }

	// Compute new total loglik if we perform the NNI "which" (0 or 1) at the nth internal edge.
	double compute_new_loglik(int edge_n, int which) const {
		std::lock_guard<std::mutex> guard(scratch_mutex);
		// Locate the nth internal edge (child >= nTips after 0-indexing)
		int cnt = 0, ind = -1;
		for (int i = 0; i < m; ++i) {
			if (E[m + i] >= nTips) {
				++cnt;
				if (cnt == edge_n) { ind = i; break; }
			}
		}
		if (ind < 0) stop("edge_n out of range in compute_new_loglik");

		const int p1 = E[ind];
		const int p2 = E[m + ind];
		const int c1 = other_child(p1, p2);
		const int c2 = children_of[p2][0];
		const int c3 = children_of[p2][1];

		const int cX    = (which == 0 ? c2 : c3);   // moves to p1
		const int cStay = (which == 0 ? c3 : c2);   // stays under p2

		// 1) Recompute F at p2 with new children {c1, cStay}
        compute_F_vectorized(
            F.data() + c1 * C * L,
            F.data() + cStay * C * L,
            logP_storage.data() + p2 * C * L,
            scratch_F_p2.data(),
            scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
        );

		// 2) Recompute F at p1 with new children {p2(new), cX}
        compute_F_vectorized(
            scratch_F_p2.data(),
            F.data() + cX * C * L,
            logP_storage.data() + p1 * C * L,
            scratch_F_p1.data(),
            scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
        );

		// If p1 is root, compute root logZ directly
		if (parent_of[p1] == -1) {
            double new_total = compute_root_logZ_vectorized(
                scratch_F_p2.data(),
                F.data() + cX * C * L,
                logP_storage.data() + p1 * C * L,
                scratch_temp.data(), scratch_max_t.data()
            );
            return new_total;
		}

		// Otherwise climb to root
		int child = p1;
        double* p_childF = scratch_childF.data();
        double* p_prev_childF = scratch_prev_childF.data();
        double* p_F_v = scratch_F_v.data();
        
        std::copy(scratch_F_p1.begin(), scratch_F_p1.end(), p_childF);

		int prev_child = -1;

		while (true) {
			const int v = parent_of[child];
			if (v == -1) {
                // child is root
				const int root_node = child;
				const int sib = other_child(root_node, prev_child);
                double new_total = compute_root_logZ_vectorized(
                    p_prev_childF,
                    F.data() + sib * C * L,
                    logP_storage.data() + root_node * C * L,
                    scratch_temp.data(), scratch_max_t.data()
                );
                return new_total;
			}

			const int sib = other_child(v, child);
            
            compute_F_vectorized(
                p_childF,
                F.data() + sib * C * L,
                logP_storage.data() + v * C * L,
                p_F_v,
                scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
            );

			prev_child = child;
            
            // Swap pointers
            double* temp = p_prev_childF;
            p_prev_childF = p_childF;
            p_childF = p_F_v;
            p_F_v = temp;
            
			child = v;
		}
	}

	// Commit the NNI: update topology and cached F/logZ across loci.
	void apply_nni(int edge_n, int which) {
		std::lock_guard<std::mutex> guard(scratch_mutex);
		// Identify key nodes from current topology
		int cnt = 0, ind = -1;
		for (int i = 0; i < m; ++i) {
			if (E[m + i] >= nTips) {
				++cnt;
				if (cnt == edge_n) { ind = i; break; }
			}
		}
		if (ind < 0) stop("edge_n out of range in apply_nni");

		const int p1 = E[ind];
		const int p2 = E[m + ind];
		const int c1 = other_child(p1, p2);
		const int c2 = children_of[p2][0];
		const int c3 = children_of[p2][1];

		const int cX    = (which == 0 ? c2 : c3);	// moves to p1
		const int cStay = (which == 0 ? c3 : c2);	// stays under p2

        // Update F
        // p2 new F
        compute_F_vectorized(
            F.data() + c1 * C * L,
            F.data() + cStay * C * L,
            logP_storage.data() + p2 * C * L,
            F.data() + p2 * C * L,
            scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
        );

        // p1 new F
        compute_F_vectorized(
            F.data() + p2 * C * L,
            F.data() + cX * C * L,
            logP_storage.data() + p1 * C * L,
            F.data() + p1 * C * L,
            scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
        );

        // Climb to root
        int child = p1;
        int prev_child = -1;

        while (true) {
            const int v = parent_of[child];
            if (v == -1) {
                // child is root
                if (prev_child == -1) {
                    // p1 is root
                    current_total_loglik_val = compute_root_logZ_vectorized(
                        F.data() + p2 * C * L,
                        F.data() + cX * C * L,
                        logP_storage.data() + p1 * C * L,
                        scratch_temp.data(), scratch_max_t.data()
                    );
                } else {
                    const int root_node = child;
                    const int sib = other_child(root_node, prev_child);
                    current_total_loglik_val = compute_root_logZ_vectorized(
                        F.data() + prev_child * C * L,
                        F.data() + sib * C * L,
                        logP_storage.data() + root_node * C * L,
                        scratch_temp.data(), scratch_max_t.data()
                    );
                }
                break;
            }

            const int sib = other_child(v, child);
            
            compute_F_vectorized(
                F.data() + child * C * L,
                F.data() + sib * C * L,
                logP_storage.data() + v * C * L,
                F.data() + v * C * L,
                scratch_u.data(), scratch_temp.data(), scratch_s.data(), scratch_max_t.data()
            );

            prev_child = child;
            child = v;
        }

		// Update adjacency: p1:{p2,cX}, p2:{c1,cStay}; parent_of for c1,cX
		auto &ch1 = children_of[p1];
		if (ch1[0] == c1) ch1[0] = cX; else ch1[1] = cX;
		auto &ch2 = children_of[p2];
		if (ch2[0] == cX) ch2[0] = c1; else ch2[1] = c1;
		parent_of[c1] = p2;
		parent_of[cX] = p1;

		// Update E's child column for the two affected edges
		for (int i = 0; i < m; ++i) {
			if (E[i] == p1 && E[m + i] == c1) { E[m + i] = cX; break; }
		}
		for (int i = 0; i < m; ++i) {
			if (E[i] == p2 && E[m + i] == cX) { E[m + i] = c1; break; }
		}

		// Keep E in postorder so the nth internal edge ordering matches nnin_cpp
		{
			arma::Col<int> E1 = E + 1;            // back to 1-indexed for reorderRcpp
			E1 = reorderRcpp(E1);
			E = E1 - 1;                            // return to 0-indexed
			root = E(m - 1);
		}
	}

	double total_loglik() const {
		return current_total_loglik_val;
	}
};

// ---- Rcpp exported wrappers around NNICache ----

SEXP nni_cache_create(arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA) {
	Rcpp::XPtr<NNICache> ptr(new NNICache(E, logP, logA), true);
	return ptr;
}

double nni_cache_loglik(SEXP xp) {
	Rcpp::XPtr<NNICache> ptr(xp);
	return ptr->total_loglik();
}

double nni_cache_delta(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
    double new_ll = ptr->compute_new_loglik(edge_n, which);
    double cur_ll = ptr->total_loglik();
    if (cur_ll == -std::numeric_limits<double>::infinity()) {
        if (new_ll > -std::numeric_limits<double>::infinity()) return std::numeric_limits<double>::infinity();
        else return 0.0;
    }
	return new_ll - cur_ll;
}

void nni_cache_apply(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
	ptr->apply_nni(edge_n, which);
}

arma::Col<int> nni_cache_current_E(SEXP xp) {
	Rcpp::XPtr<NNICache> ptr(xp);
	return ptr->E + 1; // back to 1-indexed
}

// --------------------------------------------------------------------------
// Worker struct that scores multiple trees in parallel.
// Each worker processes a subset of trees from the input list.
struct ScoreTreesWorker : public Worker {
    // Inputs are passed by const reference.
    const std::vector<arma::Col<int>>& trees;
    const std::vector< std::vector<double> >& logP;
    const std::vector<double>& logA;
    
    // Output: scores for each tree
    RVector<double> scores;
    
    // Constructor.
    ScoreTreesWorker(const std::vector<arma::Col<int>>& trees,
                     const std::vector< std::vector<double> >& logP,
                     const std::vector<double>& logA,
                     NumericVector scores)
      : trees(trees), logP(logP), logA(logA), scores(scores) {}
    
    // Operator() for processing tree indices [begin, end)
    void operator()(std::size_t begin, std::size_t end) {
        // Each iteration scores one tree.
        for (std::size_t i = begin; i < end; i++) {
            scores[i] = score_tree_bp_wrapper2(trees[i], logP, logA);
        }
    }
};

// --------------------------------------------------------------------------
// Parallel wrapper: scores multiple trees in parallel.
// Parameters:
//   trees: vector of edge matrices (each arma::Col<int>) in column-major order.
//   logP: list of flattened likelihood matrices (each in row-major order).
//   logA: flattened transition matrix (row-major order).
// Returns a vector of scores (one per tree).
// [[Rcpp::export]]
NumericVector score_trees_parallel(const std::vector<arma::Col<int>>& trees,
                                   const std::vector< std::vector<double> >& logP,
                                   const std::vector<double>& logA) {
    
    int n_trees = trees.size();
    
    // Prepare the output vector to hold scores for all trees.
    NumericVector scores(n_trees);
    
    // Create the worker that will score trees in parallel.
    ScoreTreesWorker worker(trees, logP, logA, scores);
    
    // Launch parallelFor over tree indices [0, n_trees).
    parallelFor(0, n_trees, worker);
    
    return scores;
}


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
double score_tree_bp_wrapper_multi(arma::Col<int> E,
                             const std::vector< std::vector<double> >& logP_list,
                             const std::vector< std::vector<double> >& logA_list) {
    // Compute number of loci.
    int L = logP_list.size();
    // Infer number of states from the transition matrix.
    int C = std::sqrt(logA_list[0].size());
    // Infer number of nodes from the first likelihood matrix.
    int n = logP_list[0].size() / C;
    // Compute m: number of edges.
    int m = E.n_elem / 2;

    // Adjust the edge matrix from 1-indexing (R) to 0-indexing (C++).
    E = E - 1;
    // Find the root
    int root = E(m - 1);

    double logZ = 0.0;
    for (int l = 0; l < L; l++) {
        logZ += score_tree_bp(E, logP_list[l], logA_list[l], n, C, m, root);
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
            scores[2*i] = score_tree_bp_wrapper2(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper2(Ep[1], logP, logA);
        }
    }
};

// Cached variant: efficiently scores NNI neighbours using NNICache deltas.
struct score_neighbours_cached: public Worker {
	const NNICache* cache; // read-only cache (shared across threads)
	double base_ll;        // current total log-likelihood
	RVector<double> scores;

	score_neighbours_cached(const NNICache* cache, NumericVector scores)
		: cache(cache), base_ll(cache->total_loglik()), scores(scores) {}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			const int edge_n = static_cast<int>(i) + 1; // 1-indexed internal-edge id
            
            double ll0 = cache->compute_new_loglik(edge_n, 0);
            double ll1 = cache->compute_new_loglik(edge_n, 1);
            
			scores[2*i]   = ll0;
			scores[2*i+1] = ll1;
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

// [[Rcpp::export]]
NumericVector nni_cpp_parallel_cached(arma::Col<int> E,
	const std::vector<std::vector<double>>& logP,
	const std::vector<double>& logA) {
	// Build a fresh cache for this tree
	NNICache cache(E, logP, logA);

	// Number of internal edges (same formula as the original function)
	int n = E.n_elem / 4 - 1;
	NumericVector scores(2 * n);

	// Parallel evaluation of all 2*n neighbors using cached deltas
	score_neighbours_cached worker(&cache, scores);
	parallelFor(0, n, worker);
	return scores;
}

struct score_neighbours_multi: public Worker {

    // original tree
    const arma::Col<int> E;
    const std::vector<std::vector<double>> logP;
    const std::vector<std::vector<double>> logA;
    RVector<double> scores;

    // initialize with source and destination
    score_neighbours_multi(const arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<std::vector<double>> logA, NumericVector scores): 
        E(E), logP(logP), logA(logA), scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            std::vector<arma::Col<int>> Ep = nnin_cpp(E, i+1);
            scores[2*i] = score_tree_bp_wrapper_multi(Ep[0], logP, logA);
            scores[2*i+1] = score_tree_bp_wrapper_multi(Ep[1], logP, logA);
        }
    }
};

// [[Rcpp::export]]
NumericVector nni_cpp_parallel_multi(arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<std::vector<double>> logA) {

    E = reorderRcpp(E);

    int n = E.n_elem / 4 - 1;

    NumericVector scores(2*n);

    score_neighbours_multi score_neighbours_multi(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours_multi);

    return scores;

}

/////////////////////////////////////// MCMC ////////////////////////////////////////

// [[Rcpp::export]]
std::vector<arma::Col<int>> tree_mcmc_cpp(
    arma::Col<int> E,
    const std::vector< std::vector<double> >& logP,
    const std::vector<double>& logA,
    int max_iter = 100, int seed = -1) {

    // Number of internal edges
    int n = E.n_elem / 4 - 1;

    E = reorderRcpp(E);

    std::vector<arma::Col<int>> tree_list(max_iter + 1);

    double l_0 = score_tree_bp_wrapper2(E, logP, logA);
    tree_list[0] = E;

    // random number generators
    std::mt19937 gen;
    if (seed == -1) {
        std::random_device rd;
        gen.seed(rd());
    } else {
        gen.seed(seed);
    }

    std::uniform_int_distribution<> dis1(1, n);
    std::uniform_int_distribution<> dis2(0, 1);
    std::uniform_real_distribution<> dis3(0.0, 1.0);
    
    for (int i = 0; i < max_iter; i++) {

        // propose new tree
        int r1 = dis1(gen);
        int r2 = dis2(gen);
        std::vector<arma::Col<int>> Ep = nnin_cpp(E, r1);
        arma::Col<int> E_new = std::move(Ep[r2]);
        
        double l_1 = score_tree_bp_wrapper2(E_new, logP, logA);
        
        double probab = exp(l_1 - l_0);
        double r3 = dis3(gen);

        if (r3 < probab) {
            l_0 = l_1;
            E = std::move(E_new);
        }

        tree_list[i+1] = E;

    }

    return(tree_list);
}

// [[Rcpp::export]]
std::vector<arma::Col<int>> tree_mcmc_cpp_cached(
	arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
	int max_iter = 100, int seed = -1) {

	// Number of internal edges (unchanged by reordering)
	int n = E.n_elem / 4 - 1;

	// Build message cache (handles reorder + 0-indexing internally)
	SEXP xp = nni_cache_create(E, logP, logA);

	std::vector<arma::Col<int>> tree_list(max_iter + 1);

	// Starting log-likelihood and tree
	double l_0 = nni_cache_loglik(xp);
	tree_list[0] = nni_cache_current_E(xp);

	// RNG
	std::mt19937 gen;
	if (seed == -1) {
		std::random_device rd;
		gen.seed(rd());
	} else {
		gen.seed(seed);
	}
	std::uniform_int_distribution<> dis1(1, n);
	std::uniform_int_distribution<> dis2(0, 1);
	std::uniform_real_distribution<> dis3(0.0, 1.0);

	for (int i = 0; i < max_iter; i++) {
		// propose: pick internal edge and which swap (0/1)
		int r1 = dis1(gen);
		int r2 = dis2(gen);

		// local log-likelihood delta using cached messages
		double new_ll = nni_cache_delta(xp, r1, r2); 
        
        double dl = nni_cache_delta(xp, r1, r2);

		// accept using log form (stable, avoids exp)
		double r3 = dis3(gen);
		if (std::log(r3) < dl) {
			nni_cache_apply(xp, r1, r2);
			l_0 += dl;
            l_0 = nni_cache_loglik(xp);
		}

		// store current tree (1-indexed)
		tree_list[i + 1] = nni_cache_current_E(xp);
	}

	return tree_list;
}

// Thread-safe variant: no R objects, no XPtr, safe inside RcppParallel workers.
// [[Rcpp::export]]
std::vector<arma::Col<int>> tree_mcmc_cpp_cached_threadsafe(
	arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
	int max_iter = 100, int seed = -1, bool reorder = true) {

	// Number of internal edges (unchanged by reordering)
	const int n = static_cast<int>(E.n_elem / 4) - 1;

	// Build cache locally (constructor reorders and 0-indexes internally)
	NNICache cache(E, logP, logA, reorder);

	std::vector<arma::Col<int>> tree_list(static_cast<size_t>(max_iter) + 1);

	// Starting log-likelihood and tree
	double l_0 = cache.total_loglik();
	tree_list[0] = cache.E + 1; // back to 1-indexed

	// RNG
	std::mt19937 gen;
	if (seed == -1) {
		std::random_device rd;
		gen.seed(rd());
	} else {
		gen.seed(seed);
	}
	std::uniform_int_distribution<> dis1(1, n);
	std::uniform_int_distribution<> dis2(0, 1);
	std::uniform_real_distribution<> dis3(0.0, 1.0);

	for (int i = 0; i < max_iter; ++i) {
		// propose: pick internal edge and which swap (0/1)
		const int r1 = dis1(gen);
		const int r2 = dis2(gen);

		// local log-likelihood delta using cached messages
		double new_ll = cache.compute_new_loglik(r1, r2);
        double dl;
        if (l_0 == -std::numeric_limits<double>::infinity()) {
             if (new_ll > -std::numeric_limits<double>::infinity()) dl = std::numeric_limits<double>::infinity();
             else dl = 0.0;
        } else {
             dl = new_ll - l_0;
        }

		// accept using log form (stable)
		if (std::log(dis3(gen)) < dl) {
			cache.apply_nni(r1, r2);
			l_0 = new_ll;
		}

		// store current tree (1-indexed)
		tree_list[static_cast<size_t>(i) + 1] = cache.E + 1;
	}

	return tree_list;
}

// --------------------------------------------------------------------------
// Worker struct that runs one complete MCMC chain.
// The entire chain (a vector of trees) is written to chain_results at index i.
struct TreeChainWorker : public Worker {
    // Inputs are passed by const reference.
    const arma::Col<int>& E;
    const std::vector< std::vector<double> >& logP;
    const std::vector<double>& logA;
    const int max_iter;
    
    // Output: each chain is a vector of arma::Col<int> (all iterations for that chain)
    std::vector< std::vector<arma::Col<int>> >& chain_results;
    
    // Constructor.
    TreeChainWorker(const arma::Col<int>& E,
                    const std::vector< std::vector<double> >& logP,
                    const std::vector<double>& logA,
                    int max_iter,
                    std::vector< std::vector<arma::Col<int>> >& chain_results)
      : E(E), logP(logP), logA(logA), max_iter(max_iter), chain_results(chain_results) {}
    
    // Operator() for processing chain indices [begin, end)
    void operator()(std::size_t begin, std::size_t end) {
        // Each iteration runs one full MCMC chain.
        for (std::size_t i = begin; i < end; i++) {
            // Run one chain.
            std::vector<arma::Col<int>> chain = tree_mcmc_cpp_cached_threadsafe(E, logP, logA, max_iter, static_cast<int>(i), false);
            // Store the entire chain (all iterations) in the output vector.
            chain_results[i] = std::move(chain);
        }
    }
};


// --------------------------------------------------------------------------
// Parallel wrapper: runs multiple independent MCMC chains in parallel.
// Parameters:
//   E: Edge matrix (arma::Col<int>) in column-major order (no indexing adjustments).
//   logP: list of flattened likelihood matrices (each in row-major order).
//   logA: flattened transition matrix (row-major order).
//   max_iter: number of MCMC iterations per chain.
//   nchains: number of independent chains to run in parallel.
// Returns a vector (length nchains) where each element is a vector (the chain of trees).
// [[Rcpp::export]]
std::vector< std::vector<arma::Col<int>> > tree_mcmc_parallel(arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
	int max_iter,
	int nchains) {
	// Do not adjust the indexing of E.
	// Optionally, you can reorder E if needed.
	E = reorderRcpp(E);

	// Prepare the output vector to hold nchains chains.
	std::vector< std::vector<arma::Col<int>> > chain_results(nchains);

	// Create the worker that will run chains in parallel.
	TreeChainWorker worker(E, logP, logA, max_iter, chain_results);

	// Launch parallelFor over chain indices [0, nchains).
	parallelFor(0, nchains, worker);

	return chain_results;
}


struct SeededTreeChainWorker : public Worker {
    const std::vector< arma::Col<int> >& start_edges;
    const std::vector< std::vector<double> >& logP;
    const std::vector<double>& logA;
    const std::vector<int>& max_iter_vec;
    const std::vector<int>& seeds;
    std::vector< std::vector<arma::Col<int>> >& chain_results;

    SeededTreeChainWorker(const std::vector< arma::Col<int> >& start_edges,
                          const std::vector< std::vector<double> >& logP,
                          const std::vector<double>& logA,
                          const std::vector<int>& max_iter_vec,
                          const std::vector<int>& seeds,
                          std::vector< std::vector<arma::Col<int>> >& chain_results)
        : start_edges(start_edges), logP(logP), logA(logA),
          max_iter_vec(max_iter_vec), seeds(seeds), chain_results(chain_results) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            chain_results[i] = tree_mcmc_cpp_cached_threadsafe(start_edges[i], logP, logA, max_iter_vec[i], seeds[i], false);
        }
    }
};

// [[Rcpp::export]]
std::vector< std::vector<arma::Col<int>> > tree_mcmc_parallel_seeded(std::vector< arma::Col<int> > start_edges,
    const std::vector< std::vector<double> >& logP,
    const std::vector<double>& logA,
    const std::vector<int>& max_iter_vec,
    const std::vector<int>& seeds) {

    const std::size_t nchains = start_edges.size();
    if (max_iter_vec.size() != nchains || seeds.size() != nchains) {
        stop("start_edges, max_iter_vec, and seeds must have the same length");
    }

    for (std::size_t i = 0; i < nchains; ++i) {
        start_edges[i] = reorderRcpp(start_edges[i]);
    }

    std::vector< std::vector<arma::Col<int>> > chain_results(nchains);

    SeededTreeChainWorker worker(start_edges, logP, logA, max_iter_vec, seeds, chain_results);
    parallelFor(0, static_cast<std::size_t>(nchains), worker);

    return chain_results;
}

struct LocusChunkComputeWorker : public Worker {
	const std::vector<NNICache*>& caches;
	const int edge_n;
	const int which;
	std::vector<double>& partial_totals;

	LocusChunkComputeWorker(const std::vector<NNICache*>& caches_, int edge_n_, int which_, std::vector<double>& partial_totals_)
		: caches(caches_), edge_n(edge_n_), which(which_), partial_totals(partial_totals_) {}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			partial_totals[i] = caches[i]->compute_new_loglik(edge_n, which);
		}
	}
};

struct LocusChunkApplyWorker : public Worker {
	const std::vector<NNICache*>& caches;
	const int edge_n;
	const int which;

	LocusChunkApplyWorker(const std::vector<NNICache*>& caches_, int edge_n_, int which_)
		: caches(caches_), edge_n(edge_n_), which(which_) {}

	void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; ++i) {
			caches[i]->apply_nni(edge_n, which);
		}
	}
};

std::vector<arma::Col<int>> run_locus_parallel_chain(
	arma::Col<int> start_E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
	int max_iter,
	int seed) {

	if (logP.empty()) {
		stop("logP must contain at least one locus");
	}

	arma::Col<int> ordered_E = reorderRcpp(start_E);
	const int total_loci = static_cast<int>(logP.size());
	int chunk_count = defaultNumThreads();
	if (chunk_count <= 0) {
		const unsigned hw = std::thread::hardware_concurrency();
		chunk_count = static_cast<int>(hw == 0 ? 1u : hw);
	}
	chunk_count = std::max(1, std::min(chunk_count, total_loci));

	std::vector<std::unique_ptr<NNICache>> cache_storage;
	cache_storage.reserve(chunk_count);
	std::vector<NNICache*> cache_ptrs;
	cache_ptrs.reserve(chunk_count);
	std::vector<int> offsets(chunk_count + 1);
	for (int i = 0; i <= chunk_count; ++i) {
		offsets[i] = (static_cast<long long>(i) * total_loci) / chunk_count;
	}

	for (int chunk = 0; chunk < chunk_count; ++chunk) {
		std::vector<int> locus_subset;
		const int begin = offsets[chunk];
		const int end = offsets[chunk + 1];
		locus_subset.reserve(end - begin);
		for (int idx = begin; idx < end; ++idx) locus_subset.push_back(idx);
		cache_storage.emplace_back(new NNICache(ordered_E, logP, logA, locus_subset, false));
		cache_ptrs.push_back(cache_storage.back().get());
	}

	std::vector<arma::Col<int>> tree_list(static_cast<std::size_t>(max_iter) + 1);
	if (cache_ptrs.empty()) {
		stop("Failed to initialize locus caches");
	}
	const int n_internal = static_cast<int>(ordered_E.n_elem / 4) - 1;
	std::vector<double> partial_totals(chunk_count, 0.0);

	double current_ll = 0.0;
	for (NNICache* cache : cache_ptrs) {
		current_ll += cache->total_loglik();
	}
	
	std::mt19937 gen;
	if (seed == -1) {
		std::random_device rd;
		gen.seed(rd());
	} else {
		gen.seed(seed);
	}
	std::uniform_int_distribution<> dis1(1, n_internal);
	std::uniform_int_distribution<> dis2(0, 1);
	std::uniform_real_distribution<> dis3(0.0, 1.0);

	tree_list[0] = cache_ptrs.front()->E + 1;

	auto compute_new_total = [&](int edge_n, int which) {
		if (cache_ptrs.size() == 1) {
			return cache_ptrs[0]->compute_new_loglik(edge_n, which);
		}
		LocusChunkComputeWorker worker(cache_ptrs, edge_n, which, partial_totals);
		parallelFor(0, cache_ptrs.size(), worker);
		double sum = 0.0;
		for (double val : partial_totals) sum += val;
		return sum;
	};

	auto apply_move = [&](int edge_n, int which) {
		if (cache_ptrs.size() == 1) {
			cache_ptrs[0]->apply_nni(edge_n, which);
			return;
		}
		LocusChunkApplyWorker worker(cache_ptrs, edge_n, which);
		parallelFor(0, cache_ptrs.size(), worker);
	};

	for (int iter = 0; iter < max_iter; ++iter) {
		const int edge_n = dis1(gen);
		const int which = dis2(gen);
		double new_total = compute_new_total(edge_n, which);
		double delta;
		if (current_ll == -std::numeric_limits<double>::infinity()) {
			if (new_total > -std::numeric_limits<double>::infinity()) delta = std::numeric_limits<double>::infinity();
			else delta = 0.0;
		} else {
			delta = new_total - current_ll;
		}
		if (std::log(dis3(gen)) < delta) {
			apply_move(edge_n, which);
			current_ll = new_total;
		}
		tree_list[static_cast<std::size_t>(iter) + 1] = cache_ptrs.front()->E + 1;
	}

	return tree_list;
}

// [[Rcpp::export]]
std::vector< std::vector<arma::Col<int>> > tree_mcmc_parallel_seeded_locus(
	std::vector< arma::Col<int> > start_edges,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA,
	const std::vector<int>& max_iter_vec,
	const std::vector<int>& seeds) {

	const std::size_t nchains = start_edges.size();
	if (max_iter_vec.size() != nchains || seeds.size() != nchains) {
		stop("start_edges, max_iter_vec, and seeds must have the same length");
	}

	std::vector< std::vector<arma::Col<int>> > chain_results(nchains);
	for (std::size_t i = 0; i < nchains; ++i) {
		chain_results[i] = run_locus_parallel_chain(start_edges[i], logP, logA, max_iter_vec[i], seeds[i]);
	}

	return chain_results;
}


////////////////////////  EM Helpers //////////////////////// 

struct ComputeNodeEdgeStatsWorker : public Worker {
	const int m;
	const int n;
	const int C;
	const int root;
	const std::vector<int>& tip_ids;
	const int* parent_ptr;
	const int* child_ptr;
	const std::vector< std::array<int,2> >& children_of;
	const std::vector<double>& row_maxA;
	const std::vector<double>& expA_shifted;
	const std::vector<double>& expA;
	const std::vector< std::vector<double> >& logP_list;
	std::vector< std::vector<double> >& leaf_storage;
	std::vector< std::vector<double> >& edge_counts_local;
	std::vector<double>& logZ_vals;

	ComputeNodeEdgeStatsWorker(
		int m_, int n_, int C_, int root_, const std::vector<int>& tip_ids_,
		const int* parent_ptr_, const int* child_ptr_,
		const std::vector< std::array<int,2> >& children_of_,
		const std::vector<double>& row_maxA_,
		const std::vector<double>& expA_shifted_,
		const std::vector<double>& expA_,
		const std::vector< std::vector<double> >& logP_list_,
		std::vector< std::vector<double> >& leaf_storage_,
		std::vector< std::vector<double> >& edge_counts_local_,
		std::vector<double>& logZ_vals_)
		: m(m_), n(n_), C(C_), root(root_), tip_ids(tip_ids_),
		  parent_ptr(parent_ptr_), child_ptr(child_ptr_),
		  children_of(children_of_),
		  row_maxA(row_maxA_), expA_shifted(expA_shifted_), expA(expA_),
		  logP_list(logP_list_),
		  leaf_storage(leaf_storage_), edge_counts_local(edge_counts_local_),
		  logZ_vals(logZ_vals_) {}

	void operator()(std::size_t begin, std::size_t end) {
		const double neg_inf = -std::numeric_limits<double>::infinity();
		const std::size_t nC = static_cast<std::size_t>(n) * static_cast<std::size_t>(C);

		std::vector<double> Uin(nC);
		std::vector<double> up_to_parent(nC);
		std::vector<double> down_in(nC);
		std::vector<double> temp(C);
		std::vector<double> u(C);
		std::vector<double> S_parent(C);
		std::vector<double> down_msg(C);
		std::vector<double> a(C);
		std::vector<double> b(C);
		std::vector<double> node_tmp(nC);
		std::vector<int> stack;
		stack.reserve(static_cast<std::size_t>(n));

		for (std::size_t idx_l = begin; idx_l < end; ++idx_l) {
			const double* const P = logP_list[idx_l].data();
			double logZ_root = 0.0;

			std::fill(Uin.begin(), Uin.end(), 0.0);
			std::fill(up_to_parent.begin(), up_to_parent.end(), 0.0);
			std::fill(down_in.begin(), down_in.end(), 0.0);

		for (int i = 0; i < m; ++i) {
				const int par = parent_ptr[i];
				const int node = child_ptr[i];
				double max_t = neg_inf;
				int offP = node;
				double* const Uin_node = Uin.data() + static_cast<std::size_t>(node) * C;
				double* const up_node = up_to_parent.data() + static_cast<std::size_t>(node) * C;
				double* const Uin_parent = Uin.data() + static_cast<std::size_t>(par) * C;
				for (int c = 0; c < C; ++c, offP += n) {
					const double v = P[offP] + Uin_node[c];
					temp[c] = v;
					if (v > max_t) max_t = v;
				}
				for (int c = 0; c < C; ++c) u[c] = std::exp(temp[c] - max_t);

				const double* rowA = expA_shifted.data();
				for (int r = 0; r < C; ++r, rowA += C) {
					double s = 0.0;
					for (int j = 0; j < C; ++j) s += rowA[j] * u[j];
					const double f = row_maxA[r] + max_t + std::log(s);
					up_node[r] = f;
					Uin_parent[r] += f;
				}
			}

			stack.clear();
			stack.push_back(root);
			while (!stack.empty()) {
				const int u_node = stack.back();
				stack.pop_back();
				double* const down_parent = down_in.data() + static_cast<std::size_t>(u_node) * C;
				double* const Uin_parent = Uin.data() + static_cast<std::size_t>(u_node) * C;
				const auto& ch = children_of[u_node];
				for (int t = 0; t < 2; ++t) {
					const int v = ch[t];
					if (v < 0) continue;
					double* const down_child = down_in.data() + static_cast<std::size_t>(v) * C;
					double* const up_child = up_to_parent.data() + static_cast<std::size_t>(v) * C;
					for (int r = 0; r < C; ++r) {
						const double P_ur = (u_node == root)
							? (r == 0 ? 0.0 : neg_inf)
							: P[static_cast<std::size_t>(r) * n + u_node];
						S_parent[r] = P_ur + down_parent[r] + Uin_parent[r] - up_child[r];
					}
					compute_down_msg_fast(S_parent.data(), expA.data(), C, down_msg.data(), u.data());
					for (int k = 0; k < C; ++k) down_child[k] = down_msg[k];
					stack.push_back(v);
				}
			}

			for (int v = 0; v < n; ++v) {
				double maxW = neg_inf;
				const double* const down_v = down_in.data() + static_cast<std::size_t>(v) * C;
				const double* const Uin_v = Uin.data() + static_cast<std::size_t>(v) * C;
		for (int c = 0; c < C; ++c) {
					const double P_vc = (v == root)
						? (c == 0 ? 0.0 : neg_inf)
						: P[static_cast<std::size_t>(c) * n + v];
					temp[c] = P_vc + Uin_v[c] + down_v[c];
					if (temp[c] > maxW) maxW = temp[c];
				}
				double denom = 0.0;
		for (int c = 0; c < C; ++c) {
					const double val = std::exp(temp[c] - maxW);
					node_tmp[static_cast<std::size_t>(v) + static_cast<std::size_t>(c) * n] = val;
					denom += val;
				}
				const double inv = 1.0 / denom;
				for (int c = 0; c < C; ++c) node_tmp[static_cast<std::size_t>(v) + static_cast<std::size_t>(c) * n] *= inv;
				if (v == root) logZ_root = maxW + std::log(denom);
			}

			double* const leaf_out = leaf_storage[idx_l].data();
			for (std::size_t tip_idx = 0; tip_idx < tip_ids.size(); ++tip_idx) {
				const int node_id = tip_ids[tip_idx];
				double sum = 0.0;
				for (int c = 0; c < C; ++c) {
					const double val = node_tmp[static_cast<std::size_t>(node_id) + static_cast<std::size_t>(c) * n];
					leaf_out[tip_idx + static_cast<std::size_t>(c) * tip_ids.size()] = val;
					sum += val;
				}
				sum = (sum > 0.0) ? sum : 1.0;
				for (int c = 0; c < C; ++c) {
					leaf_out[tip_idx + static_cast<std::size_t>(c) * tip_ids.size()] /= sum;
				}
			}
			logZ_vals[idx_l] = logZ_root;

			double* const edge_out = edge_counts_local[idx_l].data();
			std::fill(edge_out, edge_out + static_cast<std::size_t>(C) * static_cast<std::size_t>(C), 0.0);

		for (int i = 0; i < m; ++i) {
				const int u_node = parent_ptr[i];
				const int v = child_ptr[i];
				const double* const down_u = down_in.data() + static_cast<std::size_t>(u_node) * C;
				const double* const Uin_u = Uin.data() + static_cast<std::size_t>(u_node) * C;
				const double* const up_v = up_to_parent.data() + static_cast<std::size_t>(v) * C;
				const double* const down_v = down_in.data() + static_cast<std::size_t>(v) * C;
				const double* const Uin_v = Uin.data() + static_cast<std::size_t>(v) * C;
				for (int r = 0; r < C; ++r) {
					const double P_ur = (u_node == root)
						? (r == 0 ? 0.0 : neg_inf)
						: P[static_cast<std::size_t>(r) * n + u_node];
					S_parent[r] = P_ur + down_u[r] + Uin_u[r] - up_v[r];
				}
				for (int k = 0; k < C; ++k) down_msg[k] = P[static_cast<std::size_t>(k) * n + v] + Uin_v[k];
				double maxWr = neg_inf;
				for (int r = 0; r < C; ++r) if (S_parent[r] > maxWr) maxWr = S_parent[r];
				for (int r = 0; r < C; ++r) a[r] = std::exp(S_parent[r] - maxWr);
				double maxY = neg_inf;
				for (int k = 0; k < C; ++k) if (down_msg[k] > maxY) maxY = down_msg[k];
				for (int k = 0; k < C; ++k) b[k] = std::exp(down_msg[k] - maxY);
				double denom = 0.0;
				for (int r = 0; r < C; ++r) {
					double s = 0.0;
					const double* const Arow = expA.data() + static_cast<std::size_t>(r) * C;
					for (int k = 0; k < C; ++k) s += Arow[k] * b[k];
					denom += a[r] * s;
				}
				const double inv = 1.0 / denom;
				for (int r = 0; r < C; ++r) {
					const double* const Arow = expA.data() + static_cast<std::size_t>(r) * C;
					for (int k = 0; k < C; ++k) {
						const double prob = a[r] * Arow[k] * b[k] * inv;
						edge_out[static_cast<std::size_t>(r) * C + k] += prob;
					}
				}
			}
		}
	}
};


// [[Rcpp::export]]
Rcpp::List compute_node_edge_stats_bp2(
	arma::Col<int> E,
                             const std::vector< std::vector<double> >& logP_list,
	const std::vector<double>& logA)
{
	const int L = static_cast<int>(logP_list.size());
	const int C = static_cast<int>(std::sqrt(logA.size()));
	if (L <= 0) Rcpp::stop("logP_list must be non-empty");
	const int n = static_cast<int>(logP_list[0].size() / C);
	const int m = static_cast<int>(E.n_elem / 2);

	E = reorderRcpp(E);
	E -= 1;
    int root = E(m - 1);

	std::vector<int> parent_of(n, -1);
	std::vector< std::array<int,2> > children_of(n, std::array<int,2>{-1, -1});
	for (int i = 0; i < m; ++i) {
		const int p = E[i];
		const int c = E[m + i];
		parent_of[c] = p;
		if (children_of[p][0] == -1) children_of[p][0] = c; else children_of[p][1] = c;
	}

	std::vector<int> tip_ids;
	for (int v = 0; v < n; ++v) {
		if (children_of[v][0] == -1 && children_of[v][1] == -1) tip_ids.push_back(v);
	}
	const int nTips = static_cast<int>(tip_ids.size());

	std::vector< std::vector<double> > leaf_storage(
		static_cast<std::size_t>(L),
		std::vector<double>(static_cast<std::size_t>(nTips) * static_cast<std::size_t>(C)));
	std::vector< std::vector<double> > edge_counts_local(
		static_cast<std::size_t>(L),
		std::vector<double>(static_cast<std::size_t>(C) * static_cast<std::size_t>(C), 0.0));
	std::vector<double> logZ_vals(static_cast<std::size_t>(L), 0.0);

	std::vector<double> row_maxA(C);
	std::vector<double> expA_shifted(static_cast<size_t>(C) * C);
	for (int r = 0; r < C; ++r) {
		const double* Arow = logA.data() + static_cast<size_t>(r) * C;
		double mr = Arow[0];
		for (int j = 1; j < C; ++j) if (Arow[j] > mr) mr = Arow[j];
		row_maxA[r] = mr;
		double* out = expA_shifted.data() + static_cast<size_t>(r) * C;
		for (int j = 0; j < C; ++j) out[j] = std::exp(Arow[j] - mr);
	}
	std::vector<double> expA(static_cast<size_t>(C) * C);
	for (size_t i = 0; i < expA.size(); ++i) expA[i] = std::exp(logA[i]);

	ComputeNodeEdgeStatsWorker worker(
		m, n, C, root, tip_ids,
		E.memptr(), E.memptr() + m,
		children_of,
		row_maxA,
		expA_shifted,
		expA,
		logP_list,
		leaf_storage,
		edge_counts_local,
		logZ_vals);

	parallelFor(0, static_cast<std::size_t>(L), worker);

	double logZ_total = std::accumulate(logZ_vals.begin(), logZ_vals.end(), 0.0);

	Rcpp::List leaf_list(L);
			for (int l = 0; l < L; ++l) {
		Rcpp::NumericMatrix leaf_mat(nTips, C);
		const std::vector<double>& leaf_vec = leaf_storage[l];
		for (int v = 0; v < nTips; ++v) {
			for (int c = 0; c < C; ++c) {
				leaf_mat(v, c) = leaf_vec[static_cast<std::size_t>(v) + static_cast<std::size_t>(c) * nTips];
			}
		}
		leaf_list[l] = leaf_mat;
	}

	std::vector<double> edge_counts(static_cast<std::size_t>(C) * static_cast<std::size_t>(C), 0.0);
	for (const auto& vec : edge_counts_local) {
		for (std::size_t idx = 0; idx < edge_counts.size(); ++idx) edge_counts[idx] += vec[idx];
	}
	Rcpp::NumericMatrix edge_counts_mat(C, C);
	for (int r = 0; r < C; ++r) {
		for (int k = 0; k < C; ++k) edge_counts_mat(r, k) = edge_counts[static_cast<std::size_t>(r) * C + k];
	}

	return Rcpp::List::create(
		Rcpp::Named("leaf_beliefs") = leaf_list,
		Rcpp::Named("edge_counts") = edge_counts_mat,
		Rcpp::Named("logZ") = logZ_total
	);
}