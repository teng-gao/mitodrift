#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <limits>
#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>
#include <array>

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

	// per-locus cached child→parent contributions F (n*C) and root logZ
	std::vector< std::vector<double> > F;		// L × (n*C)
	std::vector<double> logZ;					// L

	NNICache(arma::Col<int> E_in,
		const std::vector< std::vector<double> >& logP_in,
		const std::vector<double>& logA_in) {

		L = static_cast<int>(logP_in.size());
		C = static_cast<int>(std::sqrt(logA_in.size()));
		n = static_cast<int>(logP_in[0].size() / C);

		// Bring E to postorder and 0-indexed
		E_in = reorderRcpp(E_in);
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

		// copy likelihoods
		logP_list = logP_in;

		// allocate caches
		F.assign(L, std::vector<double>(static_cast<size_t>(n) * C, 0.0));
		logZ.assign(L, 0.0);

		// initialize F and logZ per locus with one postorder pass
		std::vector<double> msg_nm(static_cast<size_t>(n) * C, 0.0);
		std::vector<double> temp(C), u(C);
		for (int l = 0; l < L; ++l) {
			std::fill(msg_nm.begin(), msg_nm.end(), 0.0);
			double* Fl = F[l].data();
			const double* P = logP_list[l].data();

			for (int i = 0; i < m; ++i) {
				const int par  = E[i];
				const int node = E[m + i];

				// temp[c] = logP[c*n + node] + sum_children(node)[c]
				double max_t = -std::numeric_limits<double>::infinity();
				int offP = node;
				const int base_node = node * C;
				for (int c = 0; c < C; ++c, offP += n) {
					const double v = P[offP] + msg_nm[base_node + c];
					temp[c] = v;
					if (v > max_t) max_t = v;
				}
				for (int c = 0; c < C; ++c) u[c] = std::exp(temp[c] - max_t);

				// F[node,·]
				const double* rowA = expA_shifted.data();
				const int base_par = par * C; // where we'll accumulate into parent
				for (int r = 0; r < C; ++r, rowA += C) {
					double s = 0.0;
					for (int j = 0; j < C; ++j) s += rowA[j] * u[j];
					const double f = row_maxA[r] + max_t + std::log(s);
					Fl[base_node + r] = f;
					msg_nm[base_par + r] += f; // accumulate into parent message
				}
			}
			// root marginal
			for (int c = 0, offP = root; c < C; ++c, offP += n) temp[c] = P[offP] + msg_nm[root * C + c];
			logZ[l] = logsumexp_array(temp.data(), C);
		}
	}

	inline int other_child(int parent, int child) const {
		const auto &ch = children_of[parent];
		return (ch[0] == child) ? ch[1] : ch[0];
	}

	// Compute F_new for a node given S[c] = sum over child's F[c] and P(node,·).
	inline void compute_F_from_sum(const double* S, const double* P_node, double* outF,
		std::vector<double>& temp, std::vector<double>& u) const {
		double max_t = -std::numeric_limits<double>::infinity();
		for (int c = 0; c < C; ++c) {
			const double v = P_node[c * n] + S[c]; // P is row-major: P[c*n + node]
			// We'll not use this form; caller will pass P offset-corrected for node.
			(void)v;
		}
		// The caller passes P_node such that P_node[c] == logP[c*n + node].
		max_t = -std::numeric_limits<double>::infinity();
		for (int c = 0; c < C; ++c) {
			if (P_node[c] + S[c] > max_t) max_t = P_node[c] + S[c];
		}
		for (int c = 0; c < C; ++c) u[c] = std::exp((P_node[c] + S[c]) - max_t);

		const double* rowA = expA_shifted.data();
		for (int r = 0; r < C; ++r, rowA += C) {
			double s = 0.0;
			for (int j = 0; j < C; ++j) s += rowA[j] * u[j];
			outF[r] = row_maxA[r] + max_t + std::log(s);
		}
	}

	// Compute new logZ for locus l if we perform the NNI "which" (0 or 1) at the nth internal edge.
	double new_logZ_for_locus(int edge_n, int which, int l) const {
		// Locate the nth internal edge (child > nTips)
		int cnt = 0, ind = -1;
		for (int i = 0; i < m; ++i) {
			if (E[m + i] > nTips) {
				++cnt;
				if (cnt == edge_n) { ind = i; break; }
			}
		}
		if (ind < 0) stop("edge_n out of range in new_logZ_for_locus");

		const int p1 = E[ind];
		const int p2 = E[m + ind];
		const int c1 = other_child(p1, p2);
		const int c2 = children_of[p2][0];
		const int c3 = children_of[p2][1];

		const int cX    = (which == 0 ? c2 : c3);   // moves to p1
		const int cStay = (which == 0 ? c3 : c2);   // stays under p2

		const double* Fl = F[l].data();
		const double* P  = logP_list[l].data();

		std::vector<double> S(C), temp(C), u(C);

		// 1) Recompute F at p2 with new children {c1, cStay}
		std::vector<double> F_p2(C);
		for (int c = 0; c < C; ++c) S[c] = Fl[c1 * C + c] + Fl[cStay * C + c];
		std::vector<double> Pnode_p2(C); for (int k = 0; k < C; ++k) Pnode_p2[k] = P[k * n + p2];
		compute_F_from_sum(S.data(), Pnode_p2.data(), F_p2.data(), temp, u);

		// 2) Recompute F at p1 with new children {p2(new), cX}
		std::vector<double> F_p1(C);
		for (int c = 0; c < C; ++c) S[c] = F_p2[c] + Fl[cX * C + c];
		std::vector<double> Pnode_p1(C); for (int k = 0; k < C; ++k) Pnode_p1[k] = P[k * n + p1];
		compute_F_from_sum(S.data(), Pnode_p1.data(), F_p1.data(), temp, u);

		// If p1 is root, compute root logZ directly from updated children {p2(new), cX}
		if (parent_of[p1] == -1) {
			for (int c = 0; c < C; ++c) S[c] = F_p2[c] + Fl[cX * C + c];
			for (int k = 0; k < C; ++k) temp[k] = P[k * n + p1] + S[k];
			return logsumexp_array(temp.data(), C);
		}

		// Otherwise climb to root, carrying the UPDATED CHILD message
		int child = p1;
		std::vector<double> childF = F_p1; // message from 'child' to its parent
		int prev_child = -1;               // the child under the next parent (root case uses this)
		std::vector<double> prev_childF;   // its updated message to that parent

		while (true) {
			const int v = parent_of[child];
			if (v == -1) {
				// child is root; combine UPDATED message from prev_child with its sibling at root
				const int root_node = child;
				const int sib = other_child(root_node, prev_child);
				for (int c = 0; c < C; ++c) S[c] = prev_childF[c] + Fl[sib * C + c];
				std::vector<double> Pnode_root(C); for (int k = 0; k < C; ++k) Pnode_root[k] = P[k * n + root_node];
				for (int c = 0; c < C; ++c) temp[c] = Pnode_root[c] + S[c];
				return logsumexp_array(temp.data(), C);
			}

			const int sib = other_child(v, child);
			for (int c = 0; c < C; ++c) S[c] = childF[c] + Fl[sib * C + c];
			std::vector<double> Pnode_v(C); for (int k = 0; k < C; ++k) Pnode_v[k] = P[k * n + v];

			std::vector<double> F_v(C);
			compute_F_from_sum(S.data(), Pnode_v.data(), F_v.data(), temp, u);

			prev_child = child;
			prev_childF = childF; // copy small C-vector
			child = v;
			childF.swap(F_v);
		}
	}

	// Commit the NNI: update topology and cached F/logZ across loci.
	void apply_nni(int edge_n, int which) {
		// Identify key nodes from current topology
		int cnt = 0, ind = -1;
		for (int i = 0; i < m; ++i) {
			if (E[m + i] > nTips) {
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

		// Update cached F along the path for every locus
		std::vector<double> S(C), curF(C), temp(C), u(C);
		for (int l = 0; l < L; ++l) {
			double* Fl = F[l].data();
			const double* P  = logP_list[l].data();

			// p2 new F
			for (int c = 0; c < C; ++c) S[c] = Fl[c1 * C + c] + Fl[cStay * C + c];
			std::vector<double> Pnode(C);
			for (int k = 0; k < C; ++k) Pnode[k] = P[k * n + p2];
			compute_F_from_sum(S.data(), Pnode.data(), curF.data(), temp, u);
			for (int r = 0; r < C; ++r) Fl[p2 * C + r] = curF[r];

			// p1 new F
			for (int c = 0; c < C; ++c) S[c] = curF[c] + Fl[cX * C + c];
			for (int k = 0; k < C; ++k) Pnode[k] = P[k * n + p1];
			compute_F_from_sum(S.data(), Pnode.data(), curF.data(), temp, u);
			for (int r = 0; r < C; ++r) Fl[p1 * C + r] = curF[r];

			// Climb to root, carrying the UPDATED child's message
			int child = p1;
			std::vector<double> childF = curF; // message from 'child' to its parent
			int prev_child = -1;
			std::vector<double> prev_childF;

			while (true) {
				const int v = parent_of[child];
				if (v == -1) {
					// If p1 was root, combine {p2(new), cX}; else combine {prev_child(updated), its sibling}
					std::vector<double> S(C);
					if (prev_child == -1) {
						for (int c = 0; c < C; ++c) S[c] = Fl[p2 * C + c] /* F_p2 stored above? see below */ + Fl[cX * C + c];
						// We must use the freshly computed F for p2 (child->root). That's `curF` when child==p2 before moving to p1,
						// but here we no longer have it. Recompute quickly:
						std::vector<double> Sroot(C), temp2(C), u2(C), Fp2_again(C);
						for (int c = 0; c < C; ++c) Sroot[c] = Fl[c1 * C + c] + Fl[cStay * C + c];
						std::vector<double> Pnode_p2b(C); for (int k = 0; k < C; ++k) Pnode_p2b[k] = P[k * n + p2];
						compute_F_from_sum(Sroot.data(), Pnode_p2b.data(), Fp2_again.data(), temp2, u2);
						for (int c = 0; c < C; ++c) S[c] = Fp2_again[c] + Fl[cX * C + c];
						std::vector<double> Pnode_root(C); for (int k = 0; k < C; ++k) Pnode_root[k] = P[k * n + v /* root == p1 */ + 0];
						for (int c = 0; c < C; ++c) temp[c] = Pnode_root[c] + S[c];
						logZ[l] = logsumexp_array(temp.data(), C);
						break;
					}

					const int root_node = child; // child == root here
					const int sib = other_child(root_node, prev_child);
					for (int c = 0; c < C; ++c) S[c] = prev_childF[c] + Fl[sib * C + c];
					std::vector<double> Pnode_root(C); for (int k = 0; k < C; ++k) Pnode_root[k] = P[k * n + root_node];
					for (int c = 0; c < C; ++c) temp[c] = Pnode_root[c] + S[c];
					logZ[l] = logsumexp_array(temp.data(), C);
					break;
				}

				const int sib = other_child(v, child);
				for (int c = 0; c < C; ++c) S[c] = childF[c] + Fl[sib * C + c];
				std::vector<double> Pnode_v(C); for (int k = 0; k < C; ++k) Pnode_v[k] = P[k * n + v];

				std::vector<double> F_v(C);
				compute_F_from_sum(S.data(), Pnode_v.data(), F_v.data(), temp, u);
				for (int r = 0; r < C; ++r) Fl[v * C + r] = F_v[r];

				prev_child = child;
				prev_childF = childF;
				child = v;
				childF.swap(F_v);
			}
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
		double s = 0.0; for (int l = 0; l < L; ++l) s += logZ[l]; return s;
	}
};

// ---- Rcpp exported wrappers around NNICache ----

// [[Rcpp::export]]
SEXP nni_cache_create(arma::Col<int> E,
	const std::vector< std::vector<double> >& logP,
	const std::vector<double>& logA) {
	Rcpp::XPtr<NNICache> ptr(new NNICache(E, logP, logA), true);
	return ptr;
}

// [[Rcpp::export]]
double nni_cache_loglik(SEXP xp) {
	Rcpp::XPtr<NNICache> ptr(xp);
	return ptr->total_loglik();
}

// [[Rcpp::export]]
double nni_cache_delta(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
	double delta = 0.0;
	for (int l = 0; l < ptr->L; ++l) {
		delta += ptr->new_logZ_for_locus(edge_n, which, l) - ptr->logZ[l];
	}
	return delta;
}

// [[Rcpp::export]]
void nni_cache_apply(SEXP xp, int edge_n, int which) {
	Rcpp::XPtr<NNICache> ptr(xp);
	ptr->apply_nni(edge_n, which);
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
NumericVector nni_cpp_parallel(arma::Col<int> E, const std::vector<std::vector<double>> logP, const std::vector<double> logA) {

    E = reorderRcpp(E);

    int n = E.n_elem / 4 - 1;

    NumericVector scores(2*n);

    score_neighbours score_neighbours(E, logP, logA, scores);

    parallelFor(0, n, score_neighbours);

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

        // double dl = l_1 - l_0;
        // Rcpp::Rcout << "dl: " << dl << std::endl;

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
		double dl = nni_cache_delta(xp, r1, r2);

		// accept using log form (stable, avoids exp)
		double r3 = dis3(gen);
		if (std::log(r3) < dl) {
			nni_cache_apply(xp, r1, r2);
			l_0 += dl;
		}

		// store current tree (1-indexed)
		tree_list[i + 1] = nni_cache_current_E(xp);
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
            std::vector<arma::Col<int>> chain = tree_mcmc_cpp(E, logP, logA, max_iter, i);
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