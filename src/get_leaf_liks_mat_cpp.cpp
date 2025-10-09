// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Leaf likelihood matrix via C++ backend
//'
//' @param amat Integer matrix of alt counts (variants x cells).
//' @param dmat Integer matrix of depths (same dimensions as `amat`).
//' @param vafs Numeric vector of VAF grid points (names optional).
//' @param eps Error rate to add to each VAF bin (default 0).
//' @param ncores Number of threads to use (default 1).
//' @param log Whether to return log-likelihoods instead of probabilities.
//' @return List of matrices, one per variant, with rows = VAF bins and columns = cells.
//' @export
// [[Rcpp::export]]
List get_leaf_liks_mat_cpp(const IntegerMatrix &amat,
						   const IntegerMatrix &dmat,
						   const NumericVector &vafs,
						   double eps = 0.0,
						   int ncores = 1,
						   bool log = false) {

	// dims check
	if (amat.nrow() != dmat.nrow() || amat.ncol() != dmat.ncol()) {
		stop("amat and dmat must have the same dimensions");
	}
	const int K      = vafs.size();     // #VAF bins
	const int nvars  = amat.nrow();     // #variants
	const int ncells = amat.ncol();     // #cells
	const int nrow   = nvars;           // for column-major indexing

	// p <- pmin(vafs + eps, 1 - eps)
	std::vector<double> p(K);
	const double cap = 1.0 - eps;
	for (int k = 0; k < K; ++k) {
		double pk = vafs[k] + eps;
		p[k] = (pk > cap) ? cap : pk;
	}

	// raw pointers (column-major)
	const int *A = amat.begin();
	const int *D = dmat.begin();

	// capture mutation (row) names to set as list names and cell names for matrix dimnames
	List dnA = amat.attr("dimnames");
	CharacterVector variants = as<CharacterVector>(dnA[0]);
	CharacterVector cell_names = as<CharacterVector>(dnA[1]);
	CharacterVector vaf_names = as<CharacterVector>(Rf_coerceVector(vafs, STRSXP));

	List dm(2);
	dm[0] = vaf_names;   // rows: VAF bins
	dm[1] = cell_names;  // cols: cells

	List out(nvars);

	for (int i = 0; i < nvars; ++i) {
		NumericMatrix m(K, ncells);      // rows=VAF bins, cols=cells
		m.attr("dimnames") = dm;
		double *M = m.begin();

		for (int j = 0; j < ncells; ++j) {
			const int idx = i + j * nrow; // (row=i, col=j)
			const int x   = A[idx];
			const int n   = D[idx];

			double *colptr = M + j * K;   // column j start
			for (int k = 0; k < K; ++k) {
				colptr[k] = R::dbinom(static_cast<double>(x),
									  static_cast<double>(n),
									  p[k],
									  log ? 1 : 0);
			}
		}

		out[i] = m;
	}

	// set list names to mutation IDs if present
	if (variants.size() == nvars) {
		out.attr("names") = variants;
	}
	return out;
}