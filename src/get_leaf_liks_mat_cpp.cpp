// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_leaf_liks_mat_cpp(const IntegerMatrix &amat,
					   	   const IntegerMatrix &dmat,
					   	   const NumericVector &vafs,
					   	   double eps = 0.0,
					   	   int ncores = 1, // kept for signature parity; not used
					   	   bool log = false) {

	// dimension checks
	if (amat.nrow() != dmat.nrow() || amat.ncol() != dmat.ncol()) {
		stop("amat and dmat must have the same dimensions");
	}

	const int K = vafs.size();
	const int nvars = amat.nrow();
	const int ncells = amat.ncol();
	const int nrow = nvars; // alias for clarity in linear indexing

	// p <- pmin(vafs + eps, 1 - eps)
	NumericVector p(K);
	const double cap = 1.0 - eps;
	for (int k = 0; k < K; ++k) {
		double pk = vafs[k] + eps;
		p[k] = (pk > cap) ? cap : pk;
	}

	// collect names
	List dnA = amat.attr("dimnames");
	CharacterVector variants = (dnA.isNULL() || dnA[0] == R_NilValue) ? CharacterVector() : as<CharacterVector>(dnA[0]);
	CharacterVector cell_names = (dnA.isNULL() || dnA.size() < 2 || dnA[1] == R_NilValue) ? CharacterVector() : as<CharacterVector>(dnA[1]);
	CharacterVector vaf_names = as<CharacterVector>(Rf_coerceVector(vafs, STRSXP));

	// raw pointers for faster access (column-major layout)
	const int *A = amat.begin();
	const int *D = dmat.begin();
	const double *P = p.begin();

	List out(nvars);

	for (int i = 0; i < nvars; ++i) {
		NumericMatrix m(K, ncells);
		double *M = m.begin();

		for (int j = 0; j < ncells; ++j) {
			const int idx = i + j * nrow; // (row=i, col=j)
			const int x = A[idx];
			const int n = D[idx];

			// fill column j: rows correspond to VAFs
			// M[k + j*K] addresses (k,j) in column-major for m (K x ncells)
			double *colptr = M + j * K;
			for (int k = 0; k < K; ++k) {
				colptr[k] = R::dbinom(static_cast<double>(x), static_cast<double>(n), P[k], log ? 1 : 0);
			}
		}

		// set dimnames to match R version
		List dn(2);
		dn[0] = vaf_names; // rownames
		dn[1] = cell_names; // colnames
		m.attr("dimnames") = dn;

		out[i] = m;
	}

	if (variants.size() == nvars) out.attr("names") = variants;
	return out;
}