// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// Coerce a numeric vector to character the same way R's as.character() would
static CharacterVector as_character_like_R(const NumericVector &x) {
	SEXP s = Rf_coerceVector(x, STRSXP);
	return CharacterVector(s);
}

// Safely fetch row/col names from a matrix
static CharacterVector get_rownames(const RObject &mat) {
	List dn = mat.attr("dimnames");
	if (dn.isNULL() || dn.size() < 1 || dn[0] == R_NilValue) return CharacterVector();
	return as<CharacterVector>(dn[0]);
}
static CharacterVector get_colnames(const RObject &mat) {
	List dn = mat.attr("dimnames");
	if (dn.isNULL() || dn.size() < 2 || dn[1] == R_NilValue) return CharacterVector();
	return as<CharacterVector>(dn[1]);
}

// [[Rcpp::export]]
List get_leaf_liks_mat_cpp(const IntegerMatrix &amat,
						   const IntegerMatrix &dmat,
						   const NumericVector &vafs,
						   double eps = 0.0,
						   int ncores = 1,			// kept for signature parity; not used
						   bool log_ = false) {

	// Basic checks to match R function expectations
	if (amat.nrow() != dmat.nrow() || amat.ncol() != dmat.ncol()) {
		stop("amat and dmat must have the same dimensions");
	}

	const int K = vafs.size();
	const int nvars = amat.nrow();
	const int ncells = amat.ncol();

	// p <- pmin(vafs + eps, 1 - eps)
	NumericVector p(K);
	for (int k = 0; k < K; ++k) {
		double pk = vafs[k] + eps;
		double cap = 1.0 - eps;
		if (pk > cap) pk = cap;
		p[k] = pk;
	}

	// Names (must mirror R output)
	CharacterVector variants = get_rownames(amat);
	CharacterVector cell_names = get_colnames(amat);
	CharacterVector vaf_names = as_character_like_R(vafs);

	List out(nvars);

	for (int i = 0; i < nvars; ++i) {
		// For variant i, compute a K x ncells matrix
		NumericMatrix m(K, ncells);

		for (int j = 0; j < ncells; ++j) {
			int x = amat(i, j);
			int n = dmat(i, j);
			// dbinom for each VAF bin (row = VAF, col = cell)
			for (int k = 0; k < K; ++k) {
				// R::dbinom(x, n, p, give_log) matches R's dbinom edge-case behavior
				m(k, j) = R::dbinom((double)x, (double)n, p[k], log_ ? 1 : 0);
			}
		}

		// Set dimnames: rows = vafs, cols = cell names
		List dn(2);
		dn[0] = vaf_names;					// rownames
		dn[1] = cell_names;					// colnames
		m.attr("dimnames") = dn;

		out[i] = m;
	}

	// Set list names = variants
	if (variants.size() == nvars) {
		out.attr("names") = variants;
	}

	return out;
}