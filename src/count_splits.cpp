#include <Rcpp.h>
#include <unordered_map>
#include <RcppParallel.h>
using namespace Rcpp;
#include <tbb/spin_mutex.h>
using namespace RcppParallel;

// --- fast hashing utils (deterministic, no RNG state) ----------------------
static inline uint64_t splitmix64(uint64_t x) {
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

struct KeyHS {
	uint64_t h;
	uint32_t s;
	bool operator==(const KeyHS& o) const noexcept { return h == o.h && s == o.s; }
};

static inline bool keyhs_lt(const KeyHS& a, const KeyHS& b) {
	if (a.s != b.s) return a.s < b.s;
	return a.h < b.h;
}

std::vector< std::vector<int> > bipartition2(Rcpp::IntegerMatrix orig, int nTips) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent), j=0;
    int nnode = m - nTips;
    // create list for results
    std::vector< std::vector<int> > out(nnode);
    std::vector<int> y;
    for(int i = 0; i<parent.size(); i++){
        j = parent[i] - nTips - 1L;
        if(children[i] > nTips){
            y = out[children[i] - nTips -1L];
            out[j].insert( out[j].end(), y.begin(), y.end() );
        }
        else out[j].push_back(children[i]);
    }
    for(int i=0; i<nnode; ++i){
        sort(out[i].begin(), out[i].end());
    }
    return out;
}

// ---- helpers for targeted counting ---------------------------------------

// compute per-tip 64-bit hashes 1..nTips
static inline std::vector<uint64_t> make_tip_hashes(int nTips) {
	std::vector<uint64_t> th((size_t)nTips + 1);
	for (int t = 1; t <= nTips; ++t) th[(size_t)t] = splitmix64((uint64_t)t);
	return th;
}

// postorder compute of (hash,size) for every node given children adjacency
static void dfs_hash_size(int u,
                          const std::vector< std::vector<int> >& ch,
                          const std::vector<uint64_t>& tipH,
                          int nTips,
                          std::vector<uint64_t>& H,
                          std::vector<uint32_t>& SZ)
{
	if (u <= nTips) { H[(size_t)u] = tipH[(size_t)u]; SZ[(size_t)u] = 1U; return; }
	uint64_t hx = 0;
	uint32_t sz = 0;
	for (int v : ch[(size_t)u]) {
		if (H[(size_t)v] == 0 && SZ[(size_t)v] == 0) dfs_hash_size(v, ch, tipH, nTips, H, SZ);
		hx ^= H[(size_t)v];
		sz += SZ[(size_t)v];
	}
	H[(size_t)u] = hx;
	SZ[(size_t)u] = sz;
}

// build children adjacency and root for a rooted tree
static inline void build_children_and_root(const Rcpp::IntegerMatrix& E,
                                           int nTips,
                                           std::vector< std::vector<int> >& ch,
                                           int& root_out)
{
	const int nEdge = E.nrow();
	int maxIdx = 0;
	for (int e = 0; e < nEdge; ++e) {
		maxIdx = std::max(maxIdx, E(e,0));
		maxIdx = std::max(maxIdx, E(e,1));
	}
	const int N = maxIdx;
	ch.assign((size_t)N + 1, {});
	std::vector<char> isChild((size_t)N + 1, 0);
	for (int e = 0; e < nEdge; ++e) {
		const int p = E(e,0), v = E(e,1);
		ch[(size_t)p].push_back(v);
		isChild[(size_t)v] = 1;
	}
	int root = nTips + 1;
	for (int v = nTips + 1; v <= N; ++v) if (!isChild[(size_t)v]) { root = v; break; }
	root_out = root;
}

struct TargetEdgesWorker : public RcppParallel::Worker {
	// inputs
	const std::vector< Rcpp::IntegerMatrix >& mats;   // edge matrices for all trees
	const std::vector<KeyHS>& keys_sorted;            // sorted (hash,size) keys
	const std::vector<int>&   pos;                    // map from sorted index -> target index
	const int K;
	const int nTips;
	const std::vector<uint64_t>& tipH;

	// shared output
	RcppParallel::RVector<double> counts;
	tbb::spin_mutex* mtx;

	TargetEdgesWorker(const std::vector<Rcpp::IntegerMatrix>& mats_,
	                  const std::vector<KeyHS>& keys_sorted_,
	                  const std::vector<int>& pos_,
	                  int K_, int nTips_,
	                  const std::vector<uint64_t>& tipH_,
	                  Rcpp::NumericVector& counts_,
	                  tbb::spin_mutex* mtx_)
		: mats(mats_), keys_sorted(keys_sorted_), pos(pos_), K(K_), nTips(nTips_), tipH(tipH_),
		  counts(counts_), mtx(mtx_) {}

	void operator()(std::size_t begin, std::size_t end) {
		std::vector<double> local((size_t)K, 0.0);

		// thread-local reusable buffers
		std::vector< std::vector<int> > ch;
		std::vector<uint64_t> H;
		std::vector<uint32_t> SZ;

		for (std::size_t k = begin; k < end; ++k) {
			const Rcpp::IntegerMatrix& E = mats[k];

			// build children and root
			int root = 0;
			build_children_and_root(E, nTips, ch, root);
			const int N = (int)ch.size() - 1;

			// per-node (hash,size)
			H.assign((size_t)N + 1, 0);
			SZ.assign((size_t)N + 1, 0);
			dfs_hash_size(root, ch, tipH, nTips, H, SZ);

			// iterate internal child edges
			const int nEdge = E.nrow();
			for (int e = 0; e < nEdge; ++e) {
				const int child = E(e,1);
				if (child <= nTips) continue;
				KeyHS key{ H[(size_t)child], SZ[(size_t)child] };

				// binary search on keys_sorted
				int lo = 0, hi = K;
				while (lo < hi) {
					int mid = (lo + hi) >> 1;
					const KeyHS& km = keys_sorted[(size_t)mid];
					if (keyhs_lt(km, key)) lo = mid + 1;
					else hi = mid;
				}
				if (lo < K) {
					const KeyHS& km = keys_sorted[(size_t)lo];
					if (km.h == key.h && km.s == key.s) {
						const int tgt = pos[(size_t)lo];
						local[(size_t)tgt] += 1.0;
					}
				}
			}
		}

		// reduce into shared counts
		tbb::spin_mutex::scoped_lock lock(*mtx);
		for (int i = 0; i < K; ++i) counts[i] += local[(size_t)i];
	}
};

static inline int infer_nTips_binary(const Rcpp::IntegerMatrix& E) {
	const int m = E.nrow();
	if ((m & 1) != 0) Rcpp::stop("infer_nTips_binary: edge count is odd (%d). Not a fully binary rooted tree.", m);
	// In a rooted full binary tree: edges m = 2 * Nnode, tips = Nnode + 1 = m/2 + 1
	return (m >> 1) + 1;
}

// [[Rcpp::export]]
Rcpp::NumericVector prop_clades_par(Rcpp::IntegerMatrix E_target,
    SEXP edges,
    bool rooted = true,
    bool normalize = true)
{
	// E_target is expected in postorder
	Rcpp::List el(edges);
	const int nbtree = el.size();
	if (nbtree <= 0) Rcpp::stop("edge_list is empty");
	if (!rooted) Rcpp::warning("prop_clades_par assumes rooted; proceeding as rooted.");

	// infer nTips from target (assumes fully binary rooted tree)
	const int nTips = infer_nTips_binary(E_target);

	// tip hashes
	std::vector<uint64_t> tipH = make_tip_hashes(nTips);

	// 1) Build target keys (hash,size) in target order
	std::vector< std::vector<int> > targetSplits = bipartition2(E_target, nTips);
	const int K = (int)targetSplits.size();
	if (K <= 0) Rcpp::stop("target tree has no internal splits");

	std::vector<KeyHS> keys((size_t)K);
	for (int i = 0; i < K; ++i) {
		uint64_t h = 0;
		for (int t : targetSplits[(size_t)i]) h ^= tipH[(size_t)t];
		keys[(size_t)i] = KeyHS{h, (uint32_t)targetSplits[(size_t)i].size()};
	}
	// sorted view for binary search
	std::vector<int> ord((size_t)K);
	for (int i = 0; i < K; ++i) ord[(size_t)i] = i;
	std::sort(ord.begin(), ord.end(), [&](int a, int b){
		return keyhs_lt(keys[(size_t)a], keys[(size_t)b]);
	});
	std::vector<KeyHS> keys_sorted((size_t)K);
	std::vector<int>   pos((size_t)K);
	for (int r = 0; r < K; ++r) {
		const int i = ord[(size_t)r];
		keys_sorted[(size_t)r] = keys[(size_t)i];
		pos[(size_t)r] = i;
	}

	// 2) Materialize IntegerMatrix views once to avoid R access in threads
	std::vector<Rcpp::IntegerMatrix> mats((size_t)nbtree);
	for (int k = 0; k < nbtree; ++k) {
		mats[(size_t)k] = Rcpp::IntegerMatrix(el[k]);
	}

	// 3) parallel scan
	Rcpp::NumericVector counts(K); // zero-initialized
	tbb::spin_mutex mtx;
	TargetEdgesWorker worker(mats, keys_sorted, pos, K, nTips, tipH, counts, &mtx);
	RcppParallel::parallelFor(0, (size_t)nbtree, worker);

	// ape parity for the first split
	if (K > 0) counts[0] = counts[0] + (double)nbtree;

	// 4) return
	if (normalize) {
		for (int i = 0; i < K; ++i) counts[i] /= (double)nbtree;
	}
	return counts;
}


struct BatchEdgesWorker : public RcppParallel::Worker {
	// inputs
	const std::vector< Rcpp::IntegerMatrix >& mats;   // edge matrices for all trees (in MCMC order)
	const std::vector<KeyHS>& keys_sorted;            // sorted (hash,size) keys
	const std::vector<int>&   pos;                    // map from sorted index -> target index
	const int K;
	const int nTips;
	const std::vector<uint64_t>& tipH;
	const int batch_size;
	const int n_use;

	// output: per-batch sums (length nbatch*K), written without locks since each batch index is exclusive
	RcppParallel::RMatrix<double> batch_sums;

	BatchEdgesWorker(const std::vector<Rcpp::IntegerMatrix>& mats_,
				  const std::vector<KeyHS>& keys_sorted_,
				  const std::vector<int>& pos_,
				  int K_, int nTips_,
				  const std::vector<uint64_t>& tipH_,
				  int batch_size_,
				  int n_use_,
				  Rcpp::NumericMatrix& batch_sums_)
		: mats(mats_), keys_sorted(keys_sorted_), pos(pos_), K(K_), nTips(nTips_), tipH(tipH_),
		  batch_size(batch_size_), n_use(n_use_), batch_sums(batch_sums_) {}

	void operator()(std::size_t begin, std::size_t end) {
		// thread-local reusable buffers
		std::vector< std::vector<int> > ch;
		std::vector<uint64_t> H;
		std::vector<uint32_t> SZ;

		for (std::size_t b = begin; b < end; ++b) {
			std::vector<double> local((size_t)K, 0.0);

			const int k0 = (int)b * batch_size;
			const int k1 = std::min(k0 + batch_size, n_use);
			for (int k = k0; k < k1; ++k) {
				const Rcpp::IntegerMatrix& E = mats[(size_t)k];

				// build children and root
				int root = 0;
				build_children_and_root(E, nTips, ch, root);
				const int N = (int)ch.size() - 1;

				// per-node (hash,size)
				H.assign((size_t)N + 1, 0);
				SZ.assign((size_t)N + 1, 0);
				dfs_hash_size(root, ch, tipH, nTips, H, SZ);

				// iterate internal child edges
				const int nEdge = E.nrow();
				for (int e = 0; e < nEdge; ++e) {
					const int child = E(e,1);
					if (child <= nTips) continue;
					KeyHS key{ H[(size_t)child], SZ[(size_t)child] };

					// binary search on keys_sorted
					int lo = 0, hi = K;
					while (lo < hi) {
						int mid = (lo + hi) >> 1;
						const KeyHS& km = keys_sorted[(size_t)mid];
						if (keyhs_lt(km, key)) lo = mid + 1;
						else hi = mid;
					}
					if (lo < K) {
						const KeyHS& km = keys_sorted[(size_t)lo];
						if (km.h == key.h && km.s == key.s) {
							const int tgt = pos[(size_t)lo];
							local[(size_t)tgt] += 1.0;
						}
					}
				}
			}

			// ape parity for the first split within this batch
			if (K > 0) local[0] += (double)(k1 - k0);

			// write local into batch_sums(b, i)
			for (int i = 0; i < K; ++i) batch_sums((int)b, i) = local[(size_t)i];
		}
	}
};


// [[Rcpp::export]]
Rcpp::List prop_clades_par_bm_se(Rcpp::IntegerMatrix E_target,
	SEXP edges,
	int nbatch = 50,
	bool rooted = true)
{
	// E_target is expected in postorder
	Rcpp::List el(edges);
	const int nbtree = el.size();
	if (nbtree <= 0) Rcpp::stop("edge_list is empty");
	if (nbatch < 2) Rcpp::stop("nbatch must be >= 2 for batch means SE");
	if (!rooted) Rcpp::warning("prop_clades_par_bm_se assumes rooted; proceeding as rooted.");

	// infer nTips from target (assumes fully binary rooted tree)
	const int nTips = infer_nTips_binary(E_target);

	// tip hashes
	std::vector<uint64_t> tipH = make_tip_hashes(nTips);

	// 1) Build target keys (hash,size) in target order
	std::vector< std::vector<int> > targetSplits = bipartition2(E_target, nTips);
	const int K = (int)targetSplits.size();
	if (K <= 0) Rcpp::stop("target tree has no internal splits");

	std::vector<KeyHS> keys((size_t)K);
	for (int i = 0; i < K; ++i) {
		uint64_t h = 0;
		for (int t : targetSplits[(size_t)i]) h ^= tipH[(size_t)t];
		keys[(size_t)i] = KeyHS{h, (uint32_t)targetSplits[(size_t)i].size()};
	}
	// sorted view for binary search
	std::vector<int> ord((size_t)K);
	for (int i = 0; i < K; ++i) ord[(size_t)i] = i;
	std::sort(ord.begin(), ord.end(), [&](int a, int b){
		return keyhs_lt(keys[(size_t)a], keys[(size_t)b]);
	});
	std::vector<KeyHS> keys_sorted((size_t)K);
	std::vector<int>   pos((size_t)K);
	for (int r = 0; r < K; ++r) {
		const int i = ord[(size_t)r];
		keys_sorted[(size_t)r] = keys[(size_t)i];
		pos[(size_t)r] = i;
	}

	// 2) Materialize IntegerMatrix views once to avoid R access in threads
	std::vector<Rcpp::IntegerMatrix> mats((size_t)nbtree);
	for (int k = 0; k < nbtree; ++k) {
		mats[(size_t)k] = Rcpp::IntegerMatrix(el[k]);
	}

	// 3) Batch setup (discard remainder)
	const int batch_size = nbtree / nbatch;
	if (batch_size < 1) Rcpp::stop("Too few trees (%d) for nbatch=%d", nbtree, nbatch);
	const int n_use = batch_size * nbatch;

	// per-batch sums: rows=nbatch, cols=K
	Rcpp::NumericMatrix batch_sums(nbatch, K);
	BatchEdgesWorker worker(mats, keys_sorted, pos, K, nTips, tipH, batch_size, n_use, batch_sums);
	RcppParallel::parallelFor(0, (size_t)nbatch, worker);

	// 4) Compute p_hat and batch means
	Rcpp::NumericVector p_hat(K);
	Rcpp::NumericVector se(K);

	for (int i = 0; i < K; ++i) {
		// mean across batches of (batch_sum / batch_size)
		double mean_bm = 0.0;
		for (int b = 0; b < nbatch; ++b) mean_bm += batch_sums(b, i) / (double)batch_size;
		mean_bm /= (double)nbatch;
		p_hat[i] = mean_bm;

		// sample variance of batch means
		double ss = 0.0;
		for (int b = 0; b < nbatch; ++b) {
			double bm = batch_sums(b, i) / (double)batch_size;
			double d = bm - mean_bm;
			ss += d * d;
		}
		// var(mean) ≈ var(batch_means) / nbatch
		double var_bm = (nbatch > 1) ? (ss / (double)(nbatch - 1)) : 0.0;
		double var_mean = var_bm / (double)nbatch;
		se[i] = std::sqrt(std::max(0.0, var_mean));
	}

	return Rcpp::List::create(
		Rcpp::Named("prop") = p_hat,
		Rcpp::Named("se") = se,
		Rcpp::Named("nbatch") = nbatch,
		Rcpp::Named("batch_size") = batch_size,
		Rcpp::Named("n_use") = n_use,
		Rcpp::Named("n_total") = nbtree
	);
}