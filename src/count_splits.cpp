#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

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
struct KeyHSHasher {
	size_t operator()(const KeyHS& k) const noexcept {
		uint64_t x = k.h ^ (uint64_t)k.s * 0x9e3779b97f4a7c15ULL;
		// one more avalanche
		x ^= (x >> 33);
		x *= 0xff51afd7ed558ccdULL;
		x ^= (x >> 33);
		return (size_t)x;
	}
};
static inline bool keyhs_lt(const KeyHS& a, const KeyHS& b) {
	if (a.s != b.s) return a.s < b.s;
	return a.h < b.h;
}

// [[Rcpp::export]]
std::vector< std::vector<int> > bipartition2(IntegerMatrix orig, int nTips) {
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

// build complement (1..nTips) \ v; assumes v is sorted ascending
static inline std::vector<int> complement_vec(const std::vector<int>& v, int nTips) {
	std::vector<int> comp;
	comp.reserve(nTips - (int)v.size());
	int j = 0;
	for (int t = 1; t <= nTips; ++t) {
		if (j < (int)v.size() && v[j] == t) {
			++j;
		} else {
			comp.push_back(t);
		}
	}
	return comp;
}

// [[Rcpp::export]]
List prop_part2(SEXP trees, int nTips){
    List tr(trees);
    int nbtree = tr.size();//, KeepPartition=1; // unused (EP 2020-05-02)
    List M = tr(0);
    IntegerMatrix E = M["edge"];
    std::vector< std::vector<int> > ans = bipartition2(E, nTips);
    std::vector<int> no;
    for(unsigned int i=0; i<ans.size();++i) no.push_back(1);
    no[0] = nbtree;
    for(int k=1; k<nbtree; ++k){
        List tmpTree = tr(k);
        IntegerMatrix tmpE = tmpTree["edge"];
        std::vector< std::vector<int> > bp = bipartition2(tmpE, nTips);
        for (unsigned int i = 1; i < bp.size(); i++) {
            unsigned int j = 1;
            next_j:
                if (bp[i] == ans[j]) {
                    no[j]++;
                    continue;
                }
                j++;
                if (j < ans.size()) goto next_j;
                else  {   //if(KeepPartition)
                    ans.push_back(bp[i]);
                    no.push_back(1);
                }
        }
    }
    List output = wrap(ans);
    output.attr("number") = no;
    output.attr("class") = "prop.part";
    return output;
}

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

// [[Rcpp::export]]
Rcpp::NumericVector prop_part2_targeted_edges(SEXP edges,
                                              Rcpp::IntegerMatrix E_target,
                                              int nTips,
                                              bool rooted = true,
                                              bool normalize = true)
{
	Rcpp::List el(edges);
	const int nbtree = el.size();
	if (nbtree <= 0) Rcpp::stop("edge_list is empty");
	if (!rooted) Rcpp::warning("prop_part2_targeted_edges assumes rooted; proceeding as rooted.");

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

	// counts aligned with target order
	std::vector<double> counts((size_t)K, 0.0);
	if (K > 0) counts[0] = (double)nbtree; // ape parity

	// reusable buffers
	std::vector< std::vector<int> > ch;
	std::vector<uint64_t> H;
	std::vector<uint32_t> SZ;

	// 2) Scan each edge matrix directly
	for (int k = 0; k < nbtree; ++k) {
		Rcpp::IntegerMatrix E = el[k];

		// children and root
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
			// binary search in keys_sorted
			int lo = 0, hi = K;
			while (lo < hi) {
				int mid = (lo + hi) >> 1;
				if (keyhs_lt(keys_sorted[(size_t)mid], key)) lo = mid + 1;
				else hi = mid;
			}
			if (lo < K && keys_sorted[(size_t)lo].h == key.h && keys_sorted[(size_t)lo].s == key.s) {
				const int tgt = pos[(size_t)lo];
				counts[(size_t)tgt] += 1.0;
			}
		}
	}

	// 3) return
	Rcpp::NumericVector out(K);
	if (normalize) {
		for (int i = 0; i < K; ++i) out[i] = counts[(size_t)i] / (double)nbtree;
	} else {
		for (int i = 0; i < K; ++i) out[i] = counts[(size_t)i];
	}
	return out;
}