#include <Rcpp.h>
#include <RcppParallel.h>
#include <tbb/spin_mutex.h>

using namespace Rcpp;
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

struct KeyIdx {
    KeyHS key;
    int idx;
};

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
    for (int v = nTips + 1; v <= N; ++v) {
        if (!isChild[(size_t)v]) { root = v; break; }
    }
    root_out = root;
}

static inline int infer_nTips_binary(const Rcpp::IntegerMatrix& E) {
    const int m = E.nrow();
    if ((m & 1) != 0) Rcpp::stop("infer_nTips_binary: edge count is odd (%d). Not a fully binary rooted tree.", m);
    // In a rooted full binary tree: edges m = 2 * Nnode, tips = Nnode + 1 = m/2 + 1
    return (m >> 1) + 1;
}

// Canonicalize a rooted clade hash/size into an unrooted split key (SHORTwise).
// Returns false for trivial splits (<= 1 on short side).
static inline bool make_unrooted_key(uint64_t h,
                                    uint32_t s,
                                    uint64_t allH,
                                    uint32_t nTips,
                                    KeyHS& out)
{
    const uint32_t comp = nTips - s;
    uint32_t s_short = s < comp ? s : comp;
    if (s_short <= 1U) return false;

    const uint64_t h_comp = allH ^ h;
    uint64_t h_short;
    if (s < comp) h_short = h;
    else if (s > comp) h_short = h_comp;
    else h_short = (h < h_comp) ? h : h_comp;

    out = KeyHS{h_short, s_short};
    return true;
}

struct UnrootedEdgesWorker : public RcppParallel::Worker {
    // inputs
    const std::vector<Rcpp::IntegerMatrix>& mats;
    const std::vector<KeyHS>& keys_unique;     // sorted unique keys
    const std::vector<int>& group_start;       // per-key start in targets_flat
    const std::vector<int>& group_len;         // per-key length
    const std::vector<int>& targets_flat;      // concatenated target indices
    const int U;
    const int K;
    const int nTips;
    const std::vector<uint64_t>& tipH;

    // shared output
    RcppParallel::RVector<double> counts;
    tbb::spin_mutex* mtx;

    UnrootedEdgesWorker(const std::vector<Rcpp::IntegerMatrix>& mats_,
                        const std::vector<KeyHS>& keys_unique_,
                        const std::vector<int>& group_start_,
                        const std::vector<int>& group_len_,
                        const std::vector<int>& targets_flat_,
                        int U_,
                        int K_,
                        int nTips_,
                        const std::vector<uint64_t>& tipH_,
                        Rcpp::NumericVector& counts_,
                        tbb::spin_mutex* mtx_)
        : mats(mats_),
          keys_unique(keys_unique_),
          group_start(group_start_),
          group_len(group_len_),
          targets_flat(targets_flat_),
          U(U_),
          K(K_),
          nTips(nTips_),
          tipH(tipH_),
          counts(counts_),
          mtx(mtx_) {}

    void operator()(std::size_t begin, std::size_t end) {
        std::vector<double> local((size_t)K, 0.0);

        std::vector< std::vector<int> > ch;
        std::vector<uint64_t> H;
        std::vector<uint32_t> SZ;

        for (std::size_t k = begin; k < end; ++k) {
            const Rcpp::IntegerMatrix& E = mats[k];

            int root = 0;
            build_children_and_root(E, nTips, ch, root);
            const int N = (int)ch.size() - 1;

            H.assign((size_t)N + 1, 0);
            SZ.assign((size_t)N + 1, 0);
            dfs_hash_size(root, ch, tipH, nTips, H, SZ);

            const uint64_t allH = H[(size_t)root];

            // scan edges and count unrooted splits
            const int nEdge = E.nrow();
            for (int e = 0; e < nEdge; ++e) {
                const int parent = E(e, 0);
                const int child = E(e, 1);
                if (child <= nTips) continue;

                const uint32_t s = SZ[(size_t)child];
                const uint64_t h = H[(size_t)child];

                // root-edge split appears twice (once per root child) when both children are internal.
                // Count it once by only taking the canonical side.
                if (parent == root) {
                    const uint32_t comp = (uint32_t)nTips - s;
                    const uint64_t h_comp = allH ^ h;
                    if (s > comp) continue;
                    if (s == comp && h > h_comp) continue;
                }

                KeyHS key;
                if (!make_unrooted_key(h, s, allH, (uint32_t)nTips, key)) continue;

                // binary search on keys_unique
                int lo = 0, hi = U;
                while (lo < hi) {
                    int mid = (lo + hi) >> 1;
                    const KeyHS& km = keys_unique[(size_t)mid];
                    if (keyhs_lt(km, key)) lo = mid + 1;
                    else hi = mid;
                }
                if (lo < U) {
                    const KeyHS& km = keys_unique[(size_t)lo];
                    if (km.h == key.h && km.s == key.s) {
                        const int start = group_start[(size_t)lo];
                        const int len = group_len[(size_t)lo];
                        for (int j = 0; j < len; ++j) {
                            const int tgt = targets_flat[(size_t)(start + j)];
                            local[(size_t)tgt] += 1.0;
                        }
                    }
                }
            }
        }

        tbb::spin_mutex::scoped_lock lock(*mtx);
        for (int i = 0; i < K; ++i) counts[i] += local[(size_t)i];
    }
};

// [[Rcpp::export]]
Rcpp::NumericVector prop_clades_unrooted_par(Rcpp::IntegerMatrix E_target,
                                             SEXP edges,
                                             bool normalize = true)
{
    // E_target is expected in postorder
    Rcpp::List el(edges);
    const int nbtree = el.size();
    if (nbtree <= 0) Rcpp::stop("edge_list is empty");

    const int nTips = infer_nTips_binary(E_target);
    std::vector<uint64_t> tipH = make_tip_hashes(nTips);

    // Build children/root for target and compute per-node hashes/sizes.
    std::vector< std::vector<int> > ch_t;
    int root_t = 0;
    build_children_and_root(E_target, nTips, ch_t, root_t);
    const int N_t = (int)ch_t.size() - 1;
    const int K = N_t - nTips; // Nnode
    if (K <= 0) Rcpp::stop("target tree has no internal nodes");

    std::vector<uint64_t> H_t((size_t)N_t + 1, 0);
    std::vector<uint32_t> SZ_t((size_t)N_t + 1, 0);
    dfs_hash_size(root_t, ch_t, tipH, nTips, H_t, SZ_t);
    const uint64_t allH_t = H_t[(size_t)root_t];

    // Collect canonical unrooted keys for each internal node index (excluding root itself).
    std::vector<KeyIdx> entries;
    entries.reserve((size_t)K);

    for (int node = nTips + 2; node <= N_t; ++node) {
        const uint32_t s = SZ_t[(size_t)node];
        const uint64_t h = H_t[(size_t)node];
        KeyHS key;
        if (!make_unrooted_key(h, s, allH_t, (uint32_t)nTips, key)) continue;
        const int idx = node - nTips - 1; // 0-based internal index (node nTips+1 -> 0)
        entries.push_back(KeyIdx{key, idx});
    }

    // Sort and group duplicates (e.g. the two root children for the root-edge split).
    std::sort(entries.begin(), entries.end(), [](const KeyIdx& a, const KeyIdx& b) {
        return keyhs_lt(a.key, b.key);
    });

    std::vector<KeyHS> keys_unique;
    std::vector<int> group_start;
    std::vector<int> group_len;
    std::vector<int> targets_flat;

    keys_unique.reserve(entries.size());
    group_start.reserve(entries.size());
    group_len.reserve(entries.size());

    int i = 0;
    while (i < (int)entries.size()) {
        const KeyHS key = entries[(size_t)i].key;
        keys_unique.push_back(key);
        group_start.push_back((int)targets_flat.size());

        int j = i;
        while (j < (int)entries.size() && entries[(size_t)j].key == key) {
            targets_flat.push_back(entries[(size_t)j].idx);
            ++j;
        }
        group_len.push_back(j - i);
        i = j;
    }

    const int U = (int)keys_unique.size();

    // Materialize IntegerMatrix views once to avoid R access in threads
    std::vector<Rcpp::IntegerMatrix> mats((size_t)nbtree);
    for (int k = 0; k < nbtree; ++k) {
        mats[(size_t)k] = Rcpp::IntegerMatrix(el[k]);
    }

    // parallel scan
    Rcpp::NumericVector counts(K); // zero-initialized
    tbb::spin_mutex mtx;
    UnrootedEdgesWorker worker(mats, keys_unique, group_start, group_len, targets_flat, U, K, nTips, tipH, counts, &mtx);
    RcppParallel::parallelFor(0, (size_t)nbtree, worker);

    // Root clade (all tips) always has support 1.
    counts[0] = counts[0] + (double)nbtree;

    if (normalize) {
        for (int k = 0; k < K; ++k) counts[k] /= (double)nbtree;
    }

    return counts;
}
