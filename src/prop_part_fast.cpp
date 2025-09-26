#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cstdint>
using namespace RcppParallel;

// --- helpers for stable keys from integer partitions ---
static inline std::string vec_to_key(const std::vector<int>& v){
	std::string s;
	// reserve rough capacity (3 chars per int on average)
	s.reserve(v.size() * 3);
	for(std::size_t i = 0; i < v.size(); ++i){
		if(i) s.push_back(',');
		s += std::to_string(v[i]);
	}
	return s;
}

static inline std::vector<int> key_to_vec(const std::string& key){
	std::vector<int> v;
	v.reserve(16);
	int num = 0;
	bool inNum = false;
	for(char c : key){
		if(c == ','){
			if(inNum){
				v.push_back(num);
				num = 0;
				inNum = false;
			}
		} else {
			inNum = true;
			num = num * 10 + (c - '0');
		}
	}
	if(inNum) v.push_back(num);
	return v;
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

// --- fast hashing & equality for integer partitions (avoid string keys) ---
struct VecHash {
	std::size_t operator()(const std::vector<int>& v) const noexcept {
		// FNV-1a style mix with size-seeding; fast & stable
		std::size_t h = v.size();
		for(int x : v){
			std::size_t k = static_cast<std::size_t>(x);
			h ^= k + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
		}
		return h;
	}
};
struct VecEq {
	bool operator()(const std::vector<int>& a, const std::vector<int>& b) const noexcept {
		return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
	}
};

// ---- Parallel reducer: aggregate counts and first-appearance order ----
struct PartReduce : public Worker {
	List trees;
	int nTips;
	// map: partition -> {count, firstpos}
	// firstpos encodes (tree_index << 32) | within_tree_index
	std::unordered_map< std::vector<int>, std::pair<int, uint64_t>, VecHash, VecEq > map;

	PartReduce(List trees_, int nTips_) : trees(trees_), nTips(nTips_) {}
	PartReduce(PartReduce& other, Split) : trees(other.trees), nTips(other.nTips) {}

	void operator()(std::size_t begin, std::size_t end){
		for(std::size_t i = begin; i < end; ++i){
			List M = trees[i];
			IntegerMatrix E = M["edge"];
			auto bp = bipartition2(E, nTips);
			for(std::size_t j = 1; j < bp.size(); ++j){ // skip root
				uint64_t first = (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
				auto it = map.find(bp[j]);
				if(it == map.end()){
					map.emplace(bp[j], std::make_pair(1, first));
				} else {
					it->second.first += 1;
					if(first < it->second.second) it->second.second = first;
				}
			}
		}
	}

	void join(const PartReduce& rhs){
		for(const auto& kv : rhs.map){
			auto it = map.find(kv.first);
			if(it == map.end()){
				map.emplace(kv.first, kv.second);
			} else {
				it->second.first += kv.second.first;
				if(kv.second.second < it->second.second) it->second.second = kv.second.second;
			}
		}
	}
};

// ---- Parallel collector: compute per-tree (non-root) partitions as integer vectors ----
struct BPCollector : public Worker {
	List trees;
	int nTips;
	std::vector< std::vector< std::vector<int> > >& perTreeParts; // per tree -> list of partitions (skip root)

	BPCollector(List trees_, int nTips_, std::vector< std::vector< std::vector<int> > >& perTreeParts_)
		: trees(trees_), nTips(nTips_), perTreeParts(perTreeParts_) {}

	void operator()(std::size_t begin, std::size_t end){
		for(std::size_t i = begin; i < end; ++i){
			List M = trees[i];
			IntegerMatrix E = M["edge"];
			std::vector< std::vector<int> > bp = bipartition2(E, nTips);
			std::vector< std::vector<int> > parts;
			if(bp.size() > 1){
				parts.reserve(bp.size() - 1);
				for(std::size_t j = 1; j < bp.size(); ++j){ // skip root at index 0
					parts.push_back(std::move(bp[j]));
				}
			}
			perTreeParts[i].swap(parts);
		}
	}
};

// [[Rcpp::export]]
List prop_part2_parallel(SEXP trees, int nTips){
	List tr(trees);
	const int nbtree = tr.size();
	if(nbtree == 0){
		List out = List::create();
		out.attr("number") = IntegerVector(0);
		out.attr("class") = "prop.part";
		return out;
	}

	// Compute first tree bipartitions to lock in canonical ordering prefix
	List M0 = tr(0);
	IntegerMatrix E0 = M0["edge"];
	auto bp0 = bipartition2(E0, nTips); // includes root at index 0

	// Parallel reduce over all trees to aggregate counts and first appearance
	PartReduce reducer(tr, nTips);
	parallelReduce(0, static_cast<std::size_t>(nbtree), reducer);

	// Seed answer with first tree's order (including root at index 0)
	std::vector< std::vector<int> > ans = bp0;
	std::vector<int> no(ans.size(), 1);
	if(!no.empty()) no[0] = nbtree; // root appears in every tree

	// Quick lookup for partitions already in ans (skip root)
	std::unordered_map< std::vector<int>, int, VecHash, VecEq > indexMap;
	indexMap.reserve(ans.size() * 2);
	for(std::size_t j = 1; j < ans.size(); ++j){
		auto it = reducer.map.find(ans[j]);
		no[j] = (it == reducer.map.end()) ? 1 : it->second.first; // if present beyond tree0, use aggregated count
		indexMap.emplace(ans[j], static_cast<int>(j));
	}

	// Collect remaining unique partitions not in the first tree, along with their first-appearance order
	struct Pending { uint64_t pos; std::vector<int> part; int cnt; };
	std::vector<Pending> pending;
	pending.reserve(reducer.map.size());
	for(const auto& kv : reducer.map){
		if(indexMap.find(kv.first) == indexMap.end()){
			pending.push_back(Pending{kv.second.second, kv.first, kv.second.first});
		}
	}
	std::sort(pending.begin(), pending.end(), [](const Pending& a, const Pending& b){ return a.pos < b.pos; });

	// Append in strict first-appearance order (tree_index, within_index) to match sequential behavior
	for(const auto& p : pending){
		ans.push_back(p.part);
		no.push_back(p.cnt);
	}

	List output = wrap(ans);
	output.attr("number") = no;
	output.attr("class") = "prop.part";
	return output;
}