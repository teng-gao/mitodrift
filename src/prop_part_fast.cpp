#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
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

// ---- Parallel wrapper for prop_part2 using RcppParallel ----
struct PartitionCollector : public Worker {
	List trees;
	int nTips;
	std::vector< std::vector<std::string> >& perTreeKeys;

	PartitionCollector(List trees_, int nTips_, std::vector< std::vector<std::string> >& perTreeKeys_)
		: trees(trees_), nTips(nTips_), perTreeKeys(perTreeKeys_) {}

	void operator()(std::size_t begin, std::size_t end){
		for(std::size_t i = begin; i < end; ++i){
			List M = trees[i];
			IntegerMatrix E = M["edge"];
			std::vector< std::vector<int> > bp = bipartition2(E, nTips);

			std::vector<std::string> keys;
			if(bp.size() > 1){
				keys.reserve(bp.size() - 1);
				for(std::size_t j = 1; j < bp.size(); ++j){ // skip root (index 0)
					keys.push_back(vec_to_key(bp[j]));
				}
			}
			perTreeKeys[i].swap(keys);
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

	// 1) Collect per-tree (non-root) partitions in parallel as stable string keys
	std::vector< std::vector<std::string> > perTreeKeys(nbtree);
	PartitionCollector worker(tr, nTips, perTreeKeys);
	parallelFor(0, static_cast<std::size_t>(nbtree), worker);

	// 2) Merge counts across trees
	std::unordered_map<std::string,int> countMap;
	countMap.reserve(static_cast<std::size_t>(nbtree) * 8);
	for(int i = 0; i < nbtree; ++i){
		for(const std::string& key : perTreeKeys[i]){
			++countMap[key];
		}
	}

	// 3) Start from the first tree to mimic original ordering for existing splits
	List M0 = tr(0);
	IntegerMatrix E0 = M0["edge"];
	std::vector< std::vector<int> > bp0 = bipartition2(E0, nTips);

	std::vector< std::vector<int> > ans = bp0;           // includes root at index 0
	std::vector<int> no(ans.size(), 0);

	// root bipartition (index 0) appears in every tree by definition
	if(!no.empty()) no[0] = nbtree;

	// existing (first-tree) splits: fill counts from the map, skipping root
	std::unordered_set<std::string> seen;
	seen.reserve(ans.size());
	for(std::size_t j = 1; j < ans.size(); ++j){
		std::string key = vec_to_key(ans[j]);
		seen.insert(key);
		auto it = countMap.find(key);
		no[j] = (it == countMap.end()) ? 0 : it->second;
	}

	// 4) Append any new splits observed in other trees, preserving original ordering:
	// iterate trees in increasing index (1..nbtree-1) and, within each tree,
	// iterate partitions in their original order; append only when first seen.
	for(int k = 1; k < nbtree; ++k){
		const std::vector<std::string>& keys = perTreeKeys[k];
		for(const std::string& key : keys){
			if(seen.find(key) == seen.end()){
				ans.push_back(key_to_vec(key));
				auto it = countMap.find(key);
				no.push_back(it == countMap.end() ? 0 : it->second);
				seen.insert(key);
			}
		}
	}

	List output = wrap(ans);
	output.attr("number") = no;
	output.attr("class") = "prop.part";
	return output;
}