#include <RcppParallel.h>
#include <RcppArmadillo.h>

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

arma::Col<int> reorderRcpp(const arma::Col<int>& E);

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
        : Uin(static_cast<size_t>(n) * static_cast<size_t>(C)),
          up_to_parent(static_cast<size_t>(n) * static_cast<size_t>(C)),
          down_in(static_cast<size_t>(n) * static_cast<size_t>(C)),
          temp(C), u(C), S_parent(C), down_msg(C), a(C), b(C) {
        stack.reserve(static_cast<size_t>(n));
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
    const size_t nC = static_cast<size_t>(n) * static_cast<size_t>(C);

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
        double* const Uin_node = Uin_data + static_cast<size_t>(node) * C;
        double* const up_node = up_data + static_cast<size_t>(node) * C;
        double* const Uin_parent = Uin_data + static_cast<size_t>(par) * C;
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
        double* const down_parent = down_data + static_cast<size_t>(u_node) * C;
        double* const Uin_parent = Uin_data + static_cast<size_t>(u_node) * C;
        const auto& ch = children_of[u_node];
        for (int t = 0; t < 2; ++t) {
            const int v = ch[t];
            if (v < 0) continue;
            double* const down_child = down_data + static_cast<size_t>(v) * C;
            double* const up_child = up_data + static_cast<size_t>(v) * C;
            for (int r = 0; r < C; ++r) {
                const double P_ur = (u_node == root)
                    ? (r == 0 ? 0.0 : neg_inf)
                    : P[static_cast<size_t>(r) * n + u_node];
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
        const double* const down_v = down_data + static_cast<size_t>(v) * C;
        const double* const Uin_v = Uin_data + static_cast<size_t>(v) * C;
        for (int c = 0; c < C; ++c) {
            const double P_vc = (v == root)
                ? (c == 0 ? 0.0 : neg_inf)
                : P[static_cast<size_t>(c) * n + v];
            temp_data[c] = P_vc + Uin_v[c] + down_v[c];
            if (temp_data[c] > maxW) maxW = temp_data[c];
        }
        double denom = 0.0;
        double* const node_col = node_out + static_cast<size_t>(v);
        for (int c = 0; c < C; ++c) {
            const double val = std::exp(temp_data[c] - maxW);
            node_col[static_cast<size_t>(c) * n] = val;
            denom += val;
        }
        const double inv = 1.0 / denom;
        for (int c = 0; c < C; ++c) node_col[static_cast<size_t>(c) * n] *= inv;
        if (v == root) logZ_root = maxW + std::log(denom);
    }

    const size_t stride_parent = static_cast<size_t>(m);
    const size_t stride_child = stride_parent * static_cast<size_t>(C);
    for (int i = 0; i < m; ++i) {
        const int u_node = parent_ptr[i];
        const int v = child_ptr[i];
        const double* const down_u = down_data + static_cast<size_t>(u_node) * C;
        const double* const Uin_u = Uin_data + static_cast<size_t>(u_node) * C;
        const double* const up_v = up_data + static_cast<size_t>(v) * C;
        const double* const down_v = down_data + static_cast<size_t>(v) * C;
        const double* const Uin_v = Uin_data + static_cast<size_t>(v) * C;
        for (int r = 0; r < C; ++r) {
            const double P_ur = (u_node == root)
                ? (r == 0 ? 0.0 : neg_inf)
                : P[static_cast<size_t>(r) * n + u_node];
            S_parent_data[r] = P_ur + down_u[r] + Uin_u[r] - up_v[r];
        }
        for (int k = 0; k < C; ++k) down_msg_data[k] = P[static_cast<size_t>(k) * n + v] + Uin_v[k];
        double maxWr = neg_inf;
        for (int r = 0; r < C; ++r) if (S_parent_data[r] > maxWr) maxWr = S_parent_data[r];
        for (int r = 0; r < C; ++r) a_data[r] = std::exp(S_parent_data[r] - maxWr);
        double maxY = neg_inf;
        for (int k = 0; k < C; ++k) if (down_msg_data[k] > maxY) maxY = down_msg_data[k];
        for (int k = 0; k < C; ++k) b_data[k] = std::exp(down_msg_data[k] - maxY);
        double denom = 0.0;
        for (int r = 0; r < C; ++r) {
            double s = 0.0;
            const double* const Arow = expA_ptr + static_cast<size_t>(r) * C;
            for (int k = 0; k < C; ++k) s += Arow[k] * b_data[k];
            denom += a_data[r] * s;
        }
        const double inv = 1.0 / denom;
        for (int r = 0; r < C; ++r) {
            const double* const Arow = expA_ptr + static_cast<size_t>(r) * C;
            const size_t base_idx_r = static_cast<size_t>(i) + static_cast<size_t>(r) * stride_parent;
            for (int k = 0; k < C; ++k) {
                edge_out[base_idx_r + static_cast<size_t>(k) * stride_child] = a_data[r] * Arow[k] * b_data[k] * inv;
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

static inline void compute_down_msg_fast(
    const double* S_parent,
    const double* expA_rowmajor,
    int C,
    double* out,
    double* su)
{
    double maxS = -std::numeric_limits<double>::infinity();
    for (int u = 0; u < C; ++u) if (S_parent[u] > maxS) maxS = S_parent[u];
    for (int u = 0; u < C; ++u) su[u] = std::exp(S_parent[u] - maxS);
    for (int k = 0; k < C; ++k) {
        double dot = 0.0;
        for (int u = 0; u < C; ++u) dot += expA_rowmajor[u * C + k] * su[u];
        out[k] = std::log(dot) + maxS;
    }
}

Rcpp::List compute_node_edge_beliefs_bp2(
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

    std::vector< std::vector<double> > node_storage(
        static_cast<size_t>(L),
        std::vector<double>(static_cast<size_t>(n) * static_cast<size_t>(C)));
    std::vector< std::vector<double> > edge_storage(
        static_cast<size_t>(L),
        std::vector<double>(static_cast<size_t>(m) * static_cast<size_t>(C) * static_cast<size_t>(C)));
    std::vector<double> logZ_vals(static_cast<size_t>(L), 0.0);

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

    parallelFor(0, static_cast<size_t>(L), worker);

    Rcpp::List node_list(L);
    Rcpp::List edge_list(L);
    Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(m, C, C);

    double logZ_total = 0.0;
    for (int l = 0; l < L; ++l) {
        Rcpp::NumericMatrix node_beliefs(n, C);
        std::copy(node_storage[static_cast<size_t>(l)].begin(),
                  node_storage[static_cast<size_t>(l)].end(),
                  node_beliefs.begin());

        Rcpp::NumericVector edge_beliefs(edge_storage[static_cast<size_t>(l)].size());
        std::copy(edge_storage[static_cast<size_t>(l)].begin(),
                  edge_storage[static_cast<size_t>(l)].end(),
                  edge_beliefs.begin());
        edge_beliefs.attr("dim") = dims;

        node_list[l] = node_beliefs;
        edge_list[l] = edge_beliefs;
        logZ_total += logZ_vals[static_cast<size_t>(l)];
    }

    return Rcpp::List::create(
        Rcpp::Named("node_beliefs") = node_list,
        Rcpp::Named("edge_beliefs") = edge_list,
        Rcpp::Named("logZ") = logZ_total
    );
}

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
        const size_t nC = static_cast<size_t>(n) * static_cast<size_t>(C);

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
        stack.reserve(static_cast<size_t>(n));

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
                double* const Uin_node = Uin.data() + static_cast<size_t>(node) * C;
                double* const up_node = up_to_parent.data() + static_cast<size_t>(node) * C;
                double* const Uin_parent = Uin.data() + static_cast<size_t>(par) * C;
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
                double* const down_parent = down_in.data() + static_cast<size_t>(u_node) * C;
                double* const Uin_parent = Uin.data() + static_cast<size_t>(u_node) * C;
                const auto& ch = children_of[u_node];
                for (int t = 0; t < 2; ++t) {
                    const int v = ch[t];
                    if (v < 0) continue;
                    double* const down_child = down_in.data() + static_cast<size_t>(v) * C;
                    double* const up_child = up_to_parent.data() + static_cast<size_t>(v) * C;
                    for (int r = 0; r < C; ++r) {
                        const double P_ur = (u_node == root)
                            ? (r == 0 ? 0.0 : neg_inf)
                            : P[static_cast<size_t>(r) * n + u_node];
                        S_parent[r] = P_ur + down_parent[r] + Uin_parent[r] - up_child[r];
                    }
                    compute_down_msg_fast(S_parent.data(), expA.data(), C, down_msg.data(), u.data());
                    for (int k = 0; k < C; ++k) down_child[k] = down_msg[k];
                    stack.push_back(v);
                }
            }

            for (int v = 0; v < n; ++v) {
                double maxW = neg_inf;
                const double* const down_v = down_in.data() + static_cast<size_t>(v) * C;
                const double* const Uin_v = Uin.data() + static_cast<size_t>(v) * C;
                for (int c = 0; c < C; ++c) {
                    const double P_vc = (v == root)
                        ? (c == 0 ? 0.0 : neg_inf)
                        : P[static_cast<size_t>(c) * n + v];
                    temp[c] = P_vc + Uin_v[c] + down_v[c];
                    if (temp[c] > maxW) maxW = temp[c];
                }
                double denom = 0.0;
                for (int c = 0; c < C; ++c) {
                    const double val = std::exp(temp[c] - maxW);
                    node_tmp[static_cast<size_t>(v) + static_cast<size_t>(c) * n] = val;
                    denom += val;
                }
                const double inv = 1.0 / denom;
                for (int c = 0; c < C; ++c) {
                    node_tmp[static_cast<size_t>(v) + static_cast<size_t>(c) * n] *= inv;
                }
                if (v == root) logZ_root = maxW + std::log(denom);
            }

            double* const leaf_out = leaf_storage[idx_l].data();
            for (size_t tip_idx = 0; tip_idx < tip_ids.size(); ++tip_idx) {
                const int node_id = tip_ids[tip_idx];
                double sum = 0.0;
                for (int c = 0; c < C; ++c) {
                    const double val = node_tmp[static_cast<size_t>(node_id) + static_cast<size_t>(c) * n];
                    leaf_out[tip_idx + static_cast<size_t>(c) * tip_ids.size()] = val;
                    sum += val;
                }
                sum = (sum > 0.0) ? sum : 1.0;
                for (int c = 0; c < C; ++c) {
                    leaf_out[tip_idx + static_cast<size_t>(c) * tip_ids.size()] /= sum;
                }
            }
            logZ_vals[idx_l] = logZ_root;

            double* const edge_out = edge_counts_local[idx_l].data();
            std::fill(edge_out, edge_out + static_cast<size_t>(C) * static_cast<size_t>(C), 0.0);

            for (int i = 0; i < m; ++i) {
                const int u_node = parent_ptr[i];
                const int v = child_ptr[i];
                const double* const down_u = down_in.data() + static_cast<size_t>(u_node) * C;
                const double* const Uin_u = Uin.data() + static_cast<size_t>(u_node) * C;
                const double* const up_v = up_to_parent.data() + static_cast<size_t>(v) * C;
                const double* const down_v = down_in.data() + static_cast<size_t>(v) * C;
                const double* const Uin_v = Uin.data() + static_cast<size_t>(v) * C;
                for (int r = 0; r < C; ++r) {
                    const double P_ur = (u_node == root)
                        ? (r == 0 ? 0.0 : neg_inf)
                        : P[static_cast<size_t>(r) * n + u_node];
                    S_parent[r] = P_ur + down_u[r] + Uin_u[r] - up_v[r];
                }
                for (int k = 0; k < C; ++k) down_msg[k] = P[static_cast<size_t>(k) * n + v] + Uin_v[k];
                double maxWr = neg_inf;
                for (int r = 0; r < C; ++r) if (S_parent[r] > maxWr) maxWr = S_parent[r];
                for (int r = 0; r < C; ++r) a[r] = std::exp(S_parent[r] - maxWr);
                double maxY = neg_inf;
                for (int k = 0; k < C; ++k) if (down_msg[k] > maxY) maxY = down_msg[k];
                for (int k = 0; k < C; ++k) b[k] = std::exp(down_msg[k] - maxY);
                double denom = 0.0;
                for (int r = 0; r < C; ++r) {
                    double s = 0.0;
                    const double* const Arow = expA.data() + static_cast<size_t>(r) * C;
                    for (int k = 0; k < C; ++k) s += Arow[k] * b[k];
                    denom += a[r] * s;
                }
                const double inv = 1.0 / denom;
                for (int r = 0; r < C; ++r) {
                    const double* const Arow = expA.data() + static_cast<size_t>(r) * C;
                    for (int k = 0; k < C; ++k) {
                        edge_out[static_cast<size_t>(r) * C + k] += a[r] * Arow[k] * b[k] * inv;
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
        static_cast<size_t>(L),
        std::vector<double>(static_cast<size_t>(nTips) * static_cast<size_t>(C)));
    std::vector< std::vector<double> > edge_counts_local(
        static_cast<size_t>(L),
        std::vector<double>(static_cast<size_t>(C) * static_cast<size_t>(C), 0.0));
    std::vector<double> logZ_vals(static_cast<size_t>(L), 0.0);

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

    parallelFor(0, static_cast<size_t>(L), worker);

    double logZ_total = std::accumulate(logZ_vals.begin(), logZ_vals.end(), 0.0);

    Rcpp::List leaf_list(L);
    for (int l = 0; l < L; ++l) {
        Rcpp::NumericMatrix leaf_mat(nTips, C);
        const std::vector<double>& leaf_vec = leaf_storage[l];
        for (int v = 0; v < nTips; ++v) {
            for (int c = 0; c < C; ++c) {
                leaf_mat(v, c) = leaf_vec[static_cast<size_t>(v) + static_cast<size_t>(c) * nTips];
            }
        }
        leaf_list[l] = leaf_mat;
    }

    std::vector<double> edge_counts(static_cast<size_t>(C) * static_cast<size_t>(C), 0.0);
    for (const auto& vec : edge_counts_local) {
        for (size_t idx = 0; idx < edge_counts.size(); ++idx) edge_counts[idx] += vec[idx];
    }
    Rcpp::NumericMatrix edge_counts_mat(C, C);
    for (int r = 0; r < C; ++r) {
        for (int k = 0; k < C; ++k) edge_counts_mat(r, k) = edge_counts[static_cast<size_t>(r) * C + k];
    }

    return Rcpp::List::create(
        Rcpp::Named("leaf_beliefs") = leaf_list,
        Rcpp::Named("edge_counts") = edge_counts_mat,
        Rcpp::Named("logZ") = logZ_total
    );
}
