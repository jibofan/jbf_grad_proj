#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <cassert>
#include <omp.h>

// =========================================================================
//  ComputeAddedEdges — 识别两帧之间新增的边
//  OpenMP 并行：每个线程维护局部 added_edges，最后合并
// =========================================================================
inline void ComputeAddedEdges(
    int M_n,
    const int* Mp,      const int* Mi,
    const int* Mp_prev, const int* Mi_prev,
    std::vector<std::pair<int,int>>& added_edges)
{
    assert(M_n >= 0);
    assert(Mp != nullptr && Mi != nullptr);
    assert(Mp_prev != nullptr && Mi_prev != nullptr);

    added_edges.clear();

    #pragma omp parallel
    {
        // 线程私有的局部结果
        std::vector<std::pair<int,int>> local_edges;

        #pragma omp for schedule(dynamic, 64) nowait
        for (int u = 0; u < M_n; ++u) {
            int cur_begin  = Mp[u],      cur_end  = Mp[u + 1];
            int prev_begin = Mp_prev[u], prev_end = Mp_prev[u + 1];
            int cur_num  = cur_end  - cur_begin;
            int prev_num = prev_end - prev_begin;

            // 快速路径
            if (cur_num == prev_num) {
                bool same = true;
                for (int k = 0; k < cur_num; ++k) {
                    if (Mi[cur_begin + k] != Mi_prev[prev_begin + k]) {
                        same = false;
                        break;
                    }
                }
                if (same) continue;
            }

            // 双指针归并
            int ci = cur_begin, pi = prev_begin;
            while (ci < cur_end) {
                if (pi < prev_end && Mi_prev[pi] == Mi[ci]) {
                    ++ci; ++pi;
                } else if (pi < prev_end && Mi_prev[pi] < Mi[ci]) {
                    ++pi;
                } else {
                    int v = Mi[ci];
                    if (u < v) local_edges.emplace_back(u, v);
                    ++ci;
                }
            }
        }

        // 合并到全局结果
        #pragma omp critical
        {
            added_edges.insert(added_edges.end(),
                               local_edges.begin(), local_edges.end());
        }
    }
}

// =========================================================================
//  LaplacianToMpMi — Laplace 稀疏矩阵 → 邻接 CSR
//  OpenMP 并行：排序去重 / Mi 填充
// =========================================================================
template <typename Scalar, int Options, typename StorageIndex>
void LaplacianToMpMi(
    const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& L,
    std::vector<int>& Mp,
    std::vector<int>& Mi)
{
    if (L.rows() != L.cols()) {
        throw std::invalid_argument("LaplacianToMpMi: input must be a square matrix.");
    }

    const int n = static_cast<int>(L.rows());
    Mp.assign(n + 1, 0);
    Mi.clear();

    std::vector<std::vector<int>> adj(n);

    // 稀疏矩阵遍历（涉及 adj[i] 和 adj[j] 的对称写入，有数据竞争，串行执行）
    for (int outer = 0; outer < L.outerSize(); ++outer) {
        for (typename Eigen::SparseMatrix<Scalar, Options, StorageIndex>::InnerIterator
                 it(L, outer); it; ++it)
        {
            int i = static_cast<int>(it.row());
            int j = static_cast<int>(it.col());
            if (i == j) continue;
            if (it.value() == Scalar(0)) continue;

            if (j >= 0 && j < n) adj[i].push_back(j);
            if (i >= 0 && i < n) adj[j].push_back(i);
        }
    }

    // 并行排序去重（每个顶点独立）
    #pragma omp parallel for schedule(dynamic, 64)
    for (int u = 0; u < n; ++u) {
        auto& nbr = adj[u];
        std::sort(nbr.begin(), nbr.end());
        nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
    }

    // 前缀和（有依赖，串行）
    for (int u = 0; u < n; ++u) {
        Mp[u + 1] = Mp[u] + static_cast<int>(adj[u].size());
    }

    // 并行填充 Mi（每个顶点写入不重叠的区间）
    Mi.resize(Mp[n]);

    #pragma omp parallel for schedule(static)
    for (int u = 0; u < n; ++u) {
        int pos = Mp[u];
        for (int v : adj[u]) {
            Mi[pos++] = v;
        }
    }
}