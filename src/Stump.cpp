#include "Stump.h"
#include "utils.h"

#include <amd.h>
#include <omp.h>
#include <chrono>
#include <cassert>
#include <algorithm>
#include <iostream>

Stump::~Stump() = default;

void Stump::clear()
{
    num_nodes = 0;
    stump_init = false;
    DOF2node.clear();
    is_cut_vertex.clear();
}

void Stump::init(
    int num_branch,
    int M_n,
    int* Mp,
    int* Mi,
    std::vector<int>& mesh_perm)
{
    assert(M_n >= 0);
    assert(Mp != nullptr);
    assert(Mi != nullptr);

    num_nodes = (num_branch > M_n) ? M_n : num_branch;
    if (num_nodes <= 0) {
        clear();
        return;
    }

    DOF2node.assign(M_n, 0);

    std::vector<int> assigned_DOF(M_n);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M_n; i++) assigned_DOF[i] = i;

    auto tk0 = std::chrono::high_resolution_clock::now();
    bool ok = this->decompose(M_n, Mp, Mi, assigned_DOF);
    auto tk1 = std::chrono::high_resolution_clock::now();
    t_kway_ms = std::chrono::duration<double, std::milli>(tk1 - tk0).count();

    mesh_perm.resize(M_n);
    auto ta0 = std::chrono::high_resolution_clock::now();
    if (M_n > 0) {
        std::vector<int> Iperm(M_n, 0);
        this->permute(M_n, Mp, Mi, mesh_perm.data(), Iperm.data());
    }
    auto ta1 = std::chrono::high_resolution_clock::now();
    t_amd_ms = std::chrono::duration<double, std::milli>(ta1 - ta0).count();

    stump_init = ok;
}

void Stump::permute(
    int M_n,
    int* Mp,
    int* Mi,
    int* perm,
    int* Iperm)
{
    assert(M_n >= 0);
    assert(Mp != nullptr);
    assert(perm != nullptr);

    if (M_n == 0) return;

    int NNZ = Mp[M_n];
    if (NNZ == 0) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M_n; i++) perm[i] = i;
    } else {
        amd_order(M_n, Mp, Mi, perm, nullptr, nullptr);
    }

    if (Iperm) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M_n; ++i) {
            Iperm[perm[i]] = i;
        }
    }
}

bool Stump::decompose(
    int M_n,
    int* Mp,
    int* Mi,
    const std::vector<int>& dofs)
{
    assert(M_n > 0);
    assert(Mp != nullptr && Mi != nullptr);
    assert((int)dofs.size() == M_n);
    assert(num_nodes > 0);

    idx_t nvtxs = static_cast<idx_t>(M_n);
    idx_t ncon = 1;
    idx_t nparts = static_cast<idx_t>(num_nodes);
    idx_t objval = 0;
    std::vector<idx_t> part(M_n, 0);

    int ret = METIS_PartGraphKway(
        &nvtxs, &ncon,
        reinterpret_cast<idx_t*>(Mp),
        reinterpret_cast<idx_t*>(Mi),
        nullptr, nullptr, nullptr,
        &nparts, nullptr, nullptr, nullptr,
        &objval, part.data());

    if (ret != METIS_OK) {
        std::cerr << "METIS_PartGraphKway failed, ret=" << ret << std::endl;
        return false;
    }

    // ---- DOF2node ----
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < M_n; ++i) {
        int global_dof = dofs[i];
        int p = static_cast<int>(part[i]);
        if (global_dof >= 0 && global_dof < (int)DOF2node.size()) {
            DOF2node[global_dof] = static_cast<idx_t>(p);
        }
    }

    // ---- is_cut_vertex ----
    is_cut_vertex.assign(M_n, 0);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int u = 0; u < M_n; ++u) {
        int pu = static_cast<int>(part[u]);
        for (int e = Mp[u]; e < Mp[u + 1]; ++e) {
            int v = Mi[e];
            if (v >= 0 && v < M_n) {
                if (static_cast<int>(part[v]) != pu) {
                    is_cut_vertex[u] = 1;
                    break;
                }
            }
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// buildSubgraph — 提取子图（计数 + 填充）
// ---------------------------------------------------------------------------
static void buildSubgraph(
    const std::vector<int>& verts,
    int /*M_n*/,
    int* Mp_new, int* Mi_new,
    std::vector<int>& g2l,
    std::vector<int>& sub_Mp,
    std::vector<int>& sub_Mi)
{
    int sub_n = static_cast<int>(verts.size());

    // 建立 global -> local 映射
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < sub_n; ++i) g2l[verts[i]] = i;

    // 计数每个局部顶点的邻居数
    sub_Mp.assign(sub_n + 1, 0);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int li = 0; li < sub_n; ++li) {
        int gi = verts[li];
        int cnt = 0;
        for (int e = Mp_new[gi]; e < Mp_new[gi + 1]; ++e) {
            if (g2l[Mi_new[e]] >= 0) ++cnt;
        }
        sub_Mp[li + 1] = cnt;
    }

    // 前缀和
    for (int i = 0; i < sub_n; ++i) sub_Mp[i + 1] += sub_Mp[i];

    // 填充
    sub_Mi.resize(sub_Mp[sub_n]);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int li = 0; li < sub_n; ++li) {
        int gi = verts[li];
        int pos = sub_Mp[li];
        for (int e = Mp_new[gi]; e < Mp_new[gi + 1]; ++e) {
            int lj = g2l[Mi_new[e]];
            if (lj >= 0) sub_Mi[pos++] = lj;
        }
    }

    // 清零映射
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < sub_n; ++i) g2l[verts[i]] = -1;
}

// ---------------------------------------------------------------------------
// amdAndUpdatePerm — 子图 AMD + 回填全局 mesh_perm
// ---------------------------------------------------------------------------
static void amdAndUpdatePerm(
    const std::vector<int>& verts,
    std::vector<int>& sub_Mp,
    std::vector<int>& sub_Mi,
    const std::vector<int>& is_in_sub,   // int 数组，1=属于子图
    int M_n,
    std::vector<int>& mesh_perm)
{
    int sub_n = static_cast<int>(verts.size());
    if (sub_n == 0) return;

    // AMD
    std::vector<int> sub_perm(sub_n);
    if (sub_Mp[sub_n] == 0) {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < sub_n; ++i) sub_perm[i] = i;
    } else {
        amd_order(sub_n, sub_Mp.data(), sub_Mi.data(),
                  sub_perm.data(), nullptr, nullptr);
    }

    // 找出 mesh_perm 中属于该子图顶点的位置
    std::vector<int> positions;
    positions.reserve(sub_n);
    for (int i = 0; i < M_n; ++i) {
        if (is_in_sub[mesh_perm[i]]) positions.push_back(i);
    }

    // 回填
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < sub_n; ++i) {
        mesh_perm[positions[i]] = verts[sub_perm[i]];
    }
}

// ---------------------------------------------------------------------------
// redecompose
// ---------------------------------------------------------------------------
void Stump::redecompose(
    int M_n,
    int* Mp_old, int* Mi_old,
    int* Mp_new, int* Mi_new,
    std::vector<int>& mesh_perm)
{
    assert(M_n >= 0);
    assert(stump_init);

    // ==== Step 1: 识别新增边 ====
    std::vector<std::pair<int,int>> added_edges;
    ComputeAddedEdges(M_n, Mp_new, Mi_new, Mp_old, Mi_old, added_edges);

    if (added_edges.empty()) return;

    // ==== Step 2: 对子图分类 ====
    std::vector<int> mark(num_nodes, 0);

    for (const auto& edge : added_edges) {
        int u = edge.first, v = edge.second;
        int pu = static_cast<int>(DOF2node[u]);
        int pv = static_cast<int>(DOF2node[v]);

        if (pu != pv || is_cut_vertex[u] || is_cut_vertex[v]) {
            mark[pu] = 1;
            mark[pv] = 1;

            // 切割顶点：将其直接关联的所有子图也标记为第一类
            if (is_cut_vertex[u]) {
                for (int e = Mp_new[u]; e < Mp_new[u + 1]; ++e) {
                    int w = Mi_new[e];
                    if (w >= 0 && w < M_n)
                        mark[static_cast<int>(DOF2node[w])] = 1;
                }
            }
            if (is_cut_vertex[v]) {
                for (int e = Mp_new[v]; e < Mp_new[v + 1]; ++e) {
                    int w = Mi_new[e];
                    if (w >= 0 && w < M_n)
                        mark[static_cast<int>(DOF2node[w])] = 1;
                }
            }
        } else {
            if (mark[pu] != 1) mark[pu] = 2;
        }
    }

    std::vector<int> type1_parts;
    for (int p = 0; p < num_nodes; ++p) {
        if (mark[p] == 1) type1_parts.push_back(p);
    }
    int k_type1 = static_cast<int>(type1_parts.size());

    std::vector<int> g2l(M_n, -1);

    // ==== Step 3: 第一类子图 —— 统一 kway + AMD ====
    if (k_type1 > 0) {
        // 并行标记哪些顶点属于第一类
        std::vector<int> is_type1(M_n, 0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M_n; ++i) {
            if (mark[static_cast<int>(DOF2node[i])] == 1)
                is_type1[i] = 1;
        }

        // 收集顶点列表
        std::vector<int> type1_verts;
        type1_verts.reserve(M_n);
        for (int i = 0; i < M_n; ++i) {
            if (is_type1[i]) type1_verts.push_back(i);
        }
        int sub_n = static_cast<int>(type1_verts.size());

        // 提取子图
        std::vector<int> sub_Mp, sub_Mi;
        buildSubgraph(type1_verts, M_n, Mp_new, Mi_new, g2l, sub_Mp, sub_Mi);

        // METIS kway 重分解
        if (k_type1 > 1 && sub_n > k_type1) {
            idx_t nvtxs = static_cast<idx_t>(sub_n);
            idx_t ncon  = 1;
            idx_t nparts = static_cast<idx_t>(k_type1);
            idx_t objval = 0;
            std::vector<idx_t> part(sub_n, 0);

            int ret = METIS_PartGraphKway(
                &nvtxs, &ncon,
                reinterpret_cast<idx_t*>(sub_Mp.data()),
                reinterpret_cast<idx_t*>(sub_Mi.data()),
                nullptr, nullptr, nullptr,
                &nparts, nullptr, nullptr, nullptr,
                &objval, part.data());

            if (ret == METIS_OK) {
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < sub_n; ++i) {
                    DOF2node[type1_verts[i]] =
                        static_cast<idx_t>(type1_parts[static_cast<int>(part[i])]);
                }
            }
        } else if (k_type1 == 1) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < sub_n; ++i) {
                DOF2node[type1_verts[i]] = static_cast<idx_t>(type1_parts[0]);
            }
        }

        amdAndUpdatePerm(type1_verts, sub_Mp, sub_Mi, is_type1, M_n, mesh_perm);
    }

    // ==== Step 4: 第二类子图 —— 逐子图 AMD ====
    for (int p = 0; p < num_nodes; ++p) {
        if (mark[p] != 2) continue;

        std::vector<int> is_part(M_n, 0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M_n; ++i) {
            if (static_cast<int>(DOF2node[i]) == p)
                is_part[i] = 1;
        }

        std::vector<int> verts;
        for (int i = 0; i < M_n; ++i) {
            if (is_part[i]) verts.push_back(i);
        }
        if (verts.empty()) continue;

        std::vector<int> sub_Mp, sub_Mi;
        buildSubgraph(verts, M_n, Mp_new, Mi_new, g2l, sub_Mp, sub_Mi);
        amdAndUpdatePerm(verts, sub_Mp, sub_Mi, is_part, M_n, mesh_perm);
    }

    is_cut_vertex.assign(M_n, 0);
    #pragma omp parallel for schedule(dynamic, 64)
    for (int u = 0; u < M_n; ++u) {
        int pu = static_cast<int>(DOF2node[u]);
        for (int e = Mp_new[u]; e < Mp_new[u + 1]; ++e) {
            int v = Mi_new[e];
            if (v >= 0 && v < M_n && static_cast<int>(DOF2node[v]) != pu) {
                is_cut_vertex[u] = 1;
                break;
            }
        }
    }
}