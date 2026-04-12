#pragma once

#include <vector>
#include <metis.h> // for idx_t

class Stump
{
public:
    int num_nodes = 0;                    // 分区数
    bool stump_init = false;
    std::vector<int> stump_node_counts;   // 各子图节点数
    std::vector<idx_t> DOF2node;          // global dof -> partition id

    // 标记顶点是否为分割边端点（大小 M_n，1=是，0=否）
    // 使用 int 而非 bool，因为 vector<bool> 按位压缩，并行写不安全
    std::vector<int> is_cut_vertex;

    // Timing (milliseconds), populated by init()
    double t_kway_ms = 0.0;
    double t_amd_ms  = 0.0;

    Stump() = default;
    ~Stump();

    void clear();

    void init(
        int num_branch,              ///<[in] Number of branches (target partitions)
        int M_n,                     ///<[in] Number of nodes
        int* Mp,                     ///<[in] CSC column pointer (size M_n+1)
        int* Mi,                     ///<[in] CSC row indices (size Mp[M_n])
        std::vector<int>& mesh_perm  ///<[out] permutation array (size M_n)
    );

    bool decompose(
        int M_n,                              ///<[in] Number of nodes
        int* Mp,                              ///<[in] CSC column pointer
        int* Mi,                              ///<[in] CSC row indices
        const std::vector<int>& assigned_DOF  ///<[in] local idx -> global dof
    );

    void permute(
        int M_n,      ///<[in] Number of DOFs
        int* Mp,      ///<[in] CSC pointer
        int* Mi,      ///<[in] CSC indices
        int* perm,    ///<[out] Permutation array
        int* Iperm    ///<[out] Inverse permutation array
    );

    /// 增量重分解：在新增边后，对受影响的子图局部重新分解和重排序
    void redecompose(
        int M_n,                     ///<[in] Number of nodes
        int* Mp_old,                 ///<[in] Old CSC column pointer
        int* Mi_old,                 ///<[in] Old CSC row indices
        int* Mp_new,                 ///<[in] New CSC column pointer (with added edges)
        int* Mi_new,                 ///<[in] New CSC row indices (with added edges)
        std::vector<int>& mesh_perm  ///<[in,out] permutation array to update
    );
};