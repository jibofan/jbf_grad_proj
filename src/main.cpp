#include "Stump.h"
#include "utils.h"
#include <matio.h>
#include <Eigen/Sparse>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <stdexcept>
#include <sys/stat.h>

/// Read a sparse matrix from a MATLAB .mat file
Eigen::SparseMatrix<double> readMAT(const std::string& filepath)
{
    mat_t* fp = Mat_Open(filepath.c_str(), MAT_ACC_RDONLY);
    if (!fp) throw std::runtime_error("readMAT: cannot open file: " + filepath);
    matvar_t* mv = nullptr;
    while ((mv = Mat_VarReadNext(fp)) != nullptr) {
        if (mv->class_type == MAT_C_SPARSE) break;
        Mat_VarFree(mv); mv = nullptr;
    }
    if (!mv) { Mat_Close(fp); throw std::runtime_error("readMAT: no sparse var"); }
    int nr = (int)mv->dims[0], nc = (int)mv->dims[1];
    mat_sparse_t* sp = (mat_sparse_t*)mv->data;
    std::vector<Eigen::Triplet<double>> T; T.reserve(sp->ndata);
    for (int j = 0; j < nc; ++j)
        for (mat_uint32_t k = sp->jc[j]; k < sp->jc[j+1]; ++k)
            T.emplace_back((int)sp->ir[k], j, ((double*)sp->data)[k]);
    Eigen::SparseMatrix<double> M(nr, nc); M.setFromTriplets(T.begin(), T.end());
    Mat_VarFree(mv); Mat_Close(fp);
    return M;
}

/// Write an integer permutation vector as a dense column vector in a .mat file
void writePermMAT(const std::string& filepath, const std::vector<int>& perm)
{
    mat_t* fp = Mat_CreateVer(filepath.c_str(), NULL, MAT_FT_MAT5);
    if (!fp) throw std::runtime_error("writePermMAT: cannot create: " + filepath);
    size_t dims[2] = {perm.size(), 1};
    std::vector<double> data(perm.begin(), perm.end());
    matvar_t* v = Mat_VarCreate("perm", MAT_C_DOUBLE, MAT_T_DOUBLE,
                                2, dims, data.data(), 0);
    Mat_VarWrite(fp, v, MAT_COMPRESSION_NONE);
    Mat_VarFree(v); Mat_Close(fp);
}

int main()
{
    using clk = std::chrono::high_resolution_clock;
    const std::string data_dir = DATA_DIR;
    const std::string result_dir = std::string(DATA_DIR) + "/../result";
    mkdir(result_dir.c_str(), 0755);

    // ---- Load matrices ----
    Eigen::SparseMatrix<double> L_old = readMAT(data_dir + "/matrix.mat");
    Eigen::SparseMatrix<double> L_new = readMAT(data_dir + "/matrix_0.mat");
    int M_n = (int)L_old.rows();

    std::vector<int> Mp_old, Mi_old, Mp_new, Mi_new;
    LaplacianToMpMi(L_old, Mp_old, Mi_old);
    LaplacianToMpMi(L_new, Mp_new, Mi_new);

    // ---- First permutation: init ----
    Stump s;
    std::vector<int> mesh_perm;
    auto t1 = clk::now();
    s.init(100, M_n, Mp_old.data(), Mi_old.data(), mesh_perm);
    auto t2 = clk::now();
    double ms1 = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "init (first perm): " << ms1 << " ms\n";
    std::cout << "  kway partition: " << s.t_kway_ms << " ms\n";
    std::cout << "  amd permutation: " << s.t_amd_ms << " ms\n";
    writePermMAT(result_dir + "/perm_vec.mat", mesh_perm);

    // ---- Second permutation: redecompose ----
    auto t3 = clk::now();
    s.redecompose(M_n, Mp_old.data(), Mi_old.data(),
                  Mp_new.data(), Mi_new.data(), mesh_perm);
    auto t4 = clk::now();
    double ms2 = std::chrono::duration<double, std::milli>(t4 - t3).count();
    std::cout << "redecompose (second perm): " << ms2 << " ms\n";
    writePermMAT(result_dir + "/perm_vec_0.mat", mesh_perm);

    std::cout << "Saved: " << result_dir << "/perm_vec.mat, perm_vec_0.mat\n";
    return 0;
}