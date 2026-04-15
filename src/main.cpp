#include "Stump.h"
#include "utils.h"
#include <matio.h>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
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

    // ---- Manual input ----
    int kway_parts = 0, N = 0;
    std::cout << "Enter branch factor (first arg of init): ";
    std::cin >> kway_parts;
    std::cout << "Enter number of runs N: ";
    std::cin >> N;

    // ---- Load matrices (once) ----
    std::cout << "Loading matrices..." << std::endl;
    Eigen::SparseMatrix<double> L_old = readMAT(data_dir + "/matrix.mat");
    Eigen::SparseMatrix<double> L_new = readMAT(data_dir + "/matrix_0.mat");
    int M_n = (int)L_old.rows();

    std::vector<int> Mp_old, Mi_old, Mp_new, Mi_new;
    LaplacianToMpMi(L_old, Mp_old, Mi_old);
    LaplacianToMpMi(L_new, Mp_new, Mi_new);
    std::cout << "Matrices loaded. n = " << M_n << std::endl;

    // ---- Storage for timing results ----
    std::vector<double> t_init(N), t_kway(N), t_amd(N), t_redecomp(N);
    std::vector<int> mesh_perm, mesh_perm_init;

    // ---- Run N independent iterations ----
    for (int i = 0; i < N; ++i) {
        std::cout << "--- Run " << (i + 1) << " / " << N << " ---" << std::endl;

        Stump s;
        mesh_perm.clear();

        // init
        auto t1 = clk::now();
        s.init(kway_parts, M_n, Mp_old.data(), Mi_old.data(), mesh_perm);
        auto t2 = clk::now();
        t_init[i] = std::chrono::duration<double, std::milli>(t2 - t1).count();
        t_kway[i] = s.t_kway_ms;
        t_amd[i]  = s.t_amd_ms;

        std::cout << "  init:  " << t_init[i] << " ms  (kway: "
                  << t_kway[i] << " ms, amd: " << t_amd[i] << " ms)" << std::endl;

        mesh_perm_init = mesh_perm; // snapshot after init, before redecompose

        // redecompose
        auto t3 = clk::now();
        s.redecompose(M_n, Mp_old.data(), Mi_old.data(),
                      Mp_new.data(), Mi_new.data(), mesh_perm);
        auto t4 = clk::now();
        t_redecomp[i] = std::chrono::duration<double, std::milli>(t4 - t3).count();

        std::cout << "  redecompose: " << t_redecomp[i] << " ms" << std::endl;

        // Save permutation vectors from the last run only
        if (i == N - 1) {
            writePermMAT(result_dir + "/perm_vec.mat", mesh_perm_init);
            writePermMAT(result_dir + "/perm_vec_0.mat", mesh_perm);
            std::cout << "Saved perm_vec.mat and perm_vec_0.mat from last run" << std::endl;
        }
    }

    // ---- Write runtime.txt ----
    std::string runtime_path = result_dir + "/runtime.txt";
    std::ofstream ofs(runtime_path);
    if (!ofs) throw std::runtime_error("Cannot open " + runtime_path + " for writing");
    ofs << N << "\n";
    for (int i = 0; i < N; ++i) {
        ofs << t_init[i] << " " << t_kway[i] << " "
            << t_amd[i]  << " " << t_redecomp[i] << "\n";
    }
    // Average line
    double avg_init = 0, avg_kway = 0, avg_amd = 0, avg_redecomp = 0;
    for (int i = 0; i < N; ++i) {
        avg_init     += t_init[i];
        avg_kway     += t_kway[i];
        avg_amd      += t_amd[i];
        avg_redecomp += t_redecomp[i];
    }
    avg_init /= N; avg_kway /= N; avg_amd /= N; avg_redecomp /= N;
    ofs << avg_init << " " << avg_kway << " " << avg_amd << " " << avg_redecomp << "\n";
    ofs.close();
    std::cout << "Timing results saved to " << runtime_path << std::endl;

    return 0;
}