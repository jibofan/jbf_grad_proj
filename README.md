# parparth

## 环境配置
目前只支持 Ubuntu，环境配置命令为：

```bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config
sudo apt install -y \
  libeigen3-dev \
  libmetis-dev \
  libsuitesparse-dev \
  libmatio-dev \
  gfortran
```

## 运行说明

### 数据生成
在 Matlab 中运行 [gen_mat](./matrices/gen_mat.m)，生成的矩阵 `matrix.mat` 和 `matrix_0.mat` 分别代表原矩阵和添加边后的矩阵。

### 求置换向量
在 [src](./src) 目录下运行：

```bash
./run_demo.sh
```

求得的置换向量保存在 [result](./result) 中，`perm_vec.mat` 和 `perm_vec_0.mat` 分别代表 `matrix.mat` 和 `matrix_0.mat` 的置换向量。

当前叉数 $k = 100$（可修改 [main](./src/main.cpp) 第 67 行，将 `init` 函数的第一个参数调整为 $k$ ）。

### 结果可视化
在 Matlab 中运行 [draw_mat](./result/draw_mat.m) 以得到置换后矩阵稀疏模式的可视化结果。

## 测试结果
已测试的数据在文件夹 [tested](./tested) 中，每个子文件夹代表一组测试，子文件夹的结构为

m

├── matrix.mat

├── matrix_0.mat

├── perm_vec.mat

├── perm_vec_0.mat

├── time.png ：运行时间

├── spy.png ：添加边后矩阵稀疏模式可视化

└── parameter.txt ：矩阵生成参数，三行分别代表 $n,\lambda,e$
