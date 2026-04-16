# Graduation Project

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
将三角网格置于[mesh](./mesh)中，在 Matlab 中运行 [gen_mat](./src/gen_mat.m)。

输入：添加的边数 $e$ 。

输出：在 [matrices](./matrices) 中，矩阵 `matrix.mat` 和 `matrix_0.mat` 分别代表原矩阵和添加边后的矩阵。

### 求解置换向量
在 [src](./src) 目录下运行：

```bash
./run.sh
```

输入：叉数 $k$ 和执行次数 N。

输出：保存在 [result](./result) 中。`perm_vec.mat` 和 `perm_vec_0.mat` 分别代表最后一次循环求得的 `matrix.mat` 和 `matrix_0.mat` 的置换向量；`runtime.txt`用于记录运行时间，其格式为：第一行一个整数 N，代表执行次数；接下来 N 行，每行 $4$ 个浮点数，代表每次执行的总初始化时间、首次 k-way 分解时间、首次 amd 时间、二次重排序时间；最后一行为上面 N 行的平均值。


### 结果可视化
在 Matlab 中运行 [draw_mat](./src/draw_mat.m) 以得到置换后矩阵稀疏模式的可视化结果。

## 测试结果
已测试的数据在文件夹 [tested](./tested) 中，每个子文件夹代表一组测试。子文件夹命名为 `网格名称_e`，子文件夹中包含测试矩阵和若干记录运行时间的txt文件，命名规则为 `runtime_k.txt`。
