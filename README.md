# parparth

# 环境配置
目前只支持Ubuntu，环境配置命令为

'''bash
sudo apt update
sudo apt install -y build-essential cmake pkg-config
sudo apt install -y \
  libeigen3-dev \
  libmetis-dev \
  libsuitesparse-dev \
  libmatio-dev \
  gfortran
'''

# 运行说明
## 数据生成
在Matlab中运行[gen_mat](./matrices/gen_mat.m)

当前数据生成参数为
