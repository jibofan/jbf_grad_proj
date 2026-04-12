clear;clc;close all;
% 手动输入参数
n      = input('请输入矩阵维度 n: ');
lambda = input('请输入稀疏程度 lambda (0-1): ');
e      = input('请输入新增非零元对数 e: ');

assert(lambda >= 0 && lambda <= 1, 'lambda 必须在 0-1 之间');

%% 1. 生成稀疏对称矩阵 A
% 仅在严格上三角按概率 lambda 采样
U = triu(rand(n) < lambda, 1);
A = U | U';              % 对称化
A = A | speye(n);        % 对角线全为 1
A = sparse(double(A));
save('../matrices/matrix.mat', 'A');
fprintf('已保存 matrix.mat\n');

%% 2. 在零位置添加 2*e 个对称非零元得到 A_0
% 找到上三角中当前为 0 的位置
[I, J] = find(triu(~A, 1));
num_zeros = length(I);
if e > num_zeros
    error('e=%d 超过上三角可用零位置数 %d', e, num_zeros);
end

idx = randperm(num_zeros, e);
addI = I(idx); addJ = J(idx);

A_0 = A;
for k = 1:e
    A_0(addI(k), addJ(k)) = 1;
    A_0(addJ(k), addI(k)) = 1;   % 对称
end
save('../matrices/matrix_0.mat', 'A_0');
fprintf('已保存 matrix_0.mat\n\n');

%% 3. 输出图参数（忽略自环）
G = graph(A_0 - diag(diag(A_0)), 'upper');

num_nodes = numnodes(G);
num_edges = numedges(G);
deg       = degree(G);
avg_deg   = mean(deg);
density   = 2*num_edges / (num_nodes*(num_nodes-1));
comps     = conncomp(G);
num_comp  = max(comps);
is_conn   = (num_comp == 1);

fprintf('===== 图参数 (A_0, 不含自环) =====\n');
fprintf('节点数       : %d\n', num_nodes);
fprintf('边数         : %d\n', num_edges);
fprintf('平均度       : %.4f\n', avg_deg);
fprintf('最大/最小度  : %d / %d\n', max(deg), min(deg));
fprintf('图密度       : %.4f\n', density);
fprintf('连通分量数   : %d\n', num_comp);
fprintf('是否连通     : %s\n', string(is_conn));
if is_conn && num_nodes <= 500
    fprintf('图直径       : %d\n', max(max(distances(G))));
end
