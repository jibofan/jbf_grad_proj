clear; clc; close all;

load('../matrices/matrix.mat', 'A');
load('../matrices/matrix_0.mat', 'A_0');

% load('../result/perm_vec.mat', 'perm');
load('../result/perm_vec_0.mat', 'perm');

perm = perm + 1;
p = amd(A_0);

A_0p = A_0;%(perm, perm);

% 计算消元树
parent = etree(A_0p);

n = length(parent);
depth = zeros(1, n);

% depth(j) 表示第 j 个节点到根节点的边数
for j = 1:n
    k = j;
    while parent(k) ~= 0
        depth(j) = depth(j) + 1;
        k = parent(k);
    end
end

tree_height_edges = max(depth);          % 按边数计的高度
tree_height_levels = tree_height_edges + 1;  % 按层数计的高度

figure;
etreeplot(A_0p);

title(sprintf('Elimination Tree, height = %d levels', ...
    tree_height_levels));

xlabel('Column Index');
ylabel('Row Index');

% 在图中添加文字说明
annotation('textbox', [0.15, 0.82, 0.35, 0.08], ...
    'String', sprintf('Height: %d levels', ...
    tree_height_levels), ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white');