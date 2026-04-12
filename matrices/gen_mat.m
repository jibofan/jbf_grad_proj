% 功能：
% 1) 随机生成 n 个顶点的无向连通图，拉普拉斯矩阵稀疏表示为 A，保存到 matrix.mat
% 2) 在图上新增 e 条无重边，得到新图拉普拉斯矩阵 A_0，保存到 matrix_0.mat
% 3) 支持手动输入参数 n, lambda(0~1), e

clc;
fprintf('--- 随机无向连通图 + 拉普拉斯矩阵生成 ---\n');

% ===== 1. 输入参数 =====
n = input('请输入顶点数 n (n>=2): ');
lambda = input('请输入稀疏程度 lambda (0~1): ');
e = input('请输入新增边数 e (e>=0): ');

% 参数检查
if isempty(n) || isempty(lambda) || isempty(e)
    error('输入不能为空。');
end
if n < 2 || floor(n) ~= n
    error('n 必须是 >=2 的整数。');
end
if lambda < 0 || lambda > 1
    error('lambda 必须在 [0,1] 区间。');
end
if e < 0 || floor(e) ~= e
    error('e 必须是 >=0 的整数。');
end

% 无向简单图最多边数
maxEdges = n*(n-1)/2;

% ===== 2. 生成“连通”的初始图 =====
% 先构造一棵随机生成树（保证连通，边数 n-1）
adj = sparse(n, n);
perm = randperm(n);
for i = 2:n
    u = perm(i);
    v = perm(randi(i-1)); % 连到前面任意一个点，保证形成树
    adj(u,v) = 1;
    adj(v,u) = 1;
end

% 当前边数
m_current = nnz(triu(adj,1));

% 根据 lambda 计算目标边数（控制稀疏程度）
% lambda 越大，图越稀疏：目标边数越接近 n-1
m_target = round((1-lambda)*maxEdges + lambda*(n-1));
m_target = max(m_target, n-1);
m_target = min(m_target, maxEdges);

% 若目标边数大于当前边数，随机补边（不重复）
edges_to_add = m_target - m_current;
if edges_to_add > 0
    % 所有上三角可能边
    [I, J] = find(triu(ones(n),1));
    allPairs = [I, J];

    % 已存在边
    [Ei, Ej] = find(triu(adj,1));
    existPairs = [Ei, Ej];

    % 找到可添加边集合
    if ~isempty(existPairs)
        [~, loc] = ismember(allPairs, existPairs, 'rows');
        candidatePairs = allPairs(loc==0, :);
    else
        candidatePairs = allPairs;
    end

    if edges_to_add > size(candidatePairs,1)
        error('无法添加足够的初始边，请调整参数。');
    end

    idx = randperm(size(candidatePairs,1), edges_to_add);
    newPairs = candidatePairs(idx, :);
    for k = 1:size(newPairs,1)
        u = newPairs(k,1); v = newPairs(k,2);
        adj(u,v) = 1;
        adj(v,u) = 1;
    end
end

% 初始图拉普拉斯矩阵 A = D - W
deg = sum(adj,2);
A = spdiags(deg, 0, n, n) - adj;
save('matrix.mat', 'A');
fprintf('已保存初始图拉普拉斯矩阵 A 到 matrix.mat\n');

% ===== 3. 在初始图上新增 e 条边（无重边） =====
m_now = nnz(triu(adj,1));
remain = maxEdges - m_now;
if e > remain
    error('新增边数 e 过大。当前最多还能加 %d 条边。', remain);
end

if e > 0
    [I, J] = find(triu(ones(n),1));
    allPairs = [I, J];

    [Ei, Ej] = find(triu(adj,1));
    existPairs = [Ei, Ej];

    [~, loc] = ismember(allPairs, existPairs, 'rows');
    candidatePairs = allPairs(loc==0, :);

    idx = randperm(size(candidatePairs,1), e);
    addPairs = candidatePairs(idx, :);

    for k = 1:size(addPairs,1)
        u = addPairs(k,1); v = addPairs(k,2);
        adj(u,v) = 1;
        adj(v,u) = 1;
    end
end

% 新图拉普拉斯矩阵 A_0
deg0 = sum(adj,2);
A_0 = spdiags(deg0, 0, n, n) - adj;
save('matrix_0.mat', 'A_0');
fprintf('已保存新增边后拉普拉斯矩阵 A_0 到 matrix_0.mat\n');

% 输出简要信息
fprintf('\n--- 结果摘要 ---\n');
fprintf('n = %d, lambda = %.4f, e = %d\n', n, lambda, e);
fprintf('初始图边数: %d\n', nnz(triu(A~=0,1))); % 仅用于参考显示，不是精确边数统计方式
fprintf('新增后边数: %d\n', nnz(triu((A_0 - spdiags(diag(A_0),0,n,n))~=0,1)));
fprintf('完成。\n');
