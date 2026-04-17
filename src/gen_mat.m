clear;clc;close all;
[x, t] = readObj('../mesh/bimba');
B = findBoundary(x, t);
nv = size(x, 1);
I = setdiff(1:nv, B);
A = sparse([t t], t(:,[2 3 1 3 1 2]), 1);
A = double( A>0 | A'>0 );
A = -spdiags(-sum(A,2), 0, A);
save('../matrices/matrix.mat', 'A');

ne = nnz(triu(A, 1));

fprintf('=============\n');
fprintf('顶点数:%d\n边数:%d\n',nv, ne);
fprintf('=============\n');

e = input('请输入整数 e: ');

added = 0;
while added < e
    i = randi(nv);
    j = randi(nv);
    if i == j, continue; end
    if i > j, [i, j] = deal(j, i); end
    if A(i, j) ~= 0, continue; end
    A(i, j) = -1;
    A(j, i) = -1;
    A(i, i) = A(i, i) + 1;
    A(j, j) = A(j, j) + 1;
    added = added + 1;
end

A_0 = A;
save('../matrices/matrix_0.mat', 'A_0');

fprintf('生成完成！\n');
fprintf('=============\n');