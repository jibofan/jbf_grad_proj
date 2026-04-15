[x, t] = readObj('../mesh/cathead');

B = findBoundary(x, t);
nv = size(x, 1);
I = setdiff(1:nv, B);

A = sparse([t t], t(:,[2 3 1 3 1 2]), 1);
A = double( A>0 | A'>0 );
A = -spdiags(-sum(A,2), 0, A);
save('../matrices/matrix.mat', 'A');

e = input('请输入整数 e: ');

U = triu(true(nv), 1);
zeroMask = (A == 0) & U;
[idx_i, idx_j] = find(zeroMask);

numZeroUpper = numel(idx_i);
if e > numZeroUpper
    error('e=%d 超过了上三角零元数量 %d。', e, numZeroUpper);
end

perm = randperm(numZeroUpper, e);
sel_i = idx_i(perm);
sel_j = idx_j(perm);

lin1 = sub2ind([nv,nv], sel_i, sel_j);
lin2 = sub2ind([nv,nv], sel_j, sel_i);
A(lin1) = -1;
A(lin2) = -1;

A(1:nv+1:end) = 0;
d = sum(abs(A), 2);
A = A + spdiags(d, 0, nv, nv);

A_0 = A;
save('../matrices/matrix_0.mat', 'A_0');