clear;clc;close all;

load('../matrices/matrix.mat', 'A');
load('../matrices/matrix_0.mat', 'A_0');

% load('../result/perm_vec.mat', 'perm');
load('../result/perm_vec_0.mat', 'perm');

perm = perm + 1;
p=amd(A_0);

% Ap = A(perm, perm);
A_0p = A_0(perm, perm);

A_0p = A_0p + 100000 * speye(size(A_0p));

L = chol(A_0p);
nnz(L)

figure;
spy(A_0p);
title('Sparsity Pattern');
xlabel('Column Index');
ylabel('Row Index');
