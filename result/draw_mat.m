clear;clc;close all;

load('matrix.mat', 'A');
load('matrix_0.mat', 'A_0');

% load('perm_vec.mat', 'perm');
load('perm_vec_0.mat', 'perm');

perm = perm + 1;

% Ap = A(perm, perm);
A0p = A_0(perm, perm);

figure;
spy(A0p);
title('Sparsity Pattern of Cholesky Factor R after AMD Permutation');
xlabel('Column Index');
ylabel('Row Index');

