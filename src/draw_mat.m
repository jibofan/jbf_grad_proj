clear;clc;close all;

load('../matrices/matrix.mat', 'A');
load('../matrices/matrix_0.mat', 'A_0');

% load('../result/perm_vec.mat', 'perm');
load('../result/perm_vec_0.mat', 'perm');

perm = perm + 1;
p=amd(A_0);

% Ap = A(perm, perm);
A0p = A_0(perm, perm);


figure;
spy(A0p);
title('Sparsity Pattern after Permutation');
xlabel('Column Index');
ylabel('Row Index');
