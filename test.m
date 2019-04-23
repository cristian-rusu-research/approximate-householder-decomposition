%% test Householder factorization algorithms

%% clear everything
close all
clear
clc

%% size of the transform
n = 32;

%% generate random orthonormal transformation
[Q, ~] = qr(randn(n));
[V, D] = eig(Q);

%% generated random symmetric (positive-definite) transformation
S = randn(n);
S = S'*S;

%% how many Householder reflectors in the decomposition, fix it now
% h = 2*round(log2(n));

%% number of pozitive and negative eigenvalues, this is optimal h
nplus = length(find(real(diag(D)) > 0));
nminus = length(find(real(diag(D)) < 0));

%% decompose the orthonormal Q into h Householder reflectors
[U1, X1_1, X2_1, theVs1, tus1, err1, theVsoriginal1] = optimizeHouseholder_decomposition(Q, nminus);
[U2, X1_2, X2_2, theVs2, tus2, err2, theVsoriginal2] = optimizeHouseholder_decomposition(-Q, nplus);

%% the Symmetric Householder Factorization (SHF) algorithm
% number of reflectors
h = 2*round(log2(n));
% allow for spectrum update?
changeSpectrum = 1;
% allow for a diagonal matrix update?
changeD = 1;
% initial reflectors
reflectors = zeros(n, h);
% initialization of the spectrum
s = diag(S);

[reflectors, s, d, val, U] = shf(S, h, changeSpectrum, changeD, s, reflectors);
