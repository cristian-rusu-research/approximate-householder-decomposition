%% test Householder factorization algorithms

%% clear everything
close all
clear
clc

%% size of the transform
n = 32;

%% generate random orthonormal transformation
[Q, ~] = qr(randn(n));

%% generated random symmetric (positive-definite) transformation
S = randn(n);
S = S'*S;

%% how many Householder reflectors in the decomposition
h = 2*round(log2(n));

%% decompose the orthonormal Q into h Householder reflectors
[U, X1, X2, theVs, tus, err, theVsoriginal] = optimizeHouseholder_decomposition(Q, h);

%% the Symmetric Householder Factorization (SHF) algorithm
% allow for spectrum update?
changeSpectrum = 1;
% allow for a diagonal matrix update?
changeD = 1;
% initial reflectors
reflectors = zeros(n, h);
% initialization of the spectrum
s = diag(S);

[reflectors, s, d, val, U] = shf(S, h, changeSpectrum, changeD, s, reflectors);
