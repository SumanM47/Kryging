%% startup file for generalized Golub-Kahan codes, 
%   including gen-LSQR and gen-HyBR
%
% These codes were used to generate the figures in the paper:
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)
% Modified for using in the paper 
%   "Kryging: Geostatistical analysis of massive spatial datasets
%     using Krylov subspaces"
%     - Majumder et al. (2020+)

directory = pwd;
path(directory, path)

% RestoreTools
path([directory, '/RestoreTools/Classes'], path)
path([directory, '/RestoreTools/TestData'], path)

% Regularization Tools
path([directory, '/REGU'], path)

% Hybrid Codes
path([directory, '/genHyBR'], path)

% Toeplitz
path([directory, '/toeplitz'], path)

% Wrappers
path([directory, '/Wrappers'], path)


close all
clear all
clc
