function k = dexpcov(r,ell,psts)
%
%  Implementation of the Matern covariance kernel
%
%   Input: r is the radial distance
%         nu is a parameter that controls smoothness of the stochastic process
%         ell controls the correlation length
%         psts is partial sill to sill ratio
%   Output: k - Matern kernel
%
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)


k = psts.*exp(-r./ell).*(r./ell.^2);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf
k(isnan(k)) = 0;

%For nu = inf, the square exponential covariance kernel
end