function k = dmatern(r,nu,ell,psts)
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

scale = sqrt(2*nu).*r./ell;
fact  = 2.^(1-nu)./gamma(nu);

k = psts.*fact.*(-(2*nu)*(scale.^nu).*besselk(nu,scale)/ell + ((scale.^(nu+1)).*besselk(nu+1,scale))./ell);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf


%For nu = inf, the square exponential covariance kernel
if nu == inf
  k = psts*exp(-((r./ell).^2)/2).*(r.^2)./(ell^3);
end
k(isnan(k)) = 0;

end