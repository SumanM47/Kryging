function k = dexpcov(r,ell,psts)
%
%  Implementation of the derivative of the exponential covariance kernel
%  wrt the spatial range
%
%   Input: r is the radial distance
%         ell controls the correlation length
%         psts is partial sill to sill ratio
%   Output: k - derivative of the exponential kernel
%


k = psts*exp(-r/ell).*(r/ell^2);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf
k(isnan(k)) = 0;
end