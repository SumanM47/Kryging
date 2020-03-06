function k = expcov(r,ell,psts)
%
%  Implementation of the Exponential covariance kernel
%
%   Input: r is the radial distance
%         ell controls the correlation length
%         psts is partial sill to sill ratio
%   Output: k - Exponential kernel
%


k = psts*exp(-r/ell);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf
k(isinf(k)) = 1;

end