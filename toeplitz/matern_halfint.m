function k = matern_halfint(r,nu,ell,psts)
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

dfn = nu-0.5;
if(~(dfn == round(dfn)))
    error('nu must be half-integer');
end
df = int8(dfn);

fact = sqrt(2*nu)*r/ell;


k = psts*(factorial(df)/factorial(2*df))*exp(-fact).*halfintpoly(fact,df);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf


%For nu = inf, the square exponential covariance kernel
if nu == inf
  k = pasts*exp(-((r/ell).^2)/2);
end
k(isinf(k)) = 1;
k(isnan(k)) = 1;

end

function poly = halfintpoly(fact,df)
if(df==0)
    poly = ones(size(fact),1);
else
   poly = zeros(size(fact),1);
   for i = 0:df
       poly = poly + (factorial(df+i)/(factorial(i)*factorial(df-i)))*(2*fact).^(df-i);             
   end
end
end