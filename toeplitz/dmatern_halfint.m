function k = dmatern_halfint(r,nu,ell,psts)
%
%  Implementation of the derivative of the Matern covariance kernel with half-integer smoothness parameters wrt
%  spatial range parameter
%
%   Input: r is the radial distance
%         nu is a parameter that controls smoothness of the stochastic process
%         ell controls the correlation length
%         psts is partial sill to sill ratio
%   Output: k - Matern kernel
%

% Checking if the input smoothness is indeed half-integer
dfn = nu-0.5;
if(~(dfn == round(dfn)))
    error('nu must be half-integer');
end
df = int8(dfn);

fact = sqrt(2*nu)*r/ell;


k = psts*(factorial(df)/factorial(2*df))*exp(-fact).*halfintpoly(fact,df).*fact/ell - psts*(factorial(df)/factorial(2*df))*exp(-fact).*dhalfintpoly(fact,df)/ell;

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf
k(isnan(k)) = 0;

%For nu = inf, the square exponential covariance kernel
if nu == inf
  k = psts*exp(-((r/ell).^2)/2).*(r.^2)/(ell^3);
end

end

% Subfunction needed for dmaternhalfint
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
% Subfunction needed for dmaternhalfint
function poly = dhalfintpoly(fact,df)
    poly = zeros(size(fact),1);
if(df~=0)
   poly = zeros(size(fact),1);
   for i = 0:(df-1)
       poly = poly + (factorial(df+i)/(factorial(i)*factorial(df-1-i)))*(2*fact).^(df-i);             
   end
end
end