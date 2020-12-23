function k = d2matern(r,nu,ell,psts)
scale = sqrt(2*nu)*r./ell;
fact = 2.^(1-nu)./gamma(nu);

k = psts*((2*nu.*matern(r,nu,ell,1)./(ell^2))-(2*nu.*dmatern(r,nu,ell,1)./ell) ...
    -(fact*(scale.^(nu+1))).*(2*nu + 3).*besselk(nu+1,scale)./(ell^2) ...
    + (fact*(scale.^(nu+2).*besselk(nu+2,scale)))./(ell^2));
k(isnan(k)) = 0;
if nu == inf
    k = psts*exp(-(r/ell).^2).*(((4*(r.^4))/(ell^6))-((6*(r.^2))/(ell^4)));
end
end