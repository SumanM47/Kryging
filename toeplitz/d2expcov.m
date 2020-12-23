function k = d2expcov(r,ell,psts)
k = psts*exp(-r./ell).*(((r.^2)/(ell^4))-((2.*r)/(ell^3)));
end