function deriv = get2numdifflogdet(ell,xmin,xmax,nvec,theta,nu,tol)
    if nu == 0.5
        k1 = @(r) expcov(r,ell+tol,1);
        k2 = @(r) expcov(r,ell-tol,1);
        k = @(r) expcov(r,ell,1);
    else
        k1 = @(r) matern(r,nu,ell+tol,1);
        k2 = @(r) matern(r,nu,ell-tol,1);
        k = @(r) matern(r,nu,ell,1);
    end
    
    Qr1 = createrow(xmin,xmax,nvec,k1,theta);
    Qr2 = createrow(xmin,xmax,nvec,k2,theta);
    Qr = createrow(xmin,xmax,nvec,k,theta);

    deriv = (logdet(Qr1,nvec) + logdet(Qr2,nvec) - 2*logdet(Qr,nvec))/(tol^2);
end