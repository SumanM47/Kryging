function deriv = getnumdifflogdet(ell,xmin,xmax,nvec,theta,nu,tol)

%     deriv = getnumdifflogdet(ell,xmin,xmax,nvec,theta,nu,tol)
%     
%     performs numerical derivative of log-determinant of Q
%     
%     Input:
%             ell - spatial range
%             xmin - min of lat-lon
%             xmax - max of lat-lon
%             nvec - gridsize
%             theta - anisotropy parameters
%             nu - smoothness parameter
%             tol - bandwidth for the numerical procedure
%     Output:
%             deriv - value of the derivative

    % create the necessary kernels
    if nu == 0.5
        k1 = @(r) expcov(r,ell+tol,1);
        k2 = @(r) expcov(r,ell-tol,1);
    else
        if round(nu-0.5) == nu-0.5 && nu~=0.5
            k1 = @(r) matern_halfint(r,nu,ell+tol,1);
            k2 = @(r) matern_halfint(r,nu,ell-tol,1);
        else
            k1 = @(r) matern(r,nu,ell+tol,1);
            k2 = @(r) matern(r,nu,ell-tol,1);
        end
    end
    
    Qr1 = createrow(xmin,xmax,nvec,k1,theta);
    Qr2 = createrow(xmin,xmax,nvec,k2,theta);
    
    % central difference
    deriv = (logdet(Qr1,nvec) - logdet(Qr2,nvec))/(2*tol);
end