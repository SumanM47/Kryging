function out = genMat(xmin,xmax,nvec,nu,range,rval,theta)
         dxy = (xmax - xmin)./(nvec-1);
         dx = dxy(1); dy = dxy(2);

         if nu == 0.5
            rhofunc = @(d) rval.*exp(-sqrt(sum((dxy.*d).^2))./range);
         else
            rhofunc = @(d) matern(sqrt(sum((dxy.*d).^2)),nu,range,rval);
         end

         [out1,out2,~,~] = stationary_Gaussian_process(nvec(1),nvec(2),rhofunc);

         out = out1(:);
    
end

