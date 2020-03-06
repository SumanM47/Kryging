function [est_x, est_y, est_theta, xl,xu,yl,yu] = estimate_ML_using_genHyBr_plain_xvaronly(muhat,sill_init,range_init,nu_opt,sigma_e_init,d,nmi,xmin,xmax,nvec,A,R,maxit,theta,W)

% [est_x, est_y, est_theta, xl,xu,yl,yu] = estimate_ML_using_genHyBr_plain_xvaronly(muhat,sill_init,range_init,nu_opt,sigma_e_init,d,nmi,xmin,xmax,nvec,A,R,maxit,theta,W)
% 
% performs profile-likelihood maximization to obtain parameter estimates and uses the said parameters for computing x,y and their pointwise variances using genHyBR
% 
% Input:
%         muhat - initial value of mean parameter
%         sill_init - initial value of sill parameter
%         range_init - initial value of range parameter
%         nu_opt - user supplied value of the smoothness parameter
%         sigma_e_init - initial value of the observation noise parameter
%         d - observed values
%         nmi - non-missing indices for d
%         xmin - min of lat-lon
%         xmax - max of lat-lon
%         nvec - gridsize
%         A - matrix
%         R - correlation matrix for observation process
%         maxit - tuning parameter k, user supplied
%         theta - anisotropy parameters
%         W - weight matrix to use
% 
% Output:
%         est_x - estimated x at all locations
%         est_y - estimated y at all locations
%         est_theta - estimated spatial parameters
%         xl - pointwise lower bound for 95% CI for x
%         xu - pointwise upper bound for 95% CI for x
%         yl - pointwise lower bound for 95% CI for y
%         yu - pointwise upper bound for 95% CI for y

%if nargout > 7
%    error('Too many output arguments');
%end

R_nm = R(nmi,nmi); A_nm = A(nmi,:); dnm = d(nmi);
ftm = @(th) obj_func_ML_prof(th,xmin,xmax,nvec,nu_opt,theta,R_nm, A_nm,dnm(:),maxit);
options = optimoptions('fmincon','Display','off','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true, 'HessianFcn','objective');

x0 = [1/sigma_e_init,range_init,1/sill_init,muhat];

theta_est = fmincon(ftm,x0,[],[],[],[],[0 0 0 -Inf],[Inf Inf Inf Inf],[],options);


if nu_opt==0.5
    k = @(r) expcov(r,theta_est(2),1);
else
    if nu_opt-0.5 == round(nu_opt-0.5) && nu_opt~=0.5
        k = @(r) matern_halfint(r,nu_opt,theta_est(2),1);
    else
        k = @(r) matern(r,nu_opt,theta_est(2),1);
    end
end

ns = prod(nvec(:));

Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
Qsub = @(x) toeplitzsub(x,Qr);
Q = funMat(Qfun,Qfun,[ns ns],Qsub);


if nargout == 3
    input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth','on','Mu',theta_est(4), 'Lambda',theta_est(3),'Sigma_e',theta_est(1),'Grad',false);    
    x_est = genHyBR_new(A_nm, dnm, Q, R_nm, input);
    est_x = W*x_est(:); est_y = A*x_est(:); est_theta = theta_est(:);
end
if nargout > 3
    input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth','on','Mu',theta_est(4), 'Lambda',theta_est(3),'Sigma_e',theta_est(1),'Grad',false);
    [x_est,output] = genHyBR_new(A_nm, dnm(:), Q, R_nm, input);
    
    est_x = W*x_est(:);
    V = output.V;
    B = output.B;
    try
      [Wd,D] = eig(B'*B);
    catch ME
      warning(ME.message);
      [Wd,D,~] = svd(B'*B);
      D = abs(D);
    end
    dD = diag(D);
    %Z = Q*(V*Wd);
    Zorg = Q*(V*Wd);
    Z = W*Zorg;
    Delt = dD./(dD+theta_est(3));
    WQWt = dAQAt(W,Q);
    
    %xsd = W*sqrt((1 - ((Z.^2)*Delt))./theta_est(3));
    xsd = sqrt((WQWt(:) - ((Z.^2)*Delt))./theta_est(3));

    xl = est_x(:) - xsd(:).*1.96;
    xu = est_x(:) + xsd(:).*1.96;
    
    est_y = A*x_est(:);
    AZ = A*Zorg;
    dAQAtv = dAQAt(A,Q);
    ysd = sqrt((dAQAtv(:) - ((AZ.^2)*Delt))./theta_est(3) + 1/theta_est(1));
    yl = est_y(:) - ysd(:).*1.96;
    yu = est_y(:) + ysd(:).*1.96;
    %yl = A*xl(:); yu = A*xu(:);

    est_theta = theta_est(:);

end
end