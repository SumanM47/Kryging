function [est_x, est_y, est_theta, numit, xl,xu,yl,yu] = Kryging_est(bhat,sill_init,range_init,nu_opt,sigma_e_init,d,X,nmi,xmin,xmax,nvec,A,R,maxit,theta,W)

bhatvec = bhat(:);

R_nm = R(nmi,nmi); A_nm = A(nmi,:); dnm = d(nmi);
ftm = @(th) obj_func_ML_prof(th,xmin,xmax,nvec,nu_opt,theta,R_nm, A_nm,dnm(:),X(nmi,:),maxit);
%ftm = @(th) obj_func_ML_prof_no_hess(th,xmin,xmax,nvec,nu_opt,theta,R_nm, A_nm,dnm(:),X(nmi,:),maxit);

options = optimoptions('fminunc','Display','off','Algorithm','trust-region','SpecifyObjectiveGradient',true, 'HessianFcn','objective');
%options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton','SpecifyObjectiveGradient',true);


x0 = [log(1/sigma_e_init),log(range_init),log(1/sill_init),bhatvec'];
[theta_est,fval,~,optout] = fminunc(ftm,x0,options);
numit = optout.iterations;
%numit = 1;
disp(theta_est)


if nu_opt==0.5
    k = @(r) expcov(r,exp(theta_est(2)),1);
else
    k = @(r) matern(r,nu_opt,exp(theta_est(2)),1);
end

ns = prod(nvec(:));

Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
Qsub = @(x) toeplitzsub(x,Qr);
Q = funMat(Qfun,Qfun,[ns ns],Qsub);

bet = theta_est(4:end);
mn = X*bet(:);

if nargout == 4
    input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth','off','Mu',0,'Lambda',exp(theta_est(3)),'Sigma_e',exp(theta_est(1)),'Grad',false);
    x_est = genHyBR_new(A_nm, dnm - mn(nmi), Q, R_nm, input);
    est_x = W*x_est(:); est_y = mn(:) + A*x_est(:);
    est_theta = [exp(theta_est(1)) exp(theta_est(2)) exp(theta_est(3)) theta_est(4:end)]';
end
if nargout > 4
    input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth','off','Mu',0,'Lambda',exp(theta_est(3)),'Sigma_e',exp(theta_est(1)),'Grad',false);
    [x_est,output] = genHyBR_new(A_nm, dnm(:)-mn(nmi), Q, R_nm, input);

    est_x = W*x_est(:);

    
    est_y = mn(:)+A*x_est(:);

    est_theta = [exp(theta_est(1)) exp(theta_est(2)) exp(theta_est(3)) theta_est(4:end)]';

    Bsize = 20;
    mi = find(isnan(d));
    n_p = size(mi,1);
    nxall = size(est_x(:),1);
    nyall = size(est_y(:),1);

    kryg_msey = zeros(n_p,1);
    kryg_msex = zeros(n_p,1);

    for b=1:Bsize
         xboot = genMat(xmin,xmax,nvec,nu_opt,exp(theta_est(2)),1,theta)./sqrt(exp(theta_est(3)));
         %xboot = genMat(xmin,xmax,nvec,nu_opt,0.28,1,theta)./sqrt(exp(theta_est(2)));
         yboot = mn(:) + A*xboot(:) + normrnd(0,1,[nyall,1])./sqrt(exp(theta_est(1)));

         xperm = W*xboot(:);
         xtest = xperm(mi);
         ytest = yboot(mi);

         yin = yboot(nmi);

         x_boot_kryg = genHyBR_new(A_nm, yin(:) - mn(nmi), Q, R_nm, input);
         xperm_kryg = W*x_boot_kryg(:);
         y_boot_kryg = mn(:) + A*x_boot_kryg(:);

         kryg_msex = kryg_msex + ((xperm_kryg(mi) - xtest(:)).^2)./Bsize;
         kryg_msey = kryg_msey + ((y_boot_kryg(mi) - ytest(:)).^2)./Bsize;
    end


    xsd = sqrt(dAQAt(W,Q)./theta_est(2));
    xsd(mi) = sqrt(kryg_msex(:));

    xl = est_x(:) - 1.96*xsd(:);
    xu = est_x(:) + 1.96*xsd(:);

    ysd = sqrt(dAQAt(A,Q)./theta_est(2) + 1/theta_est(1));
    ysd(mi) = sqrt(kryg_msey(:));

    yl = est_y(:) - 1.96*ysd(:);
    yu = est_y(:) + 1.96*ysd(:);

end

end