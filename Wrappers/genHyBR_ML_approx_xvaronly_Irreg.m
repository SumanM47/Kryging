function [x_est, y_est, theta_est, xl, xu, yl, yu] = genHyBR_ML_approx_xvaronly_Irreg(d,s,nu,A,maxit,nvec,dmaxx,dmaxy)
% [x_est, y_est, theta_est, xl, xu, yl, yu] = genHyBR_ML_approx_xvaronly(d,s,nu,A,maxit)
% 
% high level function that prepares the input and calls in estimate_ML_using_genHyBr_plain_xvaronly function
% 
% Input:
%         d - vector of data, includes missingness
%         s - matrix of coordinates, includes prediction sites
%         nu - smoothness parameter, user supplied
%         A - matrix
%         maxit - tuning parameter k
%           nvec - gridsize, user supplied
%           dmaxx - number of gridcells to extend on x-axis
%           dmaxy - number of gridcells to extend on y-axis
% Output:
%         x_est - estimated x
%         y_est - estimated y
%         theta_est - estimated spatial parameters
%         xl - lower bound for 95% CI for x
%         xu - upper bound for 95% CI for x
%         yl - lower bound for 95% CI for y
%         yu - upper bound for 95% CI for y

%if nargout > 7
%    error('Too many output arguments');
%end
ns = size(s,1);
R = speye(ns);

nmi = find(~isnan(d));

xmin = min(s);
xmax = max(s);
%nvec = [size(unique(s(:,1)),1) size(unique(s(:,2)),1)];

[W1,news,newxmin,newxmax,newnvec] = createA(s,xmin,xmax,nvec,dmaxx,dmaxy);
newns = size(news,1);
[~,sort_ind] = sortrows(news);
W = W1*sparse(sort_ind,(1:newns),ones(newns,1),newns,newns);
Anew = A*W;


theta = [1 1];
%tic
%[muhat,sill_init,range_init,sigma_e_init] = Get_initval(d(nmi),s(nmi,:),A(nmi,:),nu,25000);
%toc
muhat = mean(d(:),'omitnan')/mean(sum(Anew(nmi,:),2));
vnm = var(d(nmi));
sill_init = 0.9*vnm;
sigma_e_init = 0.1*vnm;
range_init = 0.05*sqrt(sum((xmax - xmin).^2));
if nargout > 3
    [x_est, y_est, theta_est, xl, xu, yl, yu] = estimate_ML_using_genHyBr_plain_xvaronly(muhat,sill_init,range_init,nu,sigma_e_init,d,nmi,newxmin,newxmax,newnvec,Anew,R,maxit,theta,W);
else
    [x_est, y_est, theta_est] = estimate_ML_using_genHyBr_plain_xvaronly(muhat,sill_init,range_init,nu,sigma_e_init,d,nmi,newxmin,newxmax,newnvec,Anew,R,maxit,theta,W);
end
end
