function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wocv_irreg(d,s,A,nu,X,Initvec,k,nvec)

%
%upper level function for Kryging_est for irregular data - preps user inputs for the estimation procedure
%   does not use cross-validation
%	
%Usage:
%[est_x,est_y,est_theta,numit,tim] = Kryging_wocv_irreg(d,s,A,nu,X,Initvec,k,nvec)
%[est_x, est_y, est_theta, numit, tim, xl,xu,yl,yu] = Kryging_wocv_irreg(d,s,A,nu,X,Initvec,k,nvec)
%
%Input:
%d: vector of observations
%s: observation locations
%A: matrix of relationship between Y and the latent true state values
%nu: smoothness parameter for the matern covariance to be input by the user
%X: matrix of covariates
%Initvec: vector of initial values for the optimization process - nugget, range, partial sill and mean parameters in original scale
%k: order of Krylov subspace
%nvec: size of the underlying grid
%
%Output:
%est_x: estimated values of the true state process 
%est_y: estimated and predicted values of the observed process
%est_theta: estimated mean and spatial parameters
%numit: number of iterations required for the optimization process
%tim: time taken for the estimation process
%xl: lower 95% confidence limit for the true state process estimates (optional)
%xu: upper 95% confidence limit for the true state process estimates (optional)
%yl: lower 95% confidence limit for the observed process estimates (optional)
%yu: upper 95% confidence limit for the observed process estimates (optional)
%
%
%Written for and used in "Kryging: Geostatistical analysis of large-scale datasets using Krylov subspace methods" - Majumder et al. (2020+)
%

tic
ns = size(s,1);
dmaxx = 1; dmaxy = 1;
xmin = min(s);
xmax = max(s);
[W1,news,newxmin,newxmax,newnvec] = createA(s,xmin,xmax,nvec,dmaxx,dmaxy);
newns = size(news,1);
[~,sort_ind] = sortrows(news);
W = W1*sparse(sort_ind,(1:newns),ones(newns,1),newns,newns);
Anew = A*W;R = speye(ns);

%X = [ones(ns,1) X];


%nvec = [size(unique(s(:,1)),1) size(unique(s(:,2)),1)];

theta = [1 1];

nmi = find(~isnan(d));

if nargout > 5
[x_est, y_est, theta_est, numit, xl, xu, yl, yu] = Kryging_est(Initvec(4:end),Initvec(3),Initvec(2),nu,Initvec(1),d(:),X,nmi,newxmin,newxmax,newnvec,Anew,R,k,theta,W);
else
[x_est, y_est, theta_est, numit] = Kryging_est(Initvec(4:end),Initvec(3),Initvec(2),nu,Initvec(1),d(:),X,nmi,newxmin,newxmax,newnvec,Anew,R,k,theta,W);
end
tim=toc;
end
