function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wocv(d,s,A,nu,X,Initvec,k)

%
%upper level function for Kryging_est - preps user inputs for the estimation procedure
%   does not use cross-validation
%	
%Usage:
%[est_x,est_y,est_theta,numit,tim] = Kryging_wocv(d,s,A,nu,X,Initvec,k)
%[est_x, est_y, est_theta, numit, tim, xl,xu,yl,yu] = Kryging_wocv(d,s,A,nu,X,Initvec,k)
%
%Input:
%d: vector of observations
%s: observation locations
%A: matrix of relationship between Y and the latent true state values
%nu: smoothness parameter for the matern covariance to be input by the user
%X: matrix of covariates
%Initvec: vector of initial values for the optimization process - nugget, range, partial sill and mean parameters in original scale
%k: order of Krylov subspace
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
%Written for and used in "Kryging: Geostatistical analysis of massive spatial datasets using Krylov subspaces" - Majumder et al. (2020+)
%

ns = size(s,1);
[~,sort_ind] = sortrows(s);
W = sparse(sort_ind,(1:ns),ones(ns,1),ns,ns);
Anew = A*W;
R = speye(ns);

X = [ones(ns,1) X];

xmin = min(s);
xmax = max(s);
nvec = [size(unique(s(:,1)),1) size(unique(s(:,2)),1)];

theta = [1 1];

nmi = find(~isnan(d));

tic
if nargout > 5
[x_est, y_est, theta_est, numit, xl, xu, yl, yu] = Kryging_est(Initvec(4:end),Initvec(3),Initvec(2),nu,Initvec(1),d(:),X,nmi,xmin,xmax,nvec,Anew,R,k,theta,W);
else
[x_est, y_est, theta_est, numit] = Kryging_est(Initvec(4:end),Initvec(3),Initvec(2),nu,Initvec(1),d(:),X,nmi,xmin,xmax,nvec,Anew,R,k,theta,W);
end
tim=toc;
end
