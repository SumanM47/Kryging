function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wcv(d,s,A,nu,X,Initmat,k)

%
%upper level function for Kryging_est for gridded data - preps user inputs for the estimation procedure
%   uses cross-validation -- not parallelized
%	
%Usage:
%[est_x,est_y,est_theta,numit,tim] = Kryging_wcv(d,s,A,nu,X,Initvec,k)
%[est_x, est_y, est_theta, numit, tim, xl,xu,yl,yu] = Kryging_wcv(d,s,A,nu,X,Initvec,k)
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
%Written for and used in "Kryging: Geostatistical analysis of large-scale datasets using Krylov subspace methods" - Majumder et al. (2020+)
%


ns = size(s,1);
[~,sort_ind] = sortrows(s);
W = sparse(sort_ind,(1:ns),ones(ns,1),ns,ns);
Anew = A*W;
R = speye(ns);

%X = [ones(ns,1) X];

xmin = min(s);
xmax = max(s);
nvec = [size(unique(s(:,1)),1) size(unique(s(:,2)),1)];

theta = [1 1];

nmi = find(~isnan(d));

MSEarr = zeros([size(Initmat,1),5]);

bs = [50 30]; nb = nvec./bs; nbp = prod(nb(:)); bsp = prod(bs(:));
bm = reshape(1:ns,nvec(1),nvec(2));


for fld=1:size(MSEarr,2)
	rng(1234 + fld)
	%test_ind = randsample(ns,round(0.05*ns));
	test_ind = [];
	bsamp = randsample(nbp,round(0.05*nbp));
	for bk=1:size(bsamp(:),1)
		bkval = bsamp(bk);
		bj = mod(bkval,nb(2));
		if bj == 0
			bj = nb(2);
		end
		bi = (bkval - bj)/nb(2) + 1;
		rinds = bs(1)*(bi-1) + [1:bs(1)];
		cinds = bs(2)*(bj-1) + [1:bs(2)];
		test_ind = [test_ind,reshape(bm(rinds,cinds),1,bsp)];
	end
	dtest = d(test_ind);
	din = d(:);
	din(test_ind) = NaN;

	nmiin = find(~isnan(din));

	for i=1:size(MSEarr,1)
		try
		[x_est, y_est, theta_est,numit] = Kryging_est(Initmat(i,4:end),Initmat(i,3),Initmat(i,2),nu,Initmat(i,1),din(:),X,nmiin,xmin,xmax,nvec,Anew,R,k,theta,W);

		MSEarr(i,fld) = mean((y_est(test_ind) - dtest(:)).^2,'all','omitnan');
		catch ME
		warning(ME.message)
		MSEarr(i,fld) = NaN;
		end
	end
end

MSEvec = mean(MSEarr,2,'omitnan');
disp(MSEvec)
[M,iopt] = min(MSEvec);
disp(iopt)
tic
if nargout > 5
[x_est, y_est, theta_est, numit, xl, xu, yl, yu] = Kryging_est(Initmat(iopt,4:end),Initmat(iopt,3),Initmat(iopt,2),nu,Initmat(iopt,1),d(:),X,nmi,xmin,xmax,nvec,Anew,R,k,theta,W);
else
[x_est, y_est, theta_est, numit] = Kryging_est(Initmat(iopt,4:end),Initmat(iopt,3),Initmat(iopt,2),nu,Initmat(iopt,1),d(:),X,nmi,xmin,xmax,nvec,Anew,R,k,theta,W);
end
tim=toc;
end
