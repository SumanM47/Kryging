% This is a demonstration for the estimation method described in 
%   "Kryging: Geostatistical analysis of large-scale datasets
%     using Krylov subspace methods" - Majumder et al. (2020+)

%% Add Paths

run('startup_genGK.m')


%% Create Fake Data: Gridded

% n : number of points on one side of a square grid.
n = 100;

% xmin: min of lat-lon; xmax: max of lat-lon; nvec: gridsize
xmin = [0 0]; xmax = [1 1]; nvec = [n n];

% Creating the mesh and the coordinates on the grid
[xx,yy] = meshgrid(linspace(xmin(1),xmax(1),nvec(1)),linspace(xmin(2),xmax(2),nvec(2)));
s = [xx(:) yy(:)];

%Computing the distance matrix and exponential covariance based on that
D = squareform(pdist(s));

Q = exp(-D/0.1);

% Creating x and the data
x_act = mvnrnd(zeros(n^2,1),Q,1);
dall = x_act(:) + normrnd(0,sqrt(0.1),[n^2,1]);

% Holding out 5% of the data for test sample
d = dall(:);
mis_ind = randsample(n^2,round(0.05*n^2));
d(mis_ind) = NaN;


% Start timing
tic

% Provide A. Here it is the identity
A = speye(size(s,1));

% Matern smoothness parameter; to be supplied by user
nu = 0.5;

% Tuning parameter k
maxit = 100;

% Create predictors
X = ones(size(s,1),1);

% Set initial values
vd = var(d,'omitnan');
rho = 0.05*sqrt(sum((xmax - xmin).^2,'all'));
initvec = [0.1*vd,rho,0.9*vd,mean(d,'omitnan')];

%% Call the function for Gridded data

[x,y,theta,numit,tim,xl,xu,yl,yu] = Kryging_wocv(d,s,A,nu,X,initvec,maxit); % Errors from non-positive definite embedding may occur
%[x,y,theta,numit,tim] = Kryging_wocv(d,s,A,nu,X,initvec,maxit); % This only affects the SE estimates. Uncomment and use this if only estimates are desired and not CIs.


%% Multiple initial values and Cross validation
%% Takes much much longer; not parallelized
%rhovec = [0.01 0.05 0.1 0.15].*sqrt(sum((xmax - xmin).^2,'all'));
%initmat = [(0.1*vd).*ones(4,1);rhovec;(0.9*vd).*ones(4,1);mean(d,'omitnan').*ones(4,1)]';
%[x,y,theta,numit,tim,xl,xu,yl,yu] = Kryging_wcv(d,s,A,nu,X,initvec,maxit);

tim2 = toc;

% Check performance

xrmse = sqrt(mean((x - x_act(:)).^2)); % May produce errors when Kryging produced errors due to non-positive definite embedding
xcov = mean(xl <= x_act(:) & x_act(:) <= xu); % May produce errors when Kryging produced errors due to non-positive definite embedding

yrmse = sqrt(mean((y(mis_ind) - dall(mis_ind)).^2,'all')); % May produce errors when Kryging produced errors due to non-positive definite embedding
ycov = mean(yl(mis_ind) <= dall(mis_ind) & dall(mis_ind) <= yu(mis_ind),'all'); % May produce errors when Kryging produced errors due to non-positive definite embedding

fprintf('RMSE for x is %f\n', xrmse) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Coverage for x is %f\n', xcov) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('RMSE for y is %f\n', yrmse) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Coverage for y is %f\n', ycov) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Time taken is %f seconds\n', tim) % May produce errors when Kryging produced errors due to non-positive definite embedding


%% Create Fake Data: Irregularly spaced
n = 100;
xmin = [0 0]; xmax = [1 1]; nvec = [n n];
[xx,yy] = meshgrid(linspace(xmin(1),xmax(1),nvec(1)),linspace(xmin(2),xmax(2),nvec(2)));
sall = [xx(:) yy(:)];

D = squareform(pdist(s));

Q = exp(-D/0.1);

x_act = mvnrnd(zeros(n^2,1),Q,1);
dgrid = x_act(:) + normrnd(0,sqrt(0.1),[n^2,1]);

subpoints = randsample(n^2,round(0.1*n^2));
ns = size(subpoints,1);
dall = dgrid(subpoints);

s = sall(subpoints,:);

d = dall(:);
mis_ind = randsample(ns,round(0.05*ns));
d(mis_ind) = NaN;

tic
A = speye(size(s,1));

nu = 0.5;

maxit = 100;

ngvec = [100 100];
% Create predictors
X = ones(size(s,1),1);

% Set initial values
vd = var(d,'omitnan');
rho = 0.05*sqrt(sum((xmax - xmin).^2,'all'));
initvec = [0.1*vd,rho,0.9*vd,mean(d,'omitnan')];


%% Call the function for non-gridded data

[x_irreg,y_irreg,theta_irreg,numit_irreg,tim_irreg,xl_irreg,xu_irreg,yl_irreg,yu_irreg] = Kryging_wocv_irreg(d,s,A,nu,X,initvec,maxit,ngvec); % Errors from non-positive definite embedding may occur
%[x_irreg,y_irreg,theta_irreg,numit_irreg,tim_irreg] = Kryging_wocv_irreg(d,s,A,nu,X,initvec,maxit,ngvec); % same as in regular case, no CIs produced

%% Multiple initial values and Cross validation
%% Takes much much longer; not parallelized
%rhovec = [0.01 0.05 0.1 0.15].*sqrt(sum((xmax - xmin).^2,'all'));
%initmat = [(0.1*vd).*ones(4,1);rhovec;(0.9*vd).*ones(4,1);mean(d,'omitnan').*ones(4,1)]';
%[x,y,theta,numit,tim,xl,xu,yl,yu] = Kryging_wcv_irreg(d,s,A,nu,X,initvec,maxit,ngvec);


tim3 = toc;

xrmse = sqrt(mean((x_irreg - x_act(subpoints)).^2,'all')); % May produce errors when Kryging produced errors due to non-positive definite embedding
xcov = mean(xl_irreg <= x_act(subpoints) & x_act(subpoints) <= xu_irreg,'all'); % May produce errors when Kryging produced errors due to non-positive definite embedding

yrmse = sqrt(mean((y(mis_ind) - dall(mis_ind)).^2,'all')); % May produce errors when Kryging produced errors due to non-positive definite embedding
ycov = mean(yl(mis_ind) <= dall(mis_ind) & dall(mis_ind) <= yu(mis_ind),'all'); % May produce errors when Kryging produced errors due to non-positive definite embedding


fprintf('RMSE for x is %f\n', xrmse) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Coverage for x is %f\n', xcov) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('RMSE for y is %f\n', yrmse) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Coverage for y is %f\n', ycov) % May produce errors when Kryging produced errors due to non-positive definite embedding
fprintf('Time taken is %f seconds\n', tim) % May produce errors when Kryging produced errors due to non-positive definite embedding
