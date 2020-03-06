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

Q = exp(-D/0.2);

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

%% Call the function for Gridded data

[x,y,theta,xl,xu,yl,yu] = genHyBR_ML_approx_xvaronly(d,s,nu,A,maxit);
tim = toc;

% Check performance

xrmse = sqrt(mean((x - x_act(:)).^2));
xcov = mean(xl <= x_act(:) & x_act(:) <= xu);

yrmse = sqrt(mean((y(mis_ind) - dall(mis_ind)).^2,'all'));
ycov = mean(yl(mis_ind) <= dall(mis_ind) & dall(mis_ind) <= yu(mis_ind),'all');

fprintf('RMSE for x is %f\n', xrmse)
fprintf('Coverage for x is %f\n', xcov)
fprintf('RMSE for y is %f\n', yrmse)
fprintf('Coverage for y is %f\n', ycov)
fprintf('Time taken is %f seconds\n', tim)


%% Create Fake Data: Irregularly spaced
n = 100;
xmin = [0 0]; xmax = [1 1]; nvec = [n n];
[xx,yy] = meshgrid(linspace(xmin(1),xmax(1),nvec(1)),linspace(xmin(2),xmax(2),nvec(2)));
s = [xx(:) yy(:)];

D = squareform(pdist(s));

Q = exp(-D/0.2);

x_act = mvnrnd(zeros(n^2,1),Q,1);
dgrid = x_act(:) + normrnd(0,sqrt(0.1),[n^2,1]);

subpoints = randsample(n^2,round(0.1*n^2));
ns = size(subpoints,1);
dall = dgrid(subpoints);

d = dall(:);
mis_ind = randsample(ns,round(0.05*ns));
d(mis_ind) = NaN;

tic
A = speye(size(s,1));

nu = 0.5;

maxit = 100;

ngvec = [100 100];
dmaxx = 1; dmaxy = 1; % recommended settings for maximum effect

%% Call the function for non-gridded data

[x,y,theta,xl,xu,yl,yu] = genHyBR_ML_approx_xvaronly_Irreg(d,s,nu,A,maxit,ngvec,dmaxx,dmaxy);
tim = toc;

xrmse = sqrt(mean((x - x_act(:)).^2,'all'));
xcov = mean(xl <= x_act(:) & x_act(:) <= xu,'all');

yrmse = sqrt(mean((y(mis_ind) - dall(mis_ind)).^2,'all'));
ycov = mean(yl(mis_ind) <= dall(mis_ind) & dall(mis_ind) <= yu(mis_ind),'all');


fprintf('RMSE for x is %f\n', xrmse)
fprintf('Coverage for x is %f\n', xcov)
fprintf('RMSE for y is %f\n', yrmse)
fprintf('Coverage for y is %f\n', ycov)
fprintf('Time taken is %f seconds\n', tim)