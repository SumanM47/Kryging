function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wocv(d,s,A,nu,X,Initvec,k)

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
