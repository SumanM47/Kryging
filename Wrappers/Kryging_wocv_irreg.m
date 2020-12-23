function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wocv_irreg(d,s,A,nu,X,Initvec,k,nvec)

ns = size(s,1);
dmaxx = 1; dmaxy = 1;
xmin = min(s);
xmax = max(s);
[W1,news,newxmin,newxmax,newnvec] = createA(s,xmin,xmax,nvec,dmaxx,dmaxy);
newns = size(news,1);
[~,sort_ind] = sortrows(news);
W = W1*sparse(sort_ind,(1:newns),ones(newns,1),newns,newns);
Anew = A*W;R = speye(ns);

X = [ones(ns,1) X];


%nvec = [size(unique(s(:,1)),1) size(unique(s(:,2)),1)];

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
