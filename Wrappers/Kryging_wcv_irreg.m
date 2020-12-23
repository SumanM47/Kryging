function [x_est, y_est, theta_est, numit, tim, xl, xu, yl, yu] = Kryging_wcv_irreg(d,s,A,nu,X,Initmat,k,nvec)

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
