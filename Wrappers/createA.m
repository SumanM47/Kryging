function [Abig,news,newxmin,newxmax,newnvec] = createA(s,xmin,xmax,nvec,dmaxx,dmaxy)
xseq = linspace(xmin(1),xmax(1),nvec(1));
yseq = linspace(xmin(2),xmax(2),nvec(2));

ns = size(s,1);

dmx = abs(xseq(1) - xseq(2));
dmy = abs(yseq(1) - yseq(2));

newxmin = xmin - [dmaxx*dmx dmaxy*dmy];
newxmax = xmax + [dmaxx*dmx dmaxy*dmy];
newnvec = nvec + [dmaxx dmaxy].*2;
newxseq = linspace(newxmin(1),newxmax(1),newnvec(1));
newyseq = linspace(newxmin(2),newxmax(2),newnvec(2));
[x,y] = meshgrid(newxseq,newyseq);
news = [x(:) y(:)];
newns = size(news,1);


Arowind = []; Acolind = []; Aval = [];
for i=1:ns
    dx = abs(s(i,1)-news(:,1))/(dmaxx*dmx);
    dy = abs(s(i,2)-news(:,2))/(dmaxy*dmy);
    di = max(dx,dy);
%    tentfval = (1-dx-dy).*(dx + dy <1);
%    tentfval = (1-dx).*(1-dy).*(dx < 1).*(dy < 1).*(4.*dx + 1).*(4.*dy + 1);
    tentfval = ((1-di).^4).*(di < 1).*(4.*di + 1);
    tentfval = tentfval/sum(tentfval);
    nindvec = find(tentfval);
    nnz = size(nindvec,1);
    Acolind = [Acolind nindvec'];
    Arowind = [Arowind i*ones(nnz,1)'];
    Aval = [Aval tentfval(nindvec)'];   
end
Abig = sparse(Arowind, Acolind,Aval, ns,newns);
end