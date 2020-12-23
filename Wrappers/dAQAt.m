function dd = dAQAt(A,Q)
    [nsa,nsa2] = size(A);
    nsq = size(Q,1);
    if nsa2~=nsq
        error('dimension mismatch');
    end
    dd = zeros(nsa,1);
    for i=1:nsa
        a = A(i,:);
        nzi = find(a);
        w = nonzeros(a);
        Qsubmat = Q(nzi);
        dd(i) = w'*(Qsubmat*w);
    end
end