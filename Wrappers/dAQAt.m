function dd = dAQAt(A,Q)

%
%function to compute the diagonal entries of AQA'
%	
%Usage:
%dd = dAQAt(A,Q)
%
%Input:
%A: sparse matrix
%Q: matrix of type funMat
%
%Output:
%dd: vector of diagonal values of A*Q*A'
%
%
%Written for and used in "Kryging: Geostatistical analysis of large-scale datasets using Krylov subspace methods" - Majumder et al. (2020+)
%


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
