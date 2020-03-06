function dd = dAQAt(A,Q)

% dd = dAQAt(A,Q)
% 
% computes diagonal entries of matrix of the form AQA'
% 
% Input:
%         A - matrix, sparse
%         Q - covariance matrix, funMat type
% Output:
%         dd - diagonal entries of AQA'
        
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