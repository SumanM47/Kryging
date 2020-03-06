function result = logdet(Qr, N)

% result = logdet(Qr,N)
% 
% computes the log-determinant of Q using its first row and gridsize information
% 
% Input:
%         Qr - first row of Q
%         N - gridsize
% Output:
%         result - log determinant of Q

 dim = numel(N);

 assert(numel(Qr) == prod(N(1:dim)), 'Check dimensions of Toeplitz Row');
 
 switch dim
     case 1
         circ = [Qr(end:-1:2), Qr];
         
         ff = ifft(circ)*size(circ,1);
         ff = ff(1:N);
         result = real(sum(log(ff)));
         
     case 2
         circ = reshape(Qr,N(2),N(1));
         circ = [circ(:,end:-1:2),circ];
         circ = [circ(end:-1:2,:);circ];
         [n1,n2] = size(circ);
         
         ff = ifft2(circ).*(n1*n2);
         ff = ff(1:N(2),1:N(1));
         result = reshape(ff,prod(N(1:2)),1);
         result = real(sum(log(result)));
         
     otherwise
         error('Wrong dimesnion');
 end
end
