function [U, B, V] = genGKB(A, Q, R, U, B, V, options)
%
%     [U, B, V] = genGKB(A, Q, R, U, B, V, options)
%
%  Perform one step of generalized Golub-Kahan bidiagonalization.
%
% Input:
%          A - matrix
%       Q, R - covariance matrices
%       U, V - accumulation of vectors
%          B - bidiagonal matrix
%    options - structure from HyBR (see HyBRset)
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%
%  Refs:
%   Arioli. "Generalized Golub-Kahan bidiagonalization and stopping
%       criteria. SIMAX, 2013.
%   Chung and Saibaba. "Generalized Hybrid Iterative Methods for 
%       Large-Scale Bayesian Inverse Problems", submitted 2016
%
%   J.Chung and A.Saibaba, 2017

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBR_lsmrget(options,'Reorth'), {'on'});

m = size(U,1);
k = size(B,2)+1;

if reorth % Need reorthogonalization
  if k == 1
    v = A'*(R\U(:,k));
  else
    v = A'*(R\U(:,k)) - B(k, k-1)*V(:,k-1);
  end
  
  % Reorthogonalize V
  for j = 1:k-1
    temp = (V(:,j)'*(Q*v));
    v = v - temp*V(:,j);
  end
  alpha = normM(v,Q);
  v = v / alpha;
  temp = Q*v;
  u = A*temp - alpha*U(:,k);
  
  % Reorthogonalize U
  for j = 1:k
    temp = (U(:,j)'*(R\u));
    u = u - temp*U(:,j);
  end
  beta = normM(u, @(x)R\x);
  u = u / beta;
  U = [U, u];
  
else % Do not need reorthogonalization, save on storage
  if k == 1
    v = A'*(R\U(:));
  else
    v = A'*(R\U(:,k)) - B(k, k-1)*V(:,k-1);
  end
  alpha = normM(v,Q);
  v = v / alpha;
  temp = Q*v;
  u = A*temp - alpha*U(:,k);
  
  beta = normM(u, @(x)R\x);
  u = u / beta;
  %U = u(:);
  U = [U, u];
end

V = [V, v];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];

end


function nrm = normM(v, M)
if isa(M, 'function_handle')
  Mv = M(v);
else
  Mv = M*v;
end
nrm = sqrt(v'*Mv);
end