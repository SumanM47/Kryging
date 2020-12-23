function [U, B, V,dBmu,dBsigma_e,dBrange,dVmu,dVsigma_e,dVrange,dUmu,dUsigma_e,dUrange] = genGKB_withgrad(A, Q, R, U, B, V,dQ, sigma_e, dBmu,dBsigma_e,dBrange,dVmu,dVsigma_e,dVrange,dUmu,dUsigma_e,dUrange, options)
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
  ARiu = A'*(R\U(:,k));
  if k == 1
    v = ARiu;
  else
    v = ARiu - B(k, k-1)*V(:,k-1);
  end
  
  % Reorthogonalize V
  for j = 1:k-1
    temp = (V(:,j)'*(Q*v));
    v = v - temp*V(:,j);
  end
  alpha = normM(v,Q);
  v = v / alpha;
  temp = Q*v;
  if k==1
      dalphasigma_e = 0;
      dalpharange = (normM(alpha*v,dQ)^2)/(2*alpha);
      temp1 = (A'*(R\dUmu(:,k)));
      dalphamu = temp1'*temp;
      
      dvsigma_e = zeros(size(v,1),1);
      dvrange = -v*dalpharange;
      dvmu = temp1/alpha - dalphamu*v/alpha;
  else
      temp1 = -ARiu/sigma_e + A'*(R\dUsigma_e(:,k)) - dBsigma_e(k,k-1)*V(:,k-1) - B(k,k-1)*dVsigma_e(:,k-1);
      dalphasigma_e = v'*(Q*temp1);
      dvsigma_e = temp1/alpha - dalphasigma_e*v/alpha ;
      
      temp1 = A'*(R\dUrange(:,k)) - dBrange(k,k-1)*V(:,k-1) - B(k,k-1)*dVrange(:,k-1);
      dalpharange = (normM(alpha*v,dQ)^2)/(2*alpha) + v'*(Q*temp1);
      dvrange = temp1/alpha - dalpharange*v/alpha;
      
      temp1 = A'*(R\dUmu(:,k)) - dBmu(k,k-1)*V(:,k-1) - B(k,k-1)*dVmu(:,k-1);
      dalphamu = v'*(Q*temp1);
      dvmu = temp1/alpha - dalphamu*v/alpha;
      
  end
  u = A*temp - alpha*U(:,k);
  
  % Reorthogonalize U
  for j = 1:k
    temp = (U(:,j)'*(R\u));
    u = u - temp*U(:,j);
  end
  beta = normM(u, @(x)R\x);
  u = u / beta;
  
  temp1 = A*(Q*dvsigma_e) - alphasigma_e*U(:,k) - alpha*dUsigma_e(:,k);
  dbetasigma_e = -beta/(2*sigma_e) + u'*(R\temp1);
  dusigma_e = temp1/beta - dbetasigma_e*u/beta;
  
  temp1 = A*(dQ*v) + A*(Q*dvrange) - dalpharange*U(:,k) - alpha*dUrange(:,k);
  dbetarange = u'*(R\temp1);
  durange = temp1/beta - dbetarange*u/beta;
  
  temp1 = A*(Q*dvmu) - alphamu*U(:,k) - alpha*dUmu(:,k);
  dbetamu = u'*(R\temp1);
  dumu = temp1/beta - dbetamu*u/beta;
  
  
  U = [U, u];
  dUsigma_e = [dUsigma_e, dusigma_e];
  dUrange = [dUrange, durange];
  dUmu = [dUmu, dumu];
  
else % Do not need reorthogonalization, save on storage
   ARiu = A'*(R\U(:));
  if k == 1
    v = ARiu;
  else
    v = ARiu - B(k, k-1)*V(:,k-1);
  end
  alpha = normM(v,Q);
  v = v / alpha;
  temp = Q*v;
  if k==1
      dalphasigma_e = 0;
      dalpharange = (normM(alpha*v,dQ)^2)/(2*alpha);
      temp1 = (A'*(R\dUmu(:)));
      dalphamu = temp1'*temp;
      
      dvsigma_e = zeros(size(v,1),1);
      dvrange = -v*dalpharange;
      dvmu = temp1/alpha - dalphamu*v/alpha;
  else
      temp1 = -ARiu/sigma_e + A'*(R\dUsigma_e(:)) - dBsigma_e(k,k-1)*V(:,k-1) - B(k,k-1)*dVsigma_e(:,k-1);
      dalphasigma_e = v'*(Q*temp1);
      dvsigma_e = temp1/alpha - dalphasigma_e*v/alpha ;
      
      temp1 = A'*(R\dUrange(:)) - dBrange(k,k-1)*V(:,k-1) - B(k,k-1)*dVrange(:,k-1);
      dalpharange = (normM(alpha*v,dQ)^2)/(2*alpha) + v'*(Q*temp1);
      dvrange = temp1/alpha - dalpharange*v/alpha;
      
      temp1 = A'*(R\dUmu(:)) - dBmu(k,k-1)*V(:,k-1) - B(k,k-1)*dVmu(:,k-1);
      dalphamu = v'*(Q*temp1);
      dvmu = temp1/alpha - dalphamu*v/alpha;
      
  end
  u = A*temp - alpha*U(:);
  
  beta = normM(u, @(x)R\x);
  u = u / beta;
  
  temp1 = A*(Q*dvsigma_e) - dalphasigma_e*U(:) - alpha*dUsigma_e(:);
  dbetasigma_e = -beta/(2*sigma_e) + u'*(R\temp1);
  dusigma_e = temp1/beta - dbetasigma_e*u/beta;
  
  temp1 = A*(dQ*v) + A*(Q*dvrange) - dalpharange*U(:) - alpha*dUrange(:);
  dbetarange = u'*(R\temp1);
  durange = temp1/beta - dbetarange*u/beta;
  
  temp1 = A*(Q*dvmu) - dalphamu*U(:) - alpha*dUmu(:);
  dbetamu = u'*(R\temp1);
  dumu = temp1/beta - dbetamu*u/beta;
  
  U = u(:);
  dUsigma_e = dusigma_e(:);
  dUrange = durange(:);
  dUmu = dumu(:);
  
end

V = [V, v];
dVsigma_e = [dVsigma_e, dvsigma_e];
dVrange = [dVrange, dvrange];
dVmu = [dVmu, dvmu];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
dBsigma_e = [dBsigma_e, [zeros(k-1,1);dalphasigma_e];[zeros(1,k-1), dbetasigma_e]];
dBrange = [dBrange, [zeros(k-1,1);dalpharange];[zeros(1,k-1), dbetarange]];
dBmu = [dBmu, [zeros(k-1,1);dalphamu];[zeros(1,k-1), dbetamu]];

end


function nrm = normM(v, M)
if isa(M, 'function_handle')
  Mv = M(v);
else
  Mv = M*v;
end
nrm = sqrt(v'*Mv);
end