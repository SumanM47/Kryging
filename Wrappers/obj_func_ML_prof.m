function [objval, dell, Hess] = obj_func_ML_prof(th, xmin, xmax, nvec, nu, theta, R, A, Y, maxit)

% [objval, dell, Hess] = obj_func_ML_prof(th, xmin, xmax, nvec, nu, theta, R, A, Y, maxit)
% 
% computes the approximate profile log-likelihood, approximates the gradient and Hessian
% 
% Input:
%         th - current value of the parameters
%         xmin - min of lat-lon
%         xmax - max of lat-lon
%         nvec - gridsize
%         nu - smoothness parameter, user supplied
%         theta - anisotropy parameters
%         R - covariance matrix for the observation process
%         A - matrix
%         Y - observations vector
%         maxit - tuning parameter k
% Output:
%         objval - value of the objective function, negative profile log-likelihood
%         dell - gradient of the objective function at current value of th
%         Hess - Hessian of the objective function at current value of th
 
  sigma_e = th(1);
  ell = th(2);
  lambda = th(3);
  mu = th(4);

  if nu==0.5
     k = @(r) expcov(r,ell,1);
     dk = @(r) dexpcov(r,ell,1);
  else
      if nu-0.5 == round(nu-0.5) && nu~=0.5
          k = @(r) matern_halfint(r,nu,ell,1);
          dk = @(r) dmatern_halfint(r,nu,ell,1);
      else
          k = @(r) matern(r,nu,ell,1);
          dk = @(r) dmatern(r,nu,ell,1);
      end
  end

  ns = prod(nvec(:));
 
  Qr = createrow(xmin,xmax,nvec,k,theta);
  Qfun = @(x)toeplitzproduct(x, Qr, nvec);
  Qsub = @(x) toeplitzsub(x,Qr);
  Q = funMat(Qfun,Qfun,[ns ns],Qsub);
  
  dQr = createrow(xmin,xmax,nvec,dk,theta);
  dQfun = @(x)toeplitzproduct(x, dQr, nvec);
  dQsub = @(x) toeplitzsub(x,dQr);
  dQ = funMat(dQfun,dQfun,[ns ns],dQsub);
 
  input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth', 'off', 'Mu',mu, 'Lambda',lambda,'Sigma_e', sigma_e, 'Grad',false);
  [x_genHyBR_opt,outpar] = genHyBR_new(A, Y(:), Q, R, input);

  logdetval = logdet(Qr,nvec);
  quadval = outpar.quad^2;
  res = Y(:) - A*x_genHyBR_opt(:);
  errval = normM(res(:),@(x)R\x)^2;
  
  objval = logdetval + (errval*sigma_e) + (quadval*lambda) -  size(Y,1)*log(sigma_e) - ns*log(lambda) + (log(sigma_e))^2 + (log(lambda))^2 + (log(ell))^2;   
  if nargout > 1
	  ff = outpar.f;
	  B = outpar.B;
	  V = outpar.V;
	  U = outpar.U;

  
	  try
	    [W,D] = eig(B'*B);
	  catch ME
	    warning(ME.message);
	    [W,D,~] = svd(B'*B);
	    D = abs(D);
	  end	
  
	  Delt = D./(D + lambda);
	  ImDelt = eye(size(D,1)) - Delt;
  
	  Rires = R\res;
	  AtRires = A'*Rires;
	  ZtAtRires = W'*(B'*(U'*Rires));
	  Wtf = W'*ff;
  
	  delsig = -size(Y,1)/sigma_e + errval - (2/lambda)*normM(AtRires,Q)^2 + (2/lambda)*normM(ZtAtRires,Delt)^2 + 2*Wtf'*(ImDelt*ZtAtRires) + 2*log(sigma_e)/sigma_e;
  
	  Vf = (V*ff);
	  dQVf = dQ*Vf;
	  VtdQVf = V'*dQVf;
	  WtVtdQVf = W'*VtdQVf;
  
	  DWtVtdQVf = (Delt*WtVtdQVf);
  
	  dl = getnumdifflogdet(ell,xmin,xmax,nvec,theta,nu,1e-4);
  
	  delrho = dl - 2*sigma_e*AtRires'*dQVf + 2*sigma_e*ZtAtRires'*DWtVtdQVf + lambda*ff'*VtdQVf - 2*lambda*Wtf'*DWtVtdQVf + 2*log(ell)/ell;
  
  
	  dellambda = -ns/lambda + 2*(sigma_e/lambda)*Wtf'*(ImDelt*ZtAtRires) + quadval - 2*sum(ff.^2) + 2*normM(Wtf,Delt)^2 + 2*log(lambda)/lambda;
  
	  sumVW = sum(V*W);
  
	  delmu = -2*sigma_e*sum(AtRires) + 2*sigma_e*sumVW*(Delt*ZtAtRires) - 2*lambda*sumVW*(Delt*Wtf);
  
	  dell = [delsig; delrho; dellambda; delmu];


  	  if nargout > 2
		Hess = dell*dell';
	  end
  end
end



%% Subfunction for normM

function nrm = normM(v, M)
if isa(M, 'function_handle')
  Mv = M(v);
else
  Mv = M*v;
end
nrm = sqrt(v'*Mv);
end
