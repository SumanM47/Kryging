function [objval, dell, Hess] = obj_func_ML_prof(th, xmin, xmax, nvec, nu, theta, R, A, Y, X, maxit)

%
%function to evaluate the negative profile log-likelihood function
%	also evaluates the gradient and hessian
%	
%Usage:
%[objective_val,grad,H] = obj_func_ML_prof(th,xmin,xmax,nvec,nu,theta,R,A,Y,X,maxit)
%
%Input:
%th: set of parameters-- log nugget, log spatial range, log partial sill and mean parameters
%xmin: vector of minimum of the domain in both coordinates
%xmax: vector of maximum of the domain in both coordinates
%nvec: vector of number of unique points in each coordinate
%nu: smoothness parameter for the Matern covariance function
%theta: vector of scale anisotropy parameters-- [1 1] is isotropic
%R: diagonal matrix of observational correlation, typically identity
%A: matrix of relationship between Y and the latent true state values
%Y: vector of observations
%X: matrix of covariates
%maxit: order of Krylov subspace
%
%Output:
%objval: value of the objective function
%dell: gradient of the log-likelihood with respect to the parameters in th
%Hess: Hessian-- computed using outer product of the gradient
%
%
%Written for and used in "Kryging: Geostatistical analysis of large-scale datasets using Krylov subspace methods" - Majumder et al. (2020+)
%



  sigma_e = exp(th(1));
  ell = exp(th(2));
  lambda = exp(th(3));
  beta = th(4:end);

  if nu==0.5
     k = @(r) expcov(r,ell,1);
     dk = @(r) dexpcov(r,ell,1);
     %d2k = @(r) d2expcov(r,ell,1);
  else
      k = @(r) matern(r,nu,ell,1);
      dk = @(r) dmatern(r,nu,ell,1);
      %d2k = @(r) d2matern(r,nu,ell,1);
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
 
  %d2Qr = createrow(xmin,xmax,nvec,d2k,theta);
  %d2Qfun = @(x)toeplitzproduct(x, d2Qr, nvec);
  %d2Qsub = @(x) toeplitzsub(x,d2Qr);
  %d2Q = funMat(d2Qfun,d2Qfun,[ns ns],d2Qsub);
  
  mn = X*beta(:);
  
  input = HyBR_plain_lsmrset('Iter', maxit, 'Reorth', 'off','Mu',0,'Lambda',lambda,'Sigma_e',sigma_e,'Grad',false);
  %solver = 'tikhonov';
  %input = HyBRset('InSolv', solver,'RegPar',lambda,'Iter', maxit);

  [x_genHyBR_opt,outpar] = genHyBR_new(A, Y(:) - mn(:), Q, R, input);
  %[x_genHyBR_opt,outpar] = genHyBR(A, Y(:)-mn(:), Q, R./sigma_e, input);


  [logdetval,dl] = logdet(Qr,dQr,nvec);
  quadval = sum(outpar.f.^2);
  Axhat = A*x_genHyBR_opt(:);
  res = Y(:) - Axhat(:) - mn(:);
  errval = normM(res(:),@(x)R\x)^2;
  
  objval = logdetval + (errval*sigma_e) + (quadval*lambda) -  size(Y,1)*log(sigma_e) - ns*log(lambda);   
  if nargout > 1
	  ff = outpar.f;
	  B = outpar.B;
	  V = outpar.V;
	  %U = outpar.U;

  
	  	
      
          z0 = R\res(:);
          z1 = X'*z0(:);
  
	  delsig = -size(Y,1)/sigma_e + errval;   
	  
          Vf = (V*ff);
	  dQVf = dQ*Vf;
	  VtdQVf = V'*dQVf;
	  
  
	  %dl = getnumdifflogdet(Qr,dQr,nvec);
  
	  delrho = dl - lambda*ff'*VtdQVf;
  
  
	  dellambda = -ns/lambda + quadval;

	  delbeta = -2*sigma_e*z1(:);
  
	  dell = [sigma_e*delsig; ell*delrho; lambda*dellambda; delbeta(:)];


  	  if nargout > 2
                Hess = dell*dell';
                %try
	        %  [W,D] = eig(B'*B);
	        %catch ME
	        %warning(ME.message);
	        %  [~,Dd,W] = svd(B,0);
                %  D = Dd.^2;
                %end

                %Atz0 = A'*z0(:);
                %BtUtz0 = (B'*(U'*z0(:)));
                %ZtAtz0 = W'*BtUtz0(:);

                %Delt = D./(D + lambda);
	  
  
	        %Wtf = W'*ff;

                %WtVtdQVf = W'*VtdQVf;
  
	        %DWtVtdQVf = (Delt*WtVtdQVf);

                %RiX = (R\X);
		%AtRiX = A'*RiX;
                %WtBtUtRiX = W'*(B'*(U'*RiX));

                %d2beta = 2*sigma_e*X'*RiX - 2*((sigma_e^2)/lambda)*(AtRiX'*(Q*AtRiX)) + 2*((sigma_e^2)/lambda)*(WtBtUtRiX'*(Delt*WtBtUtRiX));

                %dbetadlambda = -2*(sigma_e/lambda)*RiX'*Axhat(:) + 2*(sigma_e/lambda)*WtBtUtRiX'*(Delt*Wtf(:));
        
                %dbetadrho = 2*sigma_e*AtRiX'*dQVf(:) - 2*sigma_e*WtBtUtRiX'*DWtVtdQVf(:);

                %dbetadsigma_e = -2*z1(:) + 2*(sigma_e/lambda)*AtRiX'*(Q*Atz0(:)) - 2*(sigma_e/lambda)*WtBtUtRiX'*(Delt*ZtAtz0(:));



               %d2lambda = ns/(lambda^2) - 2*sum(ff.^2)/lambda + 2*Wtf'*(Delt*Wtf)/lambda;

               %dlambdadrho = ff'*VtdQVf - 2*Wtf'*DWtVtdQVf;

               %dlambdadsigma_e = 2*(ff'*BtUtz0)/lambda - 2*(Wtf'*(Delt*ZtAtz0))/lambda;



              %d2rho = d2l - lambda*Vf'*(d2Q*Vf) + 2*lambda*WtVtdQVf'*DWtVtdQVf;

              %drhodsigma_e = -2*dQVf'*Atz0 + 2*DWtVtdQVf'*ZtAtz0;


              %d2sigma_e = size(Y,1)/(sigma_e^2) -2*Atz0'*(Q*Atz0)/lambda + 2*ZtAtz0'*(Delt*ZtAtz0)/lambda;

              %H1 = [d2sigma_e drhodsigma_e dlambdadsigma_e dbetadsigma_e'; drhodsigma_e d2rho dlambdadrho dbetadrho'; ...
              %          dlambdadsigma_e dlambdadrho d2lambda dbetadlambda'];
              %H2 = [dbetadsigma_e(:) dbetadrho(:) dbetadlambda(:) d2beta];
              %Hess = [H1; H2];
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
