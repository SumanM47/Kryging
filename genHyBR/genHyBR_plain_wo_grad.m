%% Function for genHyBR with supplied value of lambda

function [x_out, output] = genHyBR_plain_wo_grad(A, b, Q, R, options)

%
% [x_out, output] = genHyBR(A, b, Q, R, options)
%
% genHyBR is a generalized Hybrid Bidiagonalization Regularization method used for
% solving large-scale, ill-posed inverse problems of the form:
%               b = A*x + noise
% The method combines an iterative generalize Golub-Kahan (GK) Bidiagonalization Method
% with a SVD-based regularization method to stabilize the semiconvergence
% behavior that is characteristic of many ill-posed problems.
%
% Inputs:
%     A : either (a) a full or sparse matrix
%                (b) a matrix object that performs mat*vec and mat'*vec
%         Note: If A is a function handle, create an object called funMat
%         e.g., A = funMat(Afun, Atfun) where Afun and Atfun are function
%                                       handles for A*vec and A'*vec
%     b : rhs vector
%  Q, R : covariance matrices, Q is either (a) a full or sparse matrix or
%                (b) a matrix object that performs matrix*vector operations
%         Note: If Q is a function handle, create an object called funMat
%         e.g., Q = funMat(Qfun, Qtfun) where Qfun and Qtfun are function
%                                       handles for Q*vec and Q'*vec
%
% options : structure with the following fields (optional)
%         Iter - maximum number of GK iterations:
%                       [ positive integer | {min(m,n,100)} ]
%         Reorth - reorthogonalize Lanczos subspaces: [on | {off}]
%         x_true - True solution : [ array | {off} ]
%                Returns error norms with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%        Vx - extra space needed for finding optimal reg. parameters
%        ResTol - Residual tolerance for stopping the LBD iterations,
%                    similar to the stopping criteria from [1]: [atol, btol]
%                   [non-negative scalar  | {[10^-6, 10^-6]}]
%           Mu - Mean Parameter [double | {0}]
%         Lambda - sill of the covariance process
%                   [non-negative scalar | {0}]
%          Sigma_e - Noise [positive scalar | {1}]
%          Range - Range of the spatial process [positive scalar | {1}]
%           Grad - Indicator for computing gradient [logical | {false}]
%
%       Note: options is a structure created using the function 'HyBRset'
%               (see 'HyBRset' for more details)
%
% Outputs:
%      x_out : computed solution
%     output : structure with the following fields:
%      iterations - stopping iteration (options.Iter | GCV-determined)
%            Enrm - relative error norms (requires x_true)
%            Rnrm - relative residual norms
%            Xnrm - relative solution norms
%             U,V - genGK basis vectors
%               B - bidiagonal matrix from genGK
%            flag - a flag that describes the output/stopping condition:
%                       1 - performed max number of iterations
%                       2 - achieved residual tolerance
%        comp_tol - a matrix for values computed for comparing tolerances
%
% References:
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%   [3] Chung, Nagy and O'Leary, "A Weighted-GCV Method for Lanczos-Hybrid
%       Regularization", ETNA 28 (2008), pp. 149-167.
%   [4] Chung and Saibaba. "Generalized Hybrid Iterative Methods for 
%       Large-Scale Bayesian Inverse Problems", submitted 2016
%
% J.Chung and J. Nagy 3/2007
% J. Chung and A. Saibaba, modified for gen-HyBR 2017



defaultopt = struct('Iter', [] , 'Reorth', 'off', 'x_true', 'off',...
  'Vx' , [], 'ResTol', [10^-6, 10^-6], 'Mu',0, 'Lambda', 0, 'Sigma_e',1, 'Grad', false, 'Dqrange','none');

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
  x_out = defaultopt;
  return;
end

% Check for acceptable number of input arguments
if nargin < 4
  error('genHyBR: Not Enough Inputs')
elseif nargin < 5
  options = [];
end
if isempty(options)
  options = defaultopt;
end

% Get options:
[m,n] = size(A);
defaultopt.Iter = min([m, n, 100]);
options = HyBR_plain_lsmrset(defaultopt, options);


maxiter = HyBR_plain_lsmrget(options,'Iter',[],'fast');
x_true = HyBR_plain_lsmrget(options,'x_true',[],'fast');
restol = HyBR_plain_lsmrget(options,'ResTol',[],'fast');
mu = HyBR_plain_lsmrget(options,'Mu',[],'fast');
lambda = HyBR_plain_lsmrget(options,'Lambda',[],'fast');
sigma_e = HyBR_plain_lsmrget(options,'Sigma_e',[],'fast');
grad = HyBR_plain_lsmrget(options,'Grad',[],'fast');
dqrange = HyBR_plain_lsmrget(options,'Dqrange',[],'fast');
notrue = strcmp(x_true,{'off'});

grad = false;

lastopt = struct('Iter', maxiter , 'Reorth', 'on', 'x_true', x_true,...
  'Vx' , [], 'ResTol', restol, 'Mu',mu, 'Lambda', lambda, 'Sigma_e',sigma_e, 'Grad', false, 'Dqrange','none');


if grad && strcmp(dqrange,'none')
    error("Supply a value for dQ");
end

R = sigma_e*R;



%--------------------------------------------
%  The following is needed for RestoreTools:
if isa(A, 'psfMatrix')
  bSize = size(b);
  b = b(:);
  A.imsize = bSize;
  if ~notrue
    xTrueSize = size(x_true);
    x_true = x_true(:);
  end
end
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if ~notrue
  nrmtrue = norm(x_true(:));
end

if nargout>2 && ~grad
    grad = true;
end
if nargout ==1 && grad
    grad = false;
end

grad = false;

% Set-up output parameters:
outputparams = nargout>1;
if outputparams
  output.iterations = maxiter;
  output.Enrm = ones(maxiter,1);
  output.Rnrm = ones(maxiter,1);
  output.Xnrm = ones(maxiter,1);
  output.U = [];
  output.V = [];
  output.B = [];
  output.flag = 3;
  output.comp_tol = zeros(maxiter,2);
end

A1 = sum(A,2);
if mu~=0
    b = b - mu*A1;
end

%Define GK bidiagonalization function
beta = normM(b,@(x)R\x);
U = (1 / beta)*b;
GKhandle = @genGKB;
B = []; V = []; x_out = [];
terminate = 1; norma = 0; normr = beta;
% if grad
%     GKhandle = @genGKB_withgrad;
%     dBmu = []; dBsigma_e =[]; dBrange =[];
%     dVmu = []; dVsigma_e = []; dVrange = [];
%     dbetamu = (-A1'*(R\b))/beta; dbetasigma_e = -beta/sigma_e;
%     dUmu = -A1/beta - (U/beta)*dbetamu; dUsigma_e = U/sigma_e; dUrange = zeros(size(b,1),1);
% end

for i = 1:maxiter+1 %Iteration (i=1) is just an initialization
    if i == maxiter+1
       [U, B, V] = feval(GKhandle, A, Q, R, U, B, V, lastopt);
    else
       [U, B, V] = feval(GKhandle, A, Q, R, U, B, V, options);
    end

  vector = (beta*eye(size(B,2)+1,1));
  if(i >= 2)
   if(lambda == 0)
       f = B \ vector;
  else
       f = (B'*B + lambda*speye(size(B,2)))\(B'*vector);
  end
  Vf = (V*f);
  x = Q*Vf;
    r = b(:) - A*x(:);
    normr = norm(r(:));
    
    norma = norm([norma B(i,i) B(i+1,i)]);
    if isa(A,'function_handle')
      normar = norm(A(r,'transp'));
    else
      normar = norm(A'*r);
    end
    normx = norm(x(:)+mu);
    
    quant1 = restol(1)*beta+restol(2)*norma*normx;
    quant2 = normar/(norma*normr);
    
    if outputparams
      if ~notrue
        output.Enrm(i,1) = norm(x(:)+mu-x_true(:))/nrmtrue;
      end
      output.Rnrm(i,1) = normr;
      output.Xnrm(i,1) = norm(x(:)+mu);
      output.comp_tol(i-1,:) = [quant1 quant2];
      output.quad = normM(Vf,Q);
    end
    if normr <= quant1 || quant2 <= restol(2) && terminate
      if notrue %Set all the output parameters and return
        if outputparams
          output.U = U;
          output.V = V;
          output.B = B;
          output.iterations = i-1;
          output.flag = 2;
        end
        
        if isa(A, 'psfMatrix')
          x_out = reshape(x + mu, bSize);
        else
          x_out = x + mu;
        end
        %  End of new stuff needed for RestoreTools
        %--------------------------------------------
        return
        else % Residual says stop, but continue since have x_true
        if outputparams
          output.iterations = i-1;
          output.flag = 2;
        end
      end
      terminate = 0; % Solution is already found!
    end
  else
      f = (B'*B + lambda*speye(size(B,2)))\(B'*vector);
      Vf =(V*f); 
    x = Q*Vf;
    if outputparams
      if ~notrue
        output.Enrm(i,1) = norm(x(:)+mu-x_true(:))/nrmtrue;
      end
      output.Rnrm(i,1) = normr;
      output.Xnrm(i,1) = norm(x(:)+mu);
      output.quad = normM(Vf,Q);
    end
  end
end

% if grad
%    BtBplI = (B'*B + lambda*speye(size(B,2)));
%    delx(:,1) = Q*(dVsigma_e*f) + dbetasigma_e*x/beta + Q*(V*(BtBplI\((dBsigma_e'*B + B'*dBsigma_e)*f))) + Q*(V*(BtBplI\(dBsigma_e'*vector)));
%    delx(:,2) = dqrange*(V*f) + Q*(dVrange*f) + Q*(V*(BtBplI\((dBrange'*B + B'*dBrange)*f))) + Q*(V*(BtBplI\(dBrange'*vector)));
%    delx(:,3) = Q*(V*(BtBplI\f));
%    delx(:,4) = 1 + Q*(dVmu*f) + dbetamu*x/beta + Q*(V*(BtBplI\((dBmu'*B + B'*dBmu)*f))) + Q*(V*(BtBplI\(dBmu'*vector)));    
% end

if isempty(x_out) % GCV did not stop the process, so we reached max. iterations
        x_out = x + mu;
end

%--------------------------------------------
%  The following is needed for RestoreTools:
%
    if isa(A, 'psfMatrix')
        x_out = reshape(x + mu, bSize);
    end
    
    if outputparams
        output.U = U;
        output.V = V;
        output.B = B;
    end


end



%% Subfunctions neede3d for genHyBR_plain
function nrm = normM(v, M)
if isa(M, 'function_handle')
  Mv = M(v);
else
  Mv = M*v;
end
nrm = sqrt(v'*Mv);
end