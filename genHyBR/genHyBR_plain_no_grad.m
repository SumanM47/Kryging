function [x_out, output] = genHyBR_plain_wo_grad(A, b, Q, R, options)

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

lastopt = HyBR_plain_lsmrset(options,'Reorth','on');

R = sigma_e*R;

if ~notrue
  nrmtrue = norm(x_true(:));
end


outputparams = nargout>1;
if outputparams
  output.U = [];
  output.V = [];
  output.B = [];
  output.quad = 0;
end

A1 = sum(A,2);
if mu~=0
    b = b(:) - mu*A1(:);
end

%Define GK bidiagonalization function
beta = normM(b,@(x)R\x);
U = (1 / beta)*b;
GKhandle = @genGKB;
B = []; V = []; x_out = [];

for i = 1:maxiter+1

%if i == (maxiter+1)
%[U, B, V] = feval(GKhandle, A, Q, R, U, B, V, lastopt);
%else
[U, B, V] = feval(GKhandle, A, Q, R, U, B, V, options);
%end
vector = (beta*eye(size(B,2)+1,1));

if i>=2
if lambda == 0
       f = B \ vector;
  else
       f = (B'*B + lambda*speye(size(B,2)))\(B'*vector);
  end
  Vf = (V*f);
  x = Q*Vf;
else
	f = (B'*B + lambda*speye(size(B,2)))\(B'*vector);
        Vf =(V*f); 
        x = Q*Vf;
end

end

x_out = x + mu;

if outputparams
output.U = U;
output.B = B;
output.V = V;
output.quad = normM(Vf,Q);
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