function [tenB, tenT,change] = WSNMSTIPT(tenD, lambda, mu,C,p)

% Solve the WSNM problem 
% ---------------------------------------------
% Input:
%       tenD       -    n1*n2*n3 tensor
%       lambda  -    >0, parameter
%       mu-      regelarization parameter
%       p    -    key parameter of Schattern p norm
%
% Output:
%       tenB       -    n1*n2*n3 tensor
%       tenT       -    n1*n2*n3 tensor
%       change     -    change of objective function value


% initialize
tol = 1e-7; 
max_iter = 500;
rho = 1.1;
max_mu = 1e10;
DEBUG = 1;
normD = norm(tenD(:));
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
% if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end

dim = size(tenD);
tenB = zeros(dim);
tenT = tenB;
Y = tenB;
preTnnT= 0;
NOChange_counter = 0;
weightTenT = ones(size(tenD)) ;
n=min(dim(1),dim(2));
change=zeros(1,max_iter);
for iter = 1 : max_iter
    tenBk = tenB;
    tenTk = tenT;
%% Update low-rank background tensor B
    [tenB] = prox_tnn1(-tenT+tenD-Y/mu,C/mu,p);

%% Update sparse target tensor T
    tenT = prox_non_neg_l1(-tenB+tenD-Y/mu,weightTenT*lambda/mu);
     weightTenT = 1 ./ (abs(tenT) + 0.01); %enhance sparsity
    dY = tenB+tenT-tenD;
    Y = Y + mu*dY;% Update Y
    mu = min(rho*mu,max_mu); % Update mu
%% chg+tubalrank %%运行使用
    stopCriterion = norm(tenD(:) - tenB(:) - tenT(:)) / normD;
change(iter)=(stopCriterion);
    currTnnT=tubalrank(tenT,0);
if currTnnT == preTnnT
    NOChange_counter=NOChange_counter +1;
else
    NOChange_counter = 0;
end
if (stopCriterion < tol) || (NOChange_counter >=1)
    break;
end             
    preTnnT = currTnnT; 
    disp(['#Iteration ' num2str(iter) ' trankT ' ...
        num2str(currTnnT) ...
        ' stopCriterion ' num2str(stopCriterion) ]);
end
