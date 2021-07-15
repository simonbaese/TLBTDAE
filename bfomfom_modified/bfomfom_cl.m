function param = bfomfom_cl(A,B,param)
% param = BFOMFOM_CL(A,B,param) computes the restarted classical B(FOM)^2
% approximation to the matrix function product f(A)B.
%
% This file is a part of the B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices, Electronic Transactions on Numerical Analysis,
% Vol. 44, pp. 100-126, 2017,
%
% and of the modified B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices II: Modified block FOM. Technical Report 19-04-10,
% Department of Mathematics, Temple University. 2019.

%%
mfilename = sprintf('bfomfom_%s',param.inner_product);
fprintf('%s\n',mfilename)
global delta;
global T; T = 0;                                                            % pre-set walltime
tic;
param.first_cycle = true;                                                   % store Beta
param = krylov_basis(A,B,param);
Acount = param.Acount;

% Extract quantities needed for cospatial function. We have to break
% dependence on param; otherwise the in-line functions will store all
% previous param, which includes all previous basis vectors.
Beta = param.Beta;
E1 = param.E1;
Hnext = param.Hnext;
Hshort = param.Hshort;

Xim = @(t) (Hshort + t*eye(size(Hshort)))\E1;
if strcmp(param.function,'exp')
    param.ritz = eig(Hshort);                                               % store Ritz values 
end
Cm = @(t) cospatial(Xim,Hnext,t)*Beta;

switch param.function                                                       % compute initial approximation
    case 'invSqrt'
        delta = 1;
        param.alpha = .5;
        param.Fm = param.Vshort*(sqrtm(Hshort)\E1)*Beta;
        param.quad = @(N,param) quad_invAlpha(N,.5);

    case 'invAlpha'
        delta = 1;
        param.Fm = param.Vshort*(Hshort^param.alpha\E1)*Beta;
        param.quad = @(N,param) quad_invAlpha(N,param.alpha);

    case 'logShiftInv'
        param.Fm = param.Vshort*(logm(Hshort+eye(size(Hshort))))*...
                (Hshort\E1)*Beta;
        param.quad = @(N,param) quad_logShiftInv(N);
        
    case 'exp'
        param.Fm = param.Vshort*expm(Hshort)*E1*Beta;                    
        param.quad = @(N,param) quad_exp(N,param);
        
    case 'user'
        param.Fm = param.Vshort*userf(Hshort)*E1*Beta;                
        param.quad = @(N,param) quad_userf(N,param);
end

T = T + toc;
if isfield(param,'exact')
    exact_flag = true;
    param.exact_error = block_norm(param.Fm - param.exact,param,A)/...
        param.error_scale;                                                  % compute exact error
    if param.exact_error < param.tol_err
        param.flag = 1;
        print_flag(param.flag);
        return
    end
else
    exact_flag = false;                                                     % store flag to save some time
end

%% Restart
param.Nquad_per_cycle = [];                                                 % pre-set storage for quadrature node count
param.approx_error = [];                                                    % pre-set storage for approximate error
param.first_cycle = false;                                                  % for subsequent cycles, do not redefine Beta
param.flag = 0;                                                             % initialize convergence flag
param.condH = [];
for k = 2:param.max_restarts
    T(k) = 0; tic;
    param = krylov_basis(A,param.Vnext,param);
    Acount = Acount + param.Acount;
    
    % Extract quantities needed for cospatial function update
    E1 = param.E1;
    Hnext = param.Hnext;
    Hshort = param.Hshort;
    
    Xim = @(t) (Hshort + t*eye(size(Hshort)))\E1;
    if strcmp(param.function,'exp')
        param.ritz = [param.ritz; eig(Hshort)];                             % store Ritz values
    end
    
    param = mat_err_func(param,Cm,Xim,k);                                   % error function
    if param.flag == 4
        param.Nquad_per_cycle(k-1) = param.Nquadmax;
        break
    end
    param.Nquad_per_cycle(k-1) = param.Nquad2;

    err_apx = param.Vshort*param.Delta;                                     % compute approximate error from previous cycle
    param.approx_error(k-1) = block_norm(err_apx,param,A)/...
        param.error_scale;                                                  % store approximate error from previous cycle
    
    T(k) = T(k) + toc;
    switch param.conv_check
        case 'exact'
            err_check = param.exact_error(k-1);
            if k > 2
                err_check_past = param.exact_error(1);
            end
        case 'approx'
            err_check = param.approx_error(k-1);
            if k > 2
                err_check_past = param.approx_error(1);
            end
    end
    if param.verbose
        fprintf('     %s: cycle %g, %s error = %1.4g\n',...
            mfilename,k-1,param.conv_check,err_check);
    end
    tic;
    
    if err_check < param.tol_err                                            % check for convergence
        param.flag = 1;
        break
    elseif k > 2 && param.halt_if_diverge                                     % check for stagnation
        if err_check >= err_check_past
            
            T(k) = T(k) + toc;
            if param.verbose
                fprintf('     %s: stagnation at cycle %g\n',mfilename,k-1);
            end
            tic;
            
            param.flag = 2;                                                 % stagnation
            break
        end
    end
    param.Fm = param.Fm + err_apx;                                          % update solution approximation
    
    T(k) = T(k) + toc;
    if exact_flag
        param.exact_error(k) = block_norm(param.Fm - param.exact,param,A)/...
            param.error_scale;
    end
    tic;
    
    Cm = @(t) cospatial(Xim,Hnext,t)*Cm(t);                    % next cospatial factor
    
    if param.reduce_Nquad                                                   % reduce number of quadrature nodes
        if ~param.refined && param.Nquad1 > 2
            param.Nquad2 = param.Nquad1;
            param.Nquad1 = round(param.Nquad1/sqrt(2));
            if mod(param.Nquad1,2) == 1
                param.Nquad1 = param.Nquad1-1;
            end
        end
    end
end
param.condH(end+1) = cond(Hshort);
if param.max_restarts == 1                                                  % when max_restarts = 1
    k = 1;
end
param.Acount = Acount/k;

T(k) = T(k) + toc;
param.time = T;

if k == param.max_restarts && param.flag == 0
    param.flag = 3;                                                         % maximum number of cycles reached
end

print_flag(param.flag);
end

%% Auxiliary functions
function out = cospatial(Xim,Hnext,t)
% Computes the next cospatial factor
Xi = Xim(t);
s = size(Hnext,1);
out = -Hnext*Xi(end-s+1:end,:);
end

function print_flag(flag)
switch flag
    case 0
        if param.max_restarts == 1
            fprintf('Only one cycle, huh?\n')
        else
            error('Flag unchanged.  Check code for bugs.\n')
        end
    case 1
        fprintf('Method converged to desired tolerance.\n')
    case 2
        fprintf('Method stagnated before reaching desired tolerance.\n')
    case 3
        fprintf('Maximum number of cycles reached before convergence.\n')
    case 4
        fprintf('Maximum number of quadrature nodes reached before convergence.\n')
end
end
