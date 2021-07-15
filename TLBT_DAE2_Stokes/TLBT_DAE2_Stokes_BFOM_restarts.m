function TLBT_DAE2_Stokes_BFOM_restarts()
%% Comparison of modified Arnoldi method and B(FOM)^2 with restarts
% for medium sized Stokes test system. 
% 
% Beware that generating this plot may take a while!
%
% Requires precomputed_EFt/Stokes_medium_EFt_1.mat!
% If the file is not available consult Stokes_compute_EFt.
%
% Requires forked bfomfom package [1],[2]! We introduce an additional  
% function handle to accept the projection avoiding updates of W in the  
% modified Arnoldi method (see projectionAvoidingUpdate()). 
%
% Please see LICENSE.md, MODIFICATIONS.md and README.md for bfomfom
% package in bfomfom_modified subfolder.
%  
% [1] A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods 
%     for functions of matrices, Electronic Transactions on Numerical 
%     Analysis, Vol. 47, pp. 100-126, 2017.
% [2] A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods 
%     for functions of matrices II: Modified block FOM, SIAM Journal on 
%     Matrix Analysis and Applications, 41(2), pp. 804--837, 2020.
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See LICENSE.md

%% Setup system

% These are the parameters for the medium sized test system.
m = 10;  
q = 10; 
nx = 35;
ny = 35;

% Note that B(FOM)^2 only accepts parameter t = 1!
% Further modifications to the bfomfom package are required
% to accept time parameters!
t = 1;
tol = 1e-9;

% Load model equations.
[eqn.E_,eqn.A_,eqn.Borig,eqn.Corig,nf] = stokes_ind2(m,q,nx,ny);

% Get reference solution. 
% Might take a long time for large system dimensions. One may use
% precomputed EFt (see Stokes_compute_EFt function).
fname = 'Stokes_medium_EFt_1.mat';
fname = sprintf('%s/precomputed EFt/%s', ...
                fileparts(mfilename('fullpath')),fname);
precomp = load(fname);
EFt_expm = precomp.EFt;
EFtB_expm = EFt_expm*eqn.Borig;

% You may compute EFt again using the following.
% EFt_expm = Stokes_compute_EFt(eqn,t);
% EFtB_expm = EFt_expm*eqn.Borig;

% Initialize.
full_d = nf;
len = 15;
step = ceil(full_d/len);
error_arnoldi = zeros(len,1);
error_bfomfom_2r = zeros(len,1);
error_bfomfom_3r = zeros(len,1);
error_bfomfom_4r = zeros(len,1);
cut1 = 0; cut2 = 0; cut3 = 0; cut4 = 0;

% Slowly increase Krylov subspace dimension and evaluate errors for
% different methods. Up to 4 restarts with B(FOM)^2.
for i = 1:len

    % Current Krylov subspace dimension.
    dim = i*step;

    % Arnoldi iteration as presented in Chapter 4.2.3.
    if cut1 == 0
        EFtB_arnoldi = TLBT_DAE2_Stokes_modArnoldi(eqn,t,dim,-1,'N',0);
        error = norm(EFtB_expm - EFtB_arnoldi);
        error_arnoldi(i) = error;
        if error < 10*tol      
            cut1 = i;
        end
    end

    % Lund BFOMFOM with 2 restarts.
    if cut2 == 0
        EFtB_bfomfom = bfomfom_stokes(eqn,dim,tol,2);
        error = norm(EFtB_expm - EFtB_bfomfom);
        error_bfomfom_2r(i) = error;
        if error < 10*tol      
            cut2 = i;
        end
    end

    % Lund BFOMFOM with 3 restarts.
    if cut3 == 0
        EFtB_bfomfom = bfomfom_stokes(eqn,dim,tol,3);
        error = norm(EFtB_expm - EFtB_bfomfom);
        error_bfomfom_3r(i) = error;
        if error < 10*tol       
            cut3 = i;
        end
    end

    %Lund BFOMFOM with 4 restarts.
    if cut4 == 0
        EFtB_bfomfom = bfomfom_stokes(eqn,dim,tol,4);
        error = norm(EFtB_expm - EFtB_bfomfom);
        error_bfomfom_4r(i) = error;
        if error < 10*tol      
            cut4 = i;
        end
    end
end

% Prepare data for plot.
if cut2 > 0; error_bfomfom_2r = [error_bfomfom_2r(1:cut2);zeros(len-cut2,1)]; end
if cut3 > 0; error_bfomfom_3r = [error_bfomfom_3r(1:cut3);zeros(len-cut3,1)]; end
if cut4 > 0; error_bfomfom_4r = [error_bfomfom_4r(1:cut4);zeros(len-cut4,1)]; end
figure('Name','Error EF(t)B approximates vs. Krylov subspace dimension');
dims = step*(1:len);

% Plot.
semilogy(dims,error_arnoldi,'LineWidth',1.5);
hold on
semilogy(dims,error_bfomfom_2r,'LineWidth',1.5);
semilogy(dims,error_bfomfom_3r,'LineWidth',1.5);
semilogy(dims,error_bfomfom_4r,'LineWidth',1.5);
hold off
lgd = legend('Modified Arnoldi','B(FOM)^2 with 2 restarts','B(FOM)^2 with 3 restarts','B(FOM)^2 with 4 restarts');
lgd.FontSize = 14;
title('Error EF(t)B approximates vs. Krylov subspace dimension');
xlabel('Krylov subspace dimension');
ylabel('error')

end

function EFtB = bfomfom_stokes(eqn,nf,tol,restarts)
%% Use B(FOM)^2 package to calculate EF(t)B with restarts in Arnoldi
% method. Note that the bfomfom files are modified to accept additional
% parameters.

fprintf('\nStarting Lund Arnoldi iteration for EF(t)B ...\n');

tStart = tic; 

n = size(eqn.E_,1);
nv = full(sum(diag(eqn.E_)));  
np = n - nv;
m = size(eqn.Borig,2);

% Saddle point system leading matrix.
eqn.SPM = [eqn.E_(1:nv,1:nv) eqn.A_(1:nv,nv+1:n); eqn.A_(1:nv,nv+1:n)' sparse(np,np)];

% Update matrix for arnoldi step.
eqn.OUM = eqn.A_(1:nv,1:nv);

% Set parameters for B(FOM)^2.
A = @projectionAvoidingUpdate;
B0 = eqn.SPM\[eqn.Borig(1:nv,:);sparse(np,m)];
B0 = B0(1:nv,:);

param_init = [];
param_init.SPM = eqn.SPM;
param_init.OUM = eqn.OUM;

param_init.function = 'exp';
param_init.verbose = 0;
param_init.conv_check = 'approx';
param_init.error_scale = 1;
param_init.inner_product = 'cl_full';  % choice of block inner product                                                         
param_init.restart_length = floor(nf/m); % number of basis vectors before restart
param_init.max_restarts = restarts;
param_init.tol_err = tol;

param = param_init_bfomfom(param_init); % initialize all other parameter fields
param = bfomfom(A,B0,param);

EFtB = [eqn.E_(1:nv,1:nv)*param.Fm;sparse(np,m)];

toc(tStart)

end

function W = projectionAvoidingUpdate(Vk,param)
%% Function handle for the updates of W in B(FOM)^2.

[nv,m] = size(Vk);
np = size(param.SPM,1) - nv;
Theta = param.SPM\[param.OUM*Vk;sparse(np,m)];
W = Theta;
W = W(1:nv,:);

end

