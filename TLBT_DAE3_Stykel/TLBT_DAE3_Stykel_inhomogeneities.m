function TLBT_DAE3_Stykel_inhomogeneities(problem,t)
%% Validates computation of inhomogeneities 
% for time-limited Lyapunov equations for Stykel model.
% We compare the modified Arnoldi methods with a simulation approach. We
% use EF(t)B to compute the reference solutions.
%
% Inputs:
% problem               model: 'stykel_small' or 'stykel_large'
% t                     time point
%
% Author: Simon Bäse (2021)
% Partly forked from DAE3_SO demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD

%% Read problem data
if strcmp(problem,'stykel_small')
    sys = load(sprintf('%s/../../models/ms_ind3_by_t_stykel/g600.mat',...
        fileparts(mfilename('fullpath'))));
else
    sys = load(sprintf('%s/../../models/ms_ind3_by_t_stykel/g6000.mat',...
        fileparts(mfilename('fullpath'))));
end

eqn.M_ = sys.M;
eqn.E_ = sys.D;
eqn.K_ = sys.K;
eqn.G_ = sys.G;
eqn.haveE = 1;
eqn.alpha = -0.02;
nv = size(eqn.M_,1);
np = size(eqn.G_,1);
eqn.B = full(sys.B(1:2*nv,:));
eqn.C = full(sys.C(:,1:2*nv));
clear E A B C M D K G;   

m = size(eqn.B,2);
q = size(eqn.C,1);
n = 2*nv + np;
nf = n - 3;

% Initialize.
tol = 1e-13;
        
%% Full expm used as reference.
fprintf('\nCalculating full MATLAB expm for EF(t)B ...\n');
tic;
EFt = sys.Pl*expm(t*sys.A/(sys.Pl*sys.E + (speye(n) - sys.Pl)*sys.A));
EFtB_expm = EFt*[eqn.B; sparse(np,m)];
fprintf('Full expm time EF(t)B: %.2fs\n\n',toc);
    
fprintf('\nCalculating full MATLAB expm for CF(t)E ...\n');
tic;
FtE = expm(t*speye(n)/(sys.E*sys.Pr + sys.A*(speye(n) - sys.Pr))*sys.A)*sys.Pr;
CFtE_expm = [eqn.C sparse(q,np)]*FtE;
fprintf('Full expm time CF(t)E: %.2fs\n\n',toc);   

%% Modified Arnoldi iterations
EFtB_arnoldi = TLBT_DAE3_Stykel_modArnoldi_N(eqn,t,nf,tol);
EFtCT_arnoldi = TLBT_DAE3_Stykel_modArnoldi_T(eqn,t,nf,tol);
    
%% Output difference.
fprintf('\nError with modified Arnoldi EF(t)B: %.2e\n', norm(EFtB_expm - EFtB_arnoldi));
fprintf('\nError with modified Arnoldi CF(t)E: %.2e\n', norm(CFtE_expm - EFtCT_arnoldi'));
         
end