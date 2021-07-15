function EFt = Stokes_compute_EFt(eqn,t)
%% Compute EFt for Stokes-like model
%
% This function is used to compute the EFt term with the full system
% matrices. Uncomment relevant sections to use precomputed EFt matrices.
%
% We were able to compute EFt for the medium system. With the large system
% our computer ran out of memory.
%
% Input:
% eqn           struct that contains system matrices as defined in
%               TLBT_DAE2_Stokes_mor
% t             time point
%
% Outpt:
% EFt           E times fundamental solution matrix F at t
%
% Author: Simon Bäse (2021)
%
% See License.MD

%% Total CPU time
tStart = tic;

%% Calculate full MATLAB expm for Stokes
fprintf('\nCalculating full MATLAB expm for EF(t)B ...\n');

n = size(eqn.E_,1);
nv = full(sum(diag(eqn.E_)));  
np = n - nv;

E11 = eqn.E_(1:nv,1:nv);
A11 = eqn.A_(1:nv,1:nv);
A12 = eqn.A_(1:nv,nv+1:n);

Pi1 = speye(nv,nv) - A12/(A12'/E11*A12)*A12'/E11;
Pi2 = -Pi1*A11/E11*A12/(A12'/E11*A12);
Pl = [ Pi1 Pi2; sparse(np,nv) sparse(np,np) ];

EFt = Pl*expm(t*eqn.A_/(Pl*eqn.E_ + (speye(n,n) - Pl)*eqn.A_));

%% Output CPU times
fprintf('\nTotal CPU time ...\n');
tFinish = toc(tStart);
fprintf('Expm full calculation time: %.2fs\n\n', tFinish);

%% Save a computed EFt.
% file.EFt = EFt;
% file.tFinish = tFinish;
% fname = 'Stokes_medium_EFt_1.mat';
% fname = sprintf('%s/%s',fileparts(mfilename('fullpath')),fname);
% save(fname,'-struct','file');
    
% system('shutdown -s');

%% Load a precomputed EFt.
% fname = 'Stokes_medium_EFt_1.mat';
% fname = sprintf('%s/%s', ...
%                 fileparts(mfilename('fullpath')),fname);
% precomp = load(fname);
% EFt = precomp.EFt;

end
