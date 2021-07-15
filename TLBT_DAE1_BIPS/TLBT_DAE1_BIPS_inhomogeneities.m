function TLBT_DAE1_BIPS_inhomogeneities()
%% Validates computation of inhomogeneities 
%
% Compares modified Arnoldi method and MATLAB expm with full system 
% matrices to compute inhomogeneities at different time points for
% index-1 system BIPS98_606 from 
% https://sites.google.com/site/rommes/software (compare [1]).
%
% Requires precomputed BIPS_EFt_*.mat files! Please run BIPS_compute_EFt()
% if not available.
%
% References:
% [1] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%     applied to large sparse power system descriptor models, IEEE Trans.
%     Power Syst. 23 (3) (2008) 1258-1270. 
%
% Author: Simon Bäse (2021)
% Partly forked from DAE1 BIPS demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See LICENSE.md

%% Read problem data
fname  = sprintf('%s/../../models/BIPS/bips98_606.mat',...
    fileparts(mfilename('fullpath')));
Bips = load(fname);
% from https://sites.google.com/site/rommes/software
% Note that we here use the alpha shift suggested by Rommes and coauthors
% [1], but do not perform the shift back.
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p;np];
eqn.A_ = Bips.A(pp, pp) - 0.05*Bips.E(pp,pp);
eqn.E_ = Bips.E(pp, pp);
eqn.B = Bips.b(pp, :);
eqn.C = Bips.c( : , pp);
eqn.st = length(p);
eqn.haveE = 1;
nf = length(p);
nin = length(np);

%% Turn off  close to singular warnings
% (this model is really badly conditioned)
orig_warnstate = warning('OFF','MATLAB:nearlySingularMatrix');

%% Compare modified Arnoldi to full MATLAB expm.

% Initialize tolerance and time points.
tol = 1e-12; 
tf = [1,3,10];
appendix = ["1","3","10"];

% Iterate different sample time points.
for i = 1:3 
    
    % Calculate EFtB with modified Arnoldi for semi-explicit systems
    % of index 1.
    [EFtB,~,~,~] = TLBT_DAE1_BIPS_modArnoldi(eqn,tf(i),nf,nin,tol,'N');
    
    % Load EFt precomputed with full MATLAB expm.
    fname = sprintf('BIPS_EFt_%s.mat',appendix(i));
    fprintf(['\nLoading precomputed file ',fname,' and compare ...\n']);
    fname = sprintf('%s/precomputed EFt/%s', ...
                    fileparts(mfilename('fullpath')),fname);
    precomp = load(fname);
    
    % Calculate EFtB with precomputed EFt.
    expmB = precomp.EFt*eqn.B;
    
    % Output some sample values.
    fprintf('\nSome sample values: \n');
    fprintf('%8.7f %8.7f\n', [full(expmB(1:5,1)),full(EFtB(1:5,1))]');
    
    % Output difference.
    fprintf('\nDifference in B: %.2e\n', norm(full(expmB) - full(EFtB)));
end

%% Reset warning state
warning(orig_warnstate);