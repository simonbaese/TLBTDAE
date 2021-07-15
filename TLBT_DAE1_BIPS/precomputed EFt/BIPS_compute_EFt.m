function BIPS_compute_EFt()
%% Compute EFt for BIPS model
%
% This function is used to precompute the EFt term with the full system
% matrices. Computation time can be very long! Our system at 4.4 GHz 
% took around 52 minutes to compute three time points.
%
% Author: Simon Bäse (2021)
% Partly forked from DAE1 BIPS demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD

%% Total CPU time
tStart = tic;

%% Turn off  close to singular warnings
% (this model is really badly conditioned)
orig_warnstate = warning('OFF','MATLAB:nearlySingularMatrix');

%% Read problem data
% from https://sites.google.com/site/rommes/software
fname  = sprintf('%s/../../../models/BIPS/bips98_606.mat',...
    fileparts(mfilename('fullpath')));
Bips = load(fname);

% Initialize problem data.
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p;np];
eqn.A_ = Bips.A(pp, pp)-0.05*Bips.E(pp,pp);
eqn.E_ = Bips.E(pp, pp);
nf = length(p);
nin = length(np);

% Use block matrices.
E11 = eqn.E_(1:nf,1:nf);
A12 = eqn.A_(1:nf,nf+1:end);
A21 = eqn.A_(nf+1:end,1:nf);
A22 = eqn.A_(nf+1:end,nf+1:end);

%% Precalculation.
fprintf('\nPrecalculation ...\n');
SC = A12/A22;
PL = [speye(nf) -SC; sparse(nin,nf) sparse(nin,nin)];
Omega = [E11 + SC*A21 A12; A21 A22];
AOmInv = eqn.A_/Omega;

%% Initialize time snapshots.
tf = [1,3,10];
appendix = ["1","3","10"];
cputime = zeros(3,1);

%% Calculate full MATLAB expm for BIPS
fprintf('\nCalculating full MATLAB expm ...\n');

for i = 1:3
    tic;
    fprintf('\nRun BIPS_EFt_%s ...\n',appendix(i));
    Exp = expm(tf(i)*AOmInv);
    EFt = PL*Exp;
    fname = sprintf('BIPSEFt%s.mat',appendix(i));
    fname = sprintf('%s/%s',fileparts(mfilename('fullpath')),fname);
    save(fname,'EFt');
    cputime(i) = toc;
    disp(cputime(i));
end

%% Output CPU times
fprintf('\nIndividual CPU times ...\n');
disp(cputime);

fprintf('\nTotal CPU time ...\n');
disp(toc(tStart));

%% reset warning state
warning(orig_warnstate);

% system('shutdown -s');