function hsv = TLBT_DAE1_BIPS_mor(problem,tlbt_opts)
%% Time-limited balanced truncation for semi-explicit systems of index 1
% 
% See Algorithm 1 in [1]!
% 
% Uses model data from index-1 system BIPS98_606 from 
% https://sites.google.com/site/rommes/software (compare [2]).
%
% Inputs:
% problem               problem size: 'tiny','small','medium' or 'large'
% tlbt_opts             struct with options
%          .tl          1/0 - time-limited/traditional balanced truncation
%          .tf          final time
%	       .adiplot     1/0 - plot ADI iteration residual norm
%          .analyze     1/0 - plot Hankel singular values and sigma plots
%          .simulate    1/0 - simulate system and calculate errors
%          .horizon 	time horizon for simulation, greater than tf
%          .input       input function: 'impulse' or 'step'
%
% Outputs:
% hsv           Hankel singular values
%
% References:
% [1] S. M. Bäse. "Time-Limited Balanced Truncation Model Order Reduction 
%     for Descriptor Systems". Master thesis. Technische Universität 
%     Berlin: Institut für Mathematik, 2021.
% [2] F. Freitas, J. Rommes, N. Martins, Gramian-based reduction method
%     applied to large sparse power system descriptor models, IEEE Trans.
%     Power Syst. 23 (3) (2008) 1258-1270. doi:10.1109/TPWRS.2008.926693. 
%
% Author: Simon Bäse (2021)
% Partly forked from DAE1 BIPS demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD


%% Check tlbt_opts parameter
tlbt_opts = TLBT_check_opts(tlbt_opts);

%% Read problem data
% from https://sites.google.com/site/rommes/software
% Note that we here use the alpha shift suggested by Rommes and coauthors
% [1], but do not perform the shift back.
switch lower(problem)
    case 'bips07_3078'
        fname  = sprintf('%s/models/bips07_3078.mat',...
            fileparts(mfilename('fullpath')));
        alpha = 0.08;
    otherwise
        fname  = sprintf('%s/../../models/BIPS/bips98_606.mat',...
            fileparts(mfilename('fullpath')));
        alpha = 0.05;
end
Bips = load(fname);
p = find(diag(Bips.E));
np = find(diag(Bips.E) == 0);
pp = [p;np];
eqn.A_ = Bips.A(pp,pp) - alpha*Bips.E(pp,pp);
eqn.E_ = Bips.E(pp,pp);
eqn.Borig = Bips.b(pp,:);
eqn.Corig = Bips.c(:,pp);
% eqn.D = Bips.d; % All entries are zero.
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;
eqn.st = length(p);
eqn.haveE = 1;
m = size(eqn.B,2);
q = size(eqn.C,1);
nf = length(p);
nin = length(np);

% Start timer for complete model order reduction.
tMOR = tic;

if tlbt_opts.tl
    fprintf(['\nPerforming MOR by time-limited balanced truncation ', ...
        'with end time t = %g ...\n'],tlbt_opts.tf); 
else
    fprintf('\nPerforming MOR by balanced truncation ...\n'); 
end

%% Turn off  close to singular warnings
% (this model is really badly conditioned)
orig_warnstate = warning('OFF','MATLAB:nearlySingularMatrix');

%% Calculate inhomogeneities for time-limited BT
if tlbt_opts.tl    
    
    % Tolerance and end-time for Arnoldi iteration. 
    atol = 1e-12; 
    tf = tlbt_opts.tf;
    
    % Calculate EFtB and projected initial B as ThetaB.
    [EFtB,ThetaB,~,~] = TLBT_DAE1_BIPS_modArnoldi(eqn,tf,nf,nin,atol,'N');
        
    % Calculate EFtC and projected initial C as ThetaC.
    [EFtC,ThetaC,~,~] = TLBT_DAE1_BIPS_modArnoldi(eqn,tf,nf,nin,atol,'T');  
end

%% Setup parameters for ADI iteration

fprintf('\nStarting ADI iterations ...\n\n'); 
                
% Set operation manager for the Gramian computations.
oper = operatormanager('dae_1');

% ADI tolerances and maximum iteration number.
opts.adi.maxiter = 350;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.shifts.info = 0;
opts.norm = 'fro';
opts.shifts.method = 'projection';
opts.shifts.num_desired = 30; 

%% Compute controllability Gramian Cholesky factor

fprintf('Determine controllability Gramian Cholesky factor ...\n');
eqn.type = 'N';

tADI = tic;

% Use LDL^T ADI iteration in time-limited case.
if tlbt_opts.tl
    opts.LDL_T = 1;
    % Set new eqn.B for t_s = 0 and t_f = tf.
    % Note that mess_lradi solves A*X*E' + E*X*A' + G*S*G' = 0.
    eqn.G = [ThetaB(1:nf,:), EFtB(1:nf,:)]; 
    eqn.S = [speye(m), sparse(m,m); sparse(m,m), -speye(m)];
end

% Low-rank ADI iteration for B.
tic;
outB = mess_lradi(eqn,opts,oper);
fprintf('ADI iteration for outB time: %.2fs\n', toc);

% Calculate Cholesky factors after LDL^T low-rank ADI iteration.
if tlbt_opts.tl   
    tic;    
    outB.Z = TLBT_definite_approximation(outB.Z,outB.D,eps^2);
    fprintf('Cholesky factor approximation for outB time: %.2fs\n', toc);
end

fprintf('Size outB.Z:');
disp(size(outB.Z));

%% Compute observability Gramian Cholesky factor

fprintf('Determine observability Gramian Cholesky factor ...\n');
eqn.type = 'T';

% Use LDL^T ADI iteration in time-limited case.
if tlbt_opts.tl
    % Set new eqn.C for t_s = 0 and t_f = tf.
    eqn.G = [ThetaC(1:nf,:), EFtC(1:nf,:)];
    eqn.S = [speye(q), sparse(q,q); sparse(q,q), -speye(q)];
end

% Low-rank ADI iteration for C.
tic;
outC = mess_lradi(eqn, opts, oper);
fprintf('ADI iteration for outC time: %.2fs\n', toc);

% Calculate Cholesky factors after LDL^T low-rank ADI iteration.
if tlbt_opts.tl
    tic;
    outC.Z = TLBT_definite_approximation(outC.Z,outC.D,eps^2);       
    fprintf('Cholesky factor approximation for outC time: %.2fs\n', toc);
end

fprintf('Size outC.Z:');
disp(size(outC.Z));

fprintf('Total ADI time: %.2fs\n\n', toc(tADI));

%% Compute reduced system matrices

% Balanced truncation tolerance and maximum order for the ROM.
opts.srm.tol = 1e-3;
opts.srm.max_ord = 250;
opts.srm.info = 2;

% Balanced truncation square root method.
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

% Compute ROM matrices.
B1 = TL'*(eqn.A_(1:eqn.st,1:eqn.st))*TR;
B2 = TL'*(eqn.A_(1:eqn.st,eqn.st+1:end));
A1 = eqn.A_(eqn.st+1:end,1:eqn.st)*TR;

ROM.A = B1 - B2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\A1);
ROM.B = TL'*eqn.Borig(1:eqn.st,:) - ...
     B2*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\eqn.Borig(eqn.st+1:end,:));
ROM.C = eqn.Corig(:,1:eqn.st)*TR - ...
     eqn.Corig(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\A1);
ROM.D = -eqn.Corig(:,eqn.st+1:end)*(eqn.A_(eqn.st+1:end,eqn.st+1:end)\...
    eqn.Borig(eqn.st+1:end,:)); % All entries are zero.
ROM.E = eye(size(ROM.A));

% Stop timer for model order recution.
fprintf('Total MOR time: %.2fs\n\n', toc(tMOR));

% Reset input and output matrices.
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;

% While the Gramians are computed on the hidden manifold, we need to do the
% frequency domain computations without (implicitly) using the Schur
% complement (due to the construction of the function handles).
oper = operatormanager('default');

%%  Plot normalized residual norm of ADI iterations.
if tlbt_opts.adiplot; TLBT_adi_plot(outB,outC); end

%% Output sigma plot and HSV plot

if tlbt_opts.analyze 
    
    opts.sigma.info = 2;
    opts.sigma.fmin = -3;
    opts.sigma.fmax = 4;
    mess_sigma_plot(eqn, opts, oper, ROM);
% 
%     figure;
%     mhsv = max(hsv);
%     hsv = hsv./mhsv;
%     semilogy(hsv);
%     title('Computed Hankel singular values');
%     xlabel('index');
%     ylabel('magnitude');
end

%% Evaluate ROM quality

if tlbt_opts.simulate 
        
    % Start timer for evaluation and simulation.
    tSIM = tic;
    
    fprintf(['\nComputing average & maximum output error by ', ...
            'simulating system with the implicit Euler method ...\n']);
    h = 1e-3; % step size
    tmax = tlbt_opts.horizon;
    input = tlbt_opts.input;    
    tf = tlbt_opts.tf;
    TLBT_simulate(eqn,ROM.A,ROM.B,ROM.C,h,tf,tmax,input);

    % Stop timer for evaluation and simulation.
    fprintf('Total evaluation and simulation time: %.2fs\n\n', toc(tSIM));
end

%% Reset warning state
warning(orig_warnstate);

end