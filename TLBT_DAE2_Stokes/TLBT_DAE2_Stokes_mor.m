function hsv = TLBT_DAE2_Stokes_mor(problem,tlbt_opts)
%% Time-limited balanced truncation for Stokes-like systems of index 2
%
% See Algorithm 1 in [1]!
%
% Uses model data from index-2 Semidiscretized 2D Stokes equation in 
% M.E.S.S. toolbox (compare [2]).
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
% [1] S. M. Bäse. "Time-Limited Balanced Truncation Model Order Reduction 
%     for Descriptor Systems". Master thesis. Technische Universität 
%     Berlin: Institut für Mathematik, 2021.
% [2] M.Schmidt. Systematic discretization of input/output maps and
%     other contributions to the control of distributed parameter systems.
%     Ph.D. thesis, TU Berlin, 2007.
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD

%% Check tlbt_opts parameter
tlbt_opts = TLBT_check_opts(tlbt_opts);

%% Read problem data
if strcmp(problem,'tiny')
    m = 10; % number of inputs
    q = 10; % number of outputs
    nx = 15; % number of grid points in the x-direction
    ny = 15; % number of grid points in the y-direction
elseif strcmp(problem,'small')
    m = 10; q = 10; nx = 25; ny = 25;
elseif strcmp(problem,'medium')
    m = 10; q = 10; nx = 35; ny = 35; 
elseif strcmp(problem,'large')
    m = 5; q = 5; nx = 81; ny = 81;
else
    m = 10; q = 10; nx = 35; ny = 35;
end

[eqn.E_,eqn.A_,eqn.Borig,eqn.Corig,nf] = stokes_ind2(m,q,nx,ny);
eqn.haveE = 1;
nv = full(sum(diag(eqn.E_)));
eqn.st = nv;
eqn.B = eqn.Borig(1:nv,:);
eqn.C = eqn.Corig(:,1:nv);

% Start timer for complete model order reduction.
tMOR = tic;

if tlbt_opts.tl
    fprintf(['\nPerforming MOR by time-limited balanced truncation ', ...
        'with end time t = %g ...\n'],tlbt_opts.tf); 
else
    fprintf('\nPerforming MOR by balanced truncation ...\n'); 
end

%% Calculate inhomogeneities for time-limited BT
if tlbt_opts.tl   
    % Tolerance and end-time for Arnoldi iteration. 
    atol = 1e-9; 
    tf = tlbt_opts.tf;
    
    % Calculate EFtB and projected initial B as ThetaB.
    [EFtB,ThetaB,~,~] = TLBT_DAE2_Stokes_modArnoldi(eqn,tf,nf,atol,'N',0);
    
    % Calculate EFtC and projected initial C as ThetaC.
    [EFtC,ThetaC,~,~] = TLBT_DAE2_Stokes_modArnoldi(eqn,tf,nf,atol,'T',0);
end

%% Setup parameters for ADI iteration

fprintf('\nStarting ADI iterations ...\n\n'); 

% Set operation manager for the Gramian computations.
oper = operatormanager('dae_2');

% ADI tolerance and maximum iteration number
opts.adi.maxiter = 350;
opts.adi.res_tol = 1e-8;
opts.adi.rel_diff_tol = 1e-16;
opts.adi.info = 0;
opts.shifts.info = 0;
opts.norm = 'fro';
opts.shifts.method = 'projection';
opts.shifts.num_desired = 6;

%% Compute controllability Gramian Cholesky factor
fprintf('Determine controllability Gramian Cholesky factor ...\n');
eqn.type = 'N';

% Start timer for ADI iteration.
tADI = tic;

% Use LDL^T ADI iteration in time-limited case.
if tlbt_opts.tl
    opts.LDL_T = 1;
    % Set new eqn.B for t_s = 0 and t_f = tf.
    % Note that mess_lradi solves A*X*E' + E*X*A' + G*S*G' = 0.
    eqn.G = [ThetaB(1:nv,:), EFtB(1:nv,:)];    
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
    eqn.G = [ThetaC(1:nv,:), EFtC(1:nv,:)];
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

% Stop timer for evaluation and simulation.
fprintf('Total ADI time: %.2fs\n\n', toc(tADI));

%% Compute reduced system matrices

% Balanced truncation tolerance and maximum order for the ROM.
opts.srm.tol = 1e-5;
opts.srm.max_ord = 250;
opts.srm.info = 2;

% Balanced truncation square root method.
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

% Compute ROM matrices.
ROM.A = TL'*(eqn.A_(1:nv,1:nv)*TR);
ROM.B = TL'*eqn.Borig(1:nv,:);
ROM.C = eqn.Corig(:,1:nv)*TR;

% Stop timer for model order recution.
fprintf('Total MOR time: %.2fs\n\n', toc(tMOR));

% Reset input and output matrices.
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;

%% Evaluate the ROM quality
% while the Gramians are computed exploiting the DAE structure, due to the
% construction of the function handles we can not do so for the transfer
% function. Therfore we need to extend the matrices B and C and call the
% 'default' usfs for unstructured computation:
oper = operatormanager('default');

%%  Plot normalized residual norm of ADI iterations.
if tlbt_opts.adiplot; TLBT_adi_plot(outB,outC); end

%% Output sigma plot and HSV plot

if tlbt_opts.analyze 
  
    opts.sigma.info = 2;
    opts.sigma.fmin = -3;
    opts.sigma.fmax = 4;

    mess_sigma_plot(eqn, opts, oper, ROM);
    
    figure;
    mhsv = max(hsv);
    hsv = hsv./mhsv;
    semilogy(hsv);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
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

end
