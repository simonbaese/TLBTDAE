function hsv = TLBT_DAE3_Stykel_mor(problem,tlbt_opts)
%% Time-limited balanced truncation for mechanical multibody systems 
% 
% See Algorithm 1 in [1]!
%
% Uses model data from index-3 Damped mass-spring system with a holonomic 
% constraint in M.E.S.S. toolbox (compare [2], see also [3]).
%
% Inputs:
% problem               problem size: 'tiny','small','medium' or 'large'
% tlbt_opts             struct with options
%          .tl          1/0 - time-limited/traditional balanced truncation
%          .tf          final time
%	       .adiplot     1/0 - plot ADI iteration residual norm
%          .analyze     1/0 - plot Hankel singular values and sigma plots
%
% Outputs:
% hsv           Hankel singular values
%
% [1] S. M. Bäse. "Time-Limited Balanced Truncation Model Order Reduction 
%     for Descriptor Systems". Master thesis. Technische Universität 
%     Berlin: Institut für Mathematik, 2021.
% [2] V. Mehrmann, T. Stykel. Balanced truncation model reduction for 
%     large-scale systems in descriptor form, in Dimension Reduction of 
%     Large-Scale Systems, P. Benner, V. Mehrmann, and D. Sorensen (eds.),
%     Springer-Verlag, Berlin/Heidelberg, 2005, pp.83--115
% [3] J. Saak, M. Voigt, Model reduction of constrained mechanical systems
%     in M-M.E.S.S., IFAC-PapersOnLine 9th Vienna International Conference
%     on Mathematical Modelling MATHMOD 2018, Vienna, Austria
%     February 2018 51 (2) (2018)
%
% Author: Simon Bäse (2021)
% Partly forked from DAE3_SO demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD


%% Check tlbt_opts parameter
tlbt_opts = TLBT_check_opts(tlbt_opts);

%% Problem data
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

eqn.Borig = eqn.B;
eqn.Corig = eqn.C;

m = size(eqn.B,2);
q = size(eqn.C,1);
n = 2*nv + np;
nf = n - 3;

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
    atol = 1e-13; 
    tf = tlbt_opts.tf;
    
    % Calculate EFtB and projected initial B as ThetaB.
    [EFtB,ThetaB,~,~] = TLBT_DAE3_Stykel_modArnoldi_N(eqn,tf,nf,atol);

    % Calculate EFtC and projected initial C as ThetaC.
    [EFtC,ThetaC,~,~] = TLBT_DAE3_Stykel_modArnoldi_T(eqn,tf,nf,atol);   

end

%% Setup parameters for ADI iteration

fprintf('\nStarting ADI iterations ...\n\n'); 
                
% Set operation manager for the Gramian computations.
oper = operatormanager('dae_3_so');

% ADI tolerances and maximum iteration number.
opts.adi.maxiter = 300;
opts.adi.res_tol = 1e-10;
opts.adi.rel_diff_tol = 1e-11;
opts.norm = 'fro';
opts.shifts.method = 'projection';
opts.shifts.num_desired = 10;
if ~tlbt_opts.tl; opts.shifts.p = mess_para(eqn, opts, oper); end

%% Compute controllability Gramian Cholesky factor

fprintf('Determine controllability Gramian Cholesky factor ...\n');
eqn.type = 'N';

tADI = tic;

% Use LDL^T ADI iteration in time-limited case.
if tlbt_opts.tl
    % Set new eqn.B for t_s = 0 and t_f = tf.
    % Note that mess_lradi solves A*X*E' + E*X*A' + G*S*G' = 0.
    eqn.G = [ThetaB(1:2*nv,:), EFtB(1:2*nv,:)]; 
    opts.LDL_T = 1;
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
    eqn.G = [ThetaC(1:2*nv,:), EFtC(1:2*nv,:)];
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
opts.srm.tol = 1e-6;
opts.srm.max_ord = 200;
opts.srm.info = 2;

% Balanced truncation square root method.
[TL,TR,hsv] = mess_square_root_method(eqn,opts,oper,outB.Z,outC.Z);

% Compute ROM matrices.
ROM.A = TL' * oper.mul_A(eqn, opts, 'N', TR, 'N');
ROM.B = TL' * eqn.Borig;
ROM.C = eqn.Corig * TR;
ROM.E = eye(size(ROM.A));

% Stop timer for model order recution.
fprintf('Total MOR time: %.2fs\n\n', toc(tMOR));

% Reset input and output matrices.
eqn.B = eqn.Borig;
eqn.C = eqn.Corig;

%%  Plot normalized residual norm of ADI iterations.

if tlbt_opts.adiplot; TLBT_adi_plot(outB,outC); end

%% Output sigma plot and HSV plot

if tlbt_opts.analyze 
    
    opts.sigma.nsample = 200;
    opts.sigma.info = 2;
    opts.sigma.fmin = -4;
    opts.sigma.fmax = 4;
    
    NG = sparse(np,nv);
    NS = sparse(np,np);
    eqnu.M_ = [eqn.M_ NG'; NG NS];
    eqnu.E_ = [-eqn.E_ NG'; NG NS];
    eqnu.K_ = [-eqn.K_ -eqn.G_';-eqn.G_ NS];
    eqnu.C =  [zeros(size(eqn.Corig,1),np+nv), eqn.C(:,1:nv),...
        zeros(size(eqn.Corig,1),np)];
    eqnu.B =  [eqn.Borig(nv+1:2*nv,:);zeros(2*np+nv,size(eqn.Borig,2))];
    eqnu.haveE = 1;

    operu = operatormanager('so_1');

    mess_sigma_plot(eqnu, opts, operu, ROM);

    figure;
    mhsv = max(hsv);
    hsv = hsv./mhsv;
    semilogy(hsv);
    title('Computed Hankel singular values');
    xlabel('index');
    ylabel('magnitude');
end


end