function TLBT_DAE2_Stokes_inhomogeneities(problem)
%% Validates computation of inhomogeneities 
% for time-limited Lyapunov equations for Stokes-like systems.
% We compare the modified Arnoldi methods with a simulation approach. We
% use EF(t)B to compute the reference solutions.
% 
% Inputs:
% problem               problem size: 'tiny','small','medium' or 'large'
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See License.MD

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

% Load model equations.
[eqn.E_,eqn.A_,eqn.Borig,eqn.Corig,nf] = stokes_ind2(m,q,nx,ny);

% Initialize.
t = 1;
tol = 1e-9;
        
%% Full expm used as reference.
% Can take a long time for large system dimensions!
EFt_expm = Stokes_compute_EFt(eqn,t); 
EFtB_expm = EFt_expm*eqn.Borig;
   
%% Simulate system with implicit Euler.
EFtB_sim = full(simulate_EFtB_stokes(eqn,t));
    
%% Modified Arnoldi iteration
EFtB_arnoldi = TLBT_DAE2_Stokes_modArnoldi(eqn,t,nf,tol,'N',0);
    
%% Output difference.
fprintf('\nError with modified Arnoldi: %.2e\n', norm(EFtB_expm - EFtB_arnoldi));
fprintf('\nError with simulation: %.2e\n', norm(EFtB_expm - EFtB_sim));
         
end

function EFtB = simulate_EFtB_stokes(eqn,t)
%% Simulate system with implicit Euler to obtain EF(t)B.

fprintf('\nSimulating system with implicit Euler for EF(t)B ...\n');
tic;

% Initialize projection PL.
n = size(eqn.E_,1);
nv = full(sum(diag(eqn.E_)));  
np = n - nv;

E11 = eqn.E_(1:nv,1:nv);
A11 = eqn.A_(1:nv,1:nv);
A12 = eqn.A_(1:nv,nv+1:n);

Pi1 = speye(nv,nv) - A12/(A12'/E11*A12)*A12'/E11;
Pi2 = -Pi1*A11/E11*A12/(A12'/E11*A12);
PL = [ Pi1 Pi2; sparse(np,nv) sparse(np,np) ];

% Initialize implicit Euler.
x0 = PL*eqn.Borig;
tau = 1e-3;
tmin = 0;
tmax = t;     
[L,U,P,Q] = lu(eqn.E_-tau*eqn.A_);
ntau = ceil((tmax-tmin)/tau)+1;

% Implicit Euler iteration.
for i=1:ntau
    x = Q*(U\(L\(P*(eqn.E_*x0))));
    x0 = x;
end

% Output EFtB.
EFtB = eqn.E_*x0;

fprintf('EF(t)B simulation time: %.2fs\n\n',toc);

end