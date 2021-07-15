%% RUNME for Stokes example
% Please modify the parameters and (un-)comment the desired functions.

close all;

%% Choose system size
% Computation time for large system might be very long.

% problem = 'tiny'; % n = 645, m = 10, q = 10
problem = 'small'; % n = 1825, m = 10, q = 10
% problem = 'medium'; % n = 3605, m = 10, q = 10
% problem = 'large'; % n = 19521, m = 5, q = 5

%% Setup parameters

tlbt_opts = [];
tlbt_opts.tl = 1;   
tlbt_opts.tf = 0.025;    
tlbt_opts.adiplot = 0;  
tlbt_opts.analyze = 1;
tlbt_opts.simulate = 1;
tlbt_opts.horizon = 0.1;
tlbt_opts.input = 'impulse';
% tlbt_opts.input = 'step';


%% Model order reduction with Time-Limited Balanced Truncation

TLBT_DAE2_Stokes_mor(problem,tlbt_opts);


%% Accumulated Hankel singular values plot

% TLBT_DAE2_Stokes_hankel(problem)


%% Analyze inhomogeneities computed with modified Arnoldi

% TLBT_DAE2_Stokes_inhomogeneities(problem)


%% Comparison of mod. Arnoldi and B(FOM)^2 with restarts for medium system

% TLBT_DAE2_Stokes_BFOM_restarts()

