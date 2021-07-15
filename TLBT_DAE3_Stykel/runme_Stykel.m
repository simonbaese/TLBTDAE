%% RUNME for Stykel example
% Please modify the parameters and (un-)comment the desired functions.

close all;

%% Choose problem

problem = 'stykel_small'; 
% problem = 'stykel_large';


%% Setup parameters

tlbt_opts = [];
tlbt_opts.tl = 1;   
tlbt_opts.tf = 2;    
tlbt_opts.adiplot = 0;  
tlbt_opts.analyze = 1;


%% Model order reduction with Time-Limited Balanced Truncation

TLBT_DAE3_Stykel_mor(problem,tlbt_opts);


%% Accumulated Hankel singular values plot

% TLBT_DAE3_Stykel_hankel(problem)


%% Analyze inhomogeneities computed with modified Arnoldi

% TLBT_DAE3_Stykel_inhomogeneities(problem,tlbt_opts.tf)


