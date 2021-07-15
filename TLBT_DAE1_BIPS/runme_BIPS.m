%% RUNME for BIPS model
% Please modify the parameters and (un-)comment the desired functions.

close all;

%% Choose problem

problem = 'bips98_606';
% problem = 'bips07_3078';


%% Setup parameters

tlbt_opts = [];
tlbt_opts.tl = 1;   
tlbt_opts.tf = 3;    
tlbt_opts.adiplot = 0;  
tlbt_opts.analyze = 1;
tlbt_opts.simulate = 1;
tlbt_opts.horizon = 10;
tlbt_opts.input = 'impulse';
% tlbt_opts.input = 'step';


%% Model order reduction with Time-Limited Balanced Truncation

% TLBT_DAE1_BIPS_mor(problem,tlbt_opts);


%% Script for accumulated Hankel singular values plot

% TLBT_DAE1_BIPS_hankel(problem)


%% Script for comparison computation of inhomogeneities for bips98_606

% TLBT_DAE1_BIPS_inhomogeneities()

