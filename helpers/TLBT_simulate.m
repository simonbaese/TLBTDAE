function TLBT_simulate(eqn,Ar,Br,Cr,h,tf,tmax,input)
%% Simple simulation with implicit Euler method. 
%
% Determines the average and maximum output error.
%
% Inputs:
% eqn          original system equation structure
% Ar,Br,Cr     reduced order system matrices
% h            step size
% tf           finish time
% input        system input function 'impulse' or 'step'
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 Stokes demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See LICENSE.md

%% Initialize
x0 = zeros(size(eqn.A_,1),1);
xr0 = zeros(size(Ar,1),1);
tmin = 0;
T = tmin:h:tmax;
nh = floor((tmax-tmin)/h);
hrange = 1:floor(nh/100):nh;

%% Implicit Euler
tic;
[y,yr] = TLBT_impeuler(eqn.E_,eqn.A_,eqn.B,eqn.C,eye(size(Ar)),Ar,Br,Cr,...
    h,tmin,tmax,x0,xr0,input);
fprintf('Simulation time: %.2fs\n', toc);

%% Calculate and plot maximum and average absolute error
abserr = mean(abs(y-yr));
if tmax ~= tf; nh = tf/h; end
max_abserror = max(abserr(1:nh));
avg_abserror = mean(abserr(1:nh));
fprintf('Maximum absolute output error: %.2e\n', max_abserror);
fprintf('Average absolute output error: %.2e\n\n', avg_abserror);

figure;
semilogy(T(hrange),abserr(hrange),'LineWidth',1);
xlabel('time');
ylabel('error');
if tmax ~= tf; xline(tf,':k','LineWidth',1); end
title('Absolute output error');
