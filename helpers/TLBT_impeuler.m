function [y,yr] = TLBT_impeuler(E,A,B,C,Er,Ar,Br,Cr,h,tmin,tmax,x0,xr0,input)
%% Simple implicit Euler implementation for the BIPS and Stokes model. 
%
% Inputs:
% E,A,B,C      original system matrices
% Er,Ar,Br,Cr  reduced order system matrices
% h            time step size
% tmin         start time 
% tmax         end time
% x0,x0r       initial states of the full and reduced order models
% input        system input function 'impulse' or 'step'
%
% Outputs:
% y,yr         outputs of the full and reduced systems in [tmin,tmax]
%
% Author: Simon Bäse (2021)
% Partly forked from DAE2 stokes demo in M-M.E.S.S. toolbox by
% Saak, J. and Köhler, M. and Benner, P. and others
%
% See LICENSE.md

[L,U,P,Q] = lu(E-h*A);
[Lr,Ur,Pr] = lu(Er-h*Ar);
 
nh = ceil((tmax-tmin)/h);
y = zeros(size(C,1),nh);
yr = y;

for i = 1:nh
    
    if strcmp(input,'step')
        % Implicit Euler with step like input.
        x = Q*(U\(L\(P*(E*x0+h*sum(B,2)))));
        xr = Ur\(Lr\(Pr*(Er*xr0+h*sum(Br,2))));
    else
        % Implicit Euler.
        if i > 1
            x = Q*(U\(L\(P*(E*x0))));
            xr = Ur\(Lr\(Pr*(Er*xr0)));
        % Impulse.
        else
            x = Q*(U\(L\(P*(E*x0+h*sum(B,2)))));
            xr = Ur\(Lr\(Pr*(Er*xr0+h*sum(Br,2))));
        end
    end
    
    y(:,i) = C*x;
    yr(:,i) = Cr*xr;
    x0 = x;
    xr0 = xr;
end