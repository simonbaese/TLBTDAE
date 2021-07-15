function [w, t] = quad_exp(N,param)
% [w,t] = QUAD_EXP(N,param) computes weights and nodes for approximating
% the Cauchy integral form of the exponential function, as described in [A.
% Frommer, S. Guettel, and M. Schweitzer: Efficient and stable Arnoldi
% restarts for matrix functions based on quadrature, SIAM J. Matrix Anal.
% Appl., 35:661-683, 2014].  This file has also been adapted from
% funm_quad, the body of code associated to the same paper.
%
% This file is a part of the B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices, Electronic Transactions on Numerical Analysis,
% Vol. 44, pp. 100-126, 2017,
%
% and of the modified B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices II: Modified block FOM. Technical Report 19-04-10,
% Department of Mathematics, Temple University. 2019.

%%
a = max(1,max(real(param.ritz))+1);
imtheta = imag(param.ritz);
ccs = abs((param.ritz - a - 1i*imtheta)./imtheta.^2);
cc = min(ccs)/5;                                                            % safety parameter
c = min(cc,0.25);

gamma = @(s) a + 1i*s - c*s.^2;
dgamma = @(s) 1i - 2*c*s;
strunc = sqrt(a - log(param.tol_quad)/c);
s = linspace(-strunc,strunc,N);
h = s(2)-s(1);
w = -h/(2i*pi) * exp(gamma(s)) .* dgamma(s);
t = -gamma(s);
end