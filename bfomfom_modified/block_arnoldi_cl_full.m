function param = block_arnoldi_cl_full(A,B,param)
% param = BLOCK_ARNOLDI_CL_FULL(A,B,param) computes the block Arnoldi basis
% with respect to A, B, and the classical inner product without deflation.
% One must assume that exact deflation will not occur in order to use this
% algorithm; there are no fail-safes for breakdowns.
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
m = param.restart_length;
if param.first_cycle
    [V, param.Beta] = qr(B,0);
else
    V = B;
end
H = [];

s = size(B,2);
jnew = 1:s;
param.Acount = 0;
for j = 1:m
    if isnumeric(A)
        W = A*V(:,jnew);
    else
        W = A(V(:,jnew),param);                                             % MODIFICATION for TLBT
    end
    param.Acount = param.Acount + 1;
    for jj = 1:j
        jold = (jj-1)*s+1:jj*s;
        H(jold,jnew) = V(:,jold)'*W;
        W = W - V(:,jold)*H(jold,jnew);
    end
    jold = jnew;
    jnew = jold(end)+1:jold(end)+s;
    [Q, H(jnew,jold)] = qr(W,0);    
    V = [V Q];                                                              % store new Cl basis vectors
end

Is = eye(s);
E1 = zeros(size(H,2),s); E1(1:s,1:s) = Is;
param.E1 = E1;

Em = zeros(size(H,2),s);
Em(end-s+1:end,end-s+1:end) = Is;
param.Em = Em;

param.svec = s*ones(1,m+1);
param.Hnext = H(end-s+1:end,end-s+1:end);
param.Hshort = H(1:end-s,:);
param.Hfull = H;
param.Vnext = V(:,end-s+1:end);
param.Vshort = V(:,1:end-s);
param.breakdown = m;
end