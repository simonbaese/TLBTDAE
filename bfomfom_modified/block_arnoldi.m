function param = block_arnoldi(A,B,param)
% param = BLOCK_ARNOLDI(A,B,param) computes the block orthogonal Arnoldi
% basis for K_m(A,B) via the modified Gram-Schmidt procedure, with respect
% to a given block inner product and scaling quotient.  The following
% methods are hard-coded:
%   'cl_full' - the classical method without deflation
%   'cl_defl' - the classical method with deflation
%   'gl' - the global method
%   'li' - the loop-interchange method with deflation
%   'hy' - a hybrid method without deflation
%   'ud' - a user-defined method
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
if ~isfield(param,'first_cycle')
    param.first_cycle = true;
end
switch param.inner_product
    case 'cl_full'
        param = block_arnoldi_cl_full(A,B,param);
    case 'cl_defl'
        param = block_arnoldi_cl_defl(A,B,param);
    case 'gl'
        param = block_arnoldi_gl(A,B,param);
    case 'li'
        param = block_arnoldi_li(A,B,param);
    case 'hy'
        param = block_arnoldi_hy(A,B,param);
    case 'ud'
        param = block_arnoldi_ud(A,B,param);
end
end
