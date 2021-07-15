function param = krylov_basis(A,B,param)
% param = KRYLOV_BASIS(A,B,param) is a wrapper function returning the
% Krylov basis specified by A, B, and param.
%
% This file is a part of the B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices, Electronic Transactions on Numerical Analysis,
% Vol. 44, pp. 100-126, 2017.
%
% and of the modified B(FOM)^2 code described in detail in
%
% A. Frommer, K. Lund, and D. B. Szyld: Block Krylov subspace methods for
% functions of matrices II: Modified block FOM. Technical Report 19-04-10,
% Department of Mathematics, Temple University. 2019.

%%
if param.hermitian
    if strcmp(param.inner_product,'nb')
        param = lanczos(A,B,param);
    else
        param = block_lanczos(A,B,param);
    end
else
    if strcmp(param.inner_product,'nb')
        param = arnoldi(A,B,param);
    else
        param = block_arnoldi(A,B,param);
    end
end
end