function param = bfomfom(A,B,param)
% param = BFOMFOM(A,B,param) is a wrapper function for the variations of
% (B)FOMFOM.  Note that param.mod = 'none' runs a (b)fomfom_mod with the
% modification M = 0, whereas an unspecific param.mod runs (b)fomfom.
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
if strcmp(param.function,'inv')
    param = bfom(A,B,param);
else
    if isfield(param,'mod')
        switch param.inner_product
            case {'cl_defl', 'cl_full', 'li', 'hy', 'ud'}
                param = bfomfom_cl_mod(A,B,param);
            case 'gl'
                param = bfomfom_gl_mod(A,B,param);
            case 'nb'
                param = fomfom_mod(A,B,param);
        end
    else
        switch param.inner_product
            case {'cl_defl', 'cl_full', 'li', 'hy', 'ud'}
                param = bfomfom_cl(A,B,param);
            case 'gl'
                param = bfomfom_gl(A,B,param);
            case 'nb'
                param = fomfom(A,B,param);
        end
    end
end
end