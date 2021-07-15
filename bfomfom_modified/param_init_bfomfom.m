function [param, modified] = param_init_bfomfom(param)
% [param,modified] = PARAM_INIT_BFOMFOM(param) generates and/or checks the
% input parameter struct.
%
% param = PARAM_INIT_BFOMFOM returns a default setting
% param = PARAM_INIT_BFOMFOM(param) returns a corrected parameter setting
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

modified = 0;

if ~nargin
    param = struct;
    modified = 1;
end

if ~isfield(param,'verbose')
    param.verbose = 1;
    disp('Warning: .verbose not specified, set to 1.');
    modified = 1;
end

if ~isfield(param,'div_by_fapx')
    param.div_by_fapx = 1;
    if param.verbose
        disp('Warning: .div_by_fapx not specified, set to 1.');
    end
    modified = 1;
end

if ~isfield(param,'exact')
    param.conv_check = 'approx';
    if param.verbose
        disp('Warning: .exact not specified, no error available and .conv_check set to ''approx'', which may be inaccurate.');
    end
    modified = 1;
else
    if ~isfield(param,'conv_check')
        param.conv_check = 'exact';
        if param.verbose
            disp('Warning: .conv_check not specified, set to ''exact''.');
        end
        modified = 1;
    end
end

if ~isfield(param,'error_scale') && isfield(param,'exact')
    param.error_scale = norm(param.exact,'fro');
end

if ~isfield(param,'function')
    param.function = 'invSqrt';
    if param.verbose
        disp('Warning: .function not specified, set to ''invSqrt''.');
    end
    modified = 1;
else
    if ~isfield(param,'alpha') && strcmp(param.function,'invAlpha')
        param.alpha = .25;
        if param.verbose
            disp('Warning: .alpha not specified, set to .75.');
        end
        modified = 1;
    end
end

if ~isfield(param,'hermitian')
    param.hermitian = false;
    if param.verbose
        disp('Warning: .hermitian not specified, set to false.');
    end
    modified = 1;
end
    
if ~isfield(param,'inner_product')
    param.inner_product = 'gl';
    if param.verbose
        disp('Warning: .inner_product not specified, set to ''gl''.');
    end
    modified = 1;
else
    if strcmp(param.inner_product,'cl_defl')
        if ~isfield(param,'tol_defl')
            param.tol_defl = 1e-6;
            if param.verbose
                disp('Warning: .tol_defl not specified, set to 1e-6.');
            end
            modified = 1;
        end
    end
end

if ~isfield(param,'max_restarts')
    param.max_restarts = 100;
    if param.verbose
        disp('Warning: .max_restarts not specified, set to 100.');
    end
    modified = 1;
end

if ~isfield(param,'halt_if_diverge')
    param.halt_if_diverge = 1;
    if param.verbose
        disp('Warning: .halt_if_diverge not specified, set to 1.');
    end
    modified = 1;
end

if ~isfield(param,'norm')
    if isfield(param,'hermitian')
        if param.hermitian
            param.norm = 'A-fro';
            if param.verbose
                disp('Warning: .norm not specified, set to ''A-fro''.');
            end
        else
            param.norm = 'fro';
            if param.verbose
                disp('Warning: .norm not specified, set to ''fro''.');
            end
        end
    else
        param.norm = 'fro';
        if param.verbose
            disp('Warning: .norm not specified, set to ''fro''.');
        end
    end
    
end

if ~isfield(param,'Nquad1')
    param.Nquad1 = 100;
    if param.verbose
        disp('Warning: .Nquad1 not specified, set to 32.');
    end
    modified = 1;
end

if ~isfield(param,'Nquad2')
    param.Nquad2 = round(sqrt(2)*param.Nquad1);
    if param.verbose
        disp('Warning: .Nquad2 not specified, set to round(sqrt(2)*param.Nquad1).');
    end
    modified = 1;
end

if ~isfield(param,'Nquadmax')
    param.Nquadmax = 10000;
    if param.verbose
        disp('Warning: .Nquadmax not specified, set to 10000).');
    end
    modified = 1;
end

if ~isfield(param,'reduce_Nquad')
    param.reduce_Nquad = true;
    if param.verbose
        disp('Warning: .reduce_Nquad not specified, set to 1).');
    end
end

if ~isfield(param,'restart_length')
    param.restart_length = 25;
    if param.verbose
        disp('Warning: .restart_length not specified, set to 25.');
    end
    modified = 1;
end

if ~isfield(param,'tol_err')
    param.tol_err = 1e-6;
    if param.verbose
        disp('Warning: .tol_err not specified, set to 1e-6.');
    end
    modified = 1;
end

if ~isfield(param,'tol_quad')
    param.tol_quad = param.tol_err;
    if param.verbose
        disp('Warning: .tol_quad not specified, set to 1e-16.');
    end
    modified = 1;
end

if ~isfield(param,'tol_zero')
    param.tol_zero = 1e-13;
    if param.verbose
        disp('Warning: .tol_zero not specified, set to 1e-13.');
    end
    modified = 1;
end

end