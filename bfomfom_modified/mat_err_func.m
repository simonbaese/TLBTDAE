function param = mat_err_func(param,Cm,Xim,k,i)
% param = MAT_ERR_FUNC(param,k,i) computes an approximation to the error
% update function with adaptive quadrature for all modified (FOM)^2-like
% routines.
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
global T;

Delta1 = 0;
Delta2 = 0;
acc_old = inf;                                                              % pre-set accuracy measure
accurate = false;
param.refined = false;

[w, t] = param.quad(param.Nquad1,param);
for j = 1:param.Nquad1
    Delta1 = Delta1 + w(j)*Xim(t(j))*Cm(t(j));
end

if nargin == 5  % for fomfom
    while ~accurate
        [w, t] = param.quad(param.Nquad2,param);
        for j = 1:param.Nquad2
            Delta2 = Delta2 + w(j)*Xim(t(j))*Cm(t(j));
        end

        if param.div_by_fapx
            acc = norm(Delta2-Delta1)/norm(param.Fm(:,i));
        else
            acc = norm(Delta2-Delta1);
        end

        if acc > acc_old
            accurate = true;
            T{i}(k) = T{i}(k) + toc;
            if param.verbose
                fprintf('     %s: quadrature stagnation at cycle %g, column %g, accuracy = %1.4g, Nquad = %g.\n',...
                    mfilename,k,i,acc_old,param.Nquad1)
            end
            tic;
        else
            acc_old = acc;
            accurate = acc < param.tol_quad;
        end

        if ~accurate                                                            % increase number of quadrature nodes
            param.Nquad1 = param.Nquad2;
            param.Nquad2 = ceil(sqrt(2)*param.Nquad1);
            if mod(param.Nquad2,2) == 1
                param.Nquad2 = param.Nquad2+1;
            end
            if param.Nquad2 > param.Nquadmax
                param.flag(i) = 4;
                break
            else
                T{i}(k) = T{i}(k) + toc;
                if param.verbose
                    fprintf('     %s: cycle %g, column %g, quadrature not yet accurate enough, accuracy = %1.4g, Nquad = %g.\n',...
                            mfilename,k-1,i,acc_old,param.Nquad1)
                    end
                tic;
                Delta1 = Delta2;
                Delta2 = 0;
                param.refined = true;
            end
        end
    end
elseif nargin == 4 % bfomfom-like methods
    while ~accurate
        [w, t] = param.quad(param.Nquad2,param);
        for j = 1:param.Nquad2
            Delta2 = Delta2 + w(j)*Xim(t(j))*Cm(t(j));
        end

        if param.div_by_fapx
            acc = norm(Delta2-Delta1,'fro')/norm(param.Fm,'fro');
        else
            acc = norm(Delta2-Delta1,'fro');
        end

        if acc > acc_old
            accurate = true;
            T(k) = T(k) + toc;
            if param.verbose
                fprintf('     %s: quadrature stagnation at cycle %g, accuracy = %1.4g, Nquad = %g.\n',...
                    mfilename,k,acc_old,param.Nquad1)
            end
            tic;
        else
            acc_old = acc;
            accurate = acc < param.tol_quad;
        end

        if ~accurate                                                            % increase number of quadrature nodes
            param.Nquad1 = param.Nquad2;
            param.Nquad2 = ceil(sqrt(2)*param.Nquad1);
            if mod(param.Nquad2,2) == 1
                param.Nquad2 = param.Nquad2+1;
            end
            if param.Nquad2 > param.Nquadmax
                param.flag = 4;
                break
            else
                T(k) = T(k) + toc;
                if param.verbose
                    fprintf('     %s: cycle %g, quadrature not yet accurate enough, accuracy = %1.4g, Nquad = %g.\n',...
                            mfilename,k-1,acc_old,param.Nquad1)
                end
                tic;
                Delta1 = Delta2;
                Delta2 = 0;
                param.refined = true;
            end
        end
    end
end
param.Delta = Delta2;
end