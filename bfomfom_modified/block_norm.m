function N = block_norm(X,param,A)
% N = BLOCK_NORM(X,param,A) computes a scalar norm for the block vector X.
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
switch param.norm
    case 'A-fro'                                                            % A-weighted Frobenius norm
        N = sqrt(real(trace(X'*(A*X) ) ) );
    case 'fro'                                                              % Frobenius norm; default
        N = norm(X,'fro');
    case 'radau'
        N = sqrt(real(trace(X'*...
            (A*( (param.Theta(1)*speye(size(A)) - A) \X) ) ) ) );
    case 'harmonic'
        N = norm(A*X,'fro');
    case 'max-col'                                                          % max over the 2-norms of each column
        s = size(X,2);
        N = zeros(1,s);
        for i = 1:s
            N(i) = norm(X(:,i),2);
        end
        N = max(N);
    case {2,'2'}                                                            % matrix 2-norm
        N = norm(X,2);
    case {1,'1'}                                                            % matrix 1-norm
        N = norm(X,1);
    case {inf,Inf,'inf','Inf'}                                              % matrix inf-norm
        N = norm(X,Inf);
end