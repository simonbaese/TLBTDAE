function Z = TLBT_definite_approximation(L,D,tol)
%% TLBT definite approximation
%
% Find Cholesky factor Z such that ZZ' is positive definite and 
% ZZ' approximates L*D*L' (output of LDLT ADI iteration).
% 
% Inputs:
% L,D          outputs of LDLT ADI iteration (mess_lradi)
% tol          cut off tolerance, should be very small       
%
% Outputs: 
% Z            factor of positive definite approximation
%
% Author: Simon Bäse (2021)
%
% See LICENSE.md

%% Definite approximation

% Eigenvalue factorization of out.Z.
% The first computation of T introduces rounding errors. Therefore, we
% force symmetry again. This also speeds up eig().
[Q,R] = qr(L,0); 
T = R * D * R.'; 
T = (T + T.')/2; 
[V,DIAG] = eig(T);

% Truncate eigenvalues smaller than tol.
[d,I] = sort(diag(DIAG),'descend');
for k = 1:size(d)
    if d(k) < tol; break; end
end 
I = I(1:k-1);
DIAG = DIAG(:,I);

% Set Cholesky factor.
Z = Q*V*sqrt(DIAG);

end

