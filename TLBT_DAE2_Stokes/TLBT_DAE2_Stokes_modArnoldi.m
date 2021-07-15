function [EFtB,Theta0,V,H] = TLBT_DAE2_Stokes_modArnoldi(eqn,t,nf,tol,type,symm)
%% Modified block Arnoldi iteration for Stokes-like systems of index 2
%
% Approximates the inhomogeneities of time-limited projected generalized
% Lyapunov equations for descriptor systems with Stokes-like structure.
%
%   A*X*E' + E*X*A' = -B*B' + Bt*Bt'
%   A'*X*E + E'*X*A = -C'*C - Ct*Ct'
%
% See Algorithm 7 in [1]!
% 
% Inputs:
% eqn          struct of system matrices as defined in TLBT_DAE2_Stokes_mor
% t            final time
% nf           number of finite eigenvalues of matrix pencil sE-A
% tol          error tolerance
% type         'N' - normal, 
%              'T' - transposed, 
%              for controllability or observability Lyapunov equation
% symm         0 - Arnoldi method,
%              1 - Lanczos method
%
% Outputs:
% EFtB         'N' - inhomogeneity at t: Bt = V*exp(t*H)*V'*B,
%              'T' - inhomogeneity at t: Ct' = V*exp(t*H)*V'*C',
% Theta0       inhomogeneity at time point zero
% V            projection matrix
% H            upper Hessenberg matrix
%         
% Note that the output is already projected to the correct deflating
% subspaces.
%
% [1] S. M. Bäse. "Time-Limited Balanced Truncation Model Order Reduction 
%     for Descriptor Systems". Master thesis. Technische Universität 
%     Berlin: Institut für Mathematik, 2021.
%
% Author: Simon Bäse (2021)
%
% See LICENSE.md

tic;

% Counter for evaluation in error estimate.
keval = 0;

% Use block structure of system matrices.
n = size(eqn.E_,1);
nv = full(sum(diag(eqn.E_)));    
np = n - nv; 
E11 = eqn.E_(1:nv,1:nv);
A11 = eqn.A_(1:nv,1:nv);
A12 = eqn.A_(1:nv,nv+1:n);

% Initialize right-hand side and determine relevant system dimensions.
if type == 'N'
    fprintf(['\nStarting Arnoldi iteration for EF(t)B at ', ...
                't = %g ...\n'],t); 
    m = size(eqn.Borig,2);
    B1 = eqn.Borig(1:nv,:);        
else
    fprintf(['\nStarting Arnoldi iteration for CF(t)E at ', ...
                't = %g ...\n'],t);
    m = size(eqn.Corig,1);
    B1 = eqn.Corig(:,1:nv)';      
end    
limit = floor(nf/m);

% Initialize saddle point system leading matrix.
% We use LUPQ factorization since system matrices are sparse.
SPM = [E11 A12; A12' sparse(np,np)];
[L,U,P,Q] = lu(SPM);

% Allocate memory and initialize iterator.
V = zeros(nv,limit*m);
H = zeros(limit*m,limit*m);
k = 1;

% Project initial block vector.
Theta = Q*(U\(L\(P*[B1;sparse(np,m)])));
Theta0 = [E11*Theta(1:nv,:);sparse(np,m)];
[~,R] = mess_mgs(E11*Theta(1:nv,:));  
Vk = Theta(1:nv,:)/R;
V(:,1:m) = Vk;

% Modified block Arnoldi iteration for Stokes-like systems of
% index 2.
while k <= limit

    % Projection avoiding update.
    Theta = Q*(U\(L\(P*[A11*Vk;sparse(np,m)])));
    W = Theta(1:nv,:);

    % Symmetric Lanczos method.
    if symm         
        Hkk = Vk'*W;          
        if k == 1
           H(1:m,1:m) = Hkk;
        else
           W = W - V(:,(k-2)*m+1:(k-1)*m)*...
                H((k-2)*m+1:(k-1)*m,(k-1)*m+1:k*m);
           H((k-1)*m+1:k*m,(k-1)*m+1:k*m) = Hkk;       
        end
        W = W - Vk*Hkk;

    % Modified Gram-Schmidt orthogonalization.
    else       
        for j = 1:k
           Vj = V(:,(j-1)*m+1:j*m);
           Hjk = Vj'*W;
           H((j-1)*m+1:j*m,(k-1)*m+1:k*m) = Hjk;
           W = W - Vj*Hjk;         
        end
    end

    % Set next iteratives.
    [~,Hnext] = mess_mgs(E11*W);
    Vk = W/Hnext;       
    if k < limit
       H(k*m+1:(k+1)*m,(k-1)*m+1:k*m) = Hnext;
       if symm; H((k-1)*m+1:k*m,k*m+1:(k+1)*m) = Hnext'; end
       V(:,k*m+1:(k+1)*m) = Vk; 
    end 
    k = k + 1;

    % Perform a simple error estimate in every 10 steps of the 
    % iteration. Set minimum iteration at 10 to avoid false 
    % positives.
    if k < limit && k > 11 && ~mod(k-1,10) && tol > 0
       e = eye(k-1); 
       E1 = kron(e(:,1),eye(m));
       Ek = kron(e(:,k-1),eye(m)); 
       res = Hnext*Ek'*expm(t*H(1:(k-1)*m,1:(k-1)*m))*E1*R;
       keval = keval + 1;
       if norm(res) < tol
           break;
       end
    end 
end

% Calculate full expm for low-rank H.
if tol > 0
    fprintf('Times error estimate approximated with expm: %d\n',keval);
end
fprintf('The Krylov-subspace dimension is: %d (max: %d)\n',m*(k-1),nf);
fprintf('Starting low-rank MATLAB expm ...\n');

k = k - 1;
V = [E11*V(1:nv,1:k*m);sparse(np,m*k)];
H = H(1:k*m,1:k*m);

e = eye(k);
Ehat = kron(e(:,1),eye(m));

% Return EFtB. Note V is already in the projected subspace of the
% deflating subspace corresponding to the finite eigenvalues of the
% matrix pencil sE-A.
EFtB = V*expm(t*H)*Ehat*R;

fprintf('Total Arnoldi computation time: %.2fs\n',toc);  

end