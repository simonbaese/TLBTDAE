function [EFtB,Theta0,V,H] = TLBT_DAE1_BIPS_modArnoldi(eqn,t,nf,nin,tol,type)
%% Modified block Arnoldi iteration for semi-explicit systems of index 1
%
% Approximates the inhomogeneities of time-limited projected generalized
% Lyapunov equations for descriptor systems with Stokes-like structure.
%
%   A*X*E' + E*X*A' = -B*B' + Bt*Bt'
%   A'*X*E + E'*X*A = -C'*C - Ct*Ct'
%
% See Algorithm 6 in [1]!
% 
% Inputs:
% eqn           struct of system matrices as defined in TLBT_DAE1_BIPS_mor
% t             final time
% nf            number of finite eigenvalues of matrix pencil sE-A
% nin           number of infinite eigenvalues of matrix pencil sE-A
% tol           error tolerance
% type          'N' - normal, 
%               'T' - transposed, 
%               for controllability or observability Lyapunov equation
%
% Outputs:
% EFtB          'N' - inhomogeneity at t: Bt = V*exp(t*H)*V'*B,
%               'T' - inhomogeneity at t: Ct' = V*exp(t*H)*V'*C',
% Theta0        inhomogeneity at time point zero
% V             projection matrix
% H             upper Hessenberg matrix
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
% See License.MD    
    
tic;
    
% Counter for evaluation in error estimate.
keval = 0;

% Use block structure of system matrices.
E11 = eqn.E_(1:nf,1:nf);
A11 = eqn.A_(1:nf,1:nf);
A12 = eqn.A_(1:nf,nf+1:end);
A21 = eqn.A_(nf+1:end,1:nf);
A22 = eqn.A_(nf+1:end,nf+1:end);

% Initialize right-hand side and determine relevant system dimensions.
if type == 'N'
    fprintf(['\nStarting Arnoldi iteration for EF(t)B at ', ...
                't = %g ...\n'],t); 
    m = size(eqn.B,2);
    B1 = eqn.B(1:nf,:);
    B2 = eqn.B(nf+1:end,:);
else
    fprintf(['\nStarting Arnoldi iteration for CF(t)E at ', ...
                't = %g ...\n'],t);
    m = size(eqn.C,1);
    B1 = eqn.C(:,1:nf)';  
    B2 = eqn.C(:,nf+1:end)';
end 
limit = floor(nf/m);

% Initialize saddle point system leading matrix.
% We use LUPQ factorization since system matrices are sparse.
SPM = [E11 A12;sparse(nin,nf) A22];
[L,U,P,Q] = lu(SPM);

% Allocate memory and initialize iterator.
V = zeros(nf,limit*m);
H = zeros(limit*m,limit*m);
k = 1;

% Project initial block vector.
Theta = Q*(U\(L\(P*[B1;B2])));
[~,R] = mess_mgs(E11*Theta(1:nf,:));
Theta0 = [E11*Theta(1:nf,:);sparse(nin,m)];
Vk = Theta(1:nf,:)/R;
V(:,1:m) = Vk;

% Modified block Arnoldi iteration for semi-explicit systems of 
% index 1.
while k <= limit

    % Projection avoiding update.
    Theta = Q*(U\(L\(P*[A11*Vk;A21*Vk])));
    W = Theta(1:nf,:);

    % Modified Gram-Schmidt orthogonalization.
    for j = 1:k
       Vj = V(:,(j-1)*m+1:j*m);
       Hjk = Vj'*W;
       H((j-1)*m+1:j*m,(k-1)*m+1:k*m) = Hjk;
       W = W - Vj*Hjk;         
    end

    % Set next iteratives.
    [~,Hnext] = mess_mgs(E11*W);
    Vk = W/Hnext;       
    if k < limit
       H(k*m+1:(k+1)*m,(k-1)*m+1:k*m) = Hnext;
       V(:,k*m+1:(k+1)*m) = Vk; 
    end 
    k = k + 1;

    % Perform a simple error estimate in every 10 steps of the 
    % iteration. Set minimum iteration at 11 to avoid false 
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
V = [E11*V(1:nf,1:k*m);sparse(nin,m*k)];
H = H(1:k*m,1:k*m);

e = eye(k);
Ehat = kron(e(:,1),eye(m));

% Return EFtB. Note V is already in the projected subspace of the
% deflating subspace corresponding to the finite eigenvalues of the
% matrix pencil sE-A.
EFtB = V*expm(t*H)*Ehat*R;

fprintf('Total Arnoldi computation time: %.2fs\n',toc);

end