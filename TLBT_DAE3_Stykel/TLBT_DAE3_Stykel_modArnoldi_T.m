function [EFtC,Theta0,V,H] = TLBT_DAE3_Stykel_modArnoldi_T(eqn,t,nf,tol)
%% Modified block Arnoldi iteration for mechanical multibody systems
%
% Approximates the inhomogeneities of time-limited projected generalized
% Lyapunov equations for descriptor systems with mechanical multibody
% system structure.
%
%   A'*X*E + E'*X*A = -C'*C - Ct*Ct'
%
% See Algorithm 9 in [1]!
% 
% Inputs:
% eqn          struct of system matrices as defined in TLBT_DAE3_Stykel_mor
% t            final time
% nf           number of finite eigenvalues of matrix pencil sE-A
% tol          error tolerance
%
% Outputs:
% EFtC         inhomogeneity at t: Ct' = V*exp(t*H)*V'*C',
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

fprintf(['\nStarting Arnoldi iteration for EF(t)C at ', ...
                't = %g ...\n'],t);   
tic;

% Counter for evaluation in error estimate.
keval = 0;

% Determine relevant system dimensions.
nv = size(eqn.M_,1);
np = size(eqn.G_,1);
m = size(eqn.C,1);
limit = floor(nf/m);

% Initialize saddle point system leading matrix.
% We use LUPQ factorization since system matrices are sparse.
SPM = [eqn.M_ -eqn.G_'; eqn.G_ sparse(np,np)];
[L,U,P,Q] = lu(SPM);

% Allocate memory and initialize iterator.
V = zeros(2*nv,limit*m);
H = zeros(limit*m,limit*m);
IpM = blkdiag(speye(nv),eqn.M_);
k = 1;

% Project initial block vector.
% Note that Vk carries a hidden M inverse in front of Theta2.
% Note that E_ is D.
C1 = eqn.C(:,1:nv)';
C2 = eqn.C(:,nv+1:2*nv)';

Theta2 = Q*(U\(L\(P*[C2;sparse(np,m)])));
Theta1 = Q*(U\(L\(P*[C1 + eqn.E_'*Theta2(1:nv,:);sparse(np,m)])));
Vk = [eqn.M_*Theta1(1:nv,:) - eqn.E_'*Theta2(1:nv,:);Theta2(1:nv,:)];   
Theta0 = [IpM*Vk;sparse(np,m)];

[~,RC] = mess_mgs(IpM*Vk);
Vk = Vk/RC;
Vk1 = Vk(1:nv,:);
Vk2 = Vk(nv+1:2*nv,:);
V(:,1:m) = Vk;    

% Modified block Arnoldi iteration for mechanical multi-body systems
% of index 3 for C.
while k <= limit

    % Projection avoiding update. 
    Theta2 = Q*(U\(L\(P*[Vk1 + eqn.E_'*Vk2;sparse(np,m)])));
    Theta1 = Q*(U\(L\(P*[eqn.K_'*Vk2 + eqn.E_'*Theta2(1:nv,:);...
                sparse(np,m)])));
    W = [eqn.M_*Theta1(1:nv,:) - eqn.E_'*Theta2(1:nv,:);Theta2(1:nv,:)];

    % Modified Gram-Schmidt orthogonalization.
    for j = 1:k
       Vj = V(:,(j-1)*m+1:j*m);
       Hjk = Vj'*W;
       H((j-1)*m+1:j*m,(k-1)*m+1:k*m) = Hjk;
       W = W - Vj*Hjk;         
    end

    % Set next iteratives. Here Vk holds the k+1 iterative.
    [~,Hnext] = mess_mgs(IpM*W);
    Vk = W/Hnext;         
    Vk1 = Vk(1:nv,:);
    Vk2 = Vk(nv+1:2*nv,:);
    if k < limit
       H(k*m+1:(k+1)*m,(k-1)*m+1:k*m) = Hnext;
       V(:,k*m+1:(k+1)*m) = [Vk1;Vk2]; 
    end 
    k = k + 1;

    % Perform a simple error estimate in every 3 steps of the 
    % iteration. Set minimum iteration at 5 to avoid false 
    % positives.
    if k < limit && k > 3 && ~mod(k-1,3) && tol > 0
       e = eye(k-1); 
       E1 = kron(e(:,1),eye(m));
       Ek = kron(e(:,k-1),eye(m)); 
       res = Hnext*Ek'*expm(H(1:(k-1)*m,1:(k-1)*m))*E1*RC;
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

% Reduce iterator and determine V and H. Here we eliminate the hidden 
% M inverse which was carried through the iteration in Vk2. Also, we
% add zeros to V for the evaluation of expm.
k = k - 1;
V = [IpM*V(:,1:k*m);sparse(np,m*k)];
H = H(1:k*m,1:k*m);

% Return EFtB. Note V is already in the projected subspace of the
% deflating subspace corresponding to the finite eigenvalues of the
% matrix pencil sE-A.
e = eye(k);
Ehat = kron(e(:,1),eye(m));
EFtC = V*expm(t*H)*Ehat*RC;

fprintf('Total Arnoldi computation time: %.3fs\n',toc);   

end