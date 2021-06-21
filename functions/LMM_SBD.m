function [A,Scale] = LMM_SBD(Y,bundle,groups,patch_idx,A_init,lambda1,lambda2,rho,maxiter_ADMM,tol_a)
   
   
%   This function minimizes the following cost with respect to A and S:
%  
%   J(A,S) = 1/2 * ||Y - BAS||_{F}^{2} + lambda1*||A||_{row,2,1} +
%   lambda2*||A||_{G,F,1} + I_{R+}(A) + I_{R+}(S)
%
%   Here, X is a spatial patch, B is a collection of endmember bundles, A is
%   the matrix of abundances for each spatial patch and S (or Scale) is a 
%   diognal matrix of scaling factors
%
%   
%
% Inputs:
%       Y :             (L-in-m*n) matrix of hyperspectral data
%       bundle:         LxQ matrix of the bundle dictionary
%       groups:         a vector that determines which endmember each atom of
%                       the bundle dictionary represents.
%       patch_idx:      an (m*n-in-1) vector that determines which spatial neighbor
%                       each pixel belongs to
%       A_init:         initial abundances
%       lambda1:        regularization parameter
%       lambda2:        regularization parameter
%       rho:            penalty parameter of the ADMM
%       maxiter_ADMM:   maximum number of iterations before the algorithm stops
%       tol_a:          stops the ADMM when relative variations (in norms) of the 
%                       abundance matrix goes below "tol_a"
%
% Output:
%       A:              (QxN)matrix of abundance matrix
%       Scale:          the diogonal matrix (NxN) of scaling factors 
%
% Author: Saeideh Azar
% Last edit: 2021-6-5
% 
%%

P = length(groups); %  (Q in paper): total number of endmember signature or dictionry size 
N = size(Y,2);      % total number of pixels
patch_Num = max(patch_idx);
A_init = A_init ./ repmat(sum(A_init), size(A_init, 1), 1);  % imposing sum-to-one constraint
Scale = zeros(N,N);

B = bundle;
BtBrhoI = B'*B +rho*eye(P);
A = zeros(size(A_init));

%% ADMM 
for k = 1:patch_Num  % for each patch
    
    % initialization
    Y_patch = Y(:,patch_idx==k);
    N_patch = size(Y_patch,2);
    A1 = A_init(:,patch_idx==k);
    S = diag(ones(N_patch,1));
    C = zeros(size(A1));
    
    U = A1;
    V = A1;
    W = A1;
    T = zeros(N_patch,N_patch);
    M = A1*S;
    D = C;
    E = C;
    F = C;
    G = zeros(N_patch,N_patch);
    
    for i = 1:maxiter_ADMM  % main loop of ADMM
        
        A_old = A1;
        
        %update A
        A1 = (rho*(M*S'-C*S'+U-D+V-E+W-F))/(rho*(S*S')+3*rho*eye(N_patch));
        A1 = A1./repmat(sum((A1)),P,1);
        
        %update S
        SL = rho*(A1'*A1)+rho*eye(N_patch);
        SR = rho*(A1'*M-A1'*C+T-G);
        for j=1:N_patch
            S(j,j) = SL(:,j)\SR(:,j);
        end
        
        %update M
        M = BtBrhoI\(B'*Y_patch+rho*A1*S+rho*C);
        
        %update U
        U = prox_row21(A1+D,lambda1/rho); 
        
        %update V
        V = prox_GF1(A1+E,groups,lambda2/rho);
        
        %update W
        W = max(A1+F,0);
        
        %update T
        T = max(S+G,0);
        
        % dual updates:
        C = C + A1*S - M;
        D = D + A1 - U;
        E = E + A1 - V;
        F = F + A1 - W;
        G = G + S - T;
        
        rel_A = norm(A_old-A1,'fro')/norm(A_old,'fro')^2;
      
        if i>1 && rel_A < tol_a
            break
        end
        
    end    % end of ADMM
    
A(:,patch_idx==k) = A1;   
Scale(patch_idx==k, patch_idx==k) = S;

end        % end of each patch
                    

end

