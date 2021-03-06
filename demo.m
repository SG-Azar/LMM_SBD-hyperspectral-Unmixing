%--------------------------------------------------------------------------
%        LMM-SBD for unmixing of the Houston dataset
%--------------------------------------------------------------------------

%   This demo illustrates the unmixing of a 105x128x144 subset of the real
%   Houston dataset with the LMM-SBD method introduced in the following
%   paper:
%   
%
%   S. G. Azar, S. Meshgini, S. Beheshti, and T. Y. Rezaii, "Linear Mixing
%   Model with Scaled Bundle Dictionary for Hyperspectral Unmixing with
%   Spectral Variability", Signal Processing, 2021
%   doi: https://doi.org/10.1016/j.sigpro.2021.108214
%   
%
%   Please kindly cite our paper if you find the codes useful.If you have
%   any questions don't hesitate to contact the author at:
%   sghanbariazar@tabrizu.ac.ir
%   
%
% 
%   Author: Saeideh Azar
%
%--------------------------------------------------------------------------

clear
close all

addpath('functions');
extract_ems = 0;            % flag to perform endmember bundle extraction
P = 5;                      % number of endmembers

%% load data

load Houston
[m,n,L] = size(data);
X = reshape(data,m*n,L)';

%% Extract endmember candidates and cluster them into bundles

bundle_nbr = 10;            % number of VCA runs
percent = 10;               % percentage of pixels considered in each run
clustering = 'kmeans';

if extract_ems == 1
    seed = 100;             % fix seed
    rng(seed)
    [groups, bundle] = batchvca(X, P, bundle_nbr, percent); % extract endmembers and cluster them into groups
else
    load bundles
end

%% unmixing with FCLSU used to initialize the LMM-SBD

tic
disp('unmixing with the FLCSU:')
A_FCLSU = FCLSU(X,bundle)';
toc

%% parameter setting

A_init = A_FCLSU;            % initialize abundances
rho = 10;                    % penalty parameter of the ADMM
tol_a = 10^(-6);             % stop ADMM when relative variations (in norms) of the abundance matrix goes below "tol_a" 
maxiter_ADMM = 100;          % stop ADMM after "maxiter_ADMM" iterations

% extract spatial patches:
param.patch_size = 7;        % spatial patch size for the LMM-SBD_{nol}
param.superpxNUms = 2;       % determines number of super pixels for LMM-SBD_{slic}
param.weight = 0.4;          % for LMM-SBD_{slic}
% addpath('functions');      % uncomment for slic superpixels
patch_idx = PatchExtract(data,'nonoverlapping',param); %for LMM-SBD_{slic} insert 'slic' instead of 'nonoverlapping'


% set the regularization parameters:
lambda1 = 0.1;
lambda2 = 0.2;

%% unmixing with the LMM-SBD

tic
disp('unmixing with the LMM-SBD')
[A_LMMSBD,Scale_LMMSBD] = LMM_SBD(X,bundle,groups,patch_idx,A_init,lambda1,lambda2,rho,maxiter_ADMM,tol_a);
toc

%% sum the abundances within each class 
[A_FCLSU_Global,S_FCLSU] = bundle2global(A_FCLSU,bundle,groups);
[A_LMMSBD_Global, S_LMMSBD] = bundle2global(A_LMMSBD,bundle,groups);

%% Compute reconstruction errors

H_FCLSU = bundle*A_FCLSU;      % reconstruction for FCLSU
H_LMMSBD = bundle*A_LMMSBD*Scale_LMMSBD; % reconstruction for SABLMM

%% Compute RMSE and SAM with reconstructions

RMSE_FCLSU = mean(sqrt(1/L*sum((H_FCLSU-X).^2,1)));  % bara har pixel joda mohasebe mikone (1x13440)
RMSE_LMMSBD = mean(sqrt(1/L*sum((H_LMMSBD-X).^2,1)));

N = m*n;
SAM_FCLSU = zeros(N,1);
SAM_LMMSBD = zeros(N,1);
for k = 1:N
    SAM_FCLSU(k) = 180/pi*real(acos((X(:,k)'*H_FCLSU(:,k))...
        /(norm(X(:,k))*norm(H_FCLSU(:,k)))));
end

for k = 1:N
    SAM_LMMSBD(k) = 180/pi*real(acos((X(:,k)'*H_LMMSBD(:,k))...
        /(norm(X(:,k))*norm(H_LMMSBD(:,k)))));
end

FCLSU_SAM=mean(SAM_FCLSU(:));   
LMMSBD_SAM=mean(SAM_LMMSBD(:));

%% Global abundance maps

A_FCLSU_map = reshape(A_FCLSU_Global',m,n,P);
A_EBLMM_map = reshape(A_LMMSBD_Global',m,n,P);


figure,
for p = 1:P
    
    subaxis(2,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_map(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',15)
    end
    colormap jet
    
  
    subaxis(2,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_EBLMM_map(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('LMM-SBD','fontname','times','fontsize',15)
        xlabel('Concrete','fontname','times','fontsize',15)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',15)
    elseif p == 3
        xlabel('Vegetation','fontname','times','fontsize',15)
    elseif p == 4
        xlabel('Asphalt','fontname','times','fontsize',15)
    else
        xlabel('Colored Structures','fontname','times','fontsize',15)
    end
    
end
set(gcf,'color', 'white')

