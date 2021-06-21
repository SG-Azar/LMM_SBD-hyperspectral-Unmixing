function patch_idx = PatchExtract(X,type,param)

%   This function extracts spatial neighbors from an mxnxL image and
%   returns the indices. 
%   
%   Inputs:
%           X       : (L-in-m*n) matrix contaning the hyperspectral pixels in
%                     its columns
%           type    : Determine the type of the spatial neighbors. For
%                     nonoverlapping square patches insert 'nonoverlapping'
%                     and for SLIC superpixels insert 'slic'.
%           param   : Insert the related parameters like the patch size.
%
%   Output:
%           patch_idx: an (m*n-in-1) vector that determines which spatial neighbor
%                      each pixel belongs to
%
% m and n are the spatial sizes of the HSI and L is the number of bands
% Note: The slic_HSI function is borrowed from the following paper:
% 
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, "Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing",
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, 2017.
%
% Author: Saeideh Azar
% Last edit: 2021-6-5
%%

XX = X(:,:,1);
[m,n,~]=size(X);
XX = ones(size(XX));

if strcmp(type,'nonoverlapping')
    [I,J] = find(XX);
    cord = [I,J];
    patch = ceil(cord/param.patch_size);
    patch_idx = (patch(:,1)-1)*ceil(max(cord(:,2))/param.patch_size) + patch(:,2);
elseif strcmp(type,'slic')
    Sw = param.superpxNUms;
    P = round(m*n/Sw^2);
    seg = slic_HSI(X, P, param.weight);
    patch_idx=seg.labels;
end

end

