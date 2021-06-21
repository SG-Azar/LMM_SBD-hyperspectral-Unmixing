function Y = prox_GF1(X,groups,tau)

% THis function computes the proximal operator of:
%         ||Xi||_{G,F,1}
% Author: Saeideh Azar

group_num = max(groups);
Y = zeros (size(X));

    for j=1:group_num
        
        X_group = X (groups==j,:);
        Y_group = matrix_soft(X_group,tau);
        Y(groups==j,:) = Y_group;
    end


end


function Out = matrix_soft(Z,tau)

% computes matrix soft thresholding of matrix Z

normF = norm(Z,'fro');
ZZ = max(((normF-tau)/normF),0);
Out = ZZ.*Z;

end

