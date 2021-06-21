function Y = prox_row21(X,tau)

% THis function computes the proximal operator of:
%         ||Xi||_{row,2,1}
% Author: Saeideh Azar

    Y = vector_soft_col(X',tau);
    Y = Y';

end

