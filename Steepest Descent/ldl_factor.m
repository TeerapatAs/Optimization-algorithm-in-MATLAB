function [L,D] = ldl_factor(H)

%input
% H is hessian

%output --> LDL^T = H

n = length(H) ; %number of variables
D = zeros(n,1); %Intialize for diagonal D
L = zeros(n) ; %Intialize for L


%algorithm
for j = 1:n
    % calculate D(j) 
    if j == 1
         D(j) = H(j,j)
    else
        h_jj = H(j,j);
        %summation term
        sum0 = 0;
        for k = 1:j-1
            sum0 = sum0 + D(k)*(L(j,k)^2);
        end
        D(j) = h_jj - sum0
    end
    % make the diagonal elemetns of L equal to 1
    L(j,j) = 1
    
    % compute elements of L in column j th
    for i = j+1:n
        if j == 1
            L(i,j) = (H(i,j))/D(j)
        else
            h_ij = H(i,j);
            sum1 = 0;
            %summation term
            for k_ = 1:j-1
                sum1 = sum1 + D(k_)*L(j,k_)*L(i,k_);
            end
            L(i,j) = (h_ij - sum1)/D(j)
        end
    end

end

D = diag(D);
end
            