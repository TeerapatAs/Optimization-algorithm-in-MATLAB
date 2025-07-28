function  [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG] = BFGS_test(FcnName,x0,epsilon,mu,eta,itmax)

% -------------- Function inputs --------------
%
% FcnName: function to return the value, the gradient, and the Hessian
% of the particular function.
%   Mode 1: return only f.
%   Mode 2: return f and gradient.
%
% x0: starting point of searching.
%
% epsilon: stoping criterion of the minimum search. (norm(x1-x0) < epsilon.)
%
% mu, eta: the parameters used in the stopping criterion for line search.
%
% itmax: max allowed number of iterations. 
% 
% IFRAG: success (0) or not success (-999).

% -------------- Function outputs --------------
%
% xmin, fmin: returned minimum function argument and value, respectively.
%
% Xk ,Fk, Gk, Lk: arrays to keep x, f, gradient and lambda along the search steps.
%
% nF, nG: numbers of f and gradient calculations.
%
% IFLAG: indicate the success. 0 if success, -999 otherwise.

Xk = []; % list to store x_k.
Fk = []; % list to store f_k.
Gk = []; % list to store g_k.
Lk = []; % list to store l_k.

nF = 0; % number of f calculations.
nG = 0; % number of gradient calculations.

IFLAG = -999; % IFLAG: indicate the success.

B = eye(2); % Let the first matrix B be an identity matrix.

for i = 1:itmax

    % strong backtracking.
    [f0, g0] = FcnName(x0, 2); nF = nF + 1; nG = nG + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%a = 1; % first value of lambda.
    s = B\(-g0); % set line search direction.
    [lambda,nFnew,nGnew] = AW_linesearch(FcnName, x0, s, mu, eta);
    % finding lambda that satisfied strong Wolfe's.
    nF = nF + nFnew; nG = nG + nGnew;

    % store values.
    Xk(:,i) = x0; Fk(i) = f0; Gk(:,i) = g0; Lk(i) = lambda;
    
    % update B.
    x1 = x0 + lambda*s;
    [f1,g1] = FcnName(x1, 2); nF = nF + 1; nG = nG + 1;
    delta_g = g1 - g0;
    delta_x = lambda*s;
    B = B + delta_g*delta_g'/dot(delta_g,delta_x) - B*(delta_x*delta_x')*B/(delta_x'*B*delta_x);

    % terminate
    if norm(g1) < epsilon % at local minimum, gradient converges to 0.
        xmin = x1; fmin = f1; IFLAG = 0;
        disp('search successful.');
        break
    end

    % update values
    x0 = x1; f0 = f1; g0 = g1;
end

if IFLAG == -999
    xmin = 0; fmin = 0; disp('search unsuccessful.');
end

end