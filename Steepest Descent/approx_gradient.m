function G = approx_gradient(f,x)
%inptus
% f is a function
% x is a point which we used to compute hessian
%outputs
% G- computed gradient
    h = 1e-8; % Step size 
    n = length(x); % Number of variables
    G = zeros(n, 1); % Initialize the gradient vector

    for i = 1:n
        % e_i is the unit vector for the i th variable
        e_i = zeros(n, 1);
        e_i(i) = 1;

        % Compute approximate gradient
        G(i) = (f(x + h*e_i) - f(x - h*e_i)) / (2*h);
    end
end