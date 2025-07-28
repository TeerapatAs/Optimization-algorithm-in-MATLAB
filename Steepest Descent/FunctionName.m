function [f,grad,Hess] = FunctionName(x,options)

% Define Function
f_x = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
grad_f = @(x)[400*x(1)*(x(1)^2 - x(2)) - 2*(1-x(1)); 200*(x(2)-x(1)^2)];
hess_f = @(x)[2 - 400*(x(2)-3*x(1)^2), -400*x(1); -400*x(1), 200] ;

%Find function values.
f= f_x(x);


% if-else state logics
if options == 1
    grad = [];
    Hess = [];
    return
elseif options == 2
    % Find function values, gradients
    grad = grad_f(x);
    Hess = [];
    return
elseif options == 3
    %Find function values, gradients and Hessians
    grad = grad_f(x);
    Hess = hess_f(x);
    return
else
    % Wrong Input.
    f= 0;
    grad = 0;
    Hess = 0;
    disp('error: Wrong Input.Try again.')
    return
end