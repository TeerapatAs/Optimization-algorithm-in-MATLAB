function [f,gradient] = Rosenbrock(x,options)
% Use the same Rosenbrock as HW1,HW2
f_x = @(x) 10*x(1)^2 + x(1) + 5*x(2)^2 + 2*x(2)*x(3) + 2*x(2) + 2*x(3)^2 - x(3);
grad_f = @(x) hess(x)

%Find function values.
f= f_x(x);

% if-else state logics
if options == 1
    return
elseif options == 2
    gradient = grad_f(x);
    return
else
    % Wrong Input.
    f= 0;
    gradient = 0;
    disp('error: Wrong Input.Try again.')
    return
end

end
