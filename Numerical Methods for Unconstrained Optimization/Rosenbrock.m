function [f,gradient] = Rosenbrock(x,options)
% Use the same Rosenbrock as HW1,HW2
f_x = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
grad_f = @(x)[400*x(1)*(x(1)^2 - x(2)) - 2*(1-x(1)); 200*(x(2)-x(1)^2)];

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
