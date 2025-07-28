function [f,g] = HMB(x,options)

f_x = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;

grad_f = @(x) [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7); 
            2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];

%Find function values.
f= f_x(x);

% if-else state logics
if options == 1
    return
elseif options == 2
    g = grad_f(x);
    return
else
    % Wrong Input.
    f= 0;
    g = 0;
    disp('error: Wrong Input.Try again.')
    return
end

end