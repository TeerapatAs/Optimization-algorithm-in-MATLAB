function [lambda, fmin, IFLAG, nF, nG] = golden_func2(FunctionName, x,s,gk, a, b, eta, itmax )

%inputs
% func is the function providing the valuse of f(x), gradient
% [a,b] is interval of lambda
% eps_rel and eps_abs for stopping criterion
% s is search direction
% itmax is maximum iterations

%outputs
% xmin is the golden section result. fmin is f(xmin)
% IFlag --> if success, IFlag = 0. If failed,  IFlag = -999
% nF is number of function calculations.
% nG is number of gradient calculations.

%Assign tau
tau = 0.618;

% Start values for IFLAG ,nF,nG and xmin
IFLAG = 0;
nF = 0;
nG = 0;

%{
%compute initial x1,x2 using tau
x2 = tau*(b-a)+a;
x1 = -tau*(b-a)+b;
%}

% display results
%display(['  iter    ', '    xmin   ', '        fxmin      ']);


% Calculate Wolfe's terms
Wolfe_ = -eta*dot(s,gk);

%Start Looping
for i = 1:itmax
    % compute x1,x2 using tau 
    x2 = tau*(b-a)+a;
    x1 = -tau*(b-a)+b;

    % compute function values, gradients and Hessians at x = x1,x2
    [f_x1, grad_x1] = FunctionName(x + x1*s,2);
    [f_x2, grad_x2] = FunctionName(x + x2*s,2);
    nF = nF+2 ; nG = nG+2; % update number of calculations
    
    %Start Comparing
    %if f(x2) > f(x1), set x2 as the new b
    if f_x2 > f_x1
        b = x2;
        %set x1 as xmin for now
        xmin = x1;
        fmin = f_x1;
    %if f(x2) < f(x1), set x1 as the new a
    else
        a = x1;
        xmin = x2;
        fmin = f_x2;
    end
    
    %Then, display current values (iteration, xmin, and fmin)
    %fprintf('%d      %.4f      %.4f\n',i,xmin,fmin)
    
    % check stopping criterion condition. If the condition is met, terminate.
    % If not, we continue next iteration.
    if xmin == x2
        if (abs(dot(s,grad_x2))) <= Wolfe_
            lambda = x2;
            return;
        end
    else
        if (abs(dot(s,grad_x1))) <= Wolfe_
            lambda = x1;
            return;
        end
    end

end

%if after looping the stopping criterion isn't satisfied. It failed.
IFLAG = -999;
lambda = -1;
end