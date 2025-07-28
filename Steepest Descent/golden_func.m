function [xmin, fmin, IFLAG, nF, nG] = golden_func(FunctionName, a, b, eps_rel, eps_abs,s, itmax )

%inputs
% func is the function providing the valuse of f(x), gradient, Hessian
% [a,b] is interval
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

% Start values for IFLAG ,nF and nG
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

%Due to the way the Newton method calls this function, we need to check the interval frequently.
[f_b, grad_b , null1] = FunctionName(b,2);
nF = nF+1 ; nG = nG+1;
%if slope at b <0 , this interval cannot be used.
if dot(s,grad_b) < 0 
    IFLAG = -999; % fail flag
    xmin = 0; fmin = 0;
    return;
end


%Start Looping
for i = 1:itmax
    % compute x1,x2 using tau 
    x2 = tau*(b-a)+a;
    x1 = -tau*(b-a)+b;

    % compute function values, gradients and Hessians at x = x1,x2
    [f_x1, grad_x1 , null1] = FunctionName(x1,2);
    [f_x2, grad_x2 , null1] = FunctionName(x2,2);
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
        %set x2 as xmin for now
        xmin = x2;
        fmin = f_x2;
    end
    
    %Then, display current values (iteration, xmin, and fmin)
    %fprintf('%d      %.4f      %.4f\n',i,xmin,fmin)
    
    % check stopping criterion condition. If the condition is met, terminate.
    % If not, we continue to next iteration.
    if (f_x2 < f_x1) && (abs(dot(s,grad_x2)) <= abs(dot(s,grad_x1))*eps_rel + eps_abs)
        return;
    end

end

%if after looping the stopping criterion isn't satisfied. It failed.
IFLAG = -999;
end