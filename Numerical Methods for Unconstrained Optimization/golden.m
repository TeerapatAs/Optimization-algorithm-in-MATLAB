function [xmin, fmin, IFLAG, IFunc] = golden( a, b, epsilon, itmax )

% Assign tau
tau = 0.618;

% Define the problem function
f = @(x) (x-1)^4 - 4*sin(x);
%f = @(x) x^2-2 %just for test


% Initial values for IFLAG and IFunc
IFLAG = 0;
IFunc = 0;

% display header
display(['  iter    ', '    xmin   ', '        fxmin      ']);

% Check if interval length (abs(a-b)) < epsilon in the first place? 
% If true, terminate the function, it falied to continue calculating.
if abs(a-b) < epsilon
        %Assign Values for output
        IFLAG = -999;
        xmin = a;
        fmin = f(xmin);
        IFunc = IFunc+1;
        display("enter new a,b with interval between them greater than epsilon");
        return;
end

%Start Looping
for i = 1:itmax
    %compute x1,x2 using tau 
    x2 = tau*(b-a)+a;
    x1 = -tau*(b-a)+b;

    f2 = f(x2);
    f1 = f(x1);
    IFunc = IFunc + 2;

    %if f(x2) > f(x1), set x2 as the new b
    if f2 > f1
        b = x2;
        %set x1 as xmin for now
        xmin = x1;
        fmin = f1;
    %if f(x2) < f(x1), set x1 as the new a
    elseif f2 < f1
        a = x1;
        %set x2 as xmin for now
        xmin = x2;
        fmin = f2;
    end

    %Then, display current values (iteration, xmin, and fmin)
    fprintf('   %d      %.4f      %.4f\n',i,xmin,fmin)
    
    % Finally, check if new x has at least 4 significant digit accuracy.
    % If true, stop the function
    if abs(a-b) < epsilon
        return;
    end
end

%if after looping the abs(a-b) is still not < epsilon, then it failed.
IFLAG = -999
end

     