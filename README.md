# Optimization-algorithm-in-MATLAB

The codes that I wrote in my 4th academic year. (Intro to Optimization Class).

# Example from Steepest Descent Folder:

Consider we need to find the local minimum x* of the objective function (In this case, Rosenbrock's function);
<img width="501" height="35" alt="image" src="https://github.com/user-attachments/assets/56bfdbf5-b61e-45cf-a713-eb77c0baaee8" />

By using modified Newton's Method to find the solution. The line search's search direction is <img width="221" height="31" alt="image" src="https://github.com/user-attachments/assets/92c7c7c5-6fe5-430e-8590-eb71143942fa" />. 

In case the Hessian matrix of the objective function is a singular matrix (the inverse of <img width="72" height="39" alt="image" src="https://github.com/user-attachments/assets/3e583a1a-cdcf-4ce9-ab3b-878f6de2fca4" /> cannot be found.), We change the search direction method to steepest descent <img width="128" height="29" alt="image" src="https://github.com/user-attachments/assets/32f6c1b3-ca6d-4fc2-8157-937e549c583d" />.

Quick reminder: ? What is a line search and search direction?
<img width="1364" height="544" alt="image" src="https://github.com/user-attachments/assets/915c14e4-5729-407a-915e-2a292d48f94a" />

# Result
Here is the result. (After 51 iterations, we find the correct solution x = (1,1) with 10^-4 precision).
<img width="651" height="277" alt="image" src="https://github.com/user-attachments/assets/18cb2e33-c160-4d97-9786-9c192fbb5298" />

**Code Snippet**
<sub>
function [xmin,fmin,Xk,Fk,Gk,nF,nG,nH,IFLAG] = Newton(FunctionName,x0,epsilon,e_rel,e_abs,itmax)

%inputs
% func is the function providing the valuse of f(x), gradient, Hessian
% x0 is initial point
% e_rel and e_abs for stopping criterion for line search
% epsilon for stopping criterion for Newton
% itmax is maximum iterations

%Outputs
% xmin is the golden section result. fmin is f(xmin)
% IFlag --> if success, IFlag = 0. If failed,  IFlag = -999
% Xk ,Fk, Gk are arrays that store function values, gradient and Hessians at x_k, in order.
% nF,nG,nH is number of function ,gradient and Hessian calculations, in order.


% inital values for nF,nG,nH and IFlag.
nF = 0; 
nG = 0; 
nH = 0; 
IFLAG = 0;

% array for Xk ,Fk, Gk.
Xk = {}; 
Fk = {};
Gk = {}; 

%for debug, view method
%method = []

% introduce x for computation in each iteration, let x = x0 for now.
x = x0;

display(['iter    ', 'xmin(1)     xmin(2)', '        fxmin      ',' grad(1)      grad(2)    ','nF     ','nG     ',' nH     '])

varNames = ["iter","xmin(1)","xmin2","fmin","grad(1)","grad(2)","nF","nG","nH"];
varTypes = ["double","double","double","double","double","double","double","double","double"];
sz = [1 length(varNames)];
temps = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);


%Start Looping
for i = 1:itmax
    % compute function values, gradients and Hessians at x0
    [f_x,grad,hess] = FunctionName(x, 3);
    nF = nF+1; nG = nG+1; nH = nH+1; % update number of calculations
    
    % Store values
    Xk{i} = x;
    Fk{i} = f_x;
    Gk{i} = grad;
    

    % Find search direction for Newton, In my opinion, there might be times
    % that H is singular, in such cases, we should fall back to steepest  descent
    try
        s = -hess\grad;
        %method = 'Newton'
    catch
        s = -grad;
        %method = 'Steep'
    end

    % check for descent property. The search direction is said to satisfy
    % the descent property if dot(s,gradient) < 0 switch to steepest descent
    if dot(s,grad) >= 0
        s = - grad;
        %method = 'Steep'
    end
    
    % now we identified the method we need to use (Newton or Steepest).
    % let's find the optimal step-size by using golden section method. With
    % stepsize = 2^delta searching for interval, just like HW1
    for delta = -10:10
        [x_k1,f_alp,flag,nF_gold,nG_gold] = golden_func(FunctionName,x,x+(2^delta)*s,e_rel, e_abs,s,itmax);
        % if success, break the loop out. If not, keep running till delta =
        % 100
        nF = nF + nF_gold ; nG = nG + nG_gold; %update
        if flag == 0
            break;
        end
    end
    

    if flag == -999
         disp("Fail to find optimal step-size")
         IFLAG = -999; % Fail flag from golden section after looping j
         xmin = x; fmin = f_x; % return value that failed to find step-size
         return;
    end

    % fill table
    temps(i,:) = {i,x(1),x(2),f_x,grad(1),grad(2),nF,nG,nH};


    % Check for stopping criterion, if the condition is satisfied, we can
    % return the function.
    if norm(x_k1-x) < epsilon
        [f_k1,grad1,null1] = FunctionName(x_k1, 2); % compute last f and gradient.
        %update calculations
        nF = nF + 1; nG = nG + 1;
        Xk{i} = x_k1;
        Fk{i} = f_k1;
        Gk{i} = grad1;

        %update table and display
        temps(i+1,:) = {i+1,x_k1(1),x_k1(2),f_k1,grad1(1),grad1(2),nF,nG,nH};
        display(temps);
        
        %set outputs
        xmin = x_k1; 
        fmin = f_k1;
        return;
    end
    

    x = x_k1; % update x for next iteration

end

%if after looping the norm(x_k1-x) is still not < epsilon, then it failed.
IFLAG = -999
xmin = x; fmin = f_x;

end
</sub>
