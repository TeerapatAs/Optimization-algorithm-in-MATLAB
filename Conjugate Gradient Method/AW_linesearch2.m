function [lambda,nF,nG] = AW_linesearch2(FcnName, x0, s, mu, n)

% Armijo's-Wolfe's Line Search

%inputs
% FcnName is the function providing the valuse of f(x), gradient
% x0 is initial point
% s is search direction
% n is eta for checking Wolfe's condition 
% mu is for checking Armijo's condition

%outputs
% lambda is a step-size that satisfies Armijo's and Wolfe's condition.
% nF,nG are numbers of function and gradient calculations, respectively.

% initial values
x = x0; %initial x
lambda = 1; %Initial step size
nF = 0; nG = 0;


% initial alpha and beta
a = 2; b = 0.5;

% Start Looping
% Giving States k logic----------------------------------------------------
k = 0; % Initial k
% if k = 0 ; We need to check Armijo's
% if k = 1 ; We have already checked Armijo's. We can directly go to Wolfe's.
%--------------------------------------------------------------------------

[fk,gk] = FcnName(x,2); % calculate f_values and gradient of current x
nF = nF + 1; nG = nG + 1;
while true

    xk1 = x + lambda*s; % get x_{k+1}
    [fk1,gk1] = FcnName(xk1,2);  % calculate f_values and gradient of current updated x
    nF = nF + 1; nG = nG + 1;

    if k == 0
        % If k = 0, we're checking Armijo's conditions right now.

        %For easier debugging, I define these variables.
        Arm = fk + mu*lambda*dot(s,gk);
       
        if fk1 > Arm
            % If Armijo's fails, we should reduce the step size (lambda) by update
            % lambda = beta * lambda : While lambda is 0.5.
            lambda = b*lambda;
        else
            % If Armojo's passes, we can go check Wolfe's condition. Let's
            % update k to k = 1.
            k = 1;     
        end
    end

    %elseif k == 1
    if k == 1
        % If k = 1, we're checking Wolfe's conditions right now.
        
        %For easier debugging, I declear these variables.
        Wolfe = -n*dot(s,gk);
        Grad = dot(s,gk1); % gradient 

        if abs(Grad) > Wolfe
            % If Wolfe's fails, we need to check cases of slope of function
            % at the current step size (lambda) that there is positive or
            % negative.
            if dot(s,gk1) > 0 
                % If slope of function at lambda is positive. Then, we can
                % assume that the lambda that satisfies Wolfe's condition
                % lies in between 0 and current lambda.

                % Then, we can do any type of line search to find new lambda.
                % Let's use Golden Section Method modified from HW2.
                a_g = 0; b_g = lambda; %interval
                itmax = 1000;
                % find new lambda
                [lambda,~,IFLAG, nF_g, nG_g]  = golden_func2(FcnName,x,s,gk,a_g,b_g,n,itmax);
                nF = nF + nF_g; nG = nG + nG_g;

                if IFLAG == -999 % if golden section can't find the step-size. Send error.(This 
                    %should not be the case.)
                    disp(['Golden Section Error. Please consider the new itmax parameter or' ...
                        ' new epsilon']);
                end

            else
                % If slope of function at lambda is negative. Then, we can
                % assume that the lambda that satisfy Wolfe's condition
                % lies on right side of the current lambda.

                % Let's update new lambda to far right. Note that a = 2.
                lambda = lambda * a;     
            end

        else
            % If Wolfe's passes, It's done! We can exit.
            return         
        end
    end

end
end

