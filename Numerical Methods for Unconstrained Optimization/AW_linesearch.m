function [lambda,nF,nG] = AW_linesearch(FcnName, x0, s, mu, n)

% Armijo's-Wolfe's Line Search

%inputs
% FcnName is the function providing the valuse of f(x), gradient
% x0 is initial point
% s is search direction
% n is eta for checking Wolfe's condition 
% mu is for checking Armijo's condition

% intial values
x = x0; %initial x
lambda = 1; %initial step size
a = 0; b = []; % initial interval, Let's b be unknown
nF = 0; nG = 0; %counters

% set f_0 for f(x_k) 
[f0,g0] = FcnName(x,2); nF = nF+1; nG = nG +1;

%i = 0; for debug.

%Start Looping
while true
    %i = i+1; for debug.
    x_k1 = x + lambda*s;
    [f1,g1] = FcnName(x_k1,2); % find f and gradient value of x_k + lambda*s
    nF = nF+1; nG = nG +1;
    
    % for easier debugging, I declear these variables
    %Arm = f0 + mu*lambda*dot(s,g0); % Armijo's condition
    Wolfe = (-n*dot(s,g0)); %Wolfe's condition
    SL = abs(dot(s,g1)); % slope at f1

    if f1 > f0 + mu*lambda*dot(s,g0) % Armijo's fails --> Shift b to left (to lambda)
        b = lambda;
        lambda = (a+b)/2;

    elseif SL > Wolfe % Armijo's passes, Wolfe's fails 
        if dot(s,g1) >= 0 
        % Armijo's passes, Wolfe's fails,and dot(s,g1)>=0, It means that this 
        % iteration's step size is much larger than needed. --> Shift b to left (to lambda)
            b = lambda; %new bracket
            lambda = (a+b)/2; % new step size
            %fprintf('pass arm,but not wolfe %d \n',i);

        else % Armijo's passes, Wolfe's fails,and dot(s,g1)<0, It means that this 
            % interval doesn't have local minimize step size. We should
            % find new bracket by ...
            a = lambda;
            lambda = 2*a;  
        end
    else % All conditions passes. It's done!
        return
    end

end
end