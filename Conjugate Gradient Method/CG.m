function [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(FcnName,x0,epsilon,mu,eta,itmax,option)
%inputs-----------------------
% x0 is initial point
% epsilon for stoping criterion if norm(gradient_{k+1}) < epsilon
% mu and eta for checking Armijo's and Wolfe's conditions, respectively.
% itmax is maximum iterations

%Outputs-----------------------
% xmin is the BFGS result. fmin is f(xmin)
% IFlag --> if success, IFlag = 0. If failed,  IFlag = -999
% Xk ,Fk, Gk , Lk are arrays that store x_k, function values,gradients and Lambdas at kth iteration, respectively.
% nF,nG are numbers of function and gradient calculations, respectively.
% nReset is array that store

% array for Xk ,Fk, Gk, Lk, nReset
Xk = {}; 
Fk = {};
Gk = {}; 
Lk = {};
nReset = {};

%initial value for nF,nG,nReset and iFlag
nF = 0; nG = 0;
IFLAG = 0;


% Same as HW 3, introduce x for computation in every iteration, let x = x0 for now.
x = x0;

% Step1. compute initial search direction------------------------------------
[f,g] = FcnName(x,2); nF = nF + 1; nG = nG + 1;
s = -g; % Initial search direction

% Start Looping
for i = 1:itmax
    % step2. Perform line search to obtain the step length------------------
    [lambda,nF_AW,nG_AW] = AW_linesearch2(FcnName,x,s,mu,eta); %lambda is step size
    nF = nF + nF_AW; nG = nG+nG_AW;

    % store values---------
    Xk{i} = x;
    Fk{i} = f;
    Gk{i} = g;
    Lk{i+1} = lambda;

    
    % Step 3. Set x_k+1 = x_k + lambda*s
    x1 = x + lambda*s;

    % Step 4. Check termination criterion ( norm(x1-x) < epsilon )-----------------------
    if norm(x1-x) < epsilon
        % If passes, It's done!
        xmin = x1;
        fmin = FcnName(x1,1); nF = nF + 1;
        Xk{i+1} = x1;
        Fk{i+1} = fmin;
        % There's no need find gradient and lambda of x_k+1 after termination criterion
        % so I assume that I should old valuet from x_k that came before x_k+1
        Gk{i+1} = g; Lk{i+1} = lambda; nReset{i} = 0; nReset{i+1} = [];
        return 
    else
        % Else, We should continue ....
        % First, Find the beta from different methods. to compute s_k+1 = -g_k+1 + beta_k*s_k
        [f,g1] = FcnName(x1,2); nF = nF + 1; nG = nG+1; %find g_k+1
        if option == 1
            % FR method: beta = (g_k+1)'(g_k+1)/(g_k)'(g_k)
            b = norm(g1)/norm(g);
        elseif option ==2
            % PR method: beta = [(g_k+1)'(g_k+1-g_k)]/[(g_k)'(g_k)]
            b = (norm(g1)-dot(g1,g))/(norm(g));
        else
            disp('Wrong Input of Option!')
        end

        s1 = -g1 + b*s; % completely compute s_k+1
        
        % Then, Set Watch-dog to reset the search direction 
        % either when the angle between the s1 and −g1 is too large

        cos_seta = dot(s1,-g1)/(norm(s1)*norm(-g1));

        if cos_seta < cosd(85) 
            % If angle between the s1 and −g1 is too large but still descent. 
            % set s_k+1 = -g_k+1
            s1 = -g1;
            nReset{i} = 1;
        elseif cos_seta < 0
            % If cos_seta < 0. It means that dot(s1,-g1) < 0, so
            % dot(s1,g1) > 0 , which indidcates that s1
            % can't make f(x_k+1) descend.
            s1 = -g1;
            nReset{i} = 2;
        else
            nReset{i} = 0;
        end
    end
    % Lastly, update s, x and g
    s = s1;
    x = x1;
    g = g1;
end
% if after many iterations, cannot find answer. It fails
disp('Failed. Please consider changing itmax or others inputs.')
xmin = x;
fmin = f;
IFLAG = -999;
return

end
    
    