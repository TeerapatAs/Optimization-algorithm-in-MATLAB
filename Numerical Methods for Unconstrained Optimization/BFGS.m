function  [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,nG_linesearch,IFLAG] = BFGS(FcnName,x0,epsilon,mu,eta,itmax)

%inputs-----------------------
% FcnName is the function providing the valuse of f(x) and gradient
% x0 is initial point
% epsilon for stoping criterion if norm(gradient_{k+1}) < epsilon
% mu and eta for checking Armijo's and Wolfe's conditions, respectively.
% itmax is maximum iterations

%Outputs-----------------------
% xmin is the BFGS result. fmin is f(xmin)
% IFlag --> if success, IFlag = 0. If failed,  IFlag = -999
% Xk ,Fk, Gk , Lk are arrays that store x_k, function values,gradients and Lambdas at kth iteration, respectively.
% nF,nG are numbers of function and gradient calculations, respectively.

% array for Xk ,Fk, Gk, Lk.
Xk = {}; 
Fk = {};
Gk = {}; 
Lk = {};

% initial values for nF,nG and IFlag.
nF = 0; 
nG = 0;
nG_linesearch = 0;
IFLAG = 0;

% initial values for B is an Identity Matrix.
B = eye(length(x0));

% Same as HW 2, introduce x for computation in each iteration, let x = x0 for now.
x = x0;

% Start Looping :  Algorithm : find sk -> find Lk -> update x_k1(x_{k+1})
% and B_k1 -> check criterion -> loop again
for i = 1:itmax

    % 1. find sk (Bsk = -gk)---------------------
    [f_k,g_k] = FcnName(x,2); nF = nF+1; nG = nG+1;
    s_k = -(B)\(g_k); % find sk by solving Bsk = -gk
    
    % 2. find Lk (That satisfied Armijo's and Wolfe's)---------------------
    [L_k,nF_AW,nG_AW] = AW_linesearch(FcnName,x,s_k,mu,eta); %get L_k
    nF = nF + nF_AW; nG = nG + nG_AW; % update nF,nG
    
    % store values
    Xk{i} = x;
    Fk{i} = f_k;
    Gk{i} = g_k;
    Lk{i} = L_k;


    % 3. update x_k1(x_{k+1}) and B_k1---------------------
    x_k1 = x + L_k*s_k; % compute x_k1
    
    %for B_k1 update
    [f_k1, g_k1] = FcnName(x_k1,2); nF = nF+1; nG = nG+1; %find func values and gradient at x_k1
    y_k = g_k1 - g_k; 
    d_k = L_k*s_k;

    %BFGS's update terms (Source: Wikipedia -> BFGS algorithm)
    update_B = (y_k*y_k')/(y_k'*(d_k)) - (B*(d_k)*(d_k)'*B')/((d_k)'*B*(d_k));
    B_k1 = B + update_B; % compute B_k1

    % After this, Let's update xk->xk+1 , Bk->Bk+1.
    x = x_k1; B = B_k1;

    % 4. Finally, Check Stopping Criterion (gradient_{k+1} < epsilon)----------------
    if norm(g_k1) < epsilon
        % We can stop looping. Let's add the updated values to arrays before
        % terminate.
        Xk{i+1} = x; xmin = x;
        Fk{i+1} = f_k1; fmin = f_k1;
        Gk{i+1} = g_k1;
        Lk{i+1} = L_k;
        return;
    end
end
    
%if after looping the norm(g_k1) is still not less than epsilon, then it failed.
IFLAG = -999;
xmin = x;
fmin = f_k1;
disp("This is an unsuccessful solution. Change some parameter values e.g. epsilon,x0.")
return;

end