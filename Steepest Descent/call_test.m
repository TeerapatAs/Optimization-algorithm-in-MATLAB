f = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;

x = [1,1]';

H = approx_hessian(f,x)

i = 1:0
H(i)

f = 3 - H(i)

H_1 = [-1 2 4;
    2 -3 6;
    4 6 22]

[L,D] = ldl_factor(H_1)

H_ = L*D*L'

[l,d] = ldl(H_1)

H_5 = l*d*l'

a = FunctionName(x,2)

% Define the test function
test_func = @(x) (x - 2)^2 + 3;

% Interval [a, b]
a = 0;
b = 5;

% Tolerances for stopping criteria
eps_rel = 1e-5;
eps_abs = 1e-5;

% Search direction (for 1D optimization, the direction is just 1)
s = 1;

% Maximum number of iterations
itmax = 1000;

% Call the golden section function
[xmin, fmin, IFLAG, nF, nG] = golden_func(test_func, a, b, eps_rel, eps_abs, s, itmax);

% Display the results
fprintf('xmin: %.5f\n', xmin);
fprintf('fmin: %.5f\n', fmin);
fprintf('IFLAG: %d\n', IFLAG);
fprintf('Number of function evaluations: %d\n', nF);
fprintf('Number of gradient evaluations: %d\n', nG);
