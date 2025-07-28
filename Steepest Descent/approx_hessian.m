function H = approx_hessian(f,x)

%input
%f is a function
%x is a point which we used to compute hessian

n = length(x); %number of variables
H = zeros(n); %initial H
h = 1e-8; %Step size

%Looping
for i = 1:n
    for j = 1:n
        %ei is the i th column of Identity matrix
        ei = zeros(n,1);
        ei(i) = 1;
        %ej is the j th column of Identity matrix
        ej = zeros(n,1);
        ej(j) = 1;
        
        % approximate gradient
        H(i,j) = (f(x + h*ei + h*ej) - f(x + h*ei) - f(x + h*ej) + f(x)) / (h^2);
    end
end
end