
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG] = BFGS(@Rosenbrock,[100;120],0.000002,1e-4,0.1,10000);
varNames = ["iter","xmin(1)","xmin2","fmin","grad(1)","grad(2)","nF","nG"];
fprintf('%6s %12s %12s %12s %12s %12s %5 %5',varNames)
for i = 1:length(Fk)
   fprintf('\n %6.2d %12.8f %12.8f %12.5f %12.5f %12.5f',i,Xk{i}(1),Xk{i}(2),Fk{i},Gk{i}(1),Gk{i}(2))
end

%{
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG] = BFGS(@HMB,[10;12],0.000002,1e-4,0.1,10000);

varNames = ["iter","xmin(1)","xmin2","fmin","grad(1)","grad(2)","nF","nG"];
fprintf('\n %6s %12s %12s %12s %12s %12s %5 %5',varNames)
for i = 1:length(Fk)
   fprintf('\n %6.2d %12.8f %12.8f %12.5f %12.5f %12.5f',i,Xk{i}(1),Xk{i}(2),Fk{i},Gk{i}(1),Gk{i}(2))
end
%}


