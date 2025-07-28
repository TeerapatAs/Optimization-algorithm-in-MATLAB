% FcnName,x0,epsilon,mu,eta,itmax,option
epsilon = 1e-6;
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(@Rosenbrock,[-1.2;1],epsilon,1e-4,0.25,500,1);
varNames = ["iter","xmin(1)","xmin2","fmin","grad(1)","grad(2)","nF","nG"];
fprintf('%6s %12s %12s %12s %12s %12s %5 %5',varNames)
for i = 1:length(Fk)
   fprintf('\n %6.2d %12.8f %12.8f %12.5f %12.5f %12.5f %6.1d',i,Xk{i}(1),Xk{i}(2),Fk{i},Gk{i}(1),Gk{i}(2),nReset{i})
end
%{
[xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG,nReset] = CG(@HMB,[-1.2;1],epsilon,1e-4,0.25,500,1);
varNames = ["iter","xmin(1)","xmin2","fmin","grad(1)","grad(2)","nF","nG"];
fprintf('%6s %12s %12s %12s %12s %12s %5 %5',varNames)
for i = 1:length(Fk)
   fprintf('\n %6.2d %12.8f %12.8f %12.5f %12.5f %12.5f %6.1d',i,Xk{i}(1),Xk{i}(2),Fk{i},Gk{i}(1),Gk{i}(2),nReset{i})
end
%}