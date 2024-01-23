function [X] = prox_TTNN(Y,rho,r)
dim=size(Y);
A=[1/3,1/3,1/3];
X=zeros(dim);
for k=1:3
    Yk=Unfold(Y,dim,k);
    [U,S,V] = svd(Yk,'econ');
     diagS = diag(S);
     rho1=rho*ones(size(diagS));
     rho1(1:r)=0;
     diagSS=max(diagS-rho1,0);
     Xk = U*diag(diagSS)*V';
     XK=A(k)*Fold(Xk,dim,k);
     X=X+XK;
end
    
    
