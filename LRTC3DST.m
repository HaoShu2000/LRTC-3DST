
function [RMSE,MAE,X] = LRTC3DST(missingway,missingrate,traffictensor,D1,D2,D3,alpha,nulclearnorm,r,rho,mu)
%% default paremeters setting
X0=traffictensor;
dim = size(X0);
switch missingway
      case 'Random'%%各个维度随机丢失
            Pomega=round(rand(dim(1),dim(2),dim(3)) + 0.5 - missingrate);
      case 'Non-random2'%%在每个时间点丢失情况一样（Non-random，tuba2）
            A=round(rand(dim(1), dim(3)) + 0.5 - missingrate);
            B= kron(A,ones(1,dim(2)));
            Pomega=reshape(B,[dim(1),dim(2),dim(3)]);
      case 'Non-random1'   %%在每个传感器丢失情况一样（tuba1）
            A=round(rand(dim(3), dim(2)) + 0.5 - missingrate);
            B= kron(A,ones(dim(1),1));
            Pomega=reshape(B,[dim(1),dim(2),dim(3)]);
      case 'Non-random3'%%在每天丢失数据情况完全一样（tuba3）
            A=round(rand(dim(1), dim(2)) + 0.5 - missingrate);
            Pomega=repmat(A,[1 1 dim(3)]) ;
      case 'Non-random_1'%%所有传感器在某一段连续时间全部失去数据(Blackout)
            A=round(rand(24*dim(3),1) + 0.5 - missingrate);
            A = repmat(A,1,dim(2)/24); 
            A = (reshape(A.', [], 1))';
            A= repmat(A,dim(1),1);
            Pomega=Fold(A,dim,1);

end
Pomegac=1-Pomega;
tol        = 1e-3; 
max_iter   = 100;
max_mu     = 1e10;
detail     = 1;

%% variables initialization
X            = Pomega.*X0+1/(1-missingrate)*mean(Pomega.*X0,'all')*Pomegac;

PMX          = Pomega.*X0+1/(1-missingrate)*mean(Pomega.*X0,'all')*Pomegac;

for k = 1:3
    G{k}     = porder_diff(X,k); 
end
K            = zeros(dim);
for k = 1:3
   M{k}      = zeros(dim); 
end
N            = zeros(dim);
%% FFT setting
 DD{1}=D1;
 DD{2}=D2;
 DD{3}=D3;

 
 for i=1:3
     [V{i},D{i}]=eig(DD{i}'*DD{i});
 end
 
 T=zeros(dim);
 for i=1:3
     T=T+nmodeproduct(ones(dim),D{i},i);
 end

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = K;
    %% Update X -- solve TV by FFT 
    H = zeros(dim);
    for i=1:3
        H=H+nmodeproduct(G{i}-M{i}/mu,DD{i}',i);
    end
    MAE=(1/(prod(dim)*missingrate))*sum(abs(Pomegac.*X0-Pomegac.*X),'all');
    RMSE=1/sqrt(prod(dim)*missingrate)*sqrt(sum((Pomegac.*X0-Pomegac.*X).^2,'all'));
    X= PMX-K+N/mu+H;
    for i=1:3
      X = nmodeproduct(X,V{i}',i); 
    end
    X=(X./(1+T));
    for i=1:3
      X = nmodeproduct(X,V{i},i); 
    end
    
 
    
 
    %% Updata Gi -- proximal operator of TNN
 for j=1:3
   switch nulclearnorm 
      case 'TrNN'
        [G{j}] = prox_TTNN(nmodeproduct(X,DD{j},j)+M{j}/mu,alpha(j)/mu ,r); 
   end
 end
 

      
  

 
    %% Update K 
    K          = PMX-X+N/mu;
    K          = Pomegac.*K;
    
    %% Stop criterion
    dY   =  PMX-X-K;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-K(:)));
    chg  = max([chgX chgE max(abs(dY(:)))]);
    if chg < tol
         break;
    end 
  
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            err = norm(dY(:),'fro');
            disp(['iter= ' num2str(iter) ', mu=' num2str(mu) ...
                   ', chg=' num2str(chg) ...
                     ', err=' num2str(err) ',  MAE=' num2str(MAE) ...
                     ',  RMSE=' num2str( RMSE)]); 
        end
    end   
    
    %% Update mulipliers: Lambda, Gamma, and mu
    for i=1:3
      M{i} = M{i}+mu*(nmodeproduct(X,DD{i},i)-G{i});
    end
    N = N+mu*dY;
    mu = min(rho*mu,max_mu);
end
end