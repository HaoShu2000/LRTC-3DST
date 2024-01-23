addpath(genpath('tool'));
addpath(genpath('data'));

disp("LRTC3DST-GuangZhou");
dataRoad = ['data/','GuangZhou'];
load(dataRoad);
traffictensor=GuangZhou;
dim=size(traffictensor);
D1=eye(dim(1));
D2=gammaMatrix(0.98,dim(2),dim(2));
D3=circvet2mat(num2vetT(7,dim(3)));
alpha=[0.1,0.8,0.1];
nulclearnorm='TrNN';
r=5;
rho        = 1.1;
mu         = 1e-6;
missingway='Random';
for missingrate=[0.7]
    [~,~,X7G]=LRTC3DST(missingway,missingrate,traffictensor,D1,D2,D3,alpha,nulclearnorm,r,rho,mu);
    disp(['missingrate=' num2str(missingrate)]); 
end

disp("LRTC3DST-pemsD7M");
dataRoad = ['data/','pemsD7M'];
load(dataRoad);
dataRoad = ['data/','pemsws'];
load(dataRoad);
traffictensor=pemsT;
dim=size(traffictensor);
W=exp(-pemsws/10^4);
W(W<0.5)=0;
W=W-diag(diag(W));
DW=diag(sum(W,2));
LW=(DW-W)/max(sum(W,2));
D1=LW;
D2=gammaMatrix(0.98,dim(2),dim(2));
D3=circvet2mat(num2vetT(2,dim(3)));
alpha=[0.1,0.8,0.1];
nulclearnorm='TrNN';
r=5;
rho        = 1.1;
mu         = 1e-6;
missingway='Random';
for missingrate=0.7
    [~,~,X7C1]=LRTC3DST(missingway,missingrate,traffictensor,D1,D2,D3,alpha,nulclearnorm,r,rho,mu);
    disp(['missingrate=' num2str(missingrate)]); 
end
