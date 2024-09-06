clc;
% clear all;
% close all;   Z:mean:0.074  median:0.0484             B:mean:0.0805
% median:0.0406
clear all, close all
addpath('/share/home2/ad9145/abilitionstudy/MALKBDRorl');
load /share/home2/ad9145/abilitionstudy/MALKBDRorl/logvar.mat
load /share/home2/ad9145/abilitionstudy/MALKBDRorl/logmean.mat
load /share/home2/ad9145/abilitionstudy/MALKBDRorl/ORL_32x32.mat


numclass=40;
numdata=10;
nbcluster=numclass;

%gnd=gnd(1:numclass*72,1);
%fea=fea(1:numclass*72,:);

X1=zeros(32*32,numclass*numdata);
X2=zeros(32*32,numclass*numdata);
% 
% nbcluster=20;
% gnd=gnd(1:20*72,1);
% fea=fea(1:20*72,:);
% 
% X1=zeros(32*32,20*72);
% X2=zeros(32*32,20*72);
% logmean=llmean;
% logvar=llvar;
for ii=1:numclass*numdata
    temp1=logmean(:,:,ii);
    temp2=logvar(:,:,ii);
    temp2=logm(temp2);
    temp2=abs(temp2);
    X1(:,ii)=temp1(:);
    X2(:,ii)=temp2(:);
end

[D,N] = size(X1);
t1=sqrt(sum(X1.^2,1));
t2 = repmat(t1,size(X1,1),1);
XXX1=X1./t2;

[D1,N1] = size(X2);
t11=sqrt(sum(X2.^2,1));
t22 = repmat(t11,size(X2,1),1);
XXX2=X2./t22;

X0=double(fea');
[D0,N0] = size(X0);
t10=sqrt(sum(X0.^2,1));
t20 = repmat(t10,size(X0,1),1);
XXX0=X0./t20;

lambda=100;
gamma =0.001;
tol=1e-6;

 bb1=0.1
 bb2=0.1
 bb3=0.1
        

 
q1=10000
q2=[1000,100,10,1,0.1,0.01,0,0];
max_rho=1e10;
eta=10;

KG1 = polynKernelMatrix(XXX0,0,1); 

KG2=  polynKernelMatrix(XXX1,0,1); 

KG3 = polynKernelMatrix(XXX2,0,1); 
d1=[];
d2=[];
d3=[];
d4=[];
d5=[];
d6=[];
d7=[];
d8=[];
d9=[];
d10=[];
for ooo=1:1
for ii=1:1

        %X = DataProjection(data(i).X,r);
        
        [D0,F0,Q0,phix0,E0,Y10,rho0] = BDR_solver1(XXX0,KG1,nbcluster,lambda,gamma,q1,q2(ooo));
        [D1,F1,Q1,phix1,E1,Y11,rho1] = BDR_solver2(XXX1,KG2,nbcluster,lambda,gamma,q1,q2(ooo));
        [D2,F2,Q2,phix2,E2,Y12,rho2] = BDR_solver3(XXX2,KG3,nbcluster,lambda,gamma,q1,q2(ooo));
        
        
       
        Dk=0.3333*(D0+D1+D2);
        x1=0;
        x2=0;
        x3=0;
        
   
 
        for jj=1:1500
        if(x1==0)
           [D00,F00,Q00,phix00,E00] = BDR_solver11(XXX0,KG1,nbcluster,lambda,gamma,D0,F0,Q0,Dk,bb1,q1,q2(ooo),phix0,E0,rho0,Y10);
           rho0ori=rho0;
           Y10ori=Y10;
        end
        if(x2==0)
           [D11,F11,Q11,phix11,E11] = BDR_solver11(XXX1,KG2,nbcluster,lambda,gamma,D1,F1,Q1,Dk,bb2,q1,q2(ooo),phix1,E1,rho1,Y11);
           rh10ori=rho1;
           Y11ori=Y11;
        end
        if(x3==0)
           [D22,F22,Q22,phix22,E22] = BDR_solver11(XXX2,KG3,nbcluster,lambda,gamma,D2,F2,Q2,Dk,bb3,q1,q2(ooo),phix2,E2,rho2,Y12);
           rho2ori=rho2;
           Y12ori=Y12;
        end
      
        Dktemp=(bb1*D00+bb2*D11+bb3*D22)./(bb1+bb2+bb3);
        d1=[d1,max(max(abs(Dktemp-Dk)))];
        Dk=Dktemp;
        
        diffZ0 = max(max(abs(D0-D00)));
        diffB0 = max(max(abs(F0-F00)));    
        diffC = max(max(abs(KG1-phix00'*phix00-E00)));
        stopC0 = max([diffZ0,diffB0,diffC]);
         d2=[d2,diffZ0];
         d3=[d3,diffB0];
         d4=[d4,diffC];
        if (stopC0<tol) 
            x1=1;
        else
            x1=0;
        end
        
        diffZ1 = max(max(abs(D1-D11)));
        diffB1 = max(max(abs(F1-F11)));   
        diffC = max(max(abs(KG2-phix11'*phix11-E11)));
        stopC1 = max([diffZ1,diffB1,diffC]);
        
         d5=[d5,diffZ1];
         d6=[d6,diffB1];
         d7=[d7,diffC];
         
         if (stopC1< tol) 
            x2=1;
        else
            x2=0;
         end
        
        diffZ2 = max(max(abs(D2-D22)));
        diffB2 = max(max(abs(F2-F22))); 
        diffC = max(max(abs(KG3-phix22'*phix22-E22)));
        stopC2 = max([diffZ2,diffB2,diffC]);
         d8=[d8,diffZ2];
         d9=[d9,diffB2];
         d10=[d10,diffC];
        
         if (stopC2< tol) 
            x3=1;
        else
            x3=0;
         end
         
        D0=D00;
        F0=F00;
        Q0=Q00;
        E0=E00;
        phix0=phix00;
        
        D1=D11;
        F1=F11;
        Q1=Q11;
        E1=E11;
        phix1=phix11;
        
        D2=D22;
        F2=F22;
        Q2=Q22;
        E2=E22;
        phix2=phix22;
        
        
       Y10 = Y10ori + rho0ori*(KG1-phix0'*phix0-E0);
       rho0 = min(max_rho,rho0ori*eta);
     
        
       Y11 = Y11ori + rh10ori*(KG2-phix1'*phix1-E1);
       rho1 = min(max_rho,rh10ori*eta);
     
        
       Y12 = Y12ori + rho2ori*(KG3-phix2'*phix2-E2);
       rho2 = min(max_rho,rho2ori*eta);
        end
        
        rho1=0.6;
        CKSym = BuildAdjacency(thrC(Dk,rho1));
        grps  = SpectralClustering(CKSym,nbcluster);
        rate  = 1 - evalAccuracyHungarian(grps,gnd);
        rate
        grp = bestMap(gnd,grps);
    
        [~,nmi,~]=compute_nmi(gnd',grp');
            nmi    
        missrateTot1{1}(ii) = rate;     
        missrateTot1{2}(ii) = nmi;  
        %tt=strcat('orl',num2str(ooo));
      %  save(strcat(tt,'.mat'),'Dk')
        
end


save d1.mat d1
save d2.mat d2
save d3.mat d3
save d4.mat d4
save d5.mat d5
save d6.mat d6
save d7.mat d7
save d8.mat d8
save d9.mat d9
save d10.mat d10
mean(missrateTot1{1})
median(missrateTot1{1})

mean(missrateTot1{2})
median(missrateTot1{2})
end
