clc;
% clear all;
% close all;   Z:mean:0.074  median:0.0484             B:mean:0.0805
% median:0.0406
clear all, close all
addpath('/share/home2/ad9145/newMKBDRtext');
load /share/home2/ad9145/newMKBDRtext/sublabel.mat
load /share/home2/ad9145/newMKBDRtext/oriimage.mat
load /share/home2/ad9145/newMKBDRtext/logmean.mat
load /share/home2/ad9145/newMKBDRtext/logvar.mat

numclass=11;
numdata=60;
nbcluster=numclass;
gnd=sublabel';
%gnd=gnd(1:numclass*72,1);
%fea=fea(1:numclass*72,:);

X1=zeros(48*42,numclass*numdata);
X2=zeros(48*48,numclass*numdata);
%logmean=llmean;
%logvar=llvar;

for ii=1:numclass*numdata
    temp1=logmean(:,:,ii);
    temp2=logvar(:,:,ii);
    temp2=logm(temp2);
    temp=abs(temp2);
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

X0=double(oriimage);
[D0,N0] = size(X0);
t10=sqrt(sum(X0.^2,1));
t20 = repmat(t10,size(X0,1),1);
XXX0=X0./t20;

lambda=[10];
gamma =[100,10,1,0.1,0.01,0.001,0.0001];
tol=1e-6;

 bb1=0.1
 bb2=0.1
 bb3=0.1
        
 for oo=1:1
     for qq=1:7
            lambda(oo)
            gamma(qq)
for ii=1:3
        %X = DataProjection(data(i).X,r);
        
        [D0,F0,Q0] = BDR_solver(XXX0,nbcluster,lambda(oo),gamma(qq));
        
        [D1,F1,Q1] = BDR_solver(XXX1,nbcluster,lambda(oo),gamma(qq));
        [D2,F2,Q2] = BDR_solver(XXX2,nbcluster,lambda(oo),gamma(qq));
        
        Dk=(bb1*D0+bb2*D1+bb3*D2)./(bb1+bb2+bb3);
        
        x1=0;
        x2=0;
        x3=0;
        
       
        
        for jj=1:1500
           
        if(x1==0)
            [D00,F00,Q00] = BDR_solver1(XXX0,nbcluster,lambda(oo),gamma(qq),D0,F0,Q0,Dk,bb1);
        end
        if(x2==0)
           [D11,F11,Q11] = BDR_solver1(XXX1,nbcluster,lambda(oo),gamma(qq),D1,F1,Q1,Dk,bb2);
        end
        if(x3==0)
           [D22,F22,Q22] = BDR_solver1(XXX2,nbcluster,lambda(oo),gamma(qq),D2,F2,Q2,Dk,bb3);
        end
        
       Dk=(bb1*D00+bb2*D11+bb3*D22)./(bb1+bb2+bb3);
        
        diffZ0 = max(max(abs(D0-D00)));
        diffB0 = max(max(abs(F0-F00)));    
        stopC0 = max([diffZ0,diffB0]);
        
        if (stopC0<tol) 
            x1=1;
        else
            x1=0;
        end
        
        diffZ1 = max(max(abs(D1-D11)));
        diffB1 = max(max(abs(F1-F11)));    
        stopC1 = max([diffZ1,diffB1]);
        
         if (stopC1< tol) 
            x2=1;
        else
            x2=0;
         end
        
        diffZ2 = max(max(abs(D2-D22)));
        diffB2 = max(max(abs(F2-F22))); 
        stopC2 = max([diffZ2,diffB2]);

        
         if (stopC2< tol) 
            x3=1;
        else
            x3=0;
         end
         
        D0=D00;
        F0=F00;
        Q0=Q00;
        
        D1=D11;
        F1=F11;
        Q1=Q11;
        
        D2=D22;
        F2=F22;
        Q2=Q22;
        
        end
        
        rho1=0.6;
        fprintf('fd')
        CKSym = BuildAdjacency(thrC(Dk,rho1));
        grps  = SpectralClustering(CKSym,numclass);
        rate  = 1 - evalAccuracyHungarian(grps,gnd);
        fprintf('f')
        grp = bestMap(gnd,grps);
        [~,nmi,~]=compute_nmi(gnd',grp');
                 
        missrateTot1{1}(ii) = rate;     
        missrateTot1{2}(ii) = nmi;  
end

fprintf('fdqqd')
mean(missrateTot1{1})
median(missrateTot1{1})
mean(missrateTot1{2})
median(missrateTot1{2})


     end
      
 end