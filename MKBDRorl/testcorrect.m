clc;
% clear all;
% close all;   Z:mean:0.074  median:0.0484             B:mean:0.0805
% median:0.0406
clear all, close all
addpath('/share/home2/ad9145/newMBDRcc');
load /share/home2/ad9145/newMBDRcc/sublabel.mat
load /share/home2/ad9145/newMBDRcc/oriimage.mat
load /share/home2/ad9145/newMBDRcc/logmean.mat
load /share/home2/ad9145/newMBDRcc/logvar.mat

oriimage=oriimage(:,1:6:600);
numclass=10;
numdata=10;
nbcluster=numclass;
gnd=sublabel';
gnd=[ones(1,10),2*ones(1,10),3*ones(1,10),4*ones(1,10),5*ones(1,10),6*ones(1,10),7*ones(1,10),8*ones(1,10),9*ones(1,10),10*ones(1,10)]';
%gnd=gnd(1:numclass*72,1);
%fea=fea(1:numclass*72,:);

X1=zeros(16*16,numclass*numdata);
X2=zeros(16*16,numclass*numdata);
%logmean=llmean;
%logvar=llvar;
jj=1;
for ii=1:6:600
    temp1=logmean(:,:,ii);
    temp2=logvar(:,:,ii);
    temp2=logm(temp2);
    X1(:,jj)=temp1(:);
    X2(:,jj)=temp2(:);
    jj=jj+1;
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

lambda=[100];
gamma =[0.001];
tol=1e-6;

 bb1=0.1
 bb2=0.1
 bb3=0.1
        
 for oo=1:1
     for qq=1:1
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
        CKSym = BuildAdjacency(thrC(Dk,rho1));
        grps  = SpectralClustering(CKSym,numclass);
        rate  = 1 - evalAccuracyHungarian(grps,gnd);
        
        grp = bestMap(gnd,grps);
        [~,nmi,~]=compute_nmi(gnd',grp');
                 
        missrateTot1{1}(ii) = rate;     
        missrateTot1{2}(ii) = nmi;  
end


mean(missrateTot1{1})
median(missrateTot1{1})
mean(missrateTot1{2})
median(missrateTot1{2})


     end
      
 end