clc;
% clear all;
% close all;   Z:mean:0.074  median:0.0484             B:mean:0.0805
% median:0.0406
clear all, close all
addpath('/home/lms/MKBDR');
load /home/lms/MKBDR/logvar.mat
load /home/lms/MKBDR/logmean.mat
load /home/lms/MKBDR/ORL_32x32.mat

nbcluster=40;

X1=zeros(32*32,400);
X2=zeros(32*32,400);

for ii=1:400
    temp1=logmean(:,:,ii);
    temp2=logvar(:,:,ii);
    temp2=logm(temp2);
    X1(:,ii)=temp1(:);
    X2(:,ii)=temp2(:);
end

%  [D,N] = size(X1);
%  t1=sqrt(sum(X1.^2,1));
%  t2 = repmat(t1,size(X1,1),1);
%  XXX1=X1./t2;
% 
%  [D1,N1] = size(X2);
%  t11=sqrt(sum(X2.^2,1));
%  t22 = repmat(t11,size(X2,1),1);
%  XXX2=X2./t22;
% 
 X0=fea';
%  [D0,N0] = size(X0);
%  t10=sqrt(sum(X0.^2,1));
%  t20 = repmat(t10,size(X0,1),1);
%  XXX0=X0./t20;




     if(max(X1(:)) <= 1 && min(X1(:)) >= -1)
        else
           X1 = X1 - min(X1(:));
           X1 = X1/(max(X1(:))+1);
           X1 = 2*(X1 - 0.5); 
     end
        
          if(max(X2(:)) <= 1 && min(X2(:)) >= -1)
        else
           X2 = X2 - min(X2(:));
           X2 = X2/(max(X2(:))+1);
           X2 = 2*(X2 - 0.5); 
          end
          
        if(max(X0(:)) <= 1 && min(X0(:)) >= -1)
        else
           X0 = X0 - min(X0(:));
           X0 = X0/(max(X0(:))+1);
           X0= 2*(X0- 0.5); 
       end
 
     
     
lambda=1
gamma =0.0001

for ii=1:20
        %X = DataProjection(data(i).X,r);
        
        [D0,F0,Q0] = BDR_solver00(X0,nbcluster,lambda,gamma);
        [D1,F1,Q1] = BDR_solver00(X1,nbcluster,lambda,gamma);
        [D2,F2,Q2] = BDR_solver00(X2,nbcluster,lambda,gamma);
        
        Dk=0.3333*(D0+D1+D2);
       
        rho1=0.6;
        CKSym = BuildAdjacency(thrC(Dk,rho1));
        grps  = SpectralClustering(CKSym,nbcluster);
        rate  = 1 - evalAccuracyHungarian(grps,gnd);
        missrateTot1{1}(ii) = rate;     
end

mean(missrateTot1{1})
median(missrateTot1{1})