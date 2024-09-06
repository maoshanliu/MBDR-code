clc;
% clear all;
% close all;   Z:mean:0.074  median:0.0484             B:mean:0.0805
% median:0.0406
clear all, close all
addpath('/home/lms/MKBDRbest1');
load /home/lms/MKBDRbest1/logvar.mat
load /home/lms/MKBDRbest1/logmean.mat
load /home/lms/MKBDRbest1/ORL_32x32.mat

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

 [D,N] = size(X1);
 t1=sqrt(sum(X1.^2,1));
 t2 = repmat(t1,size(X1,1),1);
 XXX1=X1./t2;
% 
 [D1,N1] = size(X2);
 t11=sqrt(sum(X2.^2,1));
 t22 = repmat(t11,size(X2,1),1);
 XXX2=X2./t22;
% 
 X0=fea';
 [D0,N0] = size(X0);
 t10=sqrt(sum(X0.^2,1));
 t20 = repmat(t10,size(X0,1),1);
 XXX0=X0./t20;




%      if(max(X1(:)) <= 1 && min(X1(:)) >= -1)
%         else
%            X1 = X1 - min(X1(:));
%            X1 = X1/(max(X1(:))+1);
%            X1 = 2*(X1 - 0.5); 
%      end
%         
%           if(max(X2(:)) <= 1 && min(X2(:)) >= -1)
%         else
%            X2 = X2 - min(X2(:));
%            X2 = X2/(max(X2(:))+1);
%            X2 = 2*(X2 - 0.5); 
%           end
%           
%         if(max(X0(:)) <= 1 && min(X0(:)) >= -1)
%         else
%            X0 = X0 - min(X0(:));
%            X0 = X0/(max(X0(:))+1);
%            X0= 2*(X0- 0.5); 
%        end
 
     
     
lambda=10
gamma =0.01
tol=1e-10;


% sig = 1;     
% KG1 = rbfKernelMatrix(XXX0, 1*sig);     
% KG2=  rbfKernelMatrix(XXX1, 1*sig);     
   

 KG1 = polynKernelMatrix(XXX0,12,2); 
% % 
 KG2=  polynKernelMatrix(XXX1,12,2); 
% % 
 KG3 = polynKernelMatrix(XXX2,12,2); 

%  for oo=1:400
%        for pp=1:400
%              KG3(oo,pp)= trace(reshape(XXX2(:,oo),[32,32]).*reshape(XXX2(:,pp),[32,32]));    
%        end
%  end
   q1=100
   q2=1
   
   
        
max_rho=1e10;
eta=20;

for ii=1:20
        %X = DataProjection(data(i).X,r);
        
        [D0,F0,Q0,phix0,E0,Y10,rho0] = BDR_solver1(XXX0,KG1,nbcluster,lambda,gamma,q1,q2);
        [D1,F1,Q1,phix1,E1,Y11,rho1] = BDR_solver2(XXX1,KG2,nbcluster,lambda,gamma,q1,q2);
        [D2,F2,Q2,phix2,E2,Y12,rho2] = BDR_solver3(XXX2,KG3,nbcluster,lambda,gamma,q1,q2);
        
        
       
        Dk=0.3333*(D0+D1+D2);
        x1=0;
        x2=0;
        x3=0;
        
        bb1=0.01
        bb2=0.01
        bb3=0.01
        

        for jj=1:500
        if(x1==0)
           [D00,F00,Q00,phix00,E00] = BDR_solver11(XXX0,KG1,nbcluster,lambda,gamma,D0,F0,Q0,Dk,bb1,q1,q2,phix0,E0,rho0,Y10);
           rho0ori=rho0;
           Y10ori=Y10;
        end
        if(x2==0)
           [D11,F11,Q11,phix11,E11] = BDR_solver11(XXX1,KG2,nbcluster,lambda,gamma,D1,F1,Q1,Dk,bb2,q1,q2,phix1,E1,rho1,Y11);
           rh10ori=rho1;
           Y11ori=Y11;
        end
        if(x3==0)
           [D22,F22,Q22,phix22,E22] = BDR_solver11(XXX2,KG3,nbcluster,lambda,gamma,D2,F2,Q2,Dk,bb3,q1,q2,phix2,E2,rho2,Y12);
           rho2ori=rho2;
           Y12ori=Y12;
        end
      
        Dk=(bb1*D00+bb2*D11+bb3*D22)./(bb1+bb2+bb3);
        
        diffZ0 = max(max(abs(D0-D00)));
        diffB0 = max(max(abs(F0-F00)));    
        diffC = max(max(abs(KG1-phix00'*phix00-E00)));
        stopC0 = max([[diffZ0,diffB0],diffC]);
        
        if (stopC0<tol) 
            x1=1;
        else
            x1=0;
        end
        
        diffZ1 = max(max(abs(D1-D11)));
        diffB1 = max(max(abs(F1-F11)));   
        diffC = max(max(abs(KG2-phix11'*phix11-E11)));
        stopC1 = max([[diffZ1,diffB1],diffC]);
        
         if (stopC1< tol) 
            x2=1;
        else
            x2=0;
         end
        
        diffZ2 = max(max(abs(D2-D22)));
        diffB2 = max(max(abs(F2-F22))); 
        diffC = max(max(abs(KG3-phix22'*phix22-E22)));
        stopC2 = max([[diffZ2,diffB2],diffC]);

        
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
        grp = bestMap(gnd,grps);
        [~,nmi,~]=compute_nmi(gnd',grp');
                 
        missrateTot1{1}(ii) = rate;     
        missrateTot1{2}(ii) = nmi;  
       
end

mean(missrateTot1{1})
mean(missrateTot1{2})