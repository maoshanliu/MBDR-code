function [B,Z,W,phix,E,Y1,rho] = BDR_solver2(X,KG,k,lambda,gamma,q1,q2) 

% min_{Z,B,W} 0.5*||X-XZ||_F^2+0.5*lambda*||Z-B||_F^2+gamma*<Diag(B1)-B,W>
% s.t. diag(B)=0, B>=0, B=B^T, 
%      0<=W<=I, Tr(W)=k.


% kType='pol';
% if(strcmpi(kType,'pol')) 
%           a = 12;
%           b = 2;
%           KG = polynKernelMatrix(X,a,b); % hopkins 2.2 3   %two frames 3.4 2  % yale face 12 2    
% end

%   a = 12;
%   b = 2;
%   KG1 = polynKernelMatrix(X,a,b);

   
%   sig = 0.8;     
%   KG = rbfKernelMatrix(X, 1*sig);         
        
  
  %KG = polynKernelMatrix(X,a,b); 
        
n = size(X,2);
tol = 1e-10;
maxIter = 1; 
one = ones(n,1);
%XtX = X'*X;

I = eye(n);
[U,S,V] = svd(KG);
phix = U*sqrt(S)*V';
K = KG;

rho=1e-8;
max_rho=1e10;
eta=10;


%invXtXI = I/(KG+lambda*I);
%invXtXI = I/(KG+lambda*I);
gammaoverlambda = gamma/lambda;
Z = zeros(n);

E=zeros(n);
Y1=zeros(n);


W = Z;
B = Z;
L = diag(B*one)-B;
iter = 0;

while iter < maxIter
    iter = iter + 1;        
      
    % update Z
    Zk = Z;
    Z = I/(phix'*phix++lambda*I)*(phix'*phix+lambda*B);
    
    % update phix
    tmp = KG - E-0.5./rho*(eye(n)-2*Z'+Z*Z')+Y1/rho;
    phix = solveB1(tmp, rho./q2);

    % update error
    
    tmp1 = KG - phix'*phix + Y1/rho;
    E = max(0,tmp1 - q1/rho)+min(0,tmp1 + q1/rho);

    
    % update B
    Bk = B;
    B = Z-gammaoverlambda*(repmat(diag(W),1,n)-W);
    B = max(0,(B+B')/2);
    B = B-diag(diag(B));
    L = diag(B*one)-B;    
    
    % update W
    [V, D] = eig(L);
    D = diag(D);
    [~, ind] = sort(D);    
    W = V(:,ind(1:k))*V(:,ind(1:k))';
    
%     diffZ = norm(Z-Zk,'fro')/norm(Zk,'fro');
%     diffB = norm(B-Bk,'fro')/norm(Bk,'fro'); 
    diffZ = max(max(abs(Z-Zk)));
    diffB = max(max(abs(B-Bk)));   
    diffC = max(max(abs(KG-phix'*phix-E)));

    stopC = max([diffZ,diffB,diffC]);

    if stopC < tol 
        break;
    end
    
     Y1 = Y1 + rho*(KG-phix'*phix-E);
     rho = min(max_rho,rho*eta);

    
end
iter


