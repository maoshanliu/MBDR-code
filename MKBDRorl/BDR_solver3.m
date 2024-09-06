function [B,Z,W] = BDR_solver3(X,KG,k,lambda,gamma,a,b) 

% min_{Z,B,W} 0.5*||X-XZ||_F^2+0.5*lambda*||Z-B||_F^2+gamma*<Diag(B1)-B,W>
% s.t. diag(B)=0, B>=0, B=B^T, 
%      0<=W<=I, Tr(W)=k.


% kType='pol';
% if(strcmpi(kType,'pol')) 
%           a = 12;
%           b = 2;
%           KG = polynKernelMatrix(X,a,b); % hopkins 2.2 3   %two frames 3.4 2  % yale face 12 2    
% end



%         if(max(X(:)) <= 1 && min(X(:)) >= -1)
%         else
%            X = X - min(X(:));
%            X = X/(max(X(:))+1);
%            X = 2*(X - 0.5); 
%         end

   
%           sig = 0.8;
%           KG = rbfKernelMatrix(X, 1*sig);         
        
     
      
n = size(X,2);
tol = 1e-10;
maxIter = 500; 
one = ones(n,1);
%XtX = X'*X;

I = eye(n);
invXtXI = I/(KG+lambda*I);
%invXtXI = I/(KG+lambda*I);
gammaoverlambda = gamma/lambda;
Z = zeros(n);
W = Z;
B = Z;
L = diag(B*one)-B;
iter = 0;

while iter < maxIter
    iter = iter + 1;        
      
    % update Z
    Zk = Z;
    Z = invXtXI*(KG+lambda*B);
    
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
    stopC = max([diffZ,diffB]);

    if stopC < tol 
        break;
    end
end
iter


