function [B,Z,W] = BDR_solver11(X,KG,k,lambda,gamma,D,F,Q,Dk,bb1) 

% min_{Z,B,W} 0.5*||X-XZ||_F^2+0.5*lambda*||Z-B||_F^2+gamma*<Diag(B1)-B,W>
% s.t. diag(B)=0, B>=0, B=B^T, 
%      0<=W<=I, Tr(W)=k.

n = size(X,2);
%maxIter = 500; 
one = ones(n,1);
%XtX = X'*X;

I = eye(n);
%invXtXI = I/(XtX+lambda*I);
invXtXI = I/(KG+lambda*I);
%gammaoverlambda = gamma/lambda;
Z=F;
W=Q;
B=D;


% update Z
Zk = Z;
Z = invXtXI*(KG+lambda*B);
    
% update B
Bk = B;
%B = (2/3*Dk-lambda*Z+gamma*(repmat(diag(W),1,n)-W))/(2/3-lambda);
B = (2*bb1*Dk-lambda*Z+gamma*(repmat(diag(W),1,n)-W))/(2*bb1-lambda);
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

end



