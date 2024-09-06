% mSpectralClustering.m
% Description: This is a modification of standard spectral clustering, with an extra step in exception detection and a
% warmstart setting for the k-means step.
%    
%  function [groups, kerN, Q, Q2, dim_kerN, num_spaU, sumD] = mSpectralClustering(CKSym, n, idx0, opt, Psi)
%
%   Inputs:  
%		CKSym - Affinity matrix
%		n - number of clusters
%		idx0 - the previous clustering result, which is used to prepare a warmstart for the next call for k-means
%		opt - if opt=0, then conduct spectral clustering with unnormalized graph Laplacian; otherwise conduct 
%			  spectral clustering with normalized graph Laplacian (i.e. Normalized Cut).
%		Psi - this option is used in CS3Cplus, in which Psi is a matrix which contains the available side-information. 
%
%  Outputs: 
%		groups - obtained membership index matrix
%		kerN - embeddings without unit L2 norm normalization
%		Q - Q = kerN (where Q was used for debug)
%		Q2 - embeddings with unit L2 norm normalization
%		dim_kerN - estimated number of zero eigen values
%		num_spaU - number of detected exceptional eigenvectors, which consist of extremely sparse vector
%		sumD - a vector consists of average distances of data points from their cluster center (returned from k-means)
%
%--------------------------------------------------------------------------
% Modified by Chunguang Li, Nov. 7, 2014.
function [groups, kerN, Q, Q2, dim_kerN, num_spaU, sumD] = mSpectralClustering(CKSym,n, idx0, opt, Psi)
if (nargin ==5)
    sideinfo =1; % use side info
else
    sideinfo =0;
end
if (nargin < 4)
    opt =0; % unnormalization for L
end
if nargin < 3
    warmstart =0;
else
    warmstart =1;
end
warning off;
N = size(CKSym,1);
MAXiter = 1000; % Maximum number of iterations for KMeans 
REPlic = 20; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
if (~opt)
    DN = diag( 1./sqrt(sum(CKSym)+eps) );
    LapN = speye(N) - DN * CKSym * DN;
    [uN,sN,vN] = svd(LapN);
    
    kerN = vN(:,N-n+1:N);
    
    %% Decide if we need to use:
    %      kerN = vN(:,N-n:N);    
    %% 1. Check the eigen-gap
    S = diag(sN);
    tmp_id =zeros(1, N);
    for i =1:N-1
        if (S(i) > eps)
            delta = S(i+1) / S(i);        
        else
            i_star = i;
        end
        if (delta < 1e-4)
            tmp_id(i) =1;
        end
    end
    [~, i_star] =find(tmp_id==1);
    i_star = min(i_star);
    
    dim_kerN =N - i_star;
    if ((N - i_star) >=n)
        disp(['dim ( kerN) is :  ', num2str(dim_kerN),' !!!'])
        kerN = vN(:,i_star:N);    
    end
    
    %% check if there is a sparse eigenvector
    thresh = 1e-3 * mean( mean( abs( kerN ) ) );
    jj_star =0;
    for jj =1:size(kerN, 2)
        tmp = (abs(kerN(:, jj)) > thresh);        
        if ( sum(tmp) < N / (n*10))
            disp('This is a sparse vector in kerN !');
            jj_star =jj_star +1;
        end
    end       
    if (jj_star >0.99)
        kerN = vN(:,N - n - jj_star+1 : N);
    end
    num_spaU =jj_star;
    
    
    Q = kerN; 
    for i = 1:N
        kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    end
    Q2 =kerNS;
    
    
else
    
    %%  Unnormalized Laplacian
    D =diag(sum(CKSym));
    Lap = D - CKSym;
    
    %% SVD
    [~,S,V] = svd(Lap);
    kerN = V(:,N-n+1:N);
    
    %     %% Re-normalization
    %     for i = 1:N
    %         kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    %     end    
    
    
    %% Decide if we need to use:
    %      kerN = vN(:,N-n:N);    
    %% 1. Check the eigen-gap
    S = diag(S);
    tmp_id =zeros(1, N);
    for i =1:N-1
        if (S(i) > eps)
            delta = S(i+1) / S(i);        
        else
            i_star = i;
        end
        if (delta < 1e-4)
            tmp_id(i) =1;
        end
    end
    [~, i_star] =find(tmp_id==1);
    i_star = min(i_star);
    
    dim_kerN =N - i_star;
    if ((N - i_star) >=n)
        disp(['dim ( kerN) is :  ', num2str(dim_kerN),' !!!'])
        kerN = V(:,i_star:N);    
    end
    
    %% check if there is a sparse eigenvector
    thresh = 1e-3 * mean( mean( abs( kerN ) ) );
    jj_star =0;
    for jj =1:size(kerN, 2)
        tmp = (abs(kerN(:, jj)) > thresh);        
        if ( sum(tmp) < N / (n*10))
            disp('This is a sparse vector in kerN !');
            jj_star =jj_star +1;
        end
    end       
    if (jj_star >0.99)
        kerN = V(:,N - n - jj_star+1 : N);
    end
    num_spaU =jj_star;
    
    kerNS =kerN;
    Q =kerN; 
    for i = 1:N
        kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    end
    Q2 =kerNS;
end

if (~sideinfo)
    %warmstart =0;
    %groups = kmeans(kerNS,n,idx0,'Start', 'initial', 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    if (~warmstart) % kmeans_version2
        %[groups, Centroid]  = kmeans2(kerNS,n,'Start', 'sample', 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
        [groups, ~, sumD]  = kmeans(kerNS,n,'Start', 'sample', 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    else
        %[groups, Centroid] = kmeans2(kerNS,n,'Start', 'warmstart', 'initialcentroid',Centroid0, 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
        [groups, ~, sumD]= kmeans(kerNS,n,'Start', 'warmstart', 'initiallabel',idx0, 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    end
    
else
    % start: sideinfo
    %[groups, ~, sumD]= kmeans(kerNS,n,'Start', 'sideinfo', 'psi',Psi, 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
    % start: sideinfo+warmstart
    [groups, ~, sumD]= kmeans(kerNS,n,'Start', 'sideinfowarmstart', 'initiallabel',idx0, 'psi',Psi, 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');    
end