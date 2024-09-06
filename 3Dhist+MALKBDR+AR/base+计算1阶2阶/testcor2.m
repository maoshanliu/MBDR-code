load ORL_32x32.mat
logvar=zeros(32,32,400);
for ii=1:400
 
    temp=fea(ii,:);
    temp1=reshape(temp,[32,32]);
    result_image = Log_gabor(temp1, 5, 8, 3,2,0.65, 1.5,0);
    temp2=(result_image{1,1}+result_image{1,2}+result_image{1,3}+result_image{1,4}+result_image{1,5}+result_image{1,6}+result_image{1,7}+result_image{1,8}+result_image{2,1}+result_image{2,2}+result_image{2,3}+result_image{2,4}+result_image{2,5}+result_image{2,6}+result_image{2,7}+result_image{2,8}+result_image{3,1}+result_image{3,2}+result_image{3,3}+result_image{3,4}+result_image{3,5}+result_image{3,6}+result_image{3,7}+result_image{3,8}+result_image{4,1}+result_image{4,2}+result_image{4,3}+result_image{4,4}+result_image{4,5}+result_image{4,6}+result_image{4,7}+result_image{4,8}+result_image{5,1}+result_image{5,2}+result_image{5,3}+result_image{5,4}+result_image{5,5}+result_image{5,6}+result_image{5,7}+result_image{5,8})/40;
    temp3=zeros(32,32);
    for kk=1:5
        for ll=1:8
          temp3=temp3+(result_image{kk,ll}-temp2)*(result_image{kk,ll}-temp2)';
        end
    end
    logvar(:,:,ii)=1/39*temp3;
end  


xx=0;