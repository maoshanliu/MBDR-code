 results=[
               39.92, 37.96,41.81,53.38;   
               37.73, 37.77,40.19,65.96;   
               37.35, 37.62,55.19,74.19;
               38.19, 39.04,64.38,75.54
               38.54, 37.04,66.58,78.62;
               38.42, 36.88,65.58,79.65;
               37.35, 37.04,65.81,80.19;
               38.19, 38.19,65.92,78.31;
                ];
%  
%       results1=[
%                0.7011, 0.7089,0.6937,0.6124;   
%                0.7097, 0.7091,0.6939,0.4328;   
%                0.7152, 0.7063,0.5565,0.3233;
%                0.7111, 0.7010,0.4461,0.3171;
%                0.7069, 0.7142,0.4320,0.2226;
%                0.7086, 0.7147,0.4366,0.2159;
%                0.7152, 0.7136,0.4341,0.2109;
%                0.7111, 0.7093,0.4345,0.2288;
%                 ];
%      
% results=zeros(8,4);
% results(:,1)=results1(:,4);
% results(:,2)=results1(:,3);
% results(:,3)=results1(:,2);
% results(:,4)=results1(:,1);
bar3(results)
set(gca,'xticklabel',{'10000','1000','100','10'})
set(gca,'yticklabel',{'1000','100','10','1','0.1','0.01','0.001','0.0001'})
 
xlabel('q_1');ylabel('q_2');zlabel('Mean-CE')


%xlabel('q_1');ylabel('q_2');zlabel('Mean-CE(%)')