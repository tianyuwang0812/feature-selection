function [covering,del] = shujuchuli(shuju,deta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d=shuju(:,end);
shuju=shuju(:,1:end-1);
[m,n]=size(shuju);
covering=shuju;
% for i=1:n
%     Max=max(shuju(:,i));
%     Min=min(shuju(:,i));
%     for j=1:m
%         covering(j,i)=(shuju(j,i)-Min)/(Max-Min);
%     end
% end

for i=1:n
    Mean=mean(covering(:,i));
    Std=std(covering(:,i));
    for j=1:m
        covering(j,i)=(atan((covering(j,i)-Mean)/Std^2)/pi)+0.5;
    end
end
covering=round(covering,2);
covering=[covering,d];
del=[];
for i=1:m
    if max(covering(i,1:end-1))<=deta
        del=[del,i];
    end
end
covering(del,:)=[];

end