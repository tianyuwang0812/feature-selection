function [lower_approximation] = lower_approximation(covering,d,beta)
[m,n]=size(covering);%m个元素，n个覆盖（属性）
% beta=0.7;
% d={[1,2],[3,4]};

z=size(d,2);

%求邻域系统
neighborhood_system=cell(1,m);
for i=1:m
    for j=1:n
        if covering(i,j)>=beta
            neighborhood_system{i}=[neighborhood_system{i},j];
        end
    end
end

%求邻域
neighborhood=[];
for i=1:m
    neighborhood(i,:)=min(covering(:,neighborhood_system{i}),[],2);%第i行是第i个元素的邻域
end

%求相似关系
relation=1-neighborhood;

%求每一类的下近似
for i=1:z
    lower(:,i)=min(relation(:,d{i}),[],2);
end
  
[~,lower_class]=min(lower,[],2);

for i=1:m
    lower_approximation_class=lower;
    lower_approximation_class(:,lower_class(i))=[];
    lower_approximation(i)=min(lower_approximation_class(i,:));
end
end