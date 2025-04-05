function [DIS,reduct,lower,d_class,neighborhood_system,neighborhood] = rdrn_fs(shuju,beta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

covering=shuju(:,1:end-1);
[m,n]=size(covering);
%决策类
d=shuju(:,end);

for i=1:m
    d_class(i,:)=(d==d(i));
end
non_class=~d_class;

%求邻域系统
neighborhood_system=(covering>=beta);

%求邻域
neighborhood=[];
for i=1:m
    neighborhood(i,:)=min(covering(:,neighborhood_system(i,:)),[],2);%第i行是第i个元素的邻域
end
 neighborhood(neighborhood>=beta)=1;

%求相似关系
relation=1-neighborhood;
%求下近似
lower=[];
for i=1:m
    lower(i)=min(relation(i,non_class(i,:)));
end

%求每个覆盖元相对分辨关系
DIS=zeros(m,m,n);

for i=1:n
    %yulinyu=find(neighborhood_system(:,i))';
    for j=find(neighborhood_system(:,i))'
        DIS(j,:,i)=(1-covering(:,i)>=lower(j))'&non_class(j,:);
    end
end

DIS=double(DIS);
dis_total=sum(DIS,3);

%判断总体相对分辨关系是否与非类向量相等
DIS_total=(dis_total>=1);
if isequal(non_class,DIS_total)==0
    disp("程序存在错误")
end

%求核
core=[];

heweizhi=dis_total==1; %寻找只出现一次的二元对，以此来判断核
[row,col]=find(heweizhi,1,"first");

while ~isempty(row)
    local=DIS(row,col,:);
    selected=find(local==1);
    if size(selected,2)>1
        disp("程序存在错误")
    else
        core=[core,selected];
        heweizhi=min(double(heweizhi),~DIS(:,:,selected));
        [row,col]=find(heweizhi==1,1);
    end
end

%根据核中的覆盖对分辨关系进行删除
feature_remain=setdiff(1:n,core);

reduct=core;
DIS_reduct=max(DIS(:,:,core),[],3);

DIS_total_remain=DIS_total-DIS_reduct;

DIS_remain=DIS;

DIS_remain(:,:,reduct)=zeros(m,m,size(reduct,2));

while ~all(DIS_total_remain(:) == 0)
    DIS_num=zeros(n,1);
    for i=feature_remain
        DIS_remain(:,:,i)=min(DIS_remain(:,:,i),DIS_total_remain);
        DIS_num(i)=sum(DIS_remain(:,:,i),"all");
    end
    [~,selected]=max(DIS_num);
    reduct=[reduct,selected];
    DIS_total_remain=DIS_total_remain-DIS_remain(:,:,selected);
    feature_remain(feature_remain==selected)=[];

end
end
