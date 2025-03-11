function [reduct_new] = Incremental_fs(shuju,beta,object_new,DIS,reduct,lower,neighborhood_system,Ker)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

[~,m_new]=size(object_new);
[~,m_old]=size(object_old);
m_total=m_old+m_new;
%更新决策类

shuju_total=shuju(([object_old,object_new]),:);
covering_total=shuju_total(:,1:end-1);
d_total=shuju_total(:,end);

for i=1:m_total
    d_class_total(i,:)=(d_total==d_total(i));
end
non_class_total=~d_class_total;

%更新邻域系统
neigh_sys_old=neighborhood_system;
neigh_sys_new=(covering_total(m_old+1:m_total,:)>=beta);
neigh_sys_total=[neigh_sys_old;neigh_sys_new];

%更新邻域
neigh_total=[];
for i=1:m_total
    neigh_total(i,:)=min(covering_total(:,neigh_sys_total(i,:)),[],2);%第i行是第i个元素的邻域
end
 neigh_total(neigh_total>=beta)=1;

 %更新相似关系
relation_total=1-neigh_total;

 %更新下近似
lower_total=[];
for i=1:m_total
    lower_total(i)=min(relation_total(i,non_class_total(i,:)));
end


%更新每个覆盖元相对分辨关系
DIS_new=zeros(m_total,m_total,size(DIS,3));
DIS_new(1:m_old,1:m_old,:)=DIS;

object_change=find(lower_total(1:m_old)<lower);

tic
DIS_add1=zeros(m_total,m_total,size(DIS,3));
%DIS_add1=DIS_new;
for i=object_change
    for j=find(neigh_sys_old(i,:))
        %DIS_add1(i,1:m_old,j)=((lower(i)>1-covering_total(1:m_old,j))&(1-covering_total(1:m_old,j)>=lower_total(i)))'&non_class_total(i,1:m_old);
        DIS_add1(i,1:m_old,j)=((lower(i)>1-covering_total(1:m_old,j))&(1-covering_total(1:m_old,j)>=lower_total(i)))'&non_class_total(i,1:m_old);
    end
end

DIS_add2=zeros(m_total,m_total,size(DIS,3));
for i=1:m_old
    for j=find(neigh_sys_old(i,:))
%         D=zeros(m_total,m_total,size(DIS,3));
%         D(i,m_old+1:m_total,j)=(1-covering_total(m_old+1:m_total,j)>=lower_total(i));
%         DIS_add2(i,m_old+1:m_total,j)=D(i,m_old+1:m_total,j)&non_class_total(i,m_old+1:m_total);
        DIS_add2(i,m_old+1:m_total,j)=(1-covering_total(m_old+1:m_total,j)>=lower_total(i))'&non_class_total(i,m_old+1:m_total);
    end
end


DIS_add3=zeros(m_total,m_total,size(DIS,3));
for i=m_old+1:m_total
    for j=find(neigh_sys_total(i,:))
%         D=zeros(m_total,m_total,size(DIS,3));
%         D(i,:,j)=(1-covering_total(:,j)>=lower_total(i));
%         DIS_add3(i,:,j)= D(i,:,j)&non_class_total(i,:);
        DIS_add3(i,:,j)= (1-covering_total(:,j)>=lower_total(i))'&non_class_total(i,:);
    end
end
DIS_add=DIS_add1|DIS_add2|DIS_add3;%DIS_add为每个覆盖元增加的相对分辨关系


DIS_new=DIS_new|DIS_add;%DIS_add为更新后的每个覆盖元总体的相对分辨关系

%判断总体相对分辨关系是否与非类向量相等
DIS_total=(sum(DIS_new,3)>=1);
if ~isequal(non_class_total,DIS_total)
    disp("程序存在错误")
end

%判断现有约简是否需要添加元素
DIS_red_old=max(DIS_new(:,:,reduct),[],3);

%根据现有约简中的覆盖对分辨关系进行添加
feature_remain=setdiff(1:n-1,reduct);
reduct_new=reduct;
DIS_total_remain=DIS_total-DIS_red_old;
DIS_remain=DIS_new;
DIS_remain(:,:,reduct)=zeros(m_total,m_total,size(reduct,2));

while ~all(DIS_total_remain(:) == 0)
    DIS_num=zeros(n,1);
    for i=feature_remain
        DIS_remain(:,:,i)=min(DIS_remain(:,:,i),DIS_total_remain);
        DIS_num(i)=sum(DIS_remain(:,:,i),"all");
    end
    [~,selected]=max(DIS_num);
    reduct_new=[reduct_new,selected];
    DIS_total_remain=DIS_total_remain-DIS_remain(:,:,selected);
    feature_remain(feature_remain==selected)=[];

end

%根据Ker对现有进行覆盖元删除
Ker_new=zeros(m_total,m_total,size(reduct_new,2));


for i=1:size(reduct,2)
    ker_del1=max(DIS_add1(1:m_old,1:m_old,setdiff(reduct,reduct(i))),[],3);
    ker_del2=max(DIS_new(1:m_old,1:m_old,setdiff(reduct_new,reduct)),[],3);
    ker_del=~(ker_del1|ker_del2);
    Ker_new(1:m_old,1:m_old,i)=Ker(:,:,i)&ker_del;
end

feature_del=[];
for i=1:size(reduct_new,2)
    Ker_new(1:m_old,m_old+1:m_total,i)=non_class_total(1:m_old,m_old+1:m_total)-max(DIS_new(1:m_old,m_old+1:m_total,setdiff(reduct_new,reduct_new(i))),[],3);
    Ker_new(m_old+1:m_total,:,i)=non_class_total(m_old+1:m_total,:)-max(DIS_new(m_old+1:m_total,:,setdiff(reduct_new,reduct_new(i))),[],3);
    if sum(Ker_new(:,:,i),"all")==0
        feature_del=[feature_del,reduct_new(i)];
    end
end

while ~isempty(feature_del)
    reduct_new=setdiff(reduct_new,feature_del(1));
    Ker_new(:,:,feature_del(1))=[];
    feature_del=[];
    for i=1:size(reduct_new,2)
        Ker_new(:,:,i)=non_class_total(:,:)-max(DIS_new(:,:,setdiff(reduct_new,reduct_new(i))),[],3);
        if sum(Ker_new(:,:,i),"all")==0
            feature_del=[feature_del,reduct_new(i)];
        end
    end
end
end