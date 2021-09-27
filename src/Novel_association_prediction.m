function Novel_association_prediction
seed = 5;                                                                  
CV= 10;                                                                     
cross_validation(seed,CV);
end

function cross_validation(seed,CV)
load PREDICT;
beta = 0.4;
lamuda = 2^(-16);
interMat = predictAdMatdgc;
[dn,dr] = size(interMat);

positive = find(interMat == 1);
negative = find(interMat == 0);

rand('state',seed);
crossval_positive_idx = crossvalind('Kfold',interMat(positive),CV);
crossval_negative_idx = crossvalind('Kfold',interMat(negative),CV);
final = zeros(dn,dr);

for fold=1:CV
    
    fprintf('round %d validation\n',fold);
    trainMat = interMat;
    one_positive_idx = find(crossval_positive_idx == fold);
    one_negative_idx = find(crossval_negative_idx == fold);
    trainMat(positive(one_positive_idx)) = 0;
    
    K1 = [];
    K1(:,:,1)=predictSimMatdg;
    K1(:,:,2)=interaction_similarity(trainMat,'1' );
    K_COM1=SKF({K1(:,:,1),K1(:,:,2)},10,10,0.7);
    
    K2 = [];
    K2(:,:,1)=predictSimMatdcChemical;
    K2(:,:,2)=predictSimMatdcGo;
    K2(:,:,3)=interaction_similarity(trainMat,'2' );
    K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},15,10,0.2);
    
    F = LapRLS(K_COM1,K_COM2,trainMat,lamuda,beta);
    final(positive(one_positive_idx)) = F(positive(one_positive_idx));
    final(negative(one_negative_idx)) = F(negative(one_negative_idx));
    
end

[~,~,~,AUC] = perfcurve(interMat(:),final(:),1);
[~,~,~,AUPR] = perfcurve(interMat(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
fprintf('CV=10 AUC: %.3f  AUPR: %.3f\n', AUC, AUPR)
end


function [LapA] = LapRLS(W1,W2,inter3, lambda, beta)
[num_1,num_2] = size(inter3);

S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;
d_tmep_1=eye(num_1)/(D_1^(1/2));
L_D_11 = d_tmep_1*L_D_1*d_tmep_1;
A_1 = W1*pinv(W1 + lambda*L_D_11*W1)*inter3;

S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;
d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*inter3';

LapA=beta*A_1+(1-beta)*A_2';
end

function [W]=SKF(Wall,K,t,ALPHA)
%This program is recoded by reference follow: 
% Wang, B., Mezlini, A. M., Demir, F., Fiume, M., Tu, Z., Brudno, M., et al. (2014). 
% Similarity network fusion for aggregating data types on a genomic scale. Nature Methods, 11, 333¨C337
C = length(Wall);
[m,n]=size(Wall{1});

for i = 1 : C
    newW1{i} = Wall{i}./repmat(sum(Wall{i},1),n,1);
end
sumW1 = zeros(m,n);
for i = 1 : C
    sumW1= sumW1 + newW1{i};
end
for i = 1 : C
    Wall{i} = Wall{i}./repmat(sum(Wall{i},1),n,1);
end
for i = 1 : C
    newW{i} = FindDominateSet(Wall{i},round(K));
end
Wsum = zeros(m,n);
for i = 1 : C
    Wsum = Wsum + Wall{i};
end
for ITER=1:t
    for i = 1 : C           
         Wall{i}=ALPHA*newW{i}*(Wsum - Wall{i})*newW{i}'/(C-1) + (1-ALPHA)*(sumW1 - newW1{i})/(C-1);
    end   
    Wsum = zeros(m,n);
    for i = 1 : C
        Wsum = Wsum + Wall{i};
    end
end
W = Wsum/C;
w = neighborhood_Com(W,K);
W= W.*w;
end

function newW = FindDominateSet(W,K)
[m,n]=size(W);
[~,IW1] = sort(W,2,'descend');
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
clear IW2;
clear temp;
end

function similarities_N = neighborhood_Com(similar_m,kk)
similarities_N=zeros(size(similar_m));
mm = size(similar_m,1);
for ii=1:mm
    for jj=ii:mm
        iu = similar_m(ii,:);
        iu_list = sort(iu,'descend');
        iu_nearest_list_end = iu_list(kk);
        ju = similar_m(:,jj);
        ju_list = sort(ju,'descend');
        ju_nearest_list_end = ju_list(kk);
        if similar_m(ii,jj)>=iu_nearest_list_end && similar_m(ii,jj)>=ju_nearest_list_end
            similarities_N(ii,jj) = 1;
            similarities_N(jj,ii) = 1;
        elseif similar_m(ii,jj)<iu_nearest_list_end && similar_m(ii,jj)<ju_nearest_list_end
            similarities_N(ii,jj) = 0;
            similarities_N(jj,ii) = 0;
        else
            similarities_N(ii,jj) = 0.5;
            similarities_N(jj,ii) = 0.5;
        end
    end
end
end

function result = interaction_similarity(inter2,type)
[rna_num, pro_num] = size(inter2);
total = sum(sum(inter2));
if type == '1'
    result=zeros(rna_num,rna_num); 
    gama=rna_num/total;
    for i=1:rna_num
        for j=i+1:rna_num
            dis=pdist2(inter2(i,:),inter2(j,:),'squaredeuclidean');
            result(i,j)=exp(-gama*sum(dis));
            result(j,i)=exp(-gama*sum(dis));
        end
    end
else
    result = zeros(pro_num,pro_num); 
    gama = pro_num/total;
    for i=1:pro_num
        for j=i+1:pro_num
            dis=pdist2(inter2(:,i)',inter2(:,j)','squaredeuclidean');
            result(i,j)=exp(-gama*sum(dis));
            result(j,i)=exp(-gama*sum(dis));
        end
    end
end
end
