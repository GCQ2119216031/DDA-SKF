function Drug_single_similarity
seed = 8;                                                                  
CV = 10;                                                                
cross_validation(seed,CV);
end

function cross_validation(seed,CV)
load PREDICT;
lamuda=2^(-16);
interMat = predictAdMatdgc;
[dn,dr] = size(interMat);
positive = find(interMat == 1);
negative = find(interMat == 0);
rand('seed', seed);
crossval_positive_idx = crossvalind('Kfold',interMat(positive),CV);
crossval_negative_idx = crossvalind('Kfold',interMat(negative),CV);
final = zeros(dn,dr);

for cv = 1:CV
    fprintf('round %d validation\n',cv);
    train_interaction_matrix = interMat;
    
    one_positive_idx = find(crossval_positive_idx == cv);
    one_negative_idx = find(crossval_negative_idx == cv);
    train_interaction_matrix(positive(one_positive_idx)) = 0;
    
    K2 = [];
    K2(:,:,1) = predictSimMatdcChemical;  
    %K2(:,:,2) = predictSimMatdcGo;
    %K2(:,:,3) = interaction_similarity(train_interaction_matrix,'2' );
    %K_COM2 = SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},15,10,0.2);    
    
    score_matrix = LapRLS_dru(K2(:,:,1),train_interaction_matrix,lamuda);
    
    final(positive(one_positive_idx)) = score_matrix(positive(one_positive_idx));
    final(negative(one_negative_idx)) = score_matrix(negative(one_negative_idx));
    
end

[~,~,~,AUC] = perfcurve(interMat(:),final(:),1);
[~,~,~,AUPR] = perfcurve(interMat(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
fprintf('AUC: %.3f  AUPR: %.3f\n',AUC, AUPR);
end

function [LapA] = LapRLS_dru(W2,inter3, lambda)
[~,num_2] = size(inter3);

S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;
d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*inter3';

LapA=A_2';
end

function [W]=SKF(Wall,K,t,ALPHA)
%This program is recoded by reference follow: 
% Wang, B., Mezlini, A. M., Demir, F., Fiume, M., Tu, Z., Brudno, M., et al. (2014). 
% Similarity network fusion for aggregating data types on a genomic scale. Nature Methods, 11, 333–337
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
        
        Wall{i} = ALPHA*newW{i}*(Wsum - Wall{i})*newW{i}'/(C-1) + (1-ALPHA)*(sumW1 - newW1{i})/(C-1);
        
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
[YW,IW1] = sort(W,2,'descend');
clear YW;
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
		if similar_m(ii,jj)>=iu_nearest_list_end & similar_m(ii,jj)>=ju_nearest_list_end
			similarities_N(ii,jj) = 1;
			similarities_N(jj,ii) = 1;
		elseif similar_m(ii,jj)<iu_nearest_list_end & similar_m(ii,jj)<ju_nearest_list_end
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
[tar_num, dru_num] = size(inter2);
total = sum(sum(inter2));                     %inter2的2范式为什么是这个？应该这样：norm(inter2,2)
if type == '1'
    result=zeros(tar_num,tar_num); 
    gama=tar_num/total;
    for i=1:tar_num
        for j=i+1:tar_num
            dis=pdist2(inter2(i,:),inter2(j,:),'squaredeuclidean');
            result(i,j)=exp(-gama*sum(dis));
            result(j,i)=exp(-gama*sum(dis));
        end
    end
else
    result = zeros(dru_num,dru_num); 
    gama = dru_num/total;
    for i=1:dru_num
        for j=i+1:dru_num
            dis=pdist2(inter2(:,i)',inter2(:,j)','squaredeuclidean');
            result(i,j)=exp(-gama*sum(dis));
            result(j,i)=exp(-gama*sum(dis));
        end
    end
end
end

function [ result ] = normFun( M )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    num = size(M,1);
    nM = zeros(num,num);
    result = zeros(num,num);
    
    for i = 1:num
        nM(i,i) = sum(M(i,:));
    end

    for i = 1:num
        rsum = nM(i,i);
        for j = 1:num
            csum = nM(j,j);
            if((rsum==0)||(csum==0))
                result(i,j) = 0;
            else
                result(i,j) = M(i,j)/sqrt(rsum*csum);
            end
        end
    end
    
end

function d = kernel_Hamads(adjmat,dim)
%tju cs for bioinformatics 
% Calculates the Hamming distance kernel from a graph adjacency 
% matrix. If tha graph is unipartite, ka = kb.
%INPUT: 
% adjmat : binary adjacency matrix
% dim    : dimension (1 - rows, 2 - cols)
%OUTPUT:
% d : kernel matrix for adjmat over dimension 'dim'
    y = adjmat;
    % Graph based kernel
	if dim == 1
        ga = squareform(pdist(y,'hamming'));
    else
        ga = squareform(pdist(y','hamming'));
    end
    d=ones(size(ga))-ga;  
end

function result = target_sim_bla(algmt_matrix)
[tar_num,~] = size(algmt_matrix);
result=zeros(tar_num,tar_num);                    
for i=1:tar_num
    for j=1:tar_num
        if i == j
            result(i,j) = 1;
        else
            result(i,j) = algmt_matrix(i,j) / (sqrt(algmt_matrix(i,i)) * sqrt(algmt_matrix(j,j)));
        end
    end
end
end