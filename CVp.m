function LPI_SKF_2
seed = 8;                                                                  %[8 7 6 5 4 ] 4(930,287)5(927,292),7(927,291)8(930,288)
CV= 10;                                                                     % 10-fold:8(935,250)7(932,243)6(933,247)5,2(934,246)4(934,245)3(932,247)42(934,247)1(934,242)0(935,242)
[AUC, AUPR] = cross_validation(seed,CV);
end

function [ALPHA_AUC , ALPHA_AUPR] = cross_validation(seed,CV)
%function cross_validation(seed,CV)
    load Enzyme_data.mat;
    beta=0.4;                                                              % 0.4 IS GOOD
    lamuda=2^(-16);
    
    
%     ALPHA_2 = 15;
     ALPHA_AUC = zeros(1,200);
     ALPHA_AUPR = zeros(1,200);
%     for ALPHA_1_num = 1: 11
%         ALPHA_1 = 6;
%     for ALPHA_num = 1 : 11
        
            positive = find(predictAdMatdgc == 1);
            negative = find(predictAdMatdgc == 0);
            
            rand('seed', seed);
            crossval_positive_idx = crossvalind('Kfold',predictAdMatdgc(positive),CV);
            
            testAUC = zeros(1,CV);
            testAUPR = zeros(1,CV);
            
            for fold=1:CV
                
                fprintf('round %d validation\n',fold);
                trainMat = predictAdMatdgc;
                final = zeros(313,593);
                F = zeros(313,593);
                one_positive_idx = find(crossval_positive_idx == fold);
                trainMat(positive(one_positive_idx)) = 0;
                
                K1 = [];
                K1(:,:,1)=predictSimMatdg;
                K1(:,:,2)=interaction_similarity(trainMat,'1' );
                K_COM1=SKF({K1(:,:,1),K1(:,:,2)},10,10,0.7);                 
                
                K2 = [];
                K2(:,:,1)=predictSimMatdcChemical;                             %simMatrix593, my simMatrix---RDKit
                %K2(:,:,3)=predictSimMatdcDomain;
                K2(:,:,2)=predictSimMatdcGo;
                K2(:,:,3)=interaction_similarity(trainMat,'2' );
                K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},15,10,0.2);   
                
                %normK1 = normFun(K_COM1);                                  %norm+BiRW一下确实管用,norm+LapRLS一下也管用（效果（926，297）不如BiRW（930，267））
                %normK2 = normFun(K_COM2);
                
                F = LapRLS(K_COM1,K_COM2,trainMat,lamuda,beta);
                %F = BiRW(normK1,normK2,trainMat,3,3,0.1);                  %ALPHA0.1,0.2最好，0.3-0.6次之,l,r相等的情况下，2以上就好(3:931,273)l,c相同效果比不同好
                final(positive(one_positive_idx)) = F(positive(one_positive_idx));
                final(negative) = F(negative);
                all = [(positive(one_positive_idx)); negative];
                
                [~,~,~,AUC] = perfcurve(predictAdMatdgc(all),final(all),1);
                %AUC
                [Recall,Precision,~,AUPR] = perfcurve(predictAdMatdgc(all),final(all),1, 'xCrit', 'reca', 'yCrit', 'prec');
                %AUPR
                fprintf('CV=10 beta=0.4 lamuda=2^(-16) K2_ALPHA = 0.9 AUC: %.3f  AUPR: %.3f \n', AUC, AUPR);
                
                [f_m,~] = size(Recall);
                F1 = [];
                for f_i = 1:f_m
                    F1(f_i) = (2*Recall(f_i)*Precision(f_i))/(Recall(f_i)+Precision(f_i));
                end
                [f1,index] = max(F1);
                fprintf('CV=10 AUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision: %.3f  F1: %.3f \n', AUC, AUPR, Recall(index),Precision(index),f1)
                fprintf('i=1,Precision: %.3f\n', Precision(2));
                
                testAUC(fold) = AUC;
                testAUPR(fold) = AUPR;
                
            end
            
            
            meanTestAUC = mean(testAUC);
            fprintf('seed: %d  meanTestAUC: %.3f\n',seed, meanTestAUC);
            
            meanTestAUPR = mean(testAUPR);
            fprintf('seed: %d  meanTestAUPR: %.3f\n',seed, meanTestAUPR);
            
        %end
        %ALPHA_AUC(ALPHA_num) = meanTestAUC;
        %ALPHA_AUPR(ALPHA_num) = meanTestAUPR;
        
        %ALPHA_1 = ALPHA_1 + 1;
    %end
        %ALPHA_2 = ALPHA_2 + 1;
%     end
end



function [LapA] = LapRLS(W1,W2,inter3, lambda, beta)
[num_1,num_2] = size(inter3);

S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;
d_tmep_1=eye(num_1)/(D_1^(1/2));                     %eye return  identity matrix,size is num_1
L_D_11 = d_tmep_1*L_D_1*d_tmep_1;
A_1 = W1*pinv(W1 + lambda*L_D_11*W1)*inter3;

S_2 = W2;
d_2 = sum(S_2);
D_2 = diag(d_2);
L_D_2 = D_2 - S_2;
d_tmep_2=eye(num_2)/(D_2^(1/2));
L_D_22 = d_tmep_2*L_D_2*d_tmep_2;
A_2 = W2*inv(W2 + lambda*L_D_22*W2)*inter3';        %pinv进行第5次交叉验证时出错，换成inv可以出来结果

LapA=beta*A_1+(1-beta)*A_2';
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
[YW,IW1] = sort(W,2,'descend');
clear YW;
newW=zeros(m,n);
temp=repmat((1:n)',1,K);             %包含列向量的水平堆栈
I1=(IW1(:,1:K)-1)*m+temp;            %前k个相似向量的索引值
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
total = sum(sum(inter2));                     %inter2的2范式为什么是这个？应该这样：norm(inter2,2)
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
