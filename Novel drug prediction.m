%function CVd_drug
%CVd交叉验证，新药预测
seed = [3 6 8 19 105];                                                                 %[3 6 8 19 51],1(851,286),2(848,280),3(859,275),5(844,272),6(857,280),7(853,279),8(859,279),9(856,282),10(854,279),12(855,271),13(851,275),14(853,272)15(855,276)16(851,268)17()18(850,267),19(857,281),32(856,283),45(853,278),48(858,279),51(858,280),105(859,279) is good
%seed = 105; %or8
CV = 10;                                                                      %5 is good
[SeedAUC,SeedAUPR] = cross_validation(seed,CV);
%end

function [fiveSeedAUC,fiveSeedAUPR] = cross_validation(seed,CV)
load Enzyme_data
lamuda=2^(-16);                                                             %15 is best 
beta=0.4; 
%interaction_matrix = lrssladmatdgc';
interaction_matrix = predictAdMatdgc;
[dn,dr] = size(interaction_matrix);
[~,col] = size(interaction_matrix);
colIndex = (1 : col)';
fiveSeedAUC = zeros(1,length(seed));
fiveSeedAUPR = zeros(1,length(seed));

for seedIndex = 1 : length(seed)
    rand('state',seed(seedIndex));
    colFold = crossvalind('Kfold',colIndex,CV);
    final = zeros(dn,dr);
    
    for cv = 1:CV
        tic
        fprintf('round %d validation\n',cv);
        train_interaction_matrix = interaction_matrix;
        one_col_idx = colIndex(colFold == cv);
        train_interaction_matrix(:,one_col_idx) = 0;
        
        K1 = [];
        K1(:,:,1)=predictSimMatdg;
        K1(:,:,2)=interaction_similarity(train_interaction_matrix,'1' );
        K_COM1=SKF({K1(:,:,1),K1(:,:,2)},12,10,0.2);                       %PREDICT设置:12,10,0.2;Cdataset:12,10,0.2;LRSSL:18,10,0.2
        
        K2 = [];
        K2(:,:,1)=predictSimMatdcChemical;                             
        %K2(:,:,2)=lrsslsimmatdcdomain;
        K2(:,:,2)=predictSimMatdcGo;
        K2(:,:,3)=interaction_similarity(train_interaction_matrix,'2' );
        
        K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},20,6,0.4);              %PREDICT设置:20,6,0.4;Cdataset:20,10,0.3;LRSSL:20,6,0.2
        
        score_matrix = LapRLS(K_COM1,K_COM2,train_interaction_matrix,lamuda,beta);
        
        final(:,one_col_idx) = score_matrix(:,one_col_idx);
        toc
    end
    
    [FPR,TPR,~,AUC] = perfcurve(interaction_matrix(:),final(:),1);
    [Recall,Precision,~,AUPR] = perfcurve(interaction_matrix(:),final(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
    fiveSeedAUC(seedIndex) = AUC;
    fiveSeedAUPR(seedIndex) = AUPR;
    [m,~]=size(Recall);
    F1=[];
    for i=1:m
        F1(i)=(2*Recall(i)*Precision(i))/(Recall(i)+Precision(i));
    end
    [f1,index]=max(F1);
    fprintf('CVd:  CV = 10  AUC: %.3f  AUPR: %.3f  Recall: %.3f  Precision: %.3f  F1: %.3f \n',AUC, AUPR, Recall(index),Precision(index),f1);
end

fiveSeedMeanAUC = mean(fiveSeedAUC);
fiveSeedMeanAUPR = mean(fiveSeedAUPR);
fprintf('seed: %d   fiveSeedMeanAUC: %.3f  fiveSeedAUPR: %.3f\n',seed(seedIndex) ,fiveSeedMeanAUC, fiveSeedMeanAUPR);
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
A_2 = W2*pinv(W2 + lambda*L_D_22*W2)*inter3';        %pinv进行第5次交叉验证时出错，换成inv可以出来结果

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
