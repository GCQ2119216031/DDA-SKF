function Novel_drug_prediction
% performing novel drug prediction
load PREDICT_divide;
CV = 10; 
Metrics = cross_validation(Novel_drug_Divide,CV);
results = Metrics_analysis(Metrics);
end

function metrics = cross_validation(Divide,CV)
load PREDICT
lamuda=2^(-16);
beta=0.4;
interaction_matrix = predictAdMatdgc;
[dn,dr] = size(interaction_matrix);
[~,col] = size(interaction_matrix);
colIndex = (1 : col)';

[num_div,~] = size(Divide);

metrics = zeros(num_div,8);

for divideIndex = 1 : num_div
    
    fprintf('No.%d division\n',divideIndex);
    colFold = Divide{divideIndex,1};
    final = zeros(dn,dr);
    
    for cv = 1:CV
        fprintf('round %d validation\n',cv);
        train_interaction_matrix = interaction_matrix;
        one_col_idx = colIndex(colFold == cv);
        train_interaction_matrix(:,one_col_idx) = 0;
        
        K1 = [];
        K1(:,:,1)=predictSimMatdg;
        K1(:,:,2)=interaction_similarity(train_interaction_matrix,'1' );
        K_COM1=SKF({K1(:,:,1),K1(:,:,2)},12,10,0.2);
        
        K2 = [];
        K2(:,:,1)=predictSimMatdcChemical;
        K2(:,:,2)=predictSimMatdcGo;
        K2(:,:,3)=interaction_similarity(train_interaction_matrix,'2' );
        
        K_COM2=SKF({K2(:,:,1),K2(:,:,2),K2(:,:,3)},20,6,0.4);
        
        score_matrix = LapRLS(K_COM1,K_COM2,train_interaction_matrix,lamuda,beta);
        
        final(:,one_col_idx) = score_matrix(:,one_col_idx);
    end
    
    metrics(divideIndex,:) = evaluation_function(interaction_matrix(:) , final(:));
    
end

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
% Similarity network fusion for aggregating data types on a genomic scale. Nature Methods, 11, 333?C337
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

function metrics = evaluation_function(interLabel , scoreLabel)
[FPR,~,~,AUC] = perfcurve(interLabel,scoreLabel,1);
[Recall,Precision,~,AUPR] = perfcurve(interLabel,scoreLabel,1, 'xCrit', 'reca', 'yCrit', 'prec');

[f_m,~] = size(Recall);
F1 = [];
for f_i = 1:f_m
    F1(f_i) = (2*Recall(f_i)*Precision(f_i))/(Recall(f_i)+Precision(f_i));
end
[f1,index] = max(F1);

% Calculate TP TN FP FN ACC MCC
num_pos = sum(interLabel == 1);
num_neg = sum(interLabel == 0);
TP = Recall(index) * num_pos;
FN = num_pos - TP;
TN = num_neg + TP - (TP / Precision(index));
FP = num_neg - TN;
PP = TP + FP;
PN = FN + TN;
ACC = (TP + TN) / (num_pos + num_neg);
MCC = ((TP * TN) - (FP * FN)) / sqrt(PP * PN * num_pos * num_neg);
Specificity = 1-FPR(index);

metrics = [AUC, AUPR, Recall(index), Specificity, Precision(index), ACC, f1, MCC];

end

function results = Metrics_analysis(Metrics)
results = zeros(8,2);
for i = 1 : 8
    results(i,1) = mean(Metrics(:,i));
    results(i,2) = std(Metrics(:,i));
end

fprintf('AUROC_mean: %.6f  AUROC_std: %.6f\n', results(1,1), results(1,2) );
fprintf('AUPR_mean: %.6f  AUPR_std: %.6f\n', results(2,1), results(2,2) );

end
