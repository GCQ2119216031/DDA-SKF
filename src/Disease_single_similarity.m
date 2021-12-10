function Disease_single_similarity
seed = 8;   
CV = 10;      
Dis_metrics = cross_validation(seed,CV);
end

function metrics = cross_validation(seed,CV)
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
    
    K1 = predictSimMatdg;                                                  %Disease_Semantic_Similarity
    %K2 = interaction_similarity(train_interaction_matrix,'1' );            %Disease_Association_Similarity
    
    %K_COM1 = SKF({K1,K2},10,10,0.7);                                       %Disease_SKF
    
    %First parameter of LapRLS_dru replace with K1, K2 and K_COM1
    score_matrix = LapRLS_dis(K1,train_interaction_matrix,lamuda);         %Change the first parameter
    
    final(positive(one_positive_idx)) = score_matrix(positive(one_positive_idx));
    final(negative(one_negative_idx)) = score_matrix(negative(one_negative_idx));
    
end

% Various evaluation metrics
metrics = evaluation_function(interMat(:) , final(:));

end

function [LapA] = LapRLS_dis(W1,inter3, lambda)
[num_1,~] = size(inter3);

S_1 = W1;
d_1 = sum(S_1);
D_1 = diag(d_1);
L_D_1 = D_1 - S_1;
d_tmep_1=eye(num_1)/(D_1^(1/2));
L_D_11 = d_tmep_1*L_D_1*d_tmep_1;
A_1 = W1*pinv(W1 + lambda*L_D_11*W1)*inter3;

LapA=A_1;
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
[~,IW1] = sort(W,2,'descend');
newW=zeros(m,n);
temp=repmat((1:n)',1,K);
I1=(IW1(:,1:K)-1)*m+temp;
newW(I1(:))=W(I1(:));
newW=newW./repmat(sum(newW,2),1,n);
clear IW1;
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
total = sum(sum(inter2));  
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

function metrics = evaluation_function(interLabel , scoreLabel)
[FPR,TPR,~,AUC] = perfcurve(interLabel,scoreLabel,1);
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

fprintf('AUC: %.6f\n', AUC);
fprintf('AUPR: %.6f\n', AUPR);
fprintf('Recall: %.6f \n', Recall(index));
fprintf('Specificity: %.6f \n', Specificity);
fprintf('Precision: %.6f \n', Precision(index));
fprintf('ACC: %.6f \n', ACC);
fprintf('F1: %.6f \n', f1);
fprintf('MCC: %.6f \n', MCC);

metrics = [AUC, AUPR, Recall(index), Specificity, Precision(index), ACC, f1, MCC];

end
