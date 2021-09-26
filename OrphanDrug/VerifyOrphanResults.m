load OrphanDrug;
Mat = Doxorubicin_score_matrix';
drug_i = 427;                                                              %Doxorubicin对应索引为427   
[row,col] = size(Mat);
[~,index] = sort(Mat,2,'descend');                                         
top = 5;
predictInter = index(:,1:top);                                             
results = zeros(top,2);
results = string(results);
results(:,1) = drugname(drug_i);

for dis_i = 1 : top
    fprintf('%s\t%s\n',drugname(drug_i),disname(predictInter(drug_i,dis_i)));
    results(dis_i, 2 ) = disname(predictInter(drug_i,dis_i));
end

