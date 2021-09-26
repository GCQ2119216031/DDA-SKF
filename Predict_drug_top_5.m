function Predict_drug_top_5
load PREDICT;
inter = predictAdMatdgc;                                                   %去除已知关联后，对预测结果中的未知关联进行排序
interIndex = find(inter);
Score_Mat = Score_matrix;
Score_Mat(interIndex) = 0;
Mat = Score_Mat';
[row,col] = size(Mat);
[~,index] = sort(Mat,2,'descend');
top = 5;                                                                   %每个药物的前5(可改20)个预测结果
predictInter = index(:,1:top);

for drg_i = 1 : 593
    for dis_i = 1 : top
        fprintf('%s\t%s\n',drugname(drg_i),disname(predictInter(drg_i,dis_i)));
    end
end
end