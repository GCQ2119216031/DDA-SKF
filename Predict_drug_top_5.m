load PREDICT;
inter = predictAdMatdgc;                                                   %去除已知关联后，对预测结果中的未知关联进行排序
interIndex = find(inter);
Score_Mat = Score_matrix;
Score_Mat(interIndex) = 0;
Mat = Score_Mat';
[row,col] = size(Mat);
[~,index] = sort(Mat,2,'descend');                                         %Index表示每行中的值从大到小排序的列号
top = 5;                                                                   %每个药物的前5(可改20)个预测结果
new_row = row*top;
predictInter = index(:,1:top); 
results = zeros(new_row,2);
results = string(results);
i = 1;
for drg_j = 1:row
    for a = 0:top-1
        results(i+a,1) = drugname(drg_j);
    end
    i = i + a + 1;
    drugname(drg_j)
end
for drg_i = 1 : row
    for dis_i = 1 : top
        fprintf('%s\t%s\n',drugname(drg_i),disname(predictInter(drg_i,dis_i)));
        results(((drg_i-1)*5 + dis_i), 2 ) = disname(predictInter(drg_i,dis_i));
    end
end

