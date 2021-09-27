load PREDICT;
inter = predictAdMatdgc;
interIndex = find(inter);
Score_Mat = Score_matrix;
Score_Mat(interIndex) = 0;
Mat = Score_Mat';
[row,col] = size(Mat);
[~,index] = sort(Mat,2,'descend');
top = 5;
predictInter = index(:,1:top); 
results = zeros(500,2);
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
        results(((drg_i-1)*5 + dis_i), 2 ) = disname(predictInter(drg_i,dis_i));
    end
end

