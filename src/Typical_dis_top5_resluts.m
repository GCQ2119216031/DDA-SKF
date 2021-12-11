% AD, PD in PREDICT dataset
load PREDICT;

% ALS in LRSSL dataset
%load LRSSL;

inter = predictAdMatdgc;
%inter = lrssladmatdgc';
interIndex = find(inter);
% known associations were removed, the uncovered associations will be ranked
Score_matrix(interIndex) = 0;
Mat = Score_matrix;
[row,col] = size(Mat);
% Sort by row, 2:row
[~,index] = sort(Mat,2,'descend');
predictInter = index(:,1:5);

% Predict drugs for AD: dis_i = 8;
% Predict drugs for PD: dis_i = 120;
% Predict drugs for ALS: dis_i = 31 and using Score_matrix in LRSSL dataset
dis_i = 8;

% PREDICT dataset: using disname and drugname
% LRSSL dataset: using lrsslDisName and lrsslDrgName
fprintf('Disease: %s\n',disname(dis_i));
fprintf('Predict_Drug: %s\n',drugname(predictInter(dis_i,:)));


