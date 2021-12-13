# DDA-SKF
In this work, we used Similarity Kernel Fusion(SKF) and Laplacian Regularized Least Squares (LapRLS) algorithms to predict uncovered drug-disease associations. 

## Dependency

MATLAB R2019b

python 3.8

Predicting drug-disease associations need to be work on MATLAB platform. Verifying predictional associations need to configure a Python 3 environment to run 'Orphan_drug_prediction.py' and 'Verify_Drug_top_5.py'

## Dataset

Under the 'data' folder, 'PREDICT.mat' include 593 drugs and 313 diseases. Meanwhile, 'PREDICT.mat' provide two drug similarity kernels, one disease similarity kernel and the known association matrix for DDA-SKF. 'OrphanDrug.mat' provide drug similarity kernels, a association matrix of each orphan drug, drug names and disease OMIM IDs. 'LRSSL.mat' includes 763 drugs and 681 diseases. 'LRSSL.mat' provides two drug similarity kernels, one disease similarity kernel and the known association matrix for DDA-SKF.

## Main modules

### 1.Three methods to evaluate the model DDA-SKF

#### Novel association prediction

* 'PREDICT.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'PREDICT_divide.mat': Providing data-division for novel association prediction. A dataset was divided for five times with different random data partition schemes.
* 'Novel_association_prediction.m': Running this code to evaluate performance of our model on predicting potential associations. 

#### Novel drug prediction

* 'PREDICT.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'PREDICT_divide.mat': Providing data-division for novel drug prediction. A dataset was divided for five times with different random data partition schemes.
* 'Novel_drug_prediction.m': Running this code to evaluate performance of our model for predicting new drugs indications. 

#### Fivefold Cross-Validation (5-CV)

* 'PREDICT.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'LRSSL.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'PREDICT_divide.mat': Providing data-division for five cross-validation. A dataset was divided for five times with different random data partition schemes.
* 'LRSSL_divide.mat': Providing data-division for five cross-validation. A dataset was divided for five times with different random data partition schemes.
* 'Five_CV_result.m': Running this code to estimate the prediction performance of DDA-SKF model.

### 2.Comparing single similarity kernel and SKF

* 'PREDICT.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'Single_divide.mat ': Providing data-division for comparing single similarity kernel and SKF. A dataset was divided for five times with different random data partition schemes.
* 'Drug_single_similarity.m': Running this code to calculate evaluation metrics with every single drug similarity kernel or SKF.
* 'Disease_single_similarity.m': Running this code to calculate evaluation metrics with every single disease similarity kernel or SKF.
* 'Gragh_Drg_singleSim.mat': Providing data for drawing receiver operating characteristics curve and precision-recall curve in drug subspace.
* 'Gragh_Dis_singleSim.mat': Providing data for drawing receiver operating characteristics curve and precision-recall curve in disease subspace.
* 'Curve_Drg_SingleSim.m': Running this code to draw receiver operating characteristics curve and precision-recall curve with every single drug similarity kernel or SKF.
* 'Curve_Dis_SingleSim.m': Running this code to draw receiver operating characteristics curve and precision-recall curve with every single disease similarity kernel or SKF.

### 3.Predicting top-5 diseases for every drug

* 'PREDICT.mat': Providing drug similarity kernels, a disease similarity kernel and an association matrix. Besides, it provides drug names and disease OMIM IDs for outputting predictional associations. We have obtained the prediction matrix (named 'Score_matrix').
* 'PredictScore_matrix.m': Running this code to obtain scores for all candidate drug-disease associations. You can obtain the prediction matrix by running this code.
* 'drug_top5_resluts.m': Running this code to calculate top-5 indications for each drug.

### 4.Predicting top-5 diseases for every orphan drug

* 'OrphanDrug.mat': Providing drug similarity kernels and a association matrix of each orphan drug (for example, named 'predictAdmatOrphanDru187'). Besides, it provides drug names and disease OMIM IDs for outputting predictional associations. For each orphan drug, we have obtained their prediction matrix (for example, named 'Celecoxib_score_matrix').
* Each orphan drug association matrix had been removed all of their known associations before making predictions.
* 'PredictOrphandrug.m': Running this code to obtain prediction matrix of an orphan drug.
* 'VerifyOrphanResults.m': Running this code to calculate top-5 indications for an orphan drug.

### 5.Predicting top-5 drugs for several complex diseases

* 'PredictSocre_matrix_for_dis.m': Running this code to obtain prediction scores for all candidate drug-disease associations.
* 'Typical_dis_top5_resluts.m': Running this code to obtain top-5 candidate drugs for an complex disease.

### 6.Prediction results (This part of code is written in Python 3.8 )

#### Verifying prediction results for each drug

* 'Associations.csv': It is a set of the relationship between drugs and diseases from CTD.
* 'drug593-dis5.csv': It based on running the code 'drug_top5_resluts.m' on MATLAB platfrom.
* 'Verify_Drug_top_5.py': Running this code to verify the prediction relationship. The results will be saved to 'Drug_top_5_results.txt'.

#### Verifying prediction results for orphan drugs

* 'Associations.csv': It is a set of the relationship between drugs and diseases from CTD.
* 'col187DB00482-dis5.csv','col228DB00563-dis5.csv' and 'col428DB00997-dis5.csv': these tables based on running the code 'VerifyOrphanResults.m' on MATLAB platfrom.
* 'col187DB00482-dis5.csv': It is the prediction results of orphan drug Celecoxib.
* 'col228DB00563-dis5.csv': It is the prediction results of orphan drug Methotrexate.
* 'col428DB00997-dis5.csv': It is the prediction results of orphan drug Doxorubicin.
* 'Orphan_drug_prediction.py': Running this code to verify the prediction relationship. The results will be saved to 'Orphan_drug_results.txt'.

