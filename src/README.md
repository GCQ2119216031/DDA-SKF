# DDA-SKF
In this work, we used Similarity Kernel Fusion(SKF) and Laplacian regularized least squares (LapRLS) algorithms to predict uncovered drug-disease associations. 

## Dependency

MATLAB R2019b

Predicting drug-disease associations need to be work on MATLAB platform. 

## Dataset

Under the 'data' folder, 'PREDICT.mat' include 593 drug names and 313 disease OMIM IDs. Meanwhile, 'PREDICT.mat' provide two drug similarity kernels, one disease similarity kernel and the known association matrix for DDA-SKF. ‘OrphanDrug.mat’ provide drug similarity kernels, each orphan drug association matrix, drug names and disease OMIM IDs.

## Main modules

### 1.Two methods to evaluate the model DDA-SKF

#### Novel association prediction

* ‘PREDICT.mat’: providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'Novel_association_prediction.m': Running this code to evaluate performance of our model on predicting potential associations. 

#### Novel drug prediction

* ‘PREDICT.mat’: providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* 'Novel_drug_prediction.m': Running this code to evaluate performance of our model for predicting new drugs indications. 

### 2.Comparing single similarity kernel and SKF

* ‘PREDICT.mat’: providing drug similarity kernels, a disease similarity kernel and an association matrix for this part.
* ‘Drug_single_similarity.m’: Running this code to calculate the values of AUROC and AUPR with every single drug similarity kernel or SKF.
* ‘Disease_single_similarity.m’: Running this code to calculate the values of AUROC and AUPR with every single disease similarity kernel or SKF.

### 3.Predicting top-5 diseases for every drug

* ‘PREDICT.mat’: providing drug similarity kernels, a disease similarity kernel and an association matrix. Besides, it provides drug names and disease OMIM IDs for outputting predictional associations. We have obtained the prediction matrix (named ‘Score_matrix’).
* ‘PredictScore_matrix.m’: Running this code to obtain scores for all candidate drug-disease associations. You can obtain the prediction matrix by running this code.
* ‘drug_top5_resluts.m’: Running this code to calculate top-5 indications for each drug.

### 4.Predicting top-5 diseases for every orphan drug

* ‘OrphanDrug.mat’: providing drug similarity kernels and each orphan drug association matrix(for example, named ‘predictAdmatOrphanDru187’). Besides, it provides drug names and disease OMIM IDs for outputting predictional associations. For each orphan drug, we have obtained their prediction matrix (for example, named ‘Celecoxib_score_matrix’).
* Each orphan drug association matrix had been removed all of their known associations before making predictions.
* ‘PredictOrphandrug.m’: Running this code to obtain prediction matrix of an orphan drug.
* ‘VerifyOrphanResults.m’: Running this code to calculate top-5 indications for an orphan drug.
