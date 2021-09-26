# DDA-SKF
Created by ChuQiao-Gao
=======

Environment: Python3

1. The relationship between each drug and their top-5 indications was verified in CTD.

(1)'关联汇总.csv' is a set of the relationship between drugs and diseases from CTD.

(2)'drug593-dis5.csv' was based on the result of 'DrugTop5_Prediction' folder.

(3)Running 'Verify_Drug_top_5.py' to verify the prediction relationship. The results will be saved to 'Drug_top_5_results.txt'.

2. The relationship between each orphan drug and their top-5 indications was verified in CTD.

(1)'关联汇总.csv' is a set of the relationship between drugs and diseases from CTD.

(2)'col187DB00482-dis5.csv','col228DB00563-dis5.csv' and 'col428DB00997-dis5.csv'was based on the result of 'Prediction_OrphanDrug' folder.

(3)'col187DB00482-dis5.csv' is the prediction results of orphan drug Celecoxib.

(4)'col228DB00563-dis5.csv' is the prediction results of orphan drug Methotrexate.

(5)'col428DB00997-dis5.csv' is the prediction results of orphan drug Doxorubicin.

(6)Running 'Orphan_drug_prediction.py' to verify the prediction relationship. The results will be saved to 'Orphan_drug_results.txt'.
# DDA-SKF
