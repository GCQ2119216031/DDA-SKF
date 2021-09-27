# DDA-SKF
This part of code could verify predictional associations. We need to configure a Python3 environment to run 'Orphan_drug_prediction.py' and 'Verify_Drug_top_5.py

Putting the data and code in the same folder and running the code.

Commands in command line mode:

    python Verify_Drug_top_5.py

and

    python Orphan_drug_prediction.py

## 1.Verifying prediction results for each drug.

* 'Associations.csv' is a set of the relationship between drugs and diseases from CTD.

* 'drug593-dis5.csv' was based on running the code 'drug_top5_resluts.m' on MATLAB platform.

* Running 'Verify_Drug_top_5.py' to verify the prediction relationship. The results will be saved to 'Drug_top_5_results.txt'.

## 2.Verifying prediction results for orphan drugs.

* 'Associations.csv' is a set of the relationship between drugs and diseases from CTD.

* 'col187DB00482-dis5.csv','col228DB00563-dis5.csv' and 'col428DB00997-dis5.csv' were based on running the code 'VerifyOrphanResults.m' on MATLAB platform.

* 'col187DB00482-dis5.csv' is the prediction results of orphan drug Celecoxib.

* 'col228DB00563-dis5.csv' is the prediction results of orphan drug Methotrexate.

* 'col428DB00997-dis5.csv' is the prediction results of orphan drug Doxorubicin.

* Running 'Orphan_drug_prediction.py' to verify the prediction relationship. The results will be saved to 'Orphan_drug_results.txt'.



