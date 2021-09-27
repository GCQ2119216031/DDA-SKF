import pandas as pd
import numpy as np

data = pd.read_csv("drug593-dis5.csv")

drug_top5 = data[["drugName","disOMIM"]]  
drug_top5 = np.array(drug_top5)

ASSO = pd.read_csv("Associations.csv")

col_12 = ASSO[["OMIM","DrugName"]]  
inter = np.array(col_12)

count = 0
f = open('Drug_top_5_results.txt','w')
for j ,y in enumerate(drug_top5):
    for i ,x in enumerate(inter):
        if x[1] == y[0] and y[1] == x[0]:
            count = count + 1
            drug = str(x[1])
            disease = str(x[0])
            f.write(drug + '\t' + disease + '\n')
print(count)    
f.close()
            