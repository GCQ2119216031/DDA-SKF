import pandas as pd

data = pd.read_csv("col187DB00482-dis5.csv")

import numpy as np
drug = data[["drugName","disOMIM"]]
drug = np.array(drug)

ASSO = pd.read_csv("Associations.csv")

col_12 = ASSO[["OMIM","DrugName"]]
inter = np.array(col_12)

f = open('Orphan_drug_results.txt','w')
count = 0
for j ,y in enumerate(drug):
    for i ,x in enumerate(inter):
        if x[1] == y[0] and y[1] == x[0]:
            print(x[0])
            count = count + 1
            orphan_drug = str(x[1])
            disease = str(x[0])
            f.write(orphan_drug + '\t' + disease + '\n')
print(count)
f.close()