#匹配预测结果是否命中
import pandas as pd
import numpy as np

data = pd.read_csv("drug593-dis5.csv")

drug_top5 = data[["drugName","disOMIM"]]  
drug_top5 = np.array(drug_top5)

ASSO = pd.read_csv("关联汇总.csv")

col_12 = ASSO[["OMIM","药物名称"]]  
inter = np.array(col_12)

#for drug_top5 in data_1:
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
            