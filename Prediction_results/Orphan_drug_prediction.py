#匹配预测结果是否命中
import pandas as pd

#其他两个孤儿药预测请选择对应文件
data = pd.read_csv("col187DB00482-dis5.csv")

import numpy as np
drug = data[["drugName","disOMIM"]]
drug = np.array(drug)

ASSO = pd.read_csv("关联汇总.csv")

col_12 = ASSO[["OMIM","药物名称"]]  #获取两列，要用二维数据
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