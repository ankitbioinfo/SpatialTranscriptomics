
import pandas as pd
import numpy as np


df=pd.read_excel('BoneMarrowCellTypeAnnotation.xlsx',engine='openpyxl',sheet_name=0,header=1)
exp=df.to_numpy()

gene={}
pvalue={}
for i in range(len(exp)):
    CT=exp[i,1]
    gene[CT]=[]
    pvalue[CT]=[]


for i in range(len(exp)):
    CT=exp[i,1]
    gene[CT].append(exp[i,2])
    pvalue[CT].append(exp[i,3])



for key in gene:
    x=np.array(pvalue[key])
    y=np.array(gene[key])
    index=np.argsort(x)
    ind=index[0:25]
    print(key,y[ind],x[ind])
