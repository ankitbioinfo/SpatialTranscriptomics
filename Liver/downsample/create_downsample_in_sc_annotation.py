



import scanpy as sc
import pandas as pd
import numpy as np

new_adata=sc.read_h5ad("sc_liver_data_downsample.h5ad")

cellname=new_adata.obs_names.to_numpy()
genename=new_adata.var_names.to_numpy()

print('gene',len(genename))

mat=new_adata.X.toarray()
umi=np.sum(mat,axis=1)
print('umi',len(umi))
umi=np.reshape(umi,(len(umi),1))


annot=pd.read_csv('./../annot_mouseStStAll.csv')
data=annot.to_numpy()
all_barcode_cell_name=data[:,5]
print(all_barcode_cell_name[0:5])

d={}
for i in range(len(data)):
    bid=data[i,5]
    d[bid]=i



index=[]
for i in range(len(cellname)):
    name=cellname[i]
    index.append(d[name])
    #print(key,dn[key],len(d[key]),count)

header=list(annot.columns)+['nUMI']
print(len(header),data[index].shape,umi.shape)

newannot=np.hstack((data[index],umi))
print(newannot.shape)

df=pd.DataFrame(data=newannot,columns=header)
df.to_csv('downsample_annotation.csv',index=False)
