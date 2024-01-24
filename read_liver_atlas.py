import scanpy as sc

import pandas as pd
import numpy as np

df=pd.read_csv('annot_mouseStStAll.csv')
data=df.to_numpy()
print(df.columns)

cluster=df['annot'].to_numpy()
umap=df[['UMAP_1', 'UMAP_2']].to_numpy()
cellid=df['cell'].to_numpy()
#cluster=df.to_numpy


df=pd.read_csv('barcodes.tsv',sep='\t',header=None)
barcode=df.to_numpy()

df=pd.read_csv('features.tsv',sep='\t',header=None)
genes=df.to_numpy()

d={}
for i in range(len(barcode)):
    d[barcode[i,0]]=i

index=[]
for i in range(len(cellid)):
    index.append(d[cellid[i]])

print(cellid[0:5],len(index))
newbarcode=barcode[index,0]

print(np.array_equal(newbarcode,cellid))
print(newbarcode.shape,umap.shape,cellid.shape,barcode.shape,genes.shape)


adata=sc.read_mtx('matrix.mtx', dtype='float32')

temp=adata.transpose()
temp.var_names=genes[:,0]
temp.obs_names=barcode[:,0]

newad=temp[index,:].copy()
newad.obs['cluster']=cluster
newad.obsm['X_umap']=umap

newad.write_h5ad('savedata.h5ad')




