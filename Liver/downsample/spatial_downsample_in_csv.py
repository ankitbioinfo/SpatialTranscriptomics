


import scanpy as sc
import pandas as pd


new_adata=sc.read_h5ad("spatial_quadrant_87210.h5ad")

sc.pp.filter_cells(new_adata, min_counts=2)

cellname=new_adata.obs_names.to_numpy()
genename=new_adata.var_names.to_numpy()

print('gene',len(genename))


spatialdata=new_adata.copy()

sc.pp.normalize_total(spatialdata, inplace=True)
sc.pp.log1p(spatialdata)
sc.pp.highly_variable_genes(spatialdata, flavor="seurat", n_top_genes=2000)
#sc.pl.highly_variable_genes(full_ad_sc,show=False,save='.png')
#newsc2 = newsc[:, newsc.var.highly_variable]

new_adata=new_adata[:,spatialdata.var.highly_variable]
mat=new_adata.X.toarray()
print('2',len(spatialdata.var.highly_variable),mat.shape)


#df=pd.DataFrame(data=mat.transpose(), index=genename , columns=cellname)
#df.to_csv("spatial_liver_data_downsample_87210.csv")

df=pd.read_csv("tissue_positions_list_quadrant.csv",header=None)
data=df.to_numpy()
print(data.shape,len(cellname))
d={}
for i in range(len(data)):
    bid=data[i,0]
    d[bid]=i


index=[]
for i in range(len(cellname)):
    name=cellname[i]
    index.append(d[name])

newdata=data[index]

header=['barcodes','X','Y']
df=pd.DataFrame(data=newdata,columns=header)
df.to_csv('tissue_positions_list_87210.csv',index=False)

count=0
for i in range(len(cellname)):
    if cellname[i]==newdata[i,0]:
        count+=1
print(count)
