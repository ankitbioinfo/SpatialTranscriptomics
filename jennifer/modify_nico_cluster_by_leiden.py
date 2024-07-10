import scanpy as sc
import numpy as np
import pandas as pd
adata = sc.read_h5ad('spatial_BM_with_nico_annotations.h5ad')

ct_leiden = adata.obs['leiden0.5']
ct_nico  = np.array(adata.obs['nico_ct'])


nctl = np.unique(ct_leiden)

'''
['0' '1' '10' '11' '12' '13' '14' '15' '16' '17' '18' '2' '3' '4' '5' '6'
 '7' '8' '9']
['B.cell' 'Endothelial' 'Erythrocyte' 'FCGR3B.neutrophil'
 'GNLY.high.T.cell' 'HSPC' 'LTF.high.neutrophil' 'MMP9.high.neutrophil'
 'Mast.cell' 'Megakaryotcyte.Epcam+' 'Mesenchymal2' 'Mesenchymal9'
 'Monocyte' 'Myeloid.prog1' 'Myeloid.prog14' 'Myeloid.prog20'
 'Myeloid.prog7' 'T.cell' 'cDC/Macrophage' 'infl.mesenchymal'
 'infl.monocyte' 'pDC']
 '''

 #Since 12 cluster id in leiden is adipocytes first find all the cells associated with that

index = np.where(ct_leiden=='12')
# these cell index put in nico as adipocytes

ct_nico[index[0]] = 'Adipocytes'
# save in anndata object

adata.obs['modify_nico_ct']= ct_nico
adata.write_h5ad('modify_spatial_BM_with_nico_annotations.h5ad')
# also save for niche prediction
nctnico = sorted(list(np.unique(ct_nico)))

data=[]
d={}
for i in range(len(nctnico)):
    data.append([i,nctnico[i]])
    d[nctnico[i]]=i

df= pd.DataFrame(data)
df.to_csv('4_nico_annotation_ct_name.csv',index=False)

data=[]
for i in range(len(ct_nico)):
    data.append(   d[ct_nico[i]]  )

df= pd.DataFrame(data,index=adata.obs_names)
df.to_csv('4_nico_annotation_cluster.csv')
