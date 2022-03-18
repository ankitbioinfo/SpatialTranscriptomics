import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
#from mnnpy import mnn_correct
from mnn import *
#from utils import *

import warnings
import time

'''
d1=pd.read_csv('/Users/agrawal/Desktop/mat1.csv')
data1=d1.to_numpy()

d2=pd.read_csv('/Users/agrawal/Desktop/mat1.csv')
data2=d2.to_numpy()

print(data1.shape)
print(data2.shape)


corrected = mnn_correct(data1.T,data2.T, k= 20,var_adj=False,batch_categories = ["sp", "sc"])

print(corrected)
'''


datapath='./../alldata/'
ad_sc_ori=sc.read_h5ad(datapath+'sc_liver_data_downsample.h5ad')
ad_sp_ori=sc.read_h5ad(datapath+'spatial_quadrant.h5ad')
ad_sp_ori.var_names_make_unique()
ad_sc_ori.var_names_make_unique()

sc.pp.filter_cells(ad_sp_ori, min_counts=2)


ad_sp1=ad_sp_ori
ad_sc1=ad_sc_ori

sp_cellname=ad_sp1.obs_names.to_numpy()
sp_genename=ad_sp1.var_names.to_numpy()
sc_cellname=ad_sc1.obs_names.to_numpy()
sc_genename=ad_sc1.var_names.to_numpy()
#print(len(sp_cellname),len(sp_genename))
#print(len(sc_cellname),len(sc_genename))

#find common genes
index_sc=[]
index_sp=[]
d={}
for j in range(len(sc_genename)):
    name=sc_genename[j]
    d[name]=j

for i in range(len(sp_genename)):
    name=sp_genename[i]
    try:
        d[name]
        flag=1
    except KeyError:
        flag=0
    if flag==1:
        index_sc.append(d[name])
        index_sp.append(i)

#print(len(index_sc),len(index_sp))

ad_sc=ad_sc1[:,index_sc].copy()
ad_sp=ad_sp1[:,index_sp].copy()


#ad_sc=ad_sc[0:100,:].copy()
#ad_sp=ad_sp[0:100,:].copy()

print(ad_sp.X.shape,len(ad_sp.var_names),len(ad_sp.obs_names))
print(ad_sc.X.shape,len(ad_sc.var_names),len(ad_sc.obs_names))


m1=ad_sp.X.transpose()
df=pd.DataFrame(data=m1, index=ad_sp.var_names , columns=ad_sp.obs_names)
df.to_csv("for_R_spatial.csv")

m2=ad_sc.X.toarray()
m2=m2.transpose()
df2=pd.DataFrame(m2, index=ad_sc.var_names , columns=ad_sc.obs_names)
df2.to_csv("for_R_sc.csv")

print('start')
sp_gene=ad_sp.var_names.to_numpy()
sc_gene=ad_sc.var_names.to_numpy()
hvg=set(sp_gene).intersection(set(sc_gene))
hvgs=sorted(list(hvg))

corrected = mnn_correct(ad_sp,ad_sc,var_subset=hvgs, k= 20,var_adj=False,batch_categories = ["sp", "sc"])
corrected[0].write_h5ad("corrected_common.h5ad")

a=np.array(corrected[1])
a[1].to_csv("pairing_separate.csv",index=False)
