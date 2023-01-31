import scanpy as sc
import pandas as pd
import numpy as np


scdatapath='./inputSC/'


name=scdatapath+'/scdata2_most_recent/cluster_poss4_25.csv'
df=pd.read_csv(name)
cluster=df.to_numpy()
#input['sc_cluster']=data


celltypefname=scdatapath+'/scdata2_most_recent/nameOfCT_poss4_25.dat'
df=pd.read_csv(celltypefname,sep='\t',header=None)
cluname=df.to_numpy()
print(cluname)

for i in range(len(cluname)):
    #if cluname[i,1]=='KCs':
    #if cluname[i,1]=='Portain Vein EC':
    #if cluname[i,1]=='Stellate cells':
    #if cluname[i,1]=='Hep2_MP':
    if cluname[i,1]=='Central Vein EC':
        kcid=cluname[i,0]
        index=np.where(cluster[:,1]==kcid)
        kc_cells=cluster[index[0],0]

print('KC',len(kc_cells))

#full_ad_sc=sc.read_h5ad(scdatapath+'sc_liver_data.h5ad')
full_ad_sc=sc.read_h5ad(scdatapath+'HVG_counts_10000.h5ad')


genename=full_ad_sc.var_names.to_numpy()
cellname=full_ad_sc.obs_names.to_numpy()

print(len(cellname),len(genename))

d={}
for i in range(len(cellname)):
    d[cellname[i]]=i

index=[]
for i in range(len(kc_cells)):
    index.append(d[kc_cells[i]])

adata=full_ad_sc[index,:]
print(adata)

m=adata.X.toarray()
print(m.shape)

mu=np.mean(m,axis=0)
std=np.std(m,axis=0)

print(len(mu),len(std))

muindex=np.argsort(-mu)
stdindex=np.argsort(-std)

fw=open('cvec_mu_std.dat','w')
for i in range(len(genename)):
    fw.write(genename[muindex[i]]+':'+str(mu[muindex[i]])+'\t'+genename[stdindex[i]]+':'+str(std[stdindex[i]])+'\n')
