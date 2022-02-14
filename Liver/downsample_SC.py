import scanpy as sc
import pandas as pd
import random

random.seed(1001)

#dataFolder='countTable_mouseStSt/'

annot=pd.read_csv('annot_mouseStStAll.csv')
data=annot.to_numpy()
all_barcode_cell_name=data[:,5]
print(all_barcode_cell_name[0:5])

d={}
dn={}
dind={}
for i in range(len(data)):
    cluster=data[i,2]
    if cluster in d:
        d[cluster].append([data[i,0],data[i,1]])
        dind[cluster].append(i)
    else:
        d[cluster]=[[data[i,0],data[i,1]]]
        dn[cluster]=data[i,3]
        dind[cluster]=[i]


clusterid=sorted(list(d.keys()))

index=[]
for key in clusterid:
    cellid=dind[key]
    for i in range(len(cellid)):
        if random.random()<0.1:
            index.append(i)
    #print(key,dn[key],len(d[key]),count)


df=pd.DataFrame(data=data[index],columns=annot.columns)
df.to_csv('downsample_annotation.csv',index=False)

ad_sc=sc.read_h5ad('sc_liver_data.h5ad')
print(ad_sc)
cellname=ad_sc.obs_names.to_numpy()
genename=ad_sc.var_names.to_numpy()

print('dim',cellname.shape,genename.shape,data.shape)

count=0
for j in range(len(cellname)):
    if cellname[j]==data[j,5]:
        count=count+1

new_ad_sc=ad_sc[index,:]
print('c',count,new_ad_sc)
new_ad_sc.write_h5ad("sc_liver_data_downsample.h5ad")


#df=pd.DataFrame(data=adata.X.transpose(), index=adata.var_names , columns=adata.obs_names)
#df.to_csv("xxoutput_filter_"+filename+".csv")

#f=open('sc_NameOfCT.dat','r')
#for line in f:
#    l=line.split('\t')
#    print(l)
