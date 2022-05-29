
import pandas as pd
import numpy as np


spatialclusterFilename='./input_vizgen_liver/data/leiden_output.dat'
df=pd.read_csv(spatialclusterFilename)
spatial_cluster_data=df.to_numpy()
annotation_spatial_barcode_id= spatial_cluster_data[:,0]
annotation_spatial_cluster_id= spatial_cluster_data[:,1]

a=np.unique(annotation_spatial_cluster_id)

print(len(a))
cname={}
cname[10]='Cholangiocytes'
cname[22]='Bcel'

cname[3]='LSEC'
cname[12]='LSEC'
cname[21]='LSEC'


cname[31]='CentVeinEC'
cname[15]='PortVeinEC'


cname[18]='Fibro'
cname[19]='Fibro'
cname[20]='Fibro'

cname[8]='Stellate'

cname[0]='Hep'
cname[1]='Hep'
cname[2]='Hep'


cname[5]='cDC2'
cname[23]='cDC1'

cname[14]='NKT'

cname[9]='Neutro'
cname[35]='Neutro'

cname[11]='Macro'
cname[16]='Macro'
cname[24]='Macro'
cname[30]='Macro'


cname[13]='KC'
cname[7]='KC'
cname[4]='KC'

cname[6]='Mono'
cname[28]='Mono'

cname[17]='pDC'
cname[33]='HSPC'
cname[32]='MigcDC'

cname[25]='Unknown'
cname[26]='Unknown'
cname[27]='Unknown'
cname[29]='Unknown'
cname[34]='Unknown'

print(len(cname.keys()))


t=[]
for key in cname:
    if cname[key] not in t:
        t.append(cname[key])

ctname=sorted(t)
print(len(t),len(cname))

d={}
f=open('NameOfCT.dat','w')
for i in range(len(ctname)):
    f.write(str(i)+'\t'+ctname[i]+'\n')
    d[ctname[i]]=i



for i in range(len(annotation_spatial_cluster_id)):
    myname=cname[annotation_spatial_cluster_id[i]]
    spatial_cluster_data[i,1]=d[myname]

#print(annotation_spatial_cluster_id.shape,annotation_spatial_barcode_id.shape)
#data=np.concatenate((annotation_spatial_barcode_id,annotation_spatial_cluster_id),axis=1)
df=pd.DataFrame(data=spatial_cluster_data, index=None)
df.to_csv('sct_leiden_MNN.dat',index=False,header=['barcode','leiden'])
