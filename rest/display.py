import pandas as pd
import numpy as np


datapath='./caseA_res_0.5/find_MNN_poss4_res0.5_dis20_/'
datapath='./caseB_res_0.6/find_MNN_poss4_res0.6_dis5_/'
datapath='./caseBasics/find_MNN_poss4_res0.8_dis10_macsct/'
datapath='./find_MNN_poss4_res0.5_dis15_/'


df=pd.read_csv(datapath+'1_deg_annotation_spatial_50_50_ct_name.dat',header=None)
data1=df.to_numpy()

df=pd.read_csv(datapath+'2_deg_annotation_spatial_50_50_ct_name.dat',header=None)
data2=df.to_numpy()

df=pd.read_csv(datapath+'3_deg_annotation_spatial_50_50_ct_name.dat',header=None)
data3=df.to_numpy()


a=np.hstack((data1,data2,data3))

for i in range(len(a)):
    l=a[i]
    print(l[0],l[1],l[2],l[5],l[8])
