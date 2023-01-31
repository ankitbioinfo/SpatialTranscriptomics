import Spatial_Annotations as spant
import pandas as pd
from types import SimpleNamespace
import scanpy as sc
import numpy as np


scdatapath='./inputSC/'
spdatapath='./inputSP/'

#outputFolder=spdatapath+'/find_MNN/'
outputFolder='./find_MNN_poss4_res0.5_dis15_/'

#dirname='./linux_results/'
#dirname='./macbasic_results/'
dirname='./macsct_results/'

neigh=50
no_of_pc=50
clusterResolution=[0.5,0.75,1,1.25,1.5]
clusterResolution=[0.25,1]
clusterResolution=[0.5]




input={}
input['neigh']=neigh
input['no_of_pc']=no_of_pc
input['outputFolder']=outputFolder

input['scfname']=scdatapath+'scTransform_singleCell.h5ad'
input['spfname']=spdatapath+'scTransform_spatial.h5ad'

'''
fname=spdatapath+'scTransform_spatial.h5ad'
adata=sc.read_h5ad(fname)
#print(adata)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
#sc.tl.umap(adata)
'''
resolutionWiseClustering=[]
resolutionWiseClusterName=[]



for fi in range(len(clusterResolution)):
    #sc.tl.leiden(adata,resolution=clusterResolution[fi])
    #n=np.unique(adata.obs.leiden)
    value=str(int(100*clusterResolution[fi]))
    celltypefname=dirname+'nameOfCTRes'+value+'.dat'
    spatialclusterFilename=dirname+'clusterRes'+value+'.csv'

    '''
    adata.obs.leiden.to_csv(spatialclusterFilename,header=True)
    fw=open(celltypefname,'w')
    for i in range(len(n)):
        fw.write(n[i]+'\t'+'c'+n[i]+'\n')
    fw.close()
    '''


    df=pd.read_csv(celltypefname,sep='\t',header=None)
    data=df.to_numpy()

    print('number of CT', data.shape)


    df=pd.read_csv(spatialclusterFilename)
    spatial_cluster_data=df.to_numpy()
    resolutionWiseClustering.append(spatial_cluster_data)
    resolutionWiseClusterName.append(data)






#adata.uns['leiden']
#adata.obs.leiden.to_csv
#load spatial cell cluster file and cell type annotation name incase if you already have it
#spatialclusterFilename=spdatapath+'authors/louvain_cluster.dat'


#load single cell cluster file and cell type annotation name

#name=scdatapath+'cluster_SI_coarse.csv'
#name=scdatapath+'cluster_SI_cellstate.csv'
name=scdatapath+'/scdata2_most_recent/cluster_poss4_25.csv'
df=pd.read_csv(name)
data=df.to_numpy()
input['sc_cluster']=data

#celltypefname=scdatapath+'nameOfCT_SI_cellstate.dat'
celltypefname=scdatapath+'/scdata2_most_recent/nameOfCT_poss4_25.dat'
df=pd.read_csv(celltypefname,sep='\t',header=None)
data=df.to_numpy()
#print(data[i,0],data[i,1],type(data[i,0]),type(data[i,1]))
input['sc_ct_name']=data
input['MNN_across_spatial_clusters_dispersion_cutoff']=0.15
input['minkowski_order']=2
input['number_of_iteration_in_degree_based_annotations']=3

#minkowski_order=[1,1]



if True:
    print('\n\n\n')
    input['sp_cluster']=resolutionWiseClustering[fi]
    input['sp_ct_name']=resolutionWiseClusterName[fi]
    input['resolutionWiseClustering']=[resolutionWiseClustering[fi]]

    #input['spatial_annotation_output_fname']=spdatapath+'/deg_annotation_spatial_'+str(neigh)
    input['spatial_annotation_output_fname']='deg_annotation_spatial_'+str(neigh)+'_'+str(int(clusterResolution[fi]*100))

    inputdata=SimpleNamespace(**input)
    spant.built_MNN_matrix_between_sc_and_sp(inputdata)


'''
for fi in range(len(clusterResolution)):
    print('\n\n\n',clusterResolution[fi],'\n')
    input['sp_cluster']=resolutionWiseClustering[fi]
    input['sp_ct_name']=resolutionWiseClusterName[fi]
    input['resolutionWiseClustering']=[resolutionWiseClustering[fi]]
    #input['minkowski_order']=minkowski_order[fi]

    #input['spatial_annotation_output_fname']=spdatapath+'/deg_annotation_spatial_'+str(neigh)
    input['spatial_annotation_output_fname']='deg_annotation_spatial_'+str(neigh)+'_'+str(int(clusterResolution[fi]*100))

    inputdata=SimpleNamespace(**input)
    spant.built_MNN_matrix_between_sc_and_sp(inputdata)
'''
