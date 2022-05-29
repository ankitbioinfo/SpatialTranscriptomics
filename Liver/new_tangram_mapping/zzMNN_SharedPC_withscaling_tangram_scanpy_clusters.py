import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import os
#import pickle
import collections

from types import SimpleNamespace
from scipy.spatial import cKDTree
#from SCTransform import SCTransform

import warnings
import time
import seaborn as snn

warnings.filterwarnings('ignore')
#export PYTHONWARNINGS='ignore:Multiprocessing-backed parallel loops:UserWarning'
os.environ["PYTHONWARNINGS"] = "ignore::UserWarning"


def create_directory(outputFolder):
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)


def find_index(sp_genename,sc_genename):
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
    return index_sp,index_sc


def find_mutual_nn(data1, data2, sp_barcode,sc_barcode, k1, k2):
    print('I am here finding mutual nearest neighbor')
    #print(data1.shape,data2.shape)

    n_jobs=-1
    k_index_1 = cKDTree(data1).query(x=data2, k=k1, n_jobs=n_jobs)[1]
    k_index_2 = cKDTree(data2).query(x=data1, k=k2, n_jobs=n_jobs)[1]
    mutual_1 = []
    mutual_2 = []
    for index_2 in range(data2.shape[0]):
        for index_1 in k_index_1[index_2]:
            if index_2 in k_index_2[index_1]:
                mutual_1.append(index_1)
                mutual_2.append(index_2)

    a1=np.array(mutual_1)
    a2=np.array(mutual_2)
    a1=sp_barcode[a1]
    a2=sc_barcode[a2]
    a1=np.reshape(a1,(1,len(a1)))
    a2=np.reshape(a2,(1,len(a2)))

    b=np.concatenate((a1,a2)).T

    #print(a1.shape,a2.shape,b.shape)

    return b




def sct_return_sc_sp_in_shared_common_PC_space(ad_sp1,ad_sc1,no_of_pc,method):
    sct_ad_sc=ad_sc1.copy()
    sct_ad_sp=ad_sp1.copy()

    sc.pp.pca(sct_ad_sc,n_comps=no_of_pc)
    sc_com_pc=sct_ad_sc.varm['PCs']

    tp_sc=str(type(sct_ad_sc.X))
    tp_sp=str(type(sct_ad_sp.X))
    if tp_sc=="<class 'scipy.sparse.csr.csr_matrix'>":
        msc=sct_ad_sc.X.toarray()
    else:
        msc=sct_ad_sc.X

    if tp_sp=="<class 'scipy.sparse.csr.csr_matrix'>":
        msp=sct_ad_sp.X.toarray()
    else:
        msp=sct_ad_sp.X

    transfer_sp_com = np.matmul(msp, sc_com_pc)
    transfer_sc_com = np.matmul(msc, sc_com_pc)

    for i in range(transfer_sp_com.shape[1]):
        mu1=np.mean(transfer_sp_com[:,i])
        svd1=np.std(transfer_sp_com[:,i])
        transfer_sp_com[:,i]= (transfer_sp_com[:,i]-mu1)/svd1

        mu2=np.mean(transfer_sc_com[:,i])
        svd2=np.std(transfer_sc_com[:,i])
        transfer_sc_com[:,i]= (transfer_sc_com[:,i]-mu2)/svd2
        #print(i,mu1,mu2,svd1,svd2)


    sc_barcode=sct_ad_sc.obs_names.to_numpy()
    sp_barcode=sct_ad_sp.obs_names.to_numpy()

    #print('sc',transfer_sc_com.shape,sc_cellname.shape)
    #print('sp',transfer_sp_com.shape,sp_cellname.shape)

    return transfer_sp_com, transfer_sc_com, sp_barcode,sc_barcode





def find_annotation_index(annot_cellname,sct_cellname):
    d={}
    for i in range(len(annot_cellname)):
        d[annot_cellname[i]]=i

    index=[]
    for i in range(len(sct_cellname)):
        index.append(d[sct_cellname[i]])

    return index




def find_commnon_MNN(input):
    print('I am in MNN')
    df=pd.read_csv(input.fname,header=None)
    #data contains 2 column files of sct_pairing_shared_common_gene_PC.csv
    # first column is MNN pairs of spatial and
    # second column is MNN pairs of single cell
    data=df.to_numpy()

    mnn_singlecell_matchpair_barcode_id=np.unique(data[:,1])
    mnn_spatial_matchpair_barcode_id=np.unique(data[:,0])

    # find the annotated indexes
    index_annot_sc=find_annotation_index(input.annotation_singlecell_barcode_id,mnn_singlecell_matchpair_barcode_id)
    index_annot_sp=find_annotation_index(input.annotation_spatial_barcode_id,mnn_spatial_matchpair_barcode_id )



    #There are many indexes for spatial and single cell data
    # 1) MNN single cell                    data[:,1]                                       90,876
    # 2) MNN unique                          mnn_singlecell_matchpair_id                    10,089
    # 3) SC transform cell id                input.sct_singlecell_barcode_id                18,754
    # 4) original matrix cell id             input.annotation_singlecell_barcode_id         185,894
    # 5) original cell type name            input.annotation_singlecell_celltypename        185,894
    # 6) MNN unique id in sct               mnn_singlecell_matchpair_barcode_id             10,089
    # 7) common index between 6 and 4       index_mnn_sc,index_annot_sc

    # 1) MNN spatial                        data[:,0]                                       90,876
    # 2) MNN unique                         mnn_spatial_matchpair_id                        8,932
    # 3) SC transform cell id               input.sct_spatial_barcode_id                    86,880
    # 4) original matrix cell id            input.annotation_spatial_barcode_id             395,215
    # 5) original cell type name            input.annotation_spatial_celltypename           395,215
    # 55) original spatial cluster id       input.annotation_spatial_cluster_id             395,215
    # 6) MNN unique id in sct               mnn_spatial_matchpair_barcode_id                8,932
    # 7) common index between 6 and 4       index_mnn_sp,index_annot_sp

    d_single={}
    for i in range(len(input.annotation_singlecell_cluster_id)):
        d_single[input.annotation_singlecell_barcode_id[i]]=input.annotation_singlecell_cluster_id[i]

    d_spatial={}
    for i in range(len(input.annotation_spatial_cluster_id)):
        d_spatial[input.annotation_spatial_barcode_id[i]]=input.annotation_spatial_cluster_id[i]

    d_single_cluster={}
    for i in range(len(input.lsc[0])):
        singlecell_unique_clusterid=input.lsc[1][i]
        d_single_cluster[singlecell_unique_clusterid]=i

    d_spatial_cluster={}
    for i in range(len(input.lsp[0])):
        spatialcell_unique_clusterid=input.lsp[1][i]
        d_spatial_cluster[spatialcell_unique_clusterid]=i



    mat2=np.zeros((len(input.lsc[0]),len(input.lsp[0])),dtype=float)
    mat1=np.zeros((len(input.lsc[0]),2),dtype=int)
    mat3=np.zeros((1,len(input.lsp[0])),dtype=float)
    total=np.zeros((1,len(input.lsp[0])),dtype=float)


    #count how many anchor points matches to each spatial clusters
    unique_spatial_barcode_in_MNN=np.unique(data[:,0])
    for i in range(len(unique_spatial_barcode_in_MNN)):
        spatialcell_cluid=d_spatial[unique_spatial_barcode_in_MNN[i]]
        col=d_spatial_cluster[spatialcell_cluid]
        #print(i,spatialcell_cluid)
        mat3[0,col]+=1

    for i in range(len(input.annotation_spatial_cluster_id)):
        spatialcell_cluid=d_spatial[input.annotation_spatial_barcode_id[i]]
        col=d_spatial_cluster[spatialcell_cluid]
        total[0,col]+=1

    anchorFreq=mat3/total

    #print(mat3)
    print(total,np.sum(total))
    #print(anchorFreq)


    for i in range(len(data)):
            spatialcell_cluid=d_spatial[data[i,0]]
            singlecell_cluid=d_single[data[i,1]]
            #print(i,spatialcell_cluid,singlecell_cluid)
            col=d_spatial_cluster[spatialcell_cluid]
            row=d_single_cluster[singlecell_cluid]
            mat2[row,col]+=1

    #col normalization
    for i in range(len(mat2[0])):
        mat2[:,i]=mat2[:,i]/np.sum(mat2[:,i])

    mat2=np.vstack((anchorFreq,mat2))

    colname=['total # of sc', 'total # of sp']
    cname1=['anchorFreq']+input.lsc[0]
    cname2=input.lsp[0]
    #print(cname2)
    #fig,ax=plt.subplots(1,2,figsize=(20,10))
    fig=plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    ax0=plt.subplot(gs[0])
    ax1=plt.subplot(gs[1])
    fig.suptitle('MNN K = ' + str(input.KNN),fontsize=12)
    #snn.heatmap(ax=ax0,data=mat1,annot=True, fmt='d',xticklabels=colname, annot_kws={"size": 5},yticklabels=cname1)

    snn.heatmap(ax=ax1,data=mat2,annot=True, fmt='0.2f',xticklabels=cname2, annot_kws={"size": 5},yticklabels=cname1)
    plt.xlabel('Spatial cell clusters')
    plt.ylabel('Single cell clusters')

    #plt.title('R = '+str(radius)+', C='+str(lambda_c))
    #g.xaxis.set_ticks_position("top")
    plt.tight_layout()
    fig.savefig(input.savepath+'MNN_K_'+str(input.KNN)+'.png',dpi=300)







def main():

    outputFolder='input_vizgen_liver/MNN_shared50PC/'
    datapath='./input_vizgen_liver/alldata/'
    create_directory(outputFolder)

    #method='gauss'
    method='umap'

    neigh=50
    no_of_pc=50
    maindir='./'

    #load spatial cells clustering and cell type name
    spatialclusterFilename='./input_vizgen_liver/data/leiden_output.dat'
    #spatialclusterFilename='./input_vizgen_liver/data/sct_leiden_cluster_res100.dat'
    df=pd.read_csv(spatialclusterFilename)
    spatial_cluster_data=df.to_numpy()
    annotation_spatial_barcode_id= spatial_cluster_data[:,0]
    annotation_spatial_cluster_id= spatial_cluster_data[:,1]

    celltypefname='./input_vizgen_liver/data/NameOfCT.dat'
    #celltypefname='./input_vizgen_liver/data/NameOfCT_sct_res100.dat'

    df=pd.read_csv(celltypefname,sep='\t',header=None)
    data=df.to_numpy()
    spatialcell_unique_clustername=data[:,1]
    spatialcell_unique_clusterid=data[:,0]

    d={}
    for i in range(len(data)):
        d[data[i,0]]=data[i,1]

    annotation_spatial_celltypename=[]
    for i in range(len(annotation_spatial_cluster_id)):
        annotation_spatial_celltypename.append(d[annotation_spatial_cluster_id[i]])

    annotation_spatial_celltypename=np.array(annotation_spatial_celltypename)

    # spatial load sct clustering
    #sct_leiden=sc.read_h5ad(datapath+'scTransform_spatial_common.h5ad')


    #load single cell annotation
    #name=datapath+'annot_mouseStStAll.csv'
    #name=datapath+'annot_mouseStStAll_correct.csv'
    name=datapath+'annot_mouseStStAll_correct_monocytes.csv'

    df=pd.read_csv(name)
    data=df.to_numpy()
    annotation_singlecell_barcode_id=data[:,5]
    annotation_singlecell_celltypename=data[:,3]

    singlecell_unique_clustername=sorted(list(np.unique(annotation_singlecell_celltypename)))
    #print(unique_name)
    d={}
    singlecell_unique_clusterid=[]
    for i in range(len(singlecell_unique_clustername)):
        d[singlecell_unique_clustername[i]]=i
        singlecell_unique_clusterid.append(i)
    annotation_singlecell_cluster_id=[]
    for i in range(len(annotation_singlecell_celltypename)):
        annotation_singlecell_cluster_id.append(d[annotation_singlecell_celltypename[i]])
    annotation_singlecell_cluster_id=np.array(annotation_singlecell_cluster_id)


    #fmnn=outputFolder+"sct_pairing_gene_expression_space.csv"
    fmnn=outputFolder+"sct_pairing_shared_common_gene_PC_"+str(neigh)+".csv"
    #fmnn=outputFolder+"sct_pairing_shared_common_gene_PC.csv"


    flag=1
    if os.path.isfile(fmnn):
        filesize = os.path.getsize(fmnn)
        if filesize>0:
            flag=0

    if flag==1:
        '''
        ad_sc_ori=sc.read_h5ad(datapath+'sc_liver_data.h5ad')
        ad_sc_ori.var_names_make_unique()

        ad_sp_ori=sc.read_h5ad(datapath+'spatial_quadrant.h5ad')
        sc.pp.filter_cells(ad_sp_ori, min_counts=2)
        ad_sp_ori.var_names_make_unique()

        sp_genename=ad_sp_ori.var_names.to_numpy()
        sc_genename=ad_sc_ori.var_names.to_numpy()
        index_sp,index_sc=find_index(sp_genename,sc_genename)

        print('before sct normalization of single cell',len(index_sp),len(index_sc))
        ad_sc=ad_sc_ori[:,index_sc].copy()
        ad_sp=ad_sp_ori[:,index_sp].copy()
        '''

        fname=outputFolder+'HVG_basic_transform_spatial_liver_data.h5ad'
        if os.path.isfile(fname):
            filesize = os.path.getsize(fname)
            if filesize>0:
                sct_ad_sp=sc.read_h5ad(fname)
        else:
            sct_ad_sp = SCTransform(ad_sp,min_cells=5,gmean_eps=1,n_genes=500,n_cells=None, #use all cells
                        bin_size=500,bw_adjust=3,inplace=False)
            sct_ad_sp.write_h5ad(fname)

        fname=outputFolder+'HVG_basic_transform_singlecell_liver_data.h5ad'
        if os.path.isfile(fname):
            filesize = os.path.getsize(fname)
            if filesize>0:
                sct_ad_sc=sc.read_h5ad(fname)
        else:
            sct_ad_sc = SCTransform(ad_sc,min_cells=5,gmean_eps=1,n_genes=500,n_cells=None, #use all cells
                        bin_size=500,bw_adjust=3,inplace=False)
            sct_ad_sc.write_h5ad(fname)



        sp_genename=sct_ad_sp.var_names.to_numpy()
        sc_genename=sct_ad_sc.var_names.to_numpy()
        index_sp,index_sc=find_index(sp_genename,sc_genename)
        print('after sct normalization of single cell',len(index_sp),len(index_sc))

        ad_sp_ori=sct_ad_sp[:,index_sp].copy()
        ad_sc_ori=sct_ad_sc[:,index_sc].copy()

        ad_sc_ori.write_h5ad(outputFolder+'final_sct_sc.h5ad')
        ad_sp_ori.write_h5ad(outputFolder+'final_sct_sp.h5ad')

        print('\n\n 3 sct Shared Common Gene PC space')
        #ad_sp_ori=ad_sp_ori[0:100,:]
        #ad_sc_ori=ad_sc_ori[0:100,:]
        #print(ad_sc_ori)
        #print(ad_sp_ori)
        input_sp,input_sc,sp_barcode,sc_barcode=sct_return_sc_sp_in_shared_common_PC_space(ad_sp_ori,ad_sc_ori,no_of_pc,method)
        print('sp',input_sp.shape,'\nsc',input_sc.shape)
        corrected = find_mutual_nn(input_sp,input_sc,sp_barcode,sc_barcode, k1= neigh,k2= neigh)
        pd.DataFrame(corrected).to_csv(fmnn,index=False,header=None)



    if flag==0:
        ad_sc_ori=sc.read_h5ad(outputFolder+'final_sct_sc.h5ad')
        ad_sp_ori=sc.read_h5ad(outputFolder+'final_sct_sp.h5ad')
        singlecell_sct_barcode_id=ad_sc_ori.obs_names.to_numpy()
        spatialcell_sct_barcode_id=ad_sp_ori.obs_names.to_numpy()
        input={}
        input['fname']=fmnn
        input['annotation_singlecell_barcode_id']=annotation_singlecell_barcode_id
        input['annotation_singlecell_celltypename']=annotation_singlecell_celltypename
        input['annotation_singlecell_cluster_id']=annotation_singlecell_cluster_id
        input['lsc']=[singlecell_unique_clustername,singlecell_unique_clusterid]
        input['sct_singlecell_barcode_id']=singlecell_sct_barcode_id
        input['sct_spatial_barcode_id']=spatialcell_sct_barcode_id
        input['annotation_spatial_barcode_id']=annotation_spatial_barcode_id
        input['annotation_spatial_celltypename']=annotation_spatial_celltypename
        input['annotation_spatial_cluster_id']=annotation_spatial_cluster_id
        input['lsp']=[spatialcell_unique_clustername, spatialcell_unique_clusterid]
        input['savepath']=outputFolder
        input['KNN']=neigh


        inputt=SimpleNamespace(**input)
        find_commnon_MNN(inputt)




    return 0





main()
