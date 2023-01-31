
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
#from scipy.spatial import Voronoi, ConvexHull,voronoi_plot_2d, Delaunay
from numpy.linalg import norm

from sklearn.datasets import make_classification
from sklearn.multioutput import MultiOutputRegressor
from sklearn.linear_model import LogisticRegression,LogisticRegressionCV, Ridge, RidgeCV, LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import RandomizedSearchCV,GridSearchCV,cross_val_predict, cross_val_score,RepeatedKFold,RepeatedStratifiedKFold,StratifiedShuffleSplit
#from sklearn.metrics import make_scorer,accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve
from sklearn.metrics import confusion_matrix,r2_score,mean_absolute_error,mean_squared_error,mean_squared_log_error,mean_absolute_percentage_error,median_absolute_error, max_error, explained_variance_score

from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import Pipeline
#from sklearn.metrics import precision_recall_fscore_support as score
#from imblearn.over_sampling import SMOTE, SMOTEN,ADASYN, KMeansSMOTE, SVMSMOTE
from sklearn.utils import class_weight
from sklearn.metrics import roc_curve, auc,consensus_score

#from sklearn.datasets import make_checkerboard
from sklearn.cluster import SpectralBiclustering
from sklearn.decomposition import NMF

#bicluster

from gseapy.plot import gseaplot, heatmap
import gseapy
from sklearn.decomposition import PCA as skPCA


#Metrics
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import log_loss
from sklearn.metrics import zero_one_loss
from sklearn.metrics import matthews_corrcoef
from scipy.stats import pearsonr

import pandas as pd
import numpy as np
import seaborn as snn
import os
import random
import warnings
import time
import scanpy as sc
import pickle
import xlsxwriter
from types import SimpleNamespace
import math
import scipy
from sklearn.utils.extmath import svd_flip
import statsmodels.api as sm
import shap




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

def read_spatial_data(clusterFilename,celltypeFilename):

    df=pd.read_csv(celltypeFilename,sep='\t',header=None)
    data=df.to_numpy()
    spatialcell_unique_clustername=data[:,1]
    spatialcell_unique_clusterid=data[:,0]
    CTname=spatialcell_unique_clustername
    CTid=spatialcell_unique_clusterid

    df=pd.read_csv(clusterFilename)
    louvainFull=df.to_numpy()

    #for i in range(len(CTname)):
    #    print(CTname[i],CTid[i])

    celltype={}
    cellsinCT={}
    index=[]
    for i in range(len(louvainFull)):
        #print(louvainFull[i],louvainFull[i][0])
        clu_id=louvainFull[i][1]
        cel_id=louvainFull[i][0]
        if clu_id in CTid:
            index.append(i)
            #celltype[cel_id]=clu_id
            if clu_id not in cellsinCT:
                cellsinCT[clu_id]=[cel_id]
            else:
                cellsinCT[clu_id].append(cel_id)

    louvain=louvainFull[index,:]
    annotation_spatial_barcode_id= louvain[:,0]
    annotation_spatial_cluster_id= louvain[:,1]



    #print(spatialcell_unique_clustername,spatialcell_unique_clusterid)
    d={}
    for i in range(len(spatialcell_unique_clustername)):
        d[spatialcell_unique_clusterid[i]]=spatialcell_unique_clustername[i]
    annotation_spatial_celltypename=[]
    for i in range(len(annotation_spatial_cluster_id)):
        annotation_spatial_celltypename.append(d[annotation_spatial_cluster_id[i]])
    annotation_spatial_celltypename=np.array(annotation_spatial_celltypename)

    return annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid


def find_correlation_bw_genes_and_PC_component_in_singlecell(KcomponentCluster,clusterExpression):
    #print('ankit',KcomponentCluster.shape,clusterExpression.shape)
    mat=np.zeros((clusterExpression.shape[1],KcomponentCluster.shape[1]),dtype=float)
    for i in range(clusterExpression.shape[1]):
        v1=clusterExpression[:,i]
        for j in range(KcomponentCluster.shape[1]):
            v2=KcomponentCluster[:,j]
            corr,_ = pearsonr(v1,v2)
            #corr=np.corrcoef(v1,v2)
            #print(corr,corr1)
            mat[i,j]=corr

    # mat shape is (# of genes x # of pc) it is a correlation between (PC and genes) of the single cell cluster
    # KcomponentCluster shape is (# of single cell in a single cell cluster x # of pc)
    # clusterExpression shape is (# of single cell in a single cell cluster x # of genes)
    #print(mat.shape)

    mat=np.nan_to_num(mat)

    return mat

'''
def mean_and_scaling_normalization(m):
    for i in range(m.shape[1]):
        mu1=np.nanmean(m[:,i])
        svd1=np.nanstd(m[:,i])
        #if (svd1==0)&(mu1==0):
        svd1=1
        #print(i,mu1,svd1)
        m[:,i]= (m[:,i]-mu1)/svd1
    return m
'''

def top_genes_in_correlation_list(genename,corr_NMFfactors_genes,n_top_words):
        top_genes_assoc_factors=[]
        for topic_idx, topic in enumerate(abs(corr_NMFfactors_genes.T)):
            top_features_ind = topic.argsort()[: -n_top_words - 1 : -1]
            #print(topic_idx,sorted(topic))
            #print(top_features_ind,topic[top_features_ind[0:5]])
            for i in top_features_ind:
                if i not in top_genes_assoc_factors:
                    top_genes_assoc_factors.append(i)


        #print(top_genes_assoc_factors)
        gname=genename[top_genes_assoc_factors]
        mat=corr_NMFfactors_genes[top_genes_assoc_factors,:]

        return gname,mat

'''
def NMF_GKL(W,H,V):

    for i in range(1):
        num=np.matmul(W.T,V)
        F=(V/np.matmul(W,H))
        print("GKL",F.shape,V.shape,W.shape,H.shape)
        #den=np.matmul(W.T,)
        H=H*(num/den)
        value=np.sqrt(np.sum((V-np.matmul(W,H))**2))

    return W,H
'''

def nndsvd_initialization(A,rank):
    u,s,v=np.linalg.svd(A,full_matrices=False)
    v=v.T
    w=np.zeros((A.shape[0],rank))
    h=np.zeros((rank,A.shape[1]))

    w[:,0]=np.sqrt(s[0])*np.abs(u[:,0])
    h[0,:]=np.sqrt(s[0])*np.abs(v[:,0].T)

    for i in range(1,rank):

        ui=u[:,i]
        vi=v[:,i]
        ui_pos=(ui>=0)*ui
        ui_neg=(ui<0)*-ui
        vi_pos=(vi>=0)*vi
        vi_neg=(vi<0)*-vi

        ui_pos_norm=np.linalg.norm(ui_pos,2)
        ui_neg_norm=np.linalg.norm(ui_neg,2)
        vi_pos_norm=np.linalg.norm(vi_pos,2)
        vi_neg_norm=np.linalg.norm(vi_neg,2)

        norm_pos=ui_pos_norm*vi_pos_norm
        norm_neg=ui_neg_norm*vi_neg_norm

        if norm_pos>=norm_neg:
            w[:,i]=np.sqrt(s[i]*norm_pos)/ui_pos_norm*ui_pos
            h[i,:]=np.sqrt(s[i]*norm_pos)/vi_pos_norm*vi_pos.T
        else:
            w[:,i]=np.sqrt(s[i]*norm_neg)/ui_neg_norm*ui_neg
            h[i,:]=np.sqrt(s[i]*norm_neg)/vi_neg_norm*vi_neg.T

    return w,h

def mu_method(A,k,max_iter,init_mode='random'):

    if init_mode == 'random':
        W ,H = random_initialization(A,k)
    elif init_mode == 'nndsvd':
        W ,H = nndsvd_initialization(A,k)
    norms = []
    e = 1.0e-10
    for n in range(max_iter):
        # Update H
        W_TA = W.T@A
        W_TWH = W.T@W@H+e
        for i in range(np.size(H, 0)):
            for j in range(np.size(H, 1)):
                H[i, j] = H[i, j] * W_TA[i, j] / W_TWH[i, j]
        # Update W
        AH_T = A@H.T
        WHH_T =  W@H@H.T+ e
        for i in range(np.size(W, 0)):
            for j in range(np.size(W, 1)):
                W[i, j] = W[i, j] * AH_T[i, j] / WHH_T[i, j]

        norm = np.linalg.norm(A - W@H, 'fro')
        norms.append(norm)
    return W ,H ,norms

def multiplicative_method(W,H,A,max_iter):
    norms = []
    e = 1.0e-10
    for n in range(max_iter):
        # Update H
        W_TA = W.T@A
        W_TWH = W.T@W@H+e
        for i in range(np.size(H, 0)):
            for j in range(np.size(H, 1)):
                H[i, j] = H[i, j] * W_TA[i, j] / W_TWH[i, j]
        # Update W
        #AH_T = A@H.T
        #WHH_T =  W@H@H.T+ e
        #for i in range(np.size(W, 0)):
        #    for j in range(np.size(W, 1)):
        #        W[i, j] = W[i, j] * AH_T[i, j] / WHH_T[i, j]

        norm = np.linalg.norm(A - W@H, 'fro')
        norms.append(norm)
    return W ,H ,norms

def find_PC_of_invidualCluster_in_SC(scbarcode,scadata,no_of_pc,spbarcode,spadata, sct_ad_sc_full,outputdir,celltype_name):
    print('\n\n\n')

    #full transcriptome single cell
    cellname=sct_ad_sc_full.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(scbarcode)):
        index.append(d[scbarcode[i]])
    full_genes_sc=sct_ad_sc_full[index,:].copy()

    #common gene single cell
    cellname=scadata.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(scbarcode)):
        index.append(d[scbarcode[i]])
    sct_ad_sc=scadata[index,:].copy()
    sc_cellname=sct_ad_sc.obs_names.to_numpy()

    #common gene spatial
    cellname=spadata.obs_names.to_numpy()
    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i
    index=[]
    for i in range(len(spbarcode)):
        index.append(d[spbarcode[i]])

    sct_ad_sp=spadata[index,:].copy()
    sp_cellname=sct_ad_sp.obs_names.to_numpy()


    tp_sc=str(type(full_genes_sc.X))
    if tp_sc=="<class 'scipy.sparse._csr.csr_matrix'>":
        CbyG=full_genes_sc.X.toarray()
    else:
        CbyG=full_genes_sc.X


    tp_sc=str(type(sct_ad_sc.X))
    if tp_sc=="<class 'scipy.sparse._csr.csr_matrix'>":
        msc=sct_ad_sc.X.toarray()
    else:
        msc=sct_ad_sc.X

    tp_sp=str(type(sct_ad_sp.X))
    if tp_sp=="<class 'scipy.sparse._csr.csr_matrix'>":
        msp=sct_ad_sp.X.toarray()
    else:
        msp=sct_ad_sp.X

    #replace nan to zero
    #msp=np.nan_to_num(msp)
    #msc=np.nan_to_num(msc)

    #print('norm1',np.sum(msc),np.sum(msp),np.sum(CbyG))
    #msc=msc/np.sum(msc)
    #msp=msp/np.sum(msp)
    #CbyG=CbyG/np.sum(CbyG)
    #print('norm2',np.sum(msc),np.sum(msp),np.sum(CbyG))
    genename_joint=sct_ad_sc.var_names.to_numpy()
    genename_full=full_genes_sc.var_names.to_numpy()
    genename_spatial=sct_ad_sp.var_names.to_numpy()





    #Fitting the NMF model (Frobenius norm) with tf-idf features,
    #model = NMF(n_components=no_of_pc, init = "nndsvda", random_state=1,beta_loss="frobenius",alpha_W=0.00005,alpha_H=0.00005,l1_ratio=1)
    #Fitting the NMF model (generalized Kullback-Leibler divergence) with tf-idf features
    #l1_ratio The regularization mixing parameter, with 0 <= l1_ratio <= 1.
    #For l1_ratio = 0 the penalty is an elementwise L2 penalty (aka Frobenius Norm).
    #For l1_ratio = 1 it is an elementwise L1 penalty. For 0 < l1_ratio < 1, the penalty is a combination of L1 and L2.
    #model = NMF(n_components=no_of_pc, init = "nndsvda", random_state=1,beta_loss="kullback-leibler",solver="mu",max_iter=1000,alpha_W=0.0,alpha_H=0.0,l1_ratio=0)

    #W = model.fit_transform(msc.T)
    #H = model.components_
    #scale_factor = np.max(H)
    #H = H/scale_factor
    #W = W*scale_factor

    #scale_factor = np.max(W)
    #H = H*scale_factor
    #W = W/scale_factor

    W ,H ,norms=mu_method(msc.T,no_of_pc,100,init_mode='nndsvd')

    value=np.sqrt(np.sum((msc.T-np.matmul(W,H))**2))
    #NMF_GKL(W,H,msc.T)



    print("W",len(np.sum(W,axis=0)), len(np.sum(W,axis=1)),     np.sum(W,axis=0),   np.sum(W,axis=1)[0:5]         )
    print("H",len(np.sum(H,axis=0)), len(np.sum(H,axis=1)),     np.sum(H,axis=1),   np.sum(H,axis=0)[0:5]         )
    print("Value",value)

    spW=W
    V=msp.T
    spH=np.ones((no_of_pc,msp.shape[0]),dtype=float)

    spW ,spH ,norms=multiplicative_method(spW,spH,V,100)
    '''
    for i in range(100):
        num=np.matmul(spW.T,V)
        den=np.matmul(spW.T,np.matmul(spW,spH))
        spH=spH*(num/den)
        spH=np.nan_to_num(spH)
        value=np.sqrt(np.sum((V-np.matmul(spW,spH))**2))
        #spH = spH * (t(W) %*% V)/(t(W) %*% W %*% H)
        #  sqrt(sum(( V - W%*%H )**2))
    #print(i,value)
    '''
    value=np.sqrt(np.sum((V-np.matmul(spW,spH))**2))
    print('nan',np.sum(np.isnan(spH)), np.sum(np.isnan(spW)), np.sum(np.isnan(msp)))
    print("SW",len(np.sum(spW,axis=0)), len(np.sum(spW,axis=1)),     np.sum(spW,axis=0),   np.sum(spW,axis=1)[0:5]         )
    print("SH",len(np.sum(spH,axis=0)), len(np.sum(spH,axis=1)),     np.sum(spH,axis=1),   np.sum(spH,axis=0)[0:5]         )
    print("value",value)

    corr_NMFfactors_genes=find_correlation_bw_genes_and_PC_component_in_singlecell(H.T,msc)
    gname,geneNMF=top_genes_in_correlation_list(genename_joint,corr_NMFfactors_genes,5)
    etcetcetc=corr_NMFfactors_genes

    fig=plt.figure(figsize=(12,5))
    xlabels=[]
    for i in range(no_of_pc):
        xlabels.append('NMF'+str(i+1))
    gs = fig.add_gridspec(ncols=3, nrows=1, wspace=0.5,width_ratios=[1, 1,2])
    ax0=fig.add_subplot(gs[0])
    ax1=fig.add_subplot(gs[1])
    ax2=fig.add_subplot(gs[2])
    b=snn.heatmap(geneNMF,yticklabels=gname,ax=ax0)#componentlabel,ax=ax
    b.set_xticklabels(xlabels,size = 8,rotation=90)

    corr_NMFfactors_genes=find_correlation_bw_genes_and_PC_component_in_singlecell(spH.T,msp)
    gname,geneNMF=top_genes_in_correlation_list(genename_spatial,corr_NMFfactors_genes,5)
    b=snn.heatmap(geneNMF,yticklabels=gname,ax=ax1)#componentlabel,ax=ax
    b.set_xticklabels(xlabels,size = 8,rotation=90)
    #b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)

    corr_NMFfactors_genes=find_correlation_bw_genes_and_PC_component_in_singlecell(H.T,CbyG)
    gname,geneNMF=top_genes_in_correlation_list(genename_full,corr_NMFfactors_genes,10)
    b=snn.heatmap(geneNMF,yticklabels=gname,ax=ax2)#componentlabel,ax=ax
    b.set_xticklabels(xlabels,size = 8,rotation=90)
    b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)

    '''
    localdir='./temp/'+celltype_name
    create_directory(localdir)
    #if celltype_name=='Portain_Vein_EC':
    if True:

        for i in range(len(genename_joint)):
            if genename_joint[i]=='Serping1':
                d1=msc[:,i]

        for i in range(len(genename_full)):
            if genename_full[i]=='Serping1':
                d2=CbyG[:,i]


        A=d1#.X.toarray()
        B=d2#.X.toarray()

        value=np.array_equal(A,B)

        print(A.shape,B.shape,value)
        print(np.mean(A),np.mean(B),np.std(A),np.std(B))



        f1=open(localdir+'etcetcetc.dat','w')
        for i in range(len(etcetcetc)):
            f1.write(genename_joint[i]+'\t'+str(etcetcetc[i])+'\n')
        f1.close()
        f1=open(localdir+'corr_NMFfactors_genes.dat','w')
        for i in range(len(corr_NMFfactors_genes)):
            f1.write(genename_full[i]+'\t'+str(corr_NMFfactors_genes[i])+'\n')
        f1.close()
    '''


    fig.savefig(outputdir+celltype_name+'.png',dpi=300)
    plt.close('all')


    #sc_com_pc=W
    #transfer_sp_com = np.matmul(msp, sc_com_pc)
    #transfer_sc_com = np.matmul(msc, sc_com_pc)
    transfer_sp_com=spH.T
    print('\n\n',celltype_name,msc.shape,msp.shape,"W",W.shape,"H",H.shape,'shape','tsp',transfer_sp_com.shape,"corr",corr_NMFfactors_genes.shape)

    sc_cluster_mean_exp=np.mean(CbyG,axis=0)
    M=corr_NMFfactors_genes
    transfer_sc_com=[]


    sc_barcode=sct_ad_sc.obs_names.to_numpy()
    sp_barcode=sct_ad_sp.obs_names.to_numpy()
    sc_genenames=full_genes_sc.var_names.to_numpy()



    return transfer_sp_com, transfer_sc_com, sp_barcode,sc_barcode, M,sc_genenames, sc_cluster_mean_exp



def makePCneighboorhoodFeatureMatrix(input):
    #print('xyz',input.annotation_spatial_barcode_id.shape)
    n=len(input.spatialcell_unique_clusterid)
    M=np.zeros((len(input.neighbors),n*input.no_of_pc),dtype=float)

    for j in range(len(input.neighbors)):
        #print(j,neighbors[j])
        CC_barcode_id=input.annotation_spatial_barcode_id[j]
        CC_cluster_id=input.annotation_spatial_cluster_id[j]
        PC_component_of_CC=input.pc_of_sp_clusterid[CC_barcode_id]
        PC_component_of_CC=PC_component_of_CC.reshape((1,input.no_of_pc))
        if j==0:
            target=PC_component_of_CC
        else:
            target=np.vstack((target,PC_component_of_CC))
        #print(target.shape)
        temp={}
        for k in range(len(input.neighbors[j])):
            id=input.neighbors[j][k]
            NC_barcode_id=input.annotation_spatial_barcode_id[id]
            NC_cluster_id=input.annotation_spatial_cluster_id[id]
            PC_component_of_NC=input.pc_of_sp_clusterid[NC_barcode_id]
            PC_component_of_NC=PC_component_of_NC.reshape((1,input.no_of_pc))
            if NC_cluster_id not in temp:
                temp[NC_cluster_id]=PC_component_of_NC
                #print(temp[NC_cluster_id].shape)
            else:
                temp[NC_cluster_id]=np.concatenate((temp[NC_cluster_id],PC_component_of_NC))
                #print(len(neighbors[j]),NC_cluster_id,temp[NC_cluster_id].shape)
            #print(PC_component_of_NC)

        for key in input.spatialcell_unique_clusterid:
            start_index=input.no_of_pc*key
            end_index=start_index+input.no_of_pc
            if key in temp:
                M[j,start_index:end_index]=np.mean(temp[key],axis=0)
        #print(j,temp)
            #d[NC_barcode_id]=1

    #cluster=input.annotation_spatial_cluster_id
    #cluster=cluster.reshape((len(cluster),1))
    #df=pd.DataFrame(np.hstack((cluster,M)))
    df=pd.DataFrame(np.hstack((target,M)))
    df.to_csv(input.outputname,index=True,header=None)





def compute_PC_space(input,sct_ad_sc_full):
    a=set(input.singlecell_unique_clustername)
    b=set(input.spatialcell_unique_clustername)
    common=a.intersection(b)


    print("\n\n Spatial and sc # of clusters",len(b),len(a))
    print('Common cell types between spatial and single cell data',len(common),common)
    print('\n Spatial cluster name not matching to single cell cluster name', b-common)
    print("Common cell types between single cell and spatial cells must be equivalent to spatial cell types otherwise method may go wrong ")
    print("If you have extra spatial cell types then you must remove it before running the spatial method itself")
    print("\n\n")

    flag=1
    if len(b-common)>0:
        flag=0

    if flag==1:
        n=len(input.spatialcell_unique_clustername)
        pc_of_sp_clusterid={}
        PCA_of_sc_cluster_accordingto_spatial_clusterid={}
        for i in range(n):
            clidsp=input.spatialcell_unique_clusterid[i]
            index=np.where(input.annotation_spatial_cluster_id==clidsp)
            spbarcode=input.annotation_spatial_barcode_id[index[0]]
            scbarcode=[]
            for j in range(len(input.singlecell_unique_clustername)):
                #print(input.singlecell_unique_clustername[j],input.spatialcell_unique_clustername[i])
                if input.singlecell_unique_clustername[j]==input.spatialcell_unique_clustername[i]:
                    clid=input.singlecell_unique_clusterid[j]
                    index=np.where(input.annotation_singlecell_cluster_id==clid)
                    scbarcode=input.annotation_singlecell_barcode_id[index[0]]
                    break

            pc_sp,pc_sc,sp_barcode,sc_barcode,sc_HVG_correlation,sc_genenames,sc_cluster_mean_exp=find_PC_of_invidualCluster_in_SC(scbarcode,input.ad_sc,input.no_of_pc,spbarcode,input.ad_sp, sct_ad_sc_full, input.nmf_output,input.spatialcell_unique_clustername[i])
            #for kkkk in range(len(sc_genenames)):
            #    if sc_genenames[kkkk]=='Mptx2':
            #        print('ankit',sc_cluster_mean_exp[kkkk])

            #check orthogonality
            #print("1 spatial orthogonality",sum(pc_sp[:,0]*pc_sp[:,1]) )
            #print("1 SC orthogonality",sum(pc_sc[:,0]*pc_sc[:,1]) )

            print(i,input.spatialcell_unique_clustername[i],pc_sp.shape,sc_HVG_correlation.shape,sc_genenames.shape)
            PCA_of_sc_cluster_accordingto_spatial_clusterid[clidsp]=[sc_HVG_correlation,pc_sp,sc_genenames,sc_cluster_mean_exp]
            for k in range(len(sp_barcode)):
                pc_of_sp_clusterid[sp_barcode[k]]=pc_sp[k]
            #pc_of_sp_clusterid[input.spatialcell_unique_clusterid[i]]=pc_sp

    return pc_of_sp_clusterid,PCA_of_sc_cluster_accordingto_spatial_clusterid



def model_linear_regression(input,savedir,logistic_predicted_interactions):
#def model_linear_regression(K_fold,n_repeats,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms,seed,n_jobs):
    shap_cluster_cutoff=input.shap_cluster_cutoff
    data1 = np.genfromtxt(open(input.outputname, "rb"), delimiter=',', skip_header=0)
    #ind=~np.isnan(data1).any(axis=1)
    #data=data1[ind,:]
    data=np.nan_to_num(data1)


    #print(data1.shape,data.shape)
    featureVector=range(input.no_of_pc+1,data.shape[1]) # #just neighborhood
    AllneighborhoodClass= data[:,featureVector]
    Alltarget= data[:,1:input.no_of_pc+1]
    #print(AllneighborhoodClass.shape,Alltarget.shape)

    count=0
    for i in range(len(input.spatialcell_unique_clusterid)):
        temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        index=temp[0]
        neighborhoodClass=AllneighborhoodClass[index,:]
        target=Alltarget[index,:]
        positive_interacted_CT= logistic_predicted_interactions[input.spatialcell_unique_clustername[i]]
        #print(i,neighborhoodClass.shape,target.shape,input.spatialcell_unique_clustername[i])
        #print(positive_interacted_CT)
        newindex=[]
        xlabel=[]
        score=[]
        for j in range(len(input.spatialcell_unique_clustername)):
            start=j*input.no_of_pc
            end=start+input.no_of_pc
            for k in range(len(positive_interacted_CT)):
                if positive_interacted_CT[k][0]==input.spatialcell_unique_clustername[j]:
                    #print(i,start,end)
                    xlabel.append(positive_interacted_CT[k][0])
                    score.append(positive_interacted_CT[k][1])
                    for kk in range(start,end):
                        newindex.append(kk)

        neighborhoodClass=neighborhoodClass[:,newindex]
        #print("\n\n\n\n",i,neighborhoodClass.shape,target.shape)
        xlabel=np.array(xlabel)
        score=np.array(score)

        ylabelname=[]
        for k in range(len(xlabel)):
            for j in range(input.no_of_pc):
                ylabelname.append(xlabel[k]+'_s'+'%0.3f'%score[k]+'_PC'+str(j+1))

        count+=neighborhoodClass.shape[0]
        saveoutname=str(input.spatialcell_unique_clusterid[i])+'_'+input.spatialcell_unique_clustername[i]
        coef,intercept,alpha,percent_variance_explained,residual_variance_explained,pv=run_ridge_regression(input,savedir,saveoutname,ylabelname,target,neighborhoodClass,shap_cluster_cutoff)
        #coef_mu,comp_score,coef_std,comp_score_std,alpha=run_ridge_regression(input.seed ,input.lambda_c,input.K_fold,input.n_repeats,target,neighborhoodClass)
        savedata=savedir+'coef'+str(input.spatialcell_unique_clusterid[i])+'.npz'
        np.savez(savedata,coef_mu=coef,intercept=intercept,alpha=alpha,xlabel=xlabel,score=score,Yreg=target,Xreg=neighborhoodClass,pvalue=pv,pve=percent_variance_explained,rve=residual_variance_explained)
        #np.savez(savedata,coef_mu=coef_mu,coef_std=coef_std,comp_score=comp_score,comp_score_std=comp_score_std,alpha=alpha,xlabel=xlabel,score=score)


    print(count)





def run_ridge_regression(input,savedir,saveoutname,ylabelname,target,neighborhoodClass,shap_cluster_cutoff):
    seed=input.seed+1
    parameters = {'alpha':input.lambda_c }
    #'''
    linear_model =Ridge(fit_intercept=True)
    flag=1
    how_many_times_repeat={}
    while(flag):
        seed=seed+1
        sss = RepeatedKFold(n_splits=input.K_fold, n_repeats=input.n_repeats ,random_state=seed)
        gs_grid = GridSearchCV(linear_model, parameters, scoring='r2',cv=sss,n_jobs=input.n_jobs)
        pipe_grid=Pipeline([ ('StandardScaler',StandardScaler()), ('linear_regression_grid',gs_grid)])
        pipe_grid.fit(neighborhoodClass,target)

        LR_grid= pipe_grid.named_steps['linear_regression_grid']
        #print('grid',LR_grid.best_params_,neighborhoodClass.shape, target.shape)
        lambda_c=LR_grid.best_params_['alpha']
        if lambda_c not in how_many_times_repeat:
            how_many_times_repeat[lambda_c]=1
        else:
            how_many_times_repeat[lambda_c]+=1

        #print('Searching hyperparameters ', 'Grid method:', lambda_c)
        for key in how_many_times_repeat:
            if how_many_times_repeat[key]>2:
                flag=0
                #print('alpha regularization found', lambda_c)
    #'''


    save_intercept=[]
    save_coef=[]

    seed=seed+1
    #sss = RepeatedStratifiedKFold(n_splits=K_fold, n_repeats=1 ,random_state=seed) for classification
    sss = RepeatedKFold(n_splits=input.K_fold, n_repeats=input.n_repeats ,random_state=seed) #for regression
    count=0
    for train_index, test_index in sss.split(neighborhoodClass,target):
        x_train,x_test=neighborhoodClass[train_index],neighborhoodClass[test_index]
        y_train,y_test=target[train_index],target[test_index]
        linear_model = Ridge(alpha=lambda_c,fit_intercept=True)
        #linear_model = LinearRegression()
        pipe=Pipeline([('StandardScaler',StandardScaler()), ('ridge_regression',linear_model)])
        pipe.fit(x_train, y_train)
        y_pred=pipe.predict(x_test)
        LR= pipe.named_steps['ridge_regression']
        #coef=_#coef.append(LR.coef_)
        save_intercept.append(LR.intercept_)
        save_coef.append(LR.coef_)

        if count==0:
            save_y_pred=y_pred
            save_y_test=y_test
            save_x_train=x_train
        else:
            save_y_pred=np.vstack((save_y_pred,y_pred))
            save_y_test=np.vstack((save_y_test,y_test))
            save_x_train=np.vstack((save_x_train,x_train))
        count+=1
        #print('AA',K_fold,y_test.shape, y_pred.shape)
        #print(save_y_test.shape,save_y_pred.shape,'\t',LR.intercept_.shape)

    save_intercept=np.array(save_intercept)
    save_coef=np.array(save_coef)
    LRI=np.mean(save_intercept,axis=0)
    LRC=np.mean(save_coef,axis=0)
    #print("BB", LRI.shape,LRC.shape)
    x_train=0
    x_test=0
    y_test=0
    pv=np.zeros(LRC.shape,dtype=float)
    EVS=[]
    for i in range(y_train.shape[1]):
        EVS.append(explained_variance_score(save_y_test[:,i], save_y_pred[:,i]))
        params = np.append(LRI[i],LRC[i,:])
        newX = np.append(np.ones((len(save_x_train),1)), save_x_train, axis=1)
        MSE = (sum((save_y_test[:,i]-save_y_pred[:,i])**2))/(len(newX)-len(newX[0]))
        try:
            var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
        except np.linalg.LinAlgError as e:
            if 'Singular matrix' in str(e):
                var_b=1# your error handling block
            else:
                raise

        sd_b = np.sqrt(var_b)
        ts_b = params/ sd_b
        df = save_x_train.shape[0] - save_x_train.shape[1]
        p_values1 =np.array([[2*(1-scipy.stats.t.cdf(np.abs(i),df-1)) for i in ts_b]])
        pv[i]=p_values1[:,1:]

    train_index=range(target.shape[0])
    test_index=[]
    if input.shap_computation:
        dir1=savedir+'Shapley_Interventional/'
        dir2=savedir+'Shapley_FullConventional/'
        create_directory(dir1)
        create_directory(dir2)
    x_test=[]
    y_test=[]
    if True:
        x_train=neighborhoodClass[train_index]
        y_train=target[train_index]
        linear_model = Ridge(alpha=lambda_c,fit_intercept=True)
        #linear_model = LinearRegression(fit_intercept=True)
        #pipe=Pipeline([('ridge_regression',linear_model)])
        pipe=Pipeline([ ('StandardScaler',StandardScaler()),('ridge_regression',linear_model)])
        #print('\n1',y_train.shape)

        pipe.fit(x_train,y_train )
        yhat=pipe.predict(x_train)
        for i in range(y_train.shape[1]):
            EVS.append(explained_variance_score(y_train[:,i], yhat[:,i]))


        pipe=Pipeline([ ('StandardScaler',StandardScaler()),('ridge_regression',linear_model)])
        y_train=y_train[:,[0]]
        #print('2',y_train.shape)
        pipe.fit(x_train,y_train )
        LR= pipe.named_steps['ridge_regression']
        #yhat=pipe.predict(x_train)
        #intercept=LR.intercept_
        #coef=LR.coef_#coef.append(LR.coef_)
        #print("intercept",intercept[3])
        #print("coef",coef[3])
        '''
        print('11',explained_variance_score(y_train, yhat))
        print('22',explained_variance_score(y_train, yhat,multioutput='raw_values'))
        print('33',explained_variance_score(y_train, yhat,multioutput='variance_weighted'))
        print('44',explained_variance_score(y_train, yhat,multioutput='uniform_average'))

        print('11',r2_score(y_train, yhat))
        print('22',r2_score(y_train, yhat,multioutput='raw_values'))
        print('33',r2_score(y_train, yhat,multioutput='variance_weighted'))
        print('44',r2_score(y_train, yhat,multioutput='uniform_average'))


        print('11',mean_absolute_error(y_test, y_pred))
        print('22',mean_absolute_error(y_test, y_pred,multioutput='raw_values'))
        print('33',mean_absolute_error(y_test, y_pred,multioutput='variance_weighted'))
        print('44',mean_absolute_error(y_test, y_pred,multioutput='uniform_average'))

        print('11',mean_squared_error(y_test, y_pred))
        print('22',mean_squared_error(y_test, y_pred,multioutput='raw_values'))
        print('33',mean_squared_error(y_test, y_pred,multioutput='variance_weighted'))
        print('44',mean_squared_error(y_test, y_pred,multioutput='uniform_average'))
        '''

        if input.shap_computation:
            #explainer = shap.LinearExplainer(LR, x_train)
            explainer = shap.explainers.Linear(LR, x_train,feature_names=ylabelname,feature_perturbation="interventional")
            #explainer = shap.Explainer(LR, x_train,feature_names=ylabelname)
            #shap_values = explainer.shap_values(x_train)
            shap_values = explainer(x_train)

            #shap.waterfall_plot(explainer.expected_value, shap_values[sample_ind], X.iloc[sample_ind], max_display=14)
            clust = shap.utils.hclust(x_train, y_train, linkage="single")
            shap.plots.bar(shap_values, clustering=clust, clustering_cutoff=shap_cluster_cutoff, show=False)
            plt.title("True to the model")
            plt.tight_layout()
            plt.savefig(dir1+saveoutname+'_.png',dpi=300, bbox_inches = "tight")
            plt.close('all')

            explainer = shap.explainers.Linear(LR, x_train,feature_names=ylabelname,feature_perturbation="correlation_dependent")
            shap_values = explainer(x_train)
            shap.plots.bar(shap_values, clustering=clust, clustering_cutoff=shap_cluster_cutoff, show=False)
            plt.title("True to the data")
            plt.tight_layout()
            plt.savefig(dir2+saveoutname+'_.png',dpi=300, bbox_inches = "tight")
            plt.close('all')



        #x_train2 = sm.add_constant(x_train)
        #print("train ",x_train.shape,y_train.shape)
        #est=sm.OLS(y_train[:,3],x_train2).fit()
        #print("summary", est.summary())

        mu=np.mean(y_train,axis=0)
        #print(y_train.shape,mu,yhat.shape)
        total_ss= np.sum((y_train-mu)**2,axis=0)
        #residual_ss=np.sum((y_train-yhat)**2,axis=0)
        #explained_ss= np.sum((yhat-mu)**2,axis=0)
        residual_ss=np.sum((save_y_test-save_y_pred)**2,axis=0)
        explained_ss= np.sum((save_y_pred-mu)**2,axis=0)

        #percent_variance_explained=explained_ss/total_ss
        #residual_variance_explained=residual_ss/total_ss

        #print('\n\n123',EVS,'\t',residual_ss,explained_ss,'\t',percent_variance_explained,residual_variance_explained)
        #cmn.append(confusion_matrix(y_test,y_pred,normalize='true'))


    #index=np.where(pv[3]<0.05)
    #print('pvalue',index)
    print(EVS)
    coef=LRC
    intercept=LRI
    residual_variance_explained=0

    return coef,intercept,lambda_c,EVS,residual_variance_explained,pv


def find_logistic_regression_interacting_score(cmn,coef,CTFeatures,nameOfCellType,logistic_coef_cutoff):
    a=np.diag(cmn)
    #b=np.diag(input.cmn_std)
    goodPredictedCellType=np.argsort(-a)
    largest=np.max(abs(coef))
    normalized_coef=coef/largest
    InteractingCTs=[]
    for k in range(len(a)):
            meanCoefficients=normalized_coef[goodPredictedCellType[k]]
            #stdCoefficients=input.coef_std[goodPredictedCellType[k]]
            highestIndex=np.argsort(-abs(meanCoefficients))
            n=len(highestIndex)
            coeff_of_CT=[]
            name_of_the_coeff=[]
            std_of_coeff=[]
            predictedCT=nameOfCellType[goodPredictedCellType[k]]
            positiveprediction=[]
            negativeprediction=[]
            score=[]
            for i in range(n):
                l=CTFeatures[highestIndex[i]].split()
                temp=''
                for j in range(len(l)):
                    temp+=nameOfCellType[int(l[j][1:])]
                    if j!=(len(l)-1):
                        temp+='--'
                if meanCoefficients[ highestIndex[i]]>logistic_coef_cutoff:
                    positiveprediction.append(temp)
                    score.append(meanCoefficients[ highestIndex[i]])
                else:
                    negativeprediction.append(temp)
            #print(predictedCT,len(score))
            InteractingCTs.append([predictedCT,positiveprediction, score   ])

    logistic_predicted_interactions={}
    for i in range(len(InteractingCTs)):
        cCT=InteractingCTs[i][0]
        nCT=InteractingCTs[i][1]
        Interacting_score=InteractingCTs[i][2]
        for j in range(len(nCT)):
            #print(i,j,cCT,nCT[j],Interacting_score[j],)
            if cCT not in logistic_predicted_interactions:
                logistic_predicted_interactions[cCT]=[[nCT[j],Interacting_score[j]]]
            else:
                logistic_predicted_interactions[cCT].append([nCT[j],Interacting_score[j]])

    return logistic_predicted_interactions


def chooseRightColorBarForPvalue(ax,M,componentlabel):
    if np.mean(M)==0:
        colors = ((1.0, 1.0, 0.0),(1.0, 1.0, 0.0))
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        b=snn.heatmap(M,cmap=cmap,yticklabels=componentlabel,ax=ax)
        colorbar=b.collections[0].colorbar
        colorbar.set_ticks([0])
        colorbar.set_ticklabels(['pv>=0.05'])
    elif np.mean(M)==1:
        colors = ((1, 0.0, 1.0),(1, 0.0, 1.0))
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        b=snn.heatmap(M,cmap=cmap,yticklabels=componentlabel,ax=ax)
        colorbar=b.collections[0].colorbar
        colorbar.set_ticks([1])
        colorbar.set_ticklabels(['pv<0.05'])
    else:
        colors = ((1.0, 1.0, 0.0), (1, 0.0, 1.0))
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        b=snn.heatmap(M,cmap=cmap,yticklabels=componentlabel,ax=ax)
        colorbar=b.collections[0].colorbar
        colorbar.set_ticks([0.25,0.75])
        colorbar.set_ticklabels(['pv>=0.05', 'pv<0.05'])

    return b



def plot_results(savedir,maindir,radius,input):

    for i in range(len(input.spatialcell_unique_clusterid)):
        filename=input.spatialcell_unique_clustername[i]
        temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        index=temp[0]
        savedata=savedir+'coef'+str(input.spatialcell_unique_clusterid[i])+'.npz'
        data=np.load(savedata,allow_pickle=True)
        coef_mu=data['coef_mu']
        intercept=data['intercept']
        pve=data['pve'] # percentage variance explanined
        rve=data['rve'] # residual variance explained
        pvalue=data['pvalue']

        #coef_std=data['coef_std']
        #comp_score=data['comp_score']
        #comp_score_std=data['comp_score_std']

        alpha=data['alpha']
        xlabel=data['xlabel']
        score=data['score']

        componentlabel=[]
        for j in range(input.no_of_pc):
            componentlabel.append('PC '+str(j+1))

        #print(pve,rve)
        percentVE=''
        percentRE=''
        for j in range(len(pve)):
            if j!=0:
                percentVE+=', '
                percentRE+=', '
            percentVE+='%0.2e'%pve[j]
            #percentRE+='%0.1f'%rve[j]

        ylabelname=[]
        for k in range(len(xlabel)):
            for j in range(input.no_of_pc):
                ylabelname.append(xlabel[k]+'_s'+'%0.3f'%score[k]+'_PC'+str(j+1))


        pvalue=pvalue<0.05
        # significant pvalue is 1  less than 0.05
        #print("anki",coef_mu.shape,rve.shape,pvalue.shape)
        fig=plt.figure(figsize=(10,5))
        gs = fig.add_gridspec(ncols=1, nrows=2, height_ratios=[3, 1])
        ax0=fig.add_subplot(gs[0])
        ax1=fig.add_subplot(gs[1])

        b=snn.heatmap(coef_mu,xticklabels=ylabelname,yticklabels=componentlabel,ax=ax0)
        xlabels= b.get_xticklabels()
        b.set_xticklabels(xlabels,size = 8,rotation=90)
        b.set_ylabel('Principal components')
        #_, ylabels= plt.yticks()
        #b.set_yticklabels(ylabels, size = 5)
        #b.set_title(filename+r',$\alpha$='+str(alpha)+',VE='+percentVE+',RE='+percentRE)
        b.set_title(filename+r',$\alpha$='+str(alpha)+',EVS='+percentVE,fontsize=6)

        b=chooseRightColorBarForPvalue(ax1,pvalue,componentlabel)
        #b.xaxis.set_label_position('top')
        b.xaxis.tick_top()
        xlabels= b.get_xticks()
        b.set_xticklabels(xlabels,size = 0)
        #b.axes.get_xaxis().set_visible(False)
        #b.set_ylabel('Principal components')
        #_, ylabels= plt.yticks()
        #b.set_yticklabels(ylabels, size = 5)
        #plt.title(filename+r',$\alpha$='+str(alpha)+',VE='+percentVE+',RE='+percentRE)
        #b.set_title('pvalue significant ')

        fig.tight_layout()
        fig.savefig(savedir+'coeff_matrix_'+str(input.spatialcell_unique_clusterid[i])+'_'+filename+'.png',dpi=300)
        plt.close()


    ylabelname=[]
    for i in range(len(input.spatialcell_unique_clustername)):
        for j in range(input.no_of_pc):
            ylabelname.append(input.spatialcell_unique_clustername[i]+'_'+'PC'+str(j+1))


    fig,axs=plt.subplots(1,1,figsize=(10,10))
    name=maindir+'Principal_component_feature_matrix'+str(input.no_of_pc)+'.csv'
    data=np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
    Feature=data[:,6:data.shape[1]]
    index=np.argsort(input.annotation_spatial_cluster_id)
    snn.heatmap(np.log(Feature[index,:]),xticklabels=ylabelname)
    fig.tight_layout()
    fig.savefig(maindir+'Feature_matrix_PC'+str(input.no_of_pc)+'.png',dpi=300)
    fig.clf()


def main(myinput):
    celltypefname=myinput.celltypefname
    clusterfname=myinput.clusterfname
    sct_ad_sp=myinput.ad_sp
    sct_ad_sc=myinput.ad_sc
    full_ad_sc=myinput.full_ad_sc
    input_spatial_analysis_dir=myinput.input_spatial_analysis_dir
    inputRadius=myinput.inputRadius
    no_of_pc=myinput.no_of_pc
    strategy=myinput.strategy
    maindir1=myinput.outdir
    gene_set_names=myinput.gene_set_names
    database=myinput.database
    pathwayorganism=myinput.pathwayorganism
    pathwayCutoff=myinput.pathwayCutoff
    shap_cluster_cutoff=myinput.shap_cluster_cutoff
    LRdb=myinput.LRdb
    pathway_plot=myinput.pathway_plot
    coeff_cutoff_for_rid_reg=myinput.coeff_cutoff_for_rid_reg
    seed=myinput.seed
    lambda_c=myinput.lambda_c
    K_fold=myinput.K_fold
    n_repeats=myinput.n_repeats
    n_jobs=myinput.n_jobs
    shap_computation=myinput.shap_computation
    logistic_coef_cutoff=myinput.logistic_coef_cutoff
    gene_correlation_plots_and_sheets=myinput.gene_correlation_plots_and_sheets



    df=pd.read_csv(clusterfname)
    data=df.to_numpy()
    annotation_singlecell_barcode_id=data[:,0]
    annotation_singlecell_cluster_id=data[:,1]

    with open(celltypefname,'r') as f:
        cont = f.read()
        lines=cont.split('\n')
        singlecell_unique_clustername=[]
        singlecell_unique_clusterid=[]
        for i in range(len(lines)):
            l=lines[i].split('\t')
            if len(l)>1:
                name=l[1].replace('/','_')
                name=name.replace(' ','_')
                name=name.replace('"','')
                name=name.replace("'",'')
                name=name.replace(')','')
                name=name.replace('(','')
                name=name.replace('+','p')
                name=name.replace('-','n')
                name=name.replace('.','')
                singlecell_unique_clustername.append(name)
                singlecell_unique_clusterid.append(int(l[0]))

    d={}
    for i in range(len(singlecell_unique_clusterid)):
        d[singlecell_unique_clusterid[i]]=singlecell_unique_clustername[i]

    annotation_singlecell_celltypename=[]
    for i in range(len(annotation_singlecell_cluster_id)):
        annotation_singlecell_celltypename.append(d[annotation_singlecell_cluster_id[i]])
    annotation_singlecell_celltypename=np.array(annotation_singlecell_celltypename)


    print('sc1 annotation_singlecell_cluster_id',len(annotation_singlecell_cluster_id))
    print('sc2 annotation_singlecell_barcode_id',len(annotation_singlecell_barcode_id))
    print('sc3 annotation_singlecell_celltypename',len(annotation_singlecell_celltypename))
    print('sc4 singlecell_unique_clustername', len(singlecell_unique_clustername))


    # load spatial dat

    sp_genename=sct_ad_sp.var_names.to_numpy()
    sc_genename=sct_ad_sc.var_names.to_numpy()
    index_sp,index_sc=find_index(sp_genename,sc_genename)
    print('common genes between sc and sp',len(index_sp),len(index_sc))
    ad_sp_ori=sct_ad_sp[:,index_sp].copy()
    ad_sc_ori=sct_ad_sc[:,index_sc].copy()



    for radius in inputRadius:
        celltypeFilename=input_spatial_analysis_dir+'BiologicalNameOfCT.dat'
        clusterFilename=input_spatial_analysis_dir+'save_clusterid_'+str(radius)+'.csv'

        annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid=read_spatial_data(clusterFilename,celltypeFilename)

        print('\nsp1 annotation_spatial_cluster_id',len(annotation_spatial_cluster_id))
        print('sp2 annotation_spatial_barcode_id',len(annotation_spatial_barcode_id))
        print('sp3 annotation_spatial_celltypename',len(annotation_spatial_celltypename))
        print('sp4 spatialcell_unique_clustername', len(spatialcell_unique_clustername))

        neighbors=pickle.load( open(input_spatial_analysis_dir+'save_neighbors_'+str(radius)+'.p', "rb" ) )
        maindir=maindir1+str(radius)+'/'
        create_directory(maindir)
        outputname=maindir+'Principal_component_feature_matrix'+str(no_of_pc)+'.csv'
        inputdata={}
        inputdata['no_of_pc']=no_of_pc
        inputdata['outputname']=outputname


        fname=input_spatial_analysis_dir+strategy+'/save_numpy_array_'+str(radius)+'.npz'
        data=np.load(fname,allow_pickle=True)
        logistic_coef=data['coef']
        logistic_cmn=data['cmn']
        logistic_cmn_std=data['cmn_std']
        logistic_coef_std=data['coef_std']
        logistic_CTFeatures=data['CTFeatures']
        #f=open(input_spatial+'BiologicalNameOfCT.dat')
        f=open(celltypeFilename)
        nameOfCellType={}
        for line in f:
            l=line[0:-1].split('\t')
            nameOfCellType[int(l[0])]=l[1]

        logistic_predicted_interactions=find_logistic_regression_interacting_score(logistic_cmn,logistic_coef,logistic_CTFeatures,nameOfCellType,logistic_coef_cutoff)

        #for key in logistic_predicted_interactions:
        #    print(key,logistic_predicted_interactions[key])



        inputdata['ad_sp']=ad_sp_ori #sct_ad_sp
        inputdata['ad_sc']=ad_sc_ori#sct_ad_sc#
        inputdata['annotation_spatial_cluster_id']=annotation_spatial_cluster_id
        inputdata['annotation_spatial_barcode_id']=annotation_spatial_barcode_id
        inputdata['annotation_spatial_celltypename']=annotation_spatial_celltypename
        inputdata['spatialcell_unique_clustername']=spatialcell_unique_clustername
        inputdata['spatialcell_unique_clusterid']=spatialcell_unique_clusterid

        inputdata['annotation_singlecell_cluster_id']=annotation_singlecell_cluster_id
        inputdata['annotation_singlecell_barcode_id']=annotation_singlecell_barcode_id
        inputdata['annotation_singlecell_celltypename']=annotation_singlecell_celltypename
        inputdata['singlecell_unique_clustername']=singlecell_unique_clustername
        inputdata['singlecell_unique_clusterid']=singlecell_unique_clusterid
        inputdata['neighbors']=neighbors
        create_directory(maindir+'NMF_output/')
        inputdata['nmf_output']=maindir+'NMF_output/'
        inputdata['seed']=seed
        inputdata['lambda_c']=lambda_c
        inputdata['K_fold']=K_fold
        inputdata['n_repeats']=n_repeats
        inputdata['n_jobs']=n_jobs
        inputdata['shap_computation']=shap_computation
        inputdata['shap_cluster_cutoff']=shap_cluster_cutoff

        inputdata['logistic_coef_cutoff']=logistic_coef_cutoff
        inputdata['coeff_cutoff_for_rid_reg']=coeff_cutoff_for_rid_reg
        inputdata['gene_set_names']=gene_set_names
        inputdata['database']=database
        inputdata['pathwayorganism']=pathwayorganism
        inputdata['pathwayCutoff']=pathwayCutoff
        inputdata['pathway_plot']=pathway_plot

        input=SimpleNamespace(**inputdata)

        flag=1
        if os.path.isfile(outputname):
            filesize = os.path.getsize(outputname)
            if filesize>0: #If file is already exist and have size greater than 0 then no need to run again. It will save some time if you want to run it again with different parameters
                flag=0

        if flag==1:
            #print(sct_ad_sc)
            pc_of_sp_clusterid,PCA_of_sc_cluster_accordingto_spatial_clusterid=compute_PC_space(input,full_ad_sc)
            # full_ad_sc use in only find_PC_of_invidualCluster_in_SC function
            # ideally it should be sctransform way of normalized matrix equivalent to sct_ad_sc but
            # if not then need to do perform scaling HVG etc
            pickle.dump(PCA_of_sc_cluster_accordingto_spatial_clusterid,open(maindir+'PCA_of_sc_cluster'+str(no_of_pc)+'.p', 'wb'))
            inputdata['pc_of_sp_clusterid']=pc_of_sp_clusterid
            input=SimpleNamespace(**inputdata)
            makePCneighboorhoodFeatureMatrix(input)

        savedir=maindir+'/save_numpy_array_PC'+str(no_of_pc)+'/'
        create_directory(savedir)
        existornot=savedir+'coef'+str(input.spatialcell_unique_clusterid[-1])+'.npz'
        if os.path.isfile(existornot):
            #pass
            plot_results(savedir,maindir,radius,input)
        else:
            model_linear_regression(input,savedir,logistic_predicted_interactions)
            plot_results(savedir,maindir,radius,input)

        totalLRpairs=read_LigRecDb(LRdb)
        find_canonical_pathways_in_interacting_cell_types(input,maindir,savedir,totalLRpairs)



        if gene_correlation_plots_and_sheets:
            made_excel_sheet_for_gene_correlation(maindir,input)



def made_excel_sheet_for_gene_correlation(maindir,input):
    workbook = xlsxwriter.Workbook(maindir+'gene_correlation'+str(input.no_of_pc)+'.xlsx')

    worksheetAvgGeneExp= workbook.add_worksheet('avg gene exp')
    worksheetFullGene=[]
    for i in range(input.no_of_pc):
        worksheetFullGene.append( workbook.add_worksheet('HVG PC'+str(i+1)))
    worksheetSpatialGene=[]
    for i in range(input.no_of_pc):
        worksheetSpatialGene.append( workbook.add_worksheet('spatial gene PC'+str(i+1)))


    PCA_of_sc_cluster_accordingto_spatial_clusterid=pickle.load(open(maindir+'PCA_of_sc_cluster'+str(input.no_of_pc)+'.p', 'rb'))

    outputFolder=maindir+'geneCorr'+str(input.no_of_pc)+'/'
    create_directory(outputFolder)


    genenames=sorted(list(input.ad_sp.var_names.to_numpy()))
    n=len(input.spatialcell_unique_clustername)

    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        #for i in range(len(cname)):
        CC_corr,CC_PCA,gene,CC_meanExpression=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        #print(CC_corr.shape)
        #print(CC_PCA.shape)
        #print(gene.shape)
        #print(CC_meanExpression.shape)
        #print('\n\n',input.spatialcell_unique_clustername[i])

        #print(cname[i][1])
        worksheetrow=0
        worksheetAvgGeneExp.write(worksheetrow,3*i,input.spatialcell_unique_clustername[i])
        for j in range(input.no_of_pc):
            worksheetFullGene[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
            worksheetSpatialGene[j].write(worksheetrow,(input.no_of_pc+2)*i,input.spatialcell_unique_clustername[i])
        worksheetrow+=1
        fixvalue=worksheetrow


        index=np.argsort(-CC_meanExpression)
        for j in range(len(index)):
            worksheetAvgGeneExp.write(j+2,3*i+1,CC_meanExpression[index[j]])
            worksheetAvgGeneExp.write(j+2,3*i,gene[index[j]])



        fig,(ax)=plt.subplots(1,1,figsize=(8,6))
        ax.plot(CC_corr[:,0],CC_corr[:,1],'.',markersize=1)

        headersave_full=[]
        headersave_common=[]
        sort_full=[]
        sort_common=[]
        for k in range(input.no_of_pc):
            sort_full.append([])
            sort_common.append([])
        #print(CC_corr)
        for j in range(len(CC_corr)):
            ind=~np.isnan(CC_corr[j]).any(axis=0)
            #print(i,ind)
            if ind==True:
                ax.text(CC_corr[j,0],CC_corr[j,1],gene[j],fontsize=5)
                header=[gene[j]]
                for k in range(input.no_of_pc):
                    sort_full[k].append(abs(CC_corr[j,k]))
                    header.append(CC_corr[j,k])
                headersave_full.append(header)
                if gene[j] in genenames:
                    headersave_common.append(header)
                    for k in range(input.no_of_pc):
                        sort_common[k].append(abs(CC_corr[j,k]))


        #print(headersave_common)

        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_title(input.spatialcell_unique_clustername[i])
        fig.tight_layout()
        fig.savefig(outputFolder+'correlation_'+input.spatialcell_unique_clustername[i]+'.png',bbox_inches='tight',dpi=300)
        plt.close('all')


        for k in range(input.no_of_pc):
            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_full[k]))
            #print('ankitankur',input.no_of_pc,len(indsort))
            for rj in range(len(indsort)):
                header=headersave_full[indsort[rj]]
                for ri in range(len(header)):
                    worksheetFullGene[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri])
                worksheetrow+=1

            worksheetrow=fixvalue
            indsort=np.argsort(-np.array(sort_common[k]))
            for rj in range(len(indsort)):
                header=headersave_common[indsort[rj]]
                for ri in range(len(header)):
                    #print("ABC",i,k,rj,ri)
                    worksheetSpatialGene[k].write(worksheetrow,(input.no_of_pc+2)*i+ri,header[ri])
                worksheetrow+=1

    workbook.close()




def find_canonical_pathways_in_interacting_cell_types(input,maindir,savedir,totalLRpairs):

    database=input.database
    pathwayorganism=input.pathwayorganism
    pathwayCutoff=input.pathwayCutoff
    pathway_plot=input.pathway_plot
    coeff_cutoff_for_log_reg=input.logistic_coef_cutoff
    coeff_cutoff_for_rid_reg=input.coeff_cutoff_for_rid_reg
    gene_set_names=input.gene_set_names

    PCA_of_sc_cluster_accordingto_spatial_clusterid=pickle.load(open(maindir+'PCA_of_sc_cluster'+str(input.no_of_pc)+'.p', 'rb'))
    n=len(input.spatialcell_unique_clustername)


    workbook = xlsxwriter.Workbook(maindir+'/Lig_and_Rec_enrichment_in_interacting_celltypes'+str(input.no_of_pc)+'.xlsx')
    fout=open(maindir+'/PC_loading_prediction_'+str(input.no_of_pc)+'.dat','w')
    worksheet = workbook.add_worksheet('LR enrichment')
    worksheetrow=0
    header=['Id','A','B','localized score','PC(A)','PC(B)', 'Coeff' ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)']
    for ri in range(len(header)):
        worksheet.write(worksheetrow,ri,header[ri])
    worksheetrow+=1

    create_directory(maindir+'/PCcomponent_wise_pathway_figures')


    d={}
    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        clname=input.spatialcell_unique_clustername[i]
        d[clname]=clid

    for i in range(n):
        clid=input.spatialcell_unique_clusterid[i]
        CC_corr,CC_PCA,CC_gene,CC_meanExpression=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
        #print(input.spatialcell_unique_clustername[i],gene[0:10],gene[-10:-1])
        CC_celltype_name=input.spatialcell_unique_clustername[i]
        #temp=np.where(input.spatialcell_unique_clusterid[i]==input.annotation_spatial_cluster_id)
        #index=temp[0]
        savedata=savedir+'coef'+str(input.spatialcell_unique_clusterid[i])+'.npz'
        data=np.load(savedata,allow_pickle=True)
        coef_mu=data['coef_mu']
        intercept=data['intercept']
        pve=data['pve'] # variance explained
        rve=data['rve'] # residual variance explained
        Xreg=data['Xreg']
        Yreg=data['Yreg']
        pvalue=data['pvalue']
        '''
        coef_std=data['coef_std']
        comp_score=data['comp_score']
        comp_score_std=data['comp_score_std']
        '''
        alpha=data['alpha']
        NC_celltype_name=data['xlabel']
        score=data['score'] # this score is from logistic regression
        largest=np.max(abs(coef_mu))
        normalized_ridge_coef=coef_mu/largest

        ylabelname=[]
        componentlabel=[]
        for j in range(input.no_of_pc):
            ylabelname.append('CC_'+CC_celltype_name+'_PC'+str(j+1))
            componentlabel.append('PC'+str(j+1))

        for k in range(len(NC_celltype_name)):
            if score[k]>coeff_cutoff_for_log_reg:
                #in ylabelname first (# of pc) is the central cell type
                #and remaining are (# of pc) from the negihborhood cell type
                if CC_celltype_name!=NC_celltype_name[k]:
                    for j in range(input.no_of_pc):
                        ylabelname.append('NC_'+NC_celltype_name[k]+'_s'+'%0.3f'%score[k]+'_PC'+str(j+1))

        pc_index_nc=[]
        for k in range(len(NC_celltype_name)):
            for j in range(input.no_of_pc):
                pc_index_nc.append(j)

        #print(cc_celltype_name,coef_mu.shape,score.shape,len(ylabelname))
        fout.write(CC_celltype_name+'\n')

        #print('god',normalized_ridge_coef.shape, ylabelname)

        #normalized_ridge_coef  noofPC x (noofPC x +ve coff in log reg)


        interaction_id=0
        pathway_already_plot={}
        for k in range(normalized_ridge_coef.shape[0]):
            #k is PC of central cell type
            for j in range(normalized_ridge_coef.shape[1]):
                interaction_id+=1
                index=math.floor(j/input.no_of_pc)
                #index is the id neighboring cell type
                #if abs(normalized_ridge_coef[k,j])>coeff_cutoff_for_rid_reg:
                if (pvalue[k,j]<0.05)&(abs(normalized_ridge_coef[k,j])>coeff_cutoff_for_rid_reg):
                    if score[index]>coeff_cutoff_for_log_reg:
                        NC_corr,NC_PCA,NC_gene,NC_meanExpression=PCA_of_sc_cluster_accordingto_spatial_clusterid[d[NC_celltype_name[index]]]
                        if str(CC_gene)==str(NC_gene):
                            print('Good')
                        #print(PCA.shape,NH_PCA.shape,k,pc_index_nc[j],nc_celltype_name[index],d[nc_celltype_name[index]])
                        top_genes_in_CC,top_genes_in_NC,genesWithUP,genesWithDown,Found1,Found2=find_fold_change(CC_corr,NC_corr,CC_gene,k,pc_index_nc[j],totalLRpairs)
                        common_genes=list(set(top_genes_in_CC).intersection(set(top_genes_in_NC)))

                        print(len(top_genes_in_CC),len(top_genes_in_NC),len(common_genes))


                        for ele in range(len(Found1)):
                            header=[str(i)+'-'+str(interaction_id),CC_celltype_name+'(cc)',NC_celltype_name[index]+'(nc)',score[index],k+1,1+pc_index_nc[j],normalized_ridge_coef[k,j] ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','Receptor(A)','Ligand(B)','GeneCor(Rec)','GeneCor(Lig)']
                            header[7]=Found1[ele][0][0]
                            header[8]=Found1[ele][1][0]
                            header[9]=Found1[ele][0][1]
                            header[10]=Found1[ele][1][1]
                            for ri in range(11):
                                worksheet.write(worksheetrow,ri,header[ri])
                            worksheetrow+=1

                        for ele in range(len(Found2)):
                            header=[str(i)+'-'+str(interaction_id),NC_celltype_name[index]+'(nc)',CC_celltype_name+'(cc)',score[index],1+pc_index_nc[j],k+1,normalized_ridge_coef[k,j] ,'Ligand(A)','Receptor(B)','GeneCor(Lig)','GeneCor(Rec)','Receptor(A)','Ligand(B)','GeneCor(Rec)','GeneCor(Lig)']
                            header[7]=Found2[ele][0][0]
                            header[8]=Found2[ele][1][0]
                            header[9]=Found2[ele][0][1]
                            header[10]=Found2[ele][1][1]
                            for ri in range(11):
                                worksheet.write(worksheetrow,ri,header[ri])
                            worksheetrow+=1

                        print('\n\nankit',CC_celltype_name,'CC-PC'+str(k+1)+'\t'+'%0.3f'%(score[index])+'\tNC-PC'+str(1+pc_index_nc[j])+'\t'+
                        NC_celltype_name[index]+'\t%0.3f'%(normalized_ridge_coef[k,j])+'\t'+str(interaction_id))
                        print('anusha',normalized_ridge_coef.shape, intercept.shape, Xreg[:,j].shape,Yreg[:,k].shape,i,k,j,index,'\t',pc_index_nc[j], CC_PCA.shape, NC_PCA.shape)
                        #save_dat_for_checking(maindir,input,Yreg[:,k],Xreg[:,j],CC_corr[:,[k]], NC_corr[:,[pc_index_nc[j]]],CC_gene,CC_celltype_name,str(k+1), NC_celltype_name[index],str(1+pc_index_nc[j]), coef_mu[k,j],normalized_ridge_coef[k,j], intercept[k], str(interaction_id))
                        #print('CC',len(top_genes_in_CC),top_genes_in_CC[0:1],genesWithUP[0:1])
                        #print('NC',len(top_genes_in_NC),top_genes_in_NC[0:1],genesWithDown[0:1])
                        savename=maindir+'PCcomponent_wise_pathway_figures/'
                        if len(genesWithUP)>0:
                            #name_cc='_CC_'+CC_celltype_name.replace('.','_')+'_PC'+str(k+1)#+'_ID'+str(interaction_id)
                            name_cc=CC_celltype_name.replace('.','_')+'_PC'+str(k+1)#+'_ID'+str(interaction_id)
                            if name_cc not in pathway_already_plot:
                                pathway_already_plot[name_cc]=1
                                flag=1
                            else:
                                flag=0
                            if (len(top_genes_in_CC)>0)&(pathway_plot)&(flag==1):
                                pathway_analysis(savename,gene_set_names,database,pathwayorganism,pathwayCutoff,name_cc,top_genes_in_CC,CC_gene,CC_meanExpression)
                        if len(genesWithDown)>0:
                            #name_cc='_NC_'+NC_celltype_name[index].replace('.','_')+ '_PC'+str(1+pc_index_nc[j])#+'_ID'+str(interaction_id)
                            name_cc=NC_celltype_name[index].replace('.','_')+ '_PC'+str(1+pc_index_nc[j])#+'_ID'+str(interaction_id)
                            if name_cc not in pathway_already_plot:
                                pathway_already_plot[name_cc]=1
                                flag=1
                            else:
                                flag=0
                            if (len(top_genes_in_NC)>0)&(pathway_plot)&(flag==1):
                                pathway_analysis(savename,gene_set_names,database,pathwayorganism,pathwayCutoff,name_cc,top_genes_in_NC,NC_gene,NC_meanExpression)


                        #if (len(common_genes)>0)&(pathway_plot):
                        #    name_cc='_Common_'+NC_celltype_name[index].replace('.','_')+'_ID'+str(interaction_id)
                        #    pathway_analysis(savename,gene_set_names,database,pathwayorganism,pathwayCutoff,name_cc,common_genes,NC_gene,NC_meanExpression)


                        fout.write('CC-PC'+str(k+1)+'\t'+'%0.3f'%(score[index])+'\tNC-PC'+str(1+pc_index_nc[j])+'\t'+NC_celltype_name[index]+'\t%0.3f'%(normalized_ridge_coef[k,j])+'\t'+str(interaction_id))
                        fout.write('\n')
                        '''
                        genesWithUP=''
                        genesWithDown=''
                        for ii in top_genes_in_CC:
                            genesWithUP+=','+ii
                        for ii in top_genes_in_NC:
                            genesWithDown+=','+ii
                        '''

                        fout.write('CC'+str(genesWithUP)+'\n')
                        fout.write('NC'+str(genesWithDown)+'\n')
                        fout.write('\n')
        fout.write('\n\n')
    workbook.close()


def pathway_analysis(savename,gene_set_names,database,pathwayorganism,pathwayCutoff,titlename,g1,background_geneName,background_expression):
    for i in range(len(database)):
        titlename1=titlename+'_'+database[i]
        finalsavename=savename+titlename1
        titlename1=titlename1+'_#G='+str(len(g1))
        #print("tt",titlename)
        #print("fin",finalsavename)
        enr_res1 = gseapy.enrichr(gene_list=g1,organism=pathwayorganism,gene_sets=database[i], description='pathway',cutoff = pathwayCutoff)
        #enr_res1 = gseapy.enrichr(gene_list=g1,organism='Mouse',gene_sets=background_model,description='pathway',cutoff = 0.5)
        finalsavename.replace(' ','_')
        gseapy.barplot(enr_res1.res2d,title=titlename1,ofname=finalsavename)#database[i]+titlename



def find_fold_change(PCA,NH_PCA,gene,CCPC,NCPC,totalLRpairs):

    #listofallLR=['Acaca', 'Acvr2a', 'Adamts13', 'Adamts3', 'Adcy2', 'Adcy5', 'Alb', 'Angpt1', 'Anxa1', 'Apoa2', 'Bcl2', 'Bmp4', 'Bmp5', 'Bmpr1b', 'C3', 'C5ar1', 'Cachd1', 'Cacna1c', 'Cadm1', 'Camp', 'Ccbe1', 'Ccl17', 'Ccl2', 'Ccl22', 'Ccl24', 'Ccl3', 'Ccl4', 'Ccl5', 'Ccl7', 'Ccl8', 'Ccr1', 'Ccr2', 'Ccr7', 'Ccr9', 'Ccrl2', 'Cd19', 'Cd79a', 'Cd83', 'Cd8b1', 'Cdh11', 'Cftr', 'Chrdl1', 'Clec3b', 'Col1a1', 'Col3a1', 'Col6a3', 'Csf1', 'Csf3r', 'Cx3cr1', 'Cxcl1', 'Cxcl10', 'Cxcl12', 'Cxcl13', 'Cxcl2', 'Cxcr2', 'Cxcr4', 'Cytl1', 'Dapk1', 'Dcc', 'Dcn', 'Ddr2', 'Ecm1', 'Eda', 'Edar', 'Edn1', 'Ednra', 'Ednrb', 'Efemp1', 'Efna3', 'Efna5', 'Efnb1', 'Egfr', 'Epcam', 'Epha3', 'Epha7', 'Ephb1', 'Erbb4', 'Ereg', 'Esr1', 'Fcgr4', 'Fgf12', 'Fgf13', 'Fgf2', 'Fgfr2', 'Flrt2', 'Fn1', 'Gdf15', 'Gdf2', 'Gfra1', 'Gpc3', 'Gpc4', 'Gpc6', 'Grin2b', 'Grm7', 'Gzma', 'Hgf', 'Hmox1', 'Hpx', 'Htr2a', 'Icam1', 'Ifitm1', 'Igf1r', 'Il10', 'Il12b', 'Il1b', 'Il1f9', 'Il1rap', 'Il1rapl1', 'Il1rn', 'Il2rb', 'Il6', 'Il7r', 'Insr', 'Itga8', 'Itgb8', 'Lama1', 'Lama2', 'Lamb1', 'Lcn2', 'Lgals3', 'Lgr5', 'Lgr6', 'Lpl', 'Lrp2', 'Lrrc4c', 'Ltbp1', 'Mmp3', 'Mmp8', 'Mmp9', 'Msmp', 'Ncam1', 'Ngf', 'Npr3', 'Nrg1', 'Nrg2', 'Nrxn1', 'Nrxn3', 'Ntn1', 'Ntn4', 'Ntng2', 'Ntrk3', 'Nts', 'Nxph1', 'Osm', 'Pard3', 'Pdgfc', 'Pdgfd', 'Pdgfrb', 'Pf4', 'Pglyrp1', 'Plcb1', 'Plscr4', 'Plxna4', 'Postn', 'Prkca', 'Prkce', 'Prkd1', 'Ptgds', 'Ptgs2', 'Pth1r', 'Pthlh', 'Ptprd', 'Ptprj', 'Ptprk', 'Ptprm', 'Rarres2', 'Reln', 'Rgs7', 'Rnf43', 'Robo1', 'Robo2', 'Ror1', 'Rspo1', 'Rspo3', 'S100a4', 'S100a6', 'Sct', 'Sema3a', 'Serpine1', 'Shank2', 'Siglece', 'Slit2', 'Slit3', 'Slpi', 'Spp1', 'Stab2', 'Stk39', 'Thbs1', 'Thbs4', 'Tnc', 'Tnf', 'Tnfrsf11b', 'Tnfrsf4', 'Unc5c', 'Vegfc', 'Vipr1', 'Vwf', 'Xcl1', 'Xcr1']

    listofallLR={}
    uniqueLRpairs={}
    for i in range(len(totalLRpairs)):
        l=totalLRpairs[i][0]
        r=totalLRpairs[i][1]
        listofallLR[l]=1
        listofallLR[r]=1
        name=l+'--'+r
        if name not in uniqueLRpairs:
            uniqueLRpairs[name]=1

    #print(PCA.shape,NH_PCA.shape)
    first=PCA[:,CCPC]
    second=NH_PCA[:,NCPC]
    ind1=np.argsort(-abs(first))
    ind2=np.argsort(-abs(second))

    cc_genes=[]
    cc_genes2=[]
    cc_genes5=[]

    nc_genes=[]
    nc_genes2=[]
    nc_genes5=[]


    no_of_gene_for_pathway_analysis=5
    for i in range(no_of_gene_for_pathway_analysis):
            cc_genes5.append([gene[ind1[i]],'%0.2f'%first[ind1[i]]])

    for i in range(no_of_gene_for_pathway_analysis):
            nc_genes5.append([gene[ ind2[i] ],'%0.2f'%second[ ind2[i] ]])


    for i in range(len(ind1)):
        if abs(first[ind1[i]])>0.1:
        #if (first[ind1[i]])<-0.4:
            cc_genes.append(gene[ind1[i]])
            if gene[ind1[i]].upper() in listofallLR:
                cc_genes2.append([gene[ind1[i]],'%0.2f'%first[ind1[i]]])


    for i in range(len(ind2)):
        if abs(second[ind2[i]])>0.1:
        #if (second[ind2[i]])<-0.4:
            nc_genes.append(gene[ind2[i]])
            if gene[ind2[i]].upper() in listofallLR:
                nc_genes2.append([gene[ ind2[i] ],'%0.2f'%second[ ind2[i] ]])

    #cc_genes=cc_genes[0:50]
    #nc_genes=nc_genes[0:50]
    #print(second[ind2])
    #print(len(nc_genes2), len(cc_genes2))

    Found1=[]
    Found2=[]
    for i in range(len(cc_genes2)):
        cc=cc_genes2[i][0].upper()
        for j in range(len(nc_genes2)):
            nc=nc_genes2[j][0].upper()
            name1=cc+'--'+nc # lig in CC and rec in NC
            name2=nc+'--'+cc # lig in NC and rec in CC
            if name1 in uniqueLRpairs:
                Found1.append([cc_genes2[i],nc_genes2[j]])  # lig in CC and rec in NC
            if name2 in uniqueLRpairs:
                Found2.append([nc_genes2[j],cc_genes2[i],])  # lig in NC and rec in CC

    #print(Found1,Found2)
    #print(len(Found1),len(Found2))

    return cc_genes, nc_genes,cc_genes5,nc_genes5,Found1,Found2



def read_LigRecDb(contdb):
    #f=open('sort_3_db_L_R_high_confident.dat','r')
    totalLRpairs=[]
    for j in range(len(contdb)):
        l=contdb[j][0:-1].split()
        if [l[0], l[1] ] not in totalLRpairs:
            totalLRpairs.append( [l[0].upper(), l[1].upper() ])

    return totalLRpairs


def save_dat_for_checking(maindir,input,yreg, xreg, CC_corr, NC_corr,g_cc,cc_celltype_name,PC_cc,nc_celltype_name, PC_nc, coeff,normalize,intercept, filename):

    savecorrelation=maindir+'/scatterPlot'+str(input.no_of_pc)+'/'
    create_directory(savecorrelation)

    myname=cc_celltype_name+'_'+PC_cc+'_nc_'+nc_celltype_name+'_'+PC_nc+'_'+filename
    yname='CC PC'+  str(PC_cc)+' '+cc_celltype_name
    xname='NC PC'+  str(PC_nc)+' '+nc_celltype_name
    #if xname!=yname:
    if True:
    #if (xname!=yname)&(cc_celltype_name==nc_celltype_name):
        '''
        f1=open(savecorrelation+myname +'.dat','w')
        for i in range(len(xreg)):
            f1.write(str(xreg[i])+','+str(yreg[i])+'\n')
        f1.close()
        f1=open(savecorrelation+'PCA'+myname + '.dat','w')
        for i in range(len(g_cc)):
            f1.write(str(pca_cc[i])+','+str(pca_nc[i])+','+g_cc[i]+','+g_nc[i]+'\n')
        f1.close()
        '''

        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12

        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

        fig,ax=plt.subplots(1,2,figsize=(5,2.4))
        #plt.scatter(first[index],second[index],c='blue')
        ax[0].scatter(xreg,yreg,s=1,marker='o',c='b')
        ax[0].set_xlabel(xname)
        ax[0].set_ylabel(yname)
        #yreg is central cell plotted in Y axis
        #xreg is neighbor cell plotted in X axis

        ax[0].set_title(  'coeff='+str('%0.2f'%coeff) +  ', normal coeff='+str('%0.2f'%normalize) +', intercept='+str('%0.2f'%intercept)  ,fontsize=7 )
        #ax[0].set_title('ortho='+str(sum(xreg*yreg)))

        #CC_corr=CC_corr.reshape((len(CC_corr),1))
        #NC_corr=NC_corr.reshape((len(NC_corr),1))
        data1=np.hstack((NC_corr,CC_corr))
        ind=~np.isnan(data1).any(axis=1)
        data=data1[ind,:]
        gene=g_cc[ind]
        #print('pavi',pca_cc.shape,pca_nc.shape,data.shape, gene.shape)

        ax[1].scatter(data[:,0],data[:,1],s=1,marker='o',c='b')
        for i in range(len(gene)):
            ax[1].text(data[i,0],data[i,1],gene[i],fontsize=5)
        ax[1].set_xlabel(xname)


        fig.tight_layout()
        #plt.matshow(np.outer(np.sort(model.row_labels_) + 1, np.sort(model.column_labels_) + 1),cmap=plt.cm.Blues)
        #plt.title("Checkerboard structure of rearranged data")
        fig.savefig(savecorrelation+myname+'.png',dpi=300)
        plt.close('all')
