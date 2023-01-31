import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.decomposition import NMF
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as snn



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

def readdata():
    scdatapath='./inputSC/'


    name=scdatapath+'/scdata2_most_recent/cluster_poss4_25.csv'
    df=pd.read_csv(name)
    cluster=df.to_numpy()
    #input['sc_cluster']=data


    celltypefname=scdatapath+'/scdata2_most_recent/nameOfCT_poss4_25.dat'
    df=pd.read_csv(celltypefname,sep='\t',header=None)
    cluname=df.to_numpy()
    #print(cluname)

    for i in range(len(cluname)):
        if cluname[i,1]=='KCs':
        #if cluname[i,1]=='Portain Vein EC':
            kcid=cluname[i,0]
            index=np.where(cluster[:,1]==kcid)
            kc_cells=cluster[index[0],0]

    print('KC',len(kc_cells))

    #full_ad_sc=sc.read_h5ad(scdatapath+'sc_liver_data.h5ad')
    #full_ad_sc=sc.read_h5ad('./countdata/HVG_on_countData_3000.h5ad')
    full_ad_sc=sc.read_h5ad(scdatapath+'HVG_counts_3000.h5ad')


    genename=full_ad_sc.var_names.to_numpy()
    cellname=full_ad_sc.obs_names.to_numpy()

    print(len(cellname),len(genename))

    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i

    cellindex=[]
    for i in range(len(kc_cells)):
        cellindex.append(d[kc_cells[i]])

    adata=full_ad_sc[cellindex,:]
    mscfull=adata.X.toarray()

    #ad_sc=sc.read_h5ad('./countdata/count_singleCell.h5ad')
    ad_sc=sc.read_h5ad(scdatapath+'common_counts_sc.h5ad')


    mygene=ad_sc.var_names.to_numpy()
    cellname=ad_sc.obs_names.to_numpy()

    d={}
    for i in range(len(cellname)):
        d[cellname[i]]=i

    cellindex=[]
    for i in range(len(kc_cells)):
        cellindex.append(d[kc_cells[i]])

    adata=ad_sc[cellindex,:]
    msc=adata.X.toarray()


    #geneIndex=[]
    #for i in range(len(genename)):
    #    if genename[i] in gn:
    #        geneIndex.append(i)

    #mygene=genename[geneIndex]
    #adata=full_ad_sc[cellindex,geneIndex]
    #msc=adata.X.toarray()


    print(msc.shape,mscfull.shape)

    return msc,mygene,mscfull,genename



def find_factors(no_of_pc,msc,genename_joint,CbyG,genename_full):
    model = NMF(n_components=no_of_pc, init = "nndsvda", random_state=1,beta_loss="kullback-leibler",solver="mu",max_iter=1000,alpha_W=0.0,alpha_H=0.0,l1_ratio=0)

    W = model.fit_transform(msc.T)
    H = model.components_
    value=np.sqrt(np.sum((msc.T-np.matmul(W,H))**2))

    w_high=np.max(W)
    w_low=np.min(W)
    h_high=np.max(H)
    h_low=np.min(H)
    print('ast',w_high,w_low,h_high,h_low)



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

    #corr_NMFfactors_genes=find_correlation_bw_genes_and_PC_component_in_singlecell(spH.T,msp)
    #gname,geneNMF=top_genes_in_correlation_list(genename_spatial,corr_NMFfactors_genes,5)
    #b=snn.heatmap(geneNMF,yticklabels=gname,ax=ax1)#componentlabel,ax=ax
    #b.set_xticklabels(xlabels,size = 8,rotation=90)
    #b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)

    corr_NMFfactors_genes=find_correlation_bw_genes_and_PC_component_in_singlecell(H.T,CbyG)
    gname,geneNMF=top_genes_in_correlation_list(genename_full,corr_NMFfactors_genes,10)
    b=snn.heatmap(geneNMF,yticklabels=gname,ax=ax2)#componentlabel,ax=ax
    b.set_xticklabels(xlabels,size = 8,rotation=90)
    b.set_yticklabels(b.get_ymajorticklabels(), fontsize = 6)


    fig.savefig('NMF_out.png',dpi=300)
    plt.close('all')



def main():
    msc,genename_joint,CbyG,genename_full=readdata()
    find_factors(5,msc,genename_joint,CbyG,genename_full)



main()
