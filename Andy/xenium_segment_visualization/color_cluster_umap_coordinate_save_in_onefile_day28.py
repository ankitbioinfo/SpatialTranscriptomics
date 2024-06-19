



import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib




def read_everything(ref_datapath,quepath,batchid=0):
        '''
        number_of_iteration_to_perform_celltype_annotations=3
        positionFilename=quepath+'./tissue_positions_list.csv'

        #fig_save_path=quepath+'MNN_based_annotations/'

        adata_all=sc.read_h5ad(quepath+'sct_spatial.h5ad')
        adata = adata_all[adata_all.obs.batch == batchid] #str
        temp=adata.obsm['X_umap']
        cellname0=adata.obs_names.to_numpy()#df.index.to_numpy()
        cellname=np.reshape(cellname0,(len(cellname0),1))
        umap_data=np.hstack((cellname,temp))
        degbased_cluster=sort_index_in_right_order(umap_data,degbased_cluster_not_order)#
        posdata=np.vstack((cellname0,adata.obs['coord_x'],adata.obs['coord_y'])).T



        #new_qpath='./../Andy_sham3/inputQuery/day7_sec1/'
        new_qpath='../day7_WA_MV_only_sec1/inputQuery/day7_sec1/'
        fid=pd.read_csv(new_qpath+'cellname.txt',sep='\t',header=None)
        viz2nico=fid.to_numpy()

        d={}
        for j in range(len(viz2nico)):
            temp=viz2nico[j,1]#+'-'+str(batchid)
            #print(temp,cellname[0:5])
            d[temp]=viz2nico[j,0]

        correct_order=[]



        nico2viz=np.array(correct_order)
        '''
        number_of_iteration_to_perform_celltype_annotations=3
        deg_annot_cluster_fname=str(number_of_iteration_to_perform_celltype_annotations)+'_nico_annotation_cluster.csv'
        deg_annot_ct_fname=str(number_of_iteration_to_perform_celltype_annotations)+'_nico_annotation_ct_name.csv'

        df_cluster=pd.read_csv(deg_annot_cluster_fname)
        cluster=df_cluster.to_numpy()
        df=pd.read_csv(deg_annot_ct_fname)
        ctname=df.to_numpy()

        nicoct =[]
        for i in range(len(cluster)):
            for j in range(len(ctname)):
                if ctname[j][0]==cluster[i][1]:
                    nicoct.append(ctname[j][1])

        cluster = np.array(nicoct)

        adata = sc.read_h5ad('spatial_tissue_with_nico_annotations.h5ad')
        print(adata)
        #cluster= adata.obs['nico_ct2']
        CT  = sorted(list(np.unique(cluster)))
        degbased_ctname=[]
        d={}
        for i in range(len(CT)):
            degbased_ctname.append([i,CT[i]])
            d[CT[i]]=i


        degbased_cluster=[]
        posdata=[]
        umap_data=[]
        nico2viz=[]

        cellname=list(adata.obs_names)
        for i in range(len(cluster)):
            degbased_cluster.append([cellname[i],d[cluster[i]]])
            u=adata.obsm['X_umap'][i]
            s=adata.obsm['spatial'][i]
            posdata.append([cellname[i],s[0],s[1]])
            umap_data.append([cellname[i],u[0],u[1]])
            nico2viz.append(cluster[i])
        degbased_cluster=np.array(degbased_cluster)
        degbased_ctname=np.array(degbased_ctname)
        posdata=np.array(posdata)
        umap_data=np.array(umap_data)
        nico2viz=np.array(nico2viz)



        #sometime if you have less number of spatial cells (due to filtering step) in the analysis than the position coordinate have
        #then need to find correct pairing.
        #print(len(cellname),len(viz2nico),len(umap_not_order))

        '''
        points=np.zeros((len(posdata),2),dtype=float)
        location_cellname2int={}
        location_int2cellname={}
        for i in range(len(posdata)):
            name=posdata[i][0]
            location_cellname2int[name]=i
            location_int2cellname[i]=name
            points[i]=[posdata[i][1], posdata[i][2]]

        cellsinCT={}
        for i in range(len(degbased_cluster)):
            id=int(degbased_cluster[i][1])
            #celltype[degbased_cluster[i][0]]=id
            if id not in cellsinCT:
                cellsinCT[id]=[ location_cellname2int[ degbased_cluster[i][0]]]
            else:
                cellsinCT[id].append(location_cellname2int[ degbased_cluster[i][0]])
        '''


        return degbased_ctname, degbased_cluster,posdata,umap_data,nico2viz

def mapping(outdir,degbased_ctname, degbased_cluster,posdata,umap_data,nico2viz,r_color_scheme):
    CTname=degbased_ctname[:,1]

    cumsum=np.linspace(0,1,len(CTname))
    #cmap=plt.cm.get_cmap('jet')
    cmap=plt.colormaps['jet']

    rgba={}
    d={}
    for i in range(len(CTname)):
        #index=cellsinCT[i]
        #labelname=str(i)+'-'+CTname[i]+'-'+str(len(index))
        value=cmap(cumsum[i])
        #rgba[degbased_ctname[i,0]]=matplotlib.colors.rgb2hex(value)
        rgba[degbased_ctname[i,0]]=r_color_scheme[degbased_ctname[i,1]]
        d[degbased_ctname[i,0]]=degbased_ctname[i,1]

    fw=open(outdir+'color_mapping.dat','w')
    fw.write('vizid,nicoid,cluid,ctname,umap_x,umap_y,pos_x,pos_y,color\n')
    for i in range(len(degbased_cluster)):
        fw.write(str(nico2viz[i])+','+degbased_cluster[i,0]+','+str(degbased_cluster[i,1])+','+d[degbased_cluster[i,1]]+','+
        str(umap_data[i,1])+','+str(umap_data[i,2])+','+str(posdata[i,1])+','+str(posdata[i,2])+','+rgba[degbased_cluster[i,1]]+'\n')



def myrun():



    #d7b0='./../MERSCOOPE_cell_segmentations/cell_boundaries.parquet (D7 batch0)'
    #d7b1='./../MERSCOOPE_cell_segmentations/cell_boundaries.parquet (D7 batch1)'
    #d28='./../MERSCOOPE_cell_segmentations/cell_boundaries.parquet (D28)'

    ref_datapath='./inputRef/sham/'
    query_datapath='./inputQuery/sham/'

    ref_datapath='./inputRef/day1/'
    query_datapath='./inputQuery/day1/'

    ref_datapath='./../Andy_sham3/inputRef/day7/'
    query_datapath='./../Andy_sham3/inputQuery/day7/'

    query_datapath='../day7_WA_MV_only_sec1/inputQuery/day7/'
    ref_datapath='../day7_WA_MV_only_sec1/inputRef/day7/'

    f=open('../r_color_scheme_detailed/celltype_color_scheme_final.dat')
    r_color_scheme={}
    for line in f:
        l=line[0:-1].split(',')
        #print(l,l[2][1:-1],l[2])
        r_color_scheme[l[0]]=l[2][1:-1]

    #print(r_color_scheme)

    CTname, degbased_cluster,posdata,umap_data,nico2viz=read_everything(ref_datapath,query_datapath)

    print(CTname[:,1])

    #outdir='./../MERSCOOPE_cell_segmentations/'
    outdir='./vis_day28/'

    mapping(outdir,CTname, degbased_cluster,posdata,umap_data,nico2viz,r_color_scheme)


def sort_index_in_right_order(correct,wrong):
    "Helper function used to visualize cell type annotations."
    d={}
    for i in range(len(wrong)):
        #print(wrong[i,0])
        d[wrong[i,0]]=i
    index=[]
    for i in range(len(correct)):
        #print(correct[i][0])
        index.append(d[correct[i,0]])
    right=wrong[index]
    return right


myrun()
