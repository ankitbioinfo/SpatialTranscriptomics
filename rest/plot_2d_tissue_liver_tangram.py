import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''
def reading_data(cont,barcodecell,cluster):
    #f=pd.read_excel('aau5324_Moffitt_Table-S7_coordinate.xlsx',engine='openpyxl',sheet_name=0,header=1)
    #data=f.to_numpy()
    #name='create_inputfile1.dat'
    #data = np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
    #print(data.shape)
    #print(data[0,7:])

    location_cellname2int={}
    location_int2cellname={}
    points=np.zeros((len(cont),2),dtype=float)
    index=[]
    for j in range(len(barcodecell)):
        index.append([])

    for i in range(len(cont)):
        name=cont[i][0]
        for j in range(len(barcodecell)):
            for k in range(len(barcodecell[j])):
                if name==barcodecell[j][k]:
                    index[j].append(i)
        location_cellname2int[name]=i
        location_int2cellname[i]=name
        points[i]=[cont[i][1], cont[i][2]]


    celltype={}
    noct=[]
    cellsinCT={}
    for i in range(len(cluster)):
        id=int(cluster[i][1])
        celltype[cluster[i][0]]=id
        if id not in noct:
            noct.append(id)

        if id not in cellsinCT:
            cellsinCT[id]=[ location_cellname2int[ cluster[i][0]]]
        else:
            cellsinCT[id].append(location_cellname2int[ cluster[i][0]])

    return points,index, cellsinCT


'''
def reading_data(cont,barcodecell,cluster,mycluster_interest_id):
    #f=pd.read_excel('aau5324_Moffitt_Table-S7_coordinate.xlsx',engine='openpyxl',sheet_name=0,header=1)
    #data=f.to_numpy()
    #name='create_inputfile1.dat'
    #data = np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
    #print(data.shape)
    #print(data[0,7:])

    location_cellname2int={}
    location_int2cellname={}
    points=np.zeros((len(cont),2),dtype=float)
    index=[]
    for j in range(len(barcodecell)):
        index.append([])


    for i in range(len(cont)):
        name=cont[i][0]
        #for j in range(len(barcodecell)):
            #for k in range(len(barcodecell[j])):
            #if name in barcodecell[j]:
            #    index[j].append(i)
        location_cellname2int[name]=i
        location_int2cellname[i]=name
        points[i]=[cont[i][1], cont[i][2]]

    for i in range(len(mycluster_interest_id)):
        ind=np.where(cluster[:,1]==mycluster_interest_id[i])
        index[i]=ind[0]


    celltype={}
    noct=[]
    cellsinCT={}
    for i in range(len(cluster)):
        id=int(cluster[i][1])
        celltype[cluster[i][0]]=id
        if id not in noct:
            noct.append(id)

        if id not in cellsinCT:
            cellsinCT[id]=[ location_cellname2int[ cluster[i][0]]]
        else:
            cellsinCT[id].append(location_cellname2int[ cluster[i][0]])

    return points,index, cellsinCT


def find_id(ctname,interest,cluster):
    barcode=[]
    for j in range(len(interest)):
        myindex='1'
        for i in range(len(ctname)):
            #print(ctname[i][1])
            if ctname[i][1]==interest[j]:
                myindex=ctname[i][0]
        if myindex=='1':
            print(ctname)

        index=np.where(cluster[:,1]==myindex)
        barcode.append(cluster[index[0],0])

    return barcode


def sort_index_in_right_order(correct,wrong):
    d={}
    for i in range(len(wrong)):
        d[wrong[i,0]]=i
    index=[]
    count=0
    for i in range(len(correct)):
        try:
            value=d[correct[i,0]]
        except KeyError:
            value=15 #monocytes
            count=count+1
        index.append(value)
    right=wrong[index]
    print('key not found',count)
    return right


#['B cells' 'Basophils' 'Central Vein EC' 'Cholangiocytes' 'Fibroblast'
# 'Hep1' 'Hep2' 'Hep3' 'Hep4' 'Hep5' 'Hep6' 'Hep7' 'Hep8' 'HsPCs' 'ILC1s'
# 'KCs' 'LSECs' 'Lymphatic EC' 'Macro' 'Mig. cDCs' 'Monocytes' 'NK cells'
# 'NM' 'Neutrophils' 'Portain Vein EC' 'Stellate cells' 'T cells' 'cDC1s'
# 'cDC2s' 'pDCs']




def main():

    res=0.5

    datapath='./Tangram/'

    df_cluster=pd.read_csv(datapath+'Tangram_derived_clusters.dat')
    degbased_cluster=df_cluster.to_numpy()
    print(degbased_cluster[0:5])
    df=pd.read_csv(datapath+'nameOfCT_poss3_24.dat',sep='\t',header=None)
    degbased_ctname=df.to_numpy()
    '''
    df_cluster=pd.read_csv('Tangram_detailed_annot_onSC/Tangram_derived_clusters.dat')
    mycluster=df_cluster.to_numpy()
    df=pd.read_csv('Tangram_detailed_annot_onSC/nameOfCT_ original.dat',sep='\t',header=None)
    mycluster_name=df.to_numpy()
    '''


    dirname='./macsct_results/'
    df=pd.read_csv(dirname+'clusterRes'+str(int(res*100))+'.csv')
    spatialcluster=df.to_numpy()
    df=pd.read_csv(dirname+'nameOfCTRes'+str(int(res*100))+'.dat',sep='\t',header=None)
    spatial_name=df.to_numpy()






    #['B cells' 'Basophils' 'Central Vein EC' 'Cholangiocytes' 'Fibroblast'
    #'Hep1245' 'Hep3678' 'HsPCs' 'ILC1s' 'KCs' 'LSECs' 'Lymphatic EC' 'Macro'
    #'Mig. cDCs' 'Monocytes' 'NK cells' 'NM' 'Neutrophils' 'Portain Vein EC'
    #'Stellate cells' 'T cells' 'cDC1s' 'cDC2s' 'pDCs'] [5, 6]
    #mycluster_interest=['Cholangiocytes','LSECs'] #,'cDC1s','cDC2s','pDCs'
    #mycluster_interest=['Hep4','Hep125_P', 'Hep3678_C']
    #mycluster_interest=['Hep1245', 'Hep3678']
    mycluster_interest=['T cells' ,'cDC2s']
    mycluster_interest=['Portain Vein EC','Central Vein EC','Lymphatic EC']
    #mycluster_interest=['Hep1', 'Hep2' ,'Hep3', 'Hep4' ,'Hep5', 'Hep6', 'Hep7', 'Hep8']
    #mycluster_interest=['Stellate cells','KCs']
    #mycluster_interest=['B cells', 'Basophils', 'Central Vein EC' ,'Cholangiocytes' ,'Fibroblast','HsPCs', 'ILC1s','LSECs', 'Lymphatic EC', 'Macro' ,'Mig. cDCs', 'NK cells','Neutrophils' ,'Portain Vein EC', 'Stellate cells', 'T cells' ,'cDC1s','cDC2s','pDCs']

    mycluster_interest_all=degbased_ctname[:,1]
    mycluster_interest_id=[]
    for i in range(len(mycluster_interest_all)):
        for j in range(len(degbased_ctname)):
            if degbased_ctname[j,1]==mycluster_interest_all[i]:
                mycluster_interest_id.append(degbased_ctname[j,0])

    print(mycluster_interest_id)

    #cl1=find_id(author_name,author_interest,authorcluster)

    df=pd.read_csv('./inputSP/tissue_positions_list.csv',header=None)
    posdata_not_order=df.to_numpy()
    posdata=sort_index_in_right_order(degbased_cluster,posdata_not_order)
    #posdata_anchors=sort_index_in_right_order(MNNcluster,posdata_not_order)

    df=pd.read_csv(dirname+'liver_umap.dat')
    umap_not_order=df.to_numpy()
    umap_data=sort_index_in_right_order(degbased_cluster,umap_not_order)
    #umap_data_anchors=sort_index_in_right_order(MNNcluster,umap_not_order)

    #spatial_interest=spatial_name[:,1]
    #cl4=find_id(spatial_name,spatial_interest,spatialcluster)
    #PP4,id4,cellsinCT4=reading_data(posdata,cl4,spatialcluster,spatial_interest)#,


    #CTname2=mycluster_name[:,1]
    CTname4=spatial_name[:,1]

    '''
    fig,(ax1)=plt.subplots(1,2,figsize=(15,6))
    plot_all_ct(CTname2,PP2,cellsinCT2,ax1[0])
    plot_all_ct(CTname2,umap_data[:,[1,2]],cellsinCT2,ax1[1])
    ax1[1].legend(loc='upper right',bbox_to_anchor=(1.30, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, shadow=True)
    fig.tight_layout()
    fig.savefig('all_cell_type_liver.png',bbox_inches='tight',dpi=300)
    '''


    for i in range(len(mycluster_interest_all)):

        mycluster_interest=[mycluster_interest_all[i]]

        cl2=find_id(degbased_ctname,mycluster_interest,degbased_cluster)
        #cl3=find_id(degbased_ctname,mycluster_interest,MNNcluster)

        PP2,id2,cellsinCT2=reading_data(posdata,cl2,degbased_cluster,[mycluster_interest_id[i]])
        #PP3,id3,cellsinCT3=reading_data(posdata_anchors,cl3,MNNcluster,[mycluster_interest_id[i]])


        fig,(ax1)=plt.subplots(1,2,figsize=(12,4))
        #plot_all_ct(CTname4,PP4,cellsinCT4,ax0[0])
        #plot_all_ct(CTname4,umap_data[:,[1,2]],cellsinCT4,ax0[1])

        plot_specific_ct(mycluster_interest,PP2,id2,ax1[0])
        plot_specific_ct(mycluster_interest,umap_data[:,[1,2]],id2,ax1[1])

        #plot_specific_ct(mycluster_interest,PP3,id3,ax2[0])
        #plot_specific_ct(mycluster_interest,umap_data_anchors[:,[1,2]],id3,ax2[1])

        ax1[0].set_xlim([3500,6000])
        ax1[0].set_ylim([2000,5000])
        #ax0[0].set_xlim([3500,6000])
        #ax0[0].set_ylim([2000,5000])
        #ax2[0].set_xlim([3500,6000])
        #ax2[0].set_ylim([2000,5000])

        xmin=-11
        xmax=17
        ymin=-11
        ymax=17
        ax1[1].set_xlim([xmin,xmax])
        ax1[1].set_ylim([ymin,ymax])
        #ax0[1].set_xlim([xmin,xmax])
        #ax0[1].set_ylim([ymin,ymax])
        #ax2[1].set_xlim([xmin,xmax])
        #ax2[1].set_ylim([ymin,ymax])

        #ax0[1].set_title('leiden clustering')
        ax1[1].set_title('degree based annotations')
        #ax2[1].set_title('Only anchored points')

        #ax0[1].legend(loc='upper right',bbox_to_anchor=(1.45, 1),ncol=2, borderaxespad=0.,prop={"size":8},fancybox=True, shadow=True)
        ax1[1].legend(loc='upper right',bbox_to_anchor=(1.25, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, shadow=True)
        #ax2[1].legend(loc='upper right',bbox_to_anchor=(1.25, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, shadow=True)

        fig.tight_layout()
        name=mycluster_interest[0].replace('/','_') + str(i)
        fig.savefig('fig/specific_cell_type_liver_'+name+'.png',bbox_inches='tight',dpi=300)
        plt.close('all')




    return 0


def plot_all_ct(CTname,PP,cellsinCT,ax):
    cmap=plt.cm.get_cmap('Spectral')
    cmap=plt.cm.get_cmap('jet')

    cumsum=np.linspace(0,1,len(CTname))

    for i in range(len(CTname)):
        index=cellsinCT[i]
        labelname=str(i)+'-'+CTname[i]
        rgba=cmap(cumsum[i])
        ax.plot(PP[index,0],PP[index,1],'o',label=labelname,color=rgba,markersize=1)
        x=np.mean(PP[index,0])
        y=np.mean(PP[index,1])
        ax.text(x,y,str(i),fontsize=12)




def plot_specific_ct(CTname,PP,index,ax):

    cmap=plt.cm.get_cmap('jet')
    #cmap=plt.cm.get_cmap('Spectral')
    cumsum=np.linspace(0,1,len(CTname))
    remaining_index=[]
    for i in range(len(PP)):
        for j in range(len(index)):
            if i not in index[j]:
                remaining_index.append(i)

    mycolor=['r','g','b','y','y','y','y','y']
    for j in range(len(index)):
        labelname=str(j)+'-'+CTname[j]+'-'+str(len(index[j]))
        #labelname=str(j)+'-'+CTname[j]
        #print(labelname)
        ms=0.5
        #if j=='0-Hep4':
        #    ms=3
        #else:
        #    ms=0.5
        rgba=cmap(cumsum[j])
        ax.plot(PP[index[j],0],PP[index[j],1],'o',label=labelname,color=rgba,markersize=ms)
        x=np.mean(PP[index[j],0])
        y=np.mean(PP[index[j],1])
        #ax.text(x,y,str(j),fontsize=12)
    ax.plot(PP[remaining_index,0],PP[remaining_index,1],'.',color="0.5",label='NA',markersize=0.05)



main()
