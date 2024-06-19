
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

#from descartes import PolygonPatch

from typing import List

import geopandas as gpd
from pyarrow import parquet
from shapely import wkb
from shapely.geometry import Polygon
import pickle
import os

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 0.1 #set the value globally

import matplotlib.pyplot as plt
plt.rc('font', family='Helvetica')

def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox,transparent=True)


def make_axis_good(ax):
    ax[0,0].set_title('region i')
    ax[0,1].set_title('region k')
    ax[1,0].set_title('region h')
    ax[1,1].set_title('region j')


    riz=[3700,4100,1200,1600]
    rbz=[2800,3200,1700,2100]
    regh=[2200,2600,1100,1500]
    regj=[700,1100,1900,2300]
    region=[200,5000,500,4800]

    ax[2,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
    ax[2,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
    ax[2,0].plot([regh[0], regh[1], regh[1], regh[0], regh[0]],   [regh[2],regh[2],regh[3],regh[3],regh[2]], 'k--')
    ax[2,0].plot([regj[0], regj[1], regj[1], regj[0], regj[0]],   [regj[2],regj[2],regj[3],regj[3],regj[2]], 'k--')

    ax[0,0].set_xlim([riz[0],riz[1]])
    ax[0,0].set_ylim([riz[2],riz[3]])

    ax[0,1].set_xlim([rbz[0],rbz[1]])
    ax[0,1].set_ylim([rbz[2],rbz[3]])

    ax[1,0].set_xlim([regh[0],regh[1]])
    ax[1,0].set_ylim([regh[2],regh[3]])

    ax[1,1].set_xlim([regj[0],regj[1]])
    ax[1,1].set_ylim([regj[2],regj[3]])


    ax[2,1].set_xticks([])
    ax[2,1].set_yticks([])
    ax[0,0].set_xticks([])
    ax[0,0].set_yticks([])
    ax[0,1].set_xticks([])
    ax[0,1].set_yticks([])

    ax[1,0].set_xticks([])
    ax[1,0].set_yticks([])
    ax[1,1].set_xticks([])
    ax[1,1].set_yticks([])

    ax[2,0].set_xlim([region[0],region[1]])
    ax[2,0].set_ylim([region[2],region[3]])

     #plt.gca().axes.get_yaxis().set_visible(False)
    ax[2,1].set_axis_off()

    return 0


grayhex='#808080'
grayhex='#B4B4B4'


f=open('../r_color_scheme_detailed/celltype_color_scheme_final.dat')
r_color_scheme={}
for line in f:
    l=line[0:-1].split(',')
    #print(l,l[2][1:-1],l[2])
    r_color_scheme[l[0]]=l[2][1:-1]


'''
And these are the colocalized cell types that we needed, using the detailed annotation scheme:
CM_Angiogenic, CM_Btk+, CM_Dedifferentiating, CM_Homeostatic, CM_Hypertrophic, CM_Prehypertrophic, CM_Slit2+

B_cells, DC_CCR7, DC_con1, DC_con2_Ifitm1, DC_con2b/moDC, DC_plasmacytoid, ILC2, ILC2_IL5, Macro_Cxcl13, Macro_MHCII_Cx3cr1, Macro_Retnla, Macro_Timd4_Lyve1, Macro_Trem2_Gpnmb, MAIT, MAIT_IL17, Mast_cells, mo_Macro_Il10_Fn1, mo_Macro_Ly6C_Isg15_Cxcl3, mo_Macro_Spp1, Neutrophil_1, Neutrophil_2, NK_Gzma, NK_Klra5, NK_T, T_CD4_effector, T_CD4_naive, T_CD8_effector, T_CD8_naive, T_gd, T_IFNg_naive, T_Isg15, T_Macro, Treg, Fibrocytes

Macro_Cxcl13, Macro_MHCII_Cx3cr1, Macro_Retnla, Macro_Timd4_Lyve1, Macro_Trem2_Gpnmb, MAIT, MAIT_IL17, Mast_cells, mo_Macro_Il10_Fn1, mo_Macro_Ly6C_Isg15_Cxcl3, mo_Macro_Spp1

ILC2, ILC2_IL5, NK_Gzma, NK_Klra5, NK_T, T_CD4_effector, T_CD4_naive, T_CD8_effector, T_CD8_naive, T_gd, T_IFNg_naive, T_Isg15, Treg

FB_Ccl2, FB_Cxcl14, FB_Duox1, FB_Fgl2, FB_Postn_Thbs4, FB_quiescent, FB_Saa3, FB_transition_Cd9, FB_transition_Il1rapl1, FB_WntX, FBmyo, Fbmyo_cycling, Fbmyo_Dkk2, Fbmyo_Hp, FBmyo_IFN, Fibrocytes

Schwann_Galectin, Schwann_IFN, Schwann_metabolic, Schwann_quiescent

vEC_angio_IFN, vEC_Areg_Dkk2_Wnt, vEC_Arterial, vEC_capillary1, vEC_capillary2, vEC_Endocardial, vEC_Immune, vEC_Lymphatic, vEC_metabolic, vEpicardial_derived, vPericyte_FB, vPericyte_INF, vPericyte_quiescent, vSMC_Ryr2, vSMC1, vSMC2
'''

outdir='./vis_day28/'

df=pd.read_csv(outdir+'color_mapping.dat')
data=df.to_numpy()



def read_segmentation_file():
    df=pq.read_table('cell_boundaries.parquet').to_pandas()
    data = df.to_numpy()
    ncell = np.unique(data[:,0])
    #ncell = ncell[0:100]
    fout='save_segmentation.p'
    #print("ncell",data.shape,len(ncell))
    flag=1
    if os.path.isfile(fout):
        filesize = os.path.getsize(fout)
        if filesize>0:
            flag=0
    if flag==1:
        lgdf=[]
        for i in range(len(ncell)):
            index=np.where(ncell[i]==data[:,0])
            vertex = data[index[0]][:,1:]
            poly = Polygon(zip(vertex[:,0],vertex[:,1]))
            lgdf.append(poly)
            #print(i,vertex)
            #print(i, poly.exterior.coords[])
        gdf = gpd.GeoSeries(lgdf)
        myfile = open(fout, 'wb')
        pickle.dump(gdf,myfile)
        pickle.dump(ncell,myfile)
        myfile.close()
    if flag==0:
        myfile = open(fout, 'rb')
        gdf=pickle.load(myfile)
        ncell=pickle.load(myfile)
        myfile.close()


    d={}
    for i in range(len(ncell)):
        d[ncell[i]]=i
    return d, gdf



d,gdf=read_segmentation_file()
#print(df.columns)
#print(df['EntityID'][1740:1750])
#print(df.loc[1740:1750,'EntityID'])
#print(df.columns)
    #value=df[df['EntityID']== cellid[i]].index.values


index=[]
mycolor=[]
for i in range(len(data)):
    goodcellid=data[i,1]
    index.append(d[goodcellid])
    mycolor.append(data[i,8])

#newdf=df.iloc[index]
gdf = gdf[index]
transparent_mode=False

#multi_polygons  = newdf["Geometry"].apply(wkb.loads)
#gdf=gpd.GeoSeries(multi_polygons)
CTname=sorted(list(np.unique(data[:,3])))

'''
andy_interest='CM_Dedifferentiating'
figsize=(8,3.5)
figsize=(12,5)
fig,(ax)=plt.subplots(1,2,figsize=figsize)
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=mycolor,ax=ax[1])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='o'
    msize1=4
    msize2=1
    rgba1=np.unique(remain[:,-1])
    rgba=rgba1[0]
    if CTname[i]==andy_interest:
        mshape='*'
        ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=10)
    ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
    ax[0].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
#plt.savefig('ankit.png',dpi=300)
leg1=ax[1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=2, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
region=[500,1000,4100,4400]
region=[1050,1150,4400,4500]
ax[1].set_xlim([region[0],region[1]])
ax[1].set_ylim([region[2],region[3]])
ax[0].set_xticks([])
ax[0].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[0].set_axis_off()
#fig.tight_layout()
fig.savefig('day_vis_'+andy_interest+'.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('day_vis_'+andy_interest+'.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
'''

andy_interest=['CM_Angiogenic', 'CM_Btk+', 'CM_Dedifferentiating', 'CM_Homeostatic1',
'CM_Homeostatic2', 'CM_Hypertrophic', 'CM_Prehypertrophic', 'CM_Slit2+',]
#andy_interest=['CM_Btk+','CM_Prehypertrophic','CM_Dedifferentiating', 'CM_Slit2+']
#cmap=plt.colormaps['jet']
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])

for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    #rgba1=np.unique(remain[:,-1])
    #rgba=rgba1[0]
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            #rgba=cmap(cumsum[j])
            rgba=r_color_scheme[CTname[i]]
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

make_axis_good(ax)
#fig.tight_layout()
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc1CM.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc1CM.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig(outdir+'coloc1CM.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')



andy_interest=['B_cells', 'DC_CCR7', 'DC_con1', 'DC_con2_Ifitm1', 'DC_con2b/moDC', 'DC_plasmacytoid', 'ILC2', 'ILC2_IL5',
 'Macro_Cxcl13', 'Macro_MHCII_Cx3cr1', 'Macro_Retnla', 'Macro_Timd4_Lyve1', 'Macro_Trem2_Gpnmb', 'MAIT', 'MAIT_IL17',
 'Mast_cells', 'mo_Macro_Il10_Fn1', 'mo_Macro_Ly6C_Isg15_Cxcl3', 'mo_Macro_Spp1', 'Neutrophil_1',
 'Neutrophil_2', 'NK_Gzma', 'NK_Klra5', 'NK_T', 'T_CD4_effector', 'T_CD4_naive', 'T_CD8_effector', 'T_CD8_naive',
 'T_gd', 'T_IFNg_naive', 'T_Isg15', 'T_Macro', 'Treg', 'Fibrocytes']
andy_interest=sorted(andy_interest)
#andy_interest=['Treg']
#cmap=plt.colormaps['jet']
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=1
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc2Immune.pdf')
ax[2,1].get_legend().remove()
#fig.tight_layout()
fig.savefig(outdir+'coloc2Immune.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc2Immune.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')


andy_interest=[ 'Macro_Retnla',  'Macro_Trem2_Gpnmb','Macro_Cxcl13', 'mo_Macro_Ly6C_Isg15_Cxcl3', 'mo_Macro_Spp1',
'MAIT', 'MAIT_IL17', 'Mast_cells', 'mo_Macro_Il10_Fn1', 'Macro_MHCII_Cx3cr1','Macro_Timd4_Lyve1']
#andy_interest=['MAIT','MAIT_IL17','Macro_Cxcl13','Macro_Retnla','Mast_cells',]
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc3Immune.pdf')
ax[2,1].get_legend().remove()
#fig.tight_layout()
fig.savefig(outdir+'coloc3Immune.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc3Immune.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')



andy_interest=[ 'ILC2_IL5', 'NK_Gzma', 'NK_Klra5', 'NK_T',  'T_CD4_naive',
'T_CD8_naive', 'T_gd',  'T_Isg15', 'ILC2','T_CD4_effector','T_IFNg_naive','Treg','T_CD8_effector']
#andy_interest=[]
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=2
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc4Lympho.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc4Lympho.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc4Lympho.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')




andy_interest=[ 'Schwann_Galectin', 'Schwann_IFN', 'Schwann_metabolic', 'Schwann_quiescent']
#andy_interest=[]
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=2
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc6Neural.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc6Neural.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc6Neural.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')





andy_interest=[ 'vEC_angio_IFN', 'vEC_Areg_Dkk2_Wnt', 'vEC_Arterial', 'vEC_capillary1', 'vEC_capillary2', 'vEC_Endocardial',
'vEC_Lymphatic', 'vSMC_Ryr2','vEC_metabolic','vPericyte_INF', 'vPericyte_FB','vEC_Immune',
 'vEpicardial_derived',  'vPericyte_quiescent',  'vSMC1', 'vSMC2']
#andy_interest=[]
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=2
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc7Vascular.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc7Vascular.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc7Vascular.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')



andy_interest=[ 'FB_Ccl2', 'FB_Cxcl14', 'FB_Duox1', 'FB_Fgl2', 'FB_Postn_Thbs4', 'FB_quiescent', 'FB_Saa3', 'FB_transition_Cd9',
'FB_transition_Il1rapl1', 'FB_WntX', 'FBmyo', 'Fbmyo_cycling', 'Fbmyo_Dkk2', 'Fbmyo_Hp', 'FBmyo_IFN', 'Fibrocytes']
#andy_interest=['FB_Ccl2','FB_Postn_Thbs4', 'FBmyo', 'Fbmyo_cycling','FB_WntX',]
#cumsum=np.linspace(0,1,len(andy_interest))
newcolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    flag=1
    for j in range(len(andy_interest)):
        if data[i,3]==andy_interest[j]:
            flag=0
            #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
            mycolor=data[i,8]
    if flag==1:
        newcolor.append(grayhex)
    else:
        newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=2
    msize2=1
    rgba=grayhex
    ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc5FB.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc5FB.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc5FB.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')


andy_interest=['B_cells', 'CM_Angiogenic', 'CM_Btk+', 'CM_Dedifferentiating', 'CM_Hypertrophic' ,'CM_Homeostatic1', 'CM_Homeostatic1',
 'CM_Prehypertrophic', 'CM_Slit2+',
 'DC_CCR7' ,'DC_con1', 'DC_con2_Ifitm1' ,'DC_con2b/moDC', 'DC_plasmacytoid','FB_WntX',
 'Erythrocyte', 'FB_Ccl2' ,'FB_Cxcl14', 'FB_Duox1', 'FB_Fgl2', 'FB_Postn_Thbs4',
 'FB_Saa3', 'FB_WntX' ,'FB_quiescent', 'FB_transition_Cd9',
 'FB_transition_Il1rapl1', 'FBmyo', 'FBmyo_IFN', 'Fbmyo_Dkk2' ,'Fbmyo_Hp',
 'Fbmyo_cycling', 'Fibrocytes' ,'ILC2' ,'ILC2_IL5', 'MAIT', 'MAIT_IL17',
 'Macro_Cxcl13', 'Macro_MHCII_Cx3cr1' ,'Macro_Retnla', 'Macro_Timd4_Lyve1',
 'Macro_Trem2_Gpnmb', 'Mast_cells' ,'NK_Gzma' 'NK_Klra5' ,'NK_T', 'NM',
 'Neutrophil_1', 'Neutrophil_2', 'Schwann_Galectin', 'Schwann_metabolic',
 'Schwann_quiescent', 'T_CD4_effector', 'T_CD4_naive', 'T_CD8_effector',
 'T_CD8_naive', 'T_IFNg_naive', 'T_Isg15', 'T_gd', 'Treg', 'mo_Macro_Il10_Fn1',
 'mo_Macro_Ly6C_Isg15_Cxcl3', 'mo_Macro_Spp1', 'vEC_Areg_Dkk2_Wnt',
 'vEC_Arterial', 'vEC_Endocardial', 'vEC_Immune', 'vEC_Lymphatic',
 'vEC_angio_IFN', 'vEC_capillary1', 'vEC_capillary2', 'vEC_metabolic',
 'vEpicardial_derived', 'vPericyte_FB', 'vPericyte_INF',
 'vPericyte_quiescent', 'vSMC1', 'vSMC2' ,'vSMC_Ryr2']

newcolor=[]
for i in range(len(data)):
     goodcellid=data[i,0]
     flag=1
     for j in range(len(andy_interest)):
         if data[i,3]==andy_interest[j]:
             flag=0
             #mycolor=matplotlib.colors.rgb2hex(cmap(cumsum[j]))
             mycolor=data[i,8]
     if flag==1:
         newcolor.append(grayhex)
     else:
         newcolor.append(mycolor)
figsize=(7,9)
fig,(ax)=plt.subplots(3,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'

gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,1])
gdf.plot(color=newcolor,ax=ax[1,0])
gdf.plot(color=newcolor,ax=ax[2,0])

for i in range(len(CTname)):
     index=np.where(CTname[i]==data[:,3])
     remain=data[index[0],:]
     labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
     mshape='.'
     msize1=2
     msize2=1
     rgba=grayhex
     ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
     for j in range(len(andy_interest)):
         if CTname[i]== andy_interest[j]:
             #mshape=markershape[j]
             mshape='o'
             rgba=r_color_scheme[CTname[i]]
             #rgba=cmap(cumsum[j])
             #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
             #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
             ax[2,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

make_axis_good(ax)
leg1=ax[2,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename=outdir+'leg_coloc_alltogether.pdf')
ax[2,1].get_legend().remove()
fig.savefig(outdir+'coloc_alltogether.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig(outdir+'coloc_alltogether.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
#'''
