
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

#from descartes import PolygonPatch

from typing import List

import geopandas as gpd
from pyarrow import parquet
from shapely import wkb

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



df=pd.read_csv('./vis_day7_simplified/color_mapping.dat')
data=df.to_numpy()

inputdir='./../MERSCOOPE_cell_segmentations/'
df=pq.read_table(inputdir+'cell_boundaries.parquet (D7 batch0)').to_pandas()

riz=[10900,11100,7450,7650]
rbz=[11750,11950,6500,6700]
region=[10800,12300,6400,7900]
grayhex='#808080'
grayhex='#B4B4B4'


f=open('r_color_scheme/celltype_color_scheme_final.dat')
r_color_scheme={}
for line in f:
    l=line[0:-1].split(',')
    #print(l,l[2][1:-1],l[2])
    r_color_scheme[l[0]]=l[2][1:-1]

#print(df['EntityID'][1740:1750])
#print(df.loc[1740:1750,'EntityID'])
#print(df.columns)
cellid=df['EntityID'].to_numpy()

d={}
for i in range(len(cellid)):
    d[cellid[i]]=i
    #value=df[df['EntityID']== cellid[i]].index.values


index=[]
mycolor=[]
for i in range(len(data)):
    goodcellid=data[i,0]
    index.append(d[goodcellid])
    mycolor.append(data[i,8])

newdf=df.iloc[index]


multi_polygons  = newdf["Geometry"].apply(wkb.loads)
gdf=gpd.GeoSeries(multi_polygons)
CTname=sorted(list(np.unique(data[:,3])))

'''
andy_interest='CM_Dedifferentiating'
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=mycolor,ax=ax[0,1])
gdf.plot(color=mycolor,ax=ax[0,0])
gdf.plot(color=mycolor,ax=ax[1,0])

for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='o'
    msize1=4
    msize2=1
    rgba1=np.unique(remain[:,-1])
    rgba=rgba1[0]
    #if CTname[i]==andy_interest:
    #    mshape='*'
        #ax[0,1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=10)
    #ax[0,1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
#plt.savefig('ankit.png',dpi=300)



ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_'+andy_interest+'.pdf')
ax[1,1].get_legend().remove()

fig.savefig('vis_day7/'+andy_interest+'.png', bbox_inches='tight',transparent=transparent_mode,dpi=600)
fig.savefig('vis_day7/'+andy_interest+'.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
plt.close('all')



andy_interest=['B_cells', 'CD4/8 T' ,'CM_BZ', 'CM_Dedifferentiating' ,'CM_General' ,'DC_CCR7',
 'DC_con1', 'DC_con2_mo', 'DC_plasmacytoid' ,'EC_merged' ,'Erythrocyte',
 'FB_Saa3', 'Fibroblasts' ,'Fibrocytes' ,'ILC2', 'MAIT', 'Macro_MHCII_Cx3cr1',
 'Macrophages_Ly6C', 'Macrophages_Timd4' ,'Mast_cells', 'NK', 'NM',
 'Neutrophil', 'Pericytes' ,'SMC' ,'SwC', 'T_gd' ,'Treg', 'vEC_Lymphatic',
 'vEpicardial_derived', 'vPericyte_quiescent']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
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
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[0,1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)


ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc_alltogether.pdf')
ax[1,1].get_legend().remove()

transparent_mode=False
fig.savefig('vis_day7/coloc_alltogether.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc_alltogether.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)




andy_interest=['CM_BZ', 'CM_General', 'CM_Dedifferentiating', 'EC_merged', 'SMC', 'Neutrophil']
#andy_interest=['CM_BZ', 'CM_Dedifferentiating','SMC', 'Neutrophil']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
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
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[0,1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)


ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc1.pdf')
ax[1,1].get_legend().remove()

transparent_mode=False
fig.savefig('vis_day7/coloc1.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc1.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)




andy_interest=['CM_General', 'CM_Dedifferentiating', 'Treg', 'ILC2']
#andy_interest=[]
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
#axs.set_aspect('equal', 'datalim'
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])

for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc2_way.pdf')
ax[1,1].get_legend().remove()

transparent_mode=False
fig.savefig('vis_day7/coloc2_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc2_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)




andy_interest=['CM_General', 'CM_Dedifferentiating','EC_merged']
#andy_interest=['CM_Dedifferentiating']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']
gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])

for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)
ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc3_way.pdf')
ax[1,1].get_legend().remove()
transparent_mode=False
fig.savefig('vis_day7/coloc3_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc3_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)





andy_interest=['FB_Saa3',  'Macro_Trem2_Gpnmb', 'Macro_MHCII_Cx3cr1','Fibroblasts']#
#andy_interest=['CM_Dedifferentiating']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']

gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc4_way.pdf')
ax[1,1].get_legend().remove()
transparent_mode=False
fig.savefig('vis_day7/coloc4_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc4_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)






andy_interest=['vEC_Lymphatic', 'DC_CCR7', 'DC_con1', 'DC_con2_mo', 'DC_plasmacytoid', 'ILC2', 'CD4/8 T', 'T_gd', 'Treg', 'MAIT']#
#andy_interest=['CM_Dedifferentiating']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']

gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7/leg_coloc5_way.pdf')
ax[1,1].get_legend().remove()
transparent_mode=False
fig.savefig('vis_day7/coloc5_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7/coloc5_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)




andy_interest=['Neutrophil' ,'vEpicardial_derived']#
#andy_interest=['CM_Dedifferentiating']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']

gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7_simplified/leg_coloc6_way.pdf')
ax[1,1].get_legend().remove()
transparent_mode=False
fig.savefig('vis_day7_simplified/coloc6_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7_simplified/coloc6_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
'''



andy_interest=['MAIT', 'MAIT_Il17', 'ILC2', 'T_gd', 'Treg', 'Macro_MHCII_Cx3cr1']#
#andy_interest=['CM_Dedifferentiating']
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
figsize=(8,7)
fig,(ax)=plt.subplots(2,2,figsize=figsize)
markershape=['*','d','s','1','2','x','p']

gdf.plot(color=newcolor,ax=ax[0,1])
gdf.plot(color=newcolor,ax=ax[0,0])
gdf.plot(color=newcolor,ax=ax[1,0])
for i in range(len(CTname)):
    index=np.where(CTname[i]==data[:,3])
    remain=data[index[0],:]
    labelname=str(i)+'-'+CTname[i]+'-'+str(len(index[0]))
    mshape='.'
    msize1=4
    msize2=1
    rgba=grayhex
    ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,markersize=msize2)
    for j in range(len(andy_interest)):
        if CTname[i]== andy_interest[j]:
            #mshape=markershape[j]
            mshape='o'
            rgba=r_color_scheme[CTname[i]]
            #rgba=cmap(cumsum[j])
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,color='k',markersize=6)
            #ax[1].plot(data[index[0],6],data[index[0],7],mshape,label=labelname,color=rgba,markersize=msize1)
            ax[1,1].plot(data[index[0],4],data[index[0],5],mshape,color=rgba,label=labelname,markersize=msize2)

ax[0,0].set_title('IZ')
ax[0,1].set_title('BZ')
ax[1,0].plot([riz[0], riz[1], riz[1], riz[0], riz[0]],   [riz[2],riz[2],riz[3],riz[3],riz[2]], 'k--')
ax[1,0].plot([rbz[0], rbz[1], rbz[1], rbz[0], rbz[0]],   [rbz[2],rbz[2],rbz[3],rbz[3],rbz[2]], 'k--')
ax[0,0].set_xlim([riz[0],riz[1]])
ax[0,0].set_ylim([riz[2],riz[3]])
ax[0,1].set_xlim([rbz[0],rbz[1]])
ax[0,1].set_ylim([rbz[2],rbz[3]])
ax[1,0].set_xlim([region[0],region[1]])
ax[1,0].set_ylim([region[2],region[3]])
ax[1,1].set_xticks([])
ax[1,1].set_yticks([])
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])
#plt.gca().axes.get_yaxis().set_visible(False)
ax[1,1].set_axis_off()
fig.tight_layout()
transparent_mode=False
leg1=ax[1,1].legend(loc='upper right',bbox_to_anchor=(2, 1),ncol=1, borderaxespad=0.,prop={"size":8},fancybox=True, frameon=False,shadow=True)
export_legend(leg1,filename='vis_day7_simplified/leg_coloc7_way.pdf')
ax[1,1].get_legend().remove()
transparent_mode=False
fig.savefig('vis_day7_simplified/coloc7_way.png', bbox_inches='tight',transparent=transparent_mode,dpi=300)
fig.savefig('vis_day7_simplified/coloc7_way.pdf', bbox_inches='tight',transparent=transparent_mode,dpi=300)
