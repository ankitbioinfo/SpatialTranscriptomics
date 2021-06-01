import numpy as np
import matplotlib.pyplot as plt
from colormap import Color
import matplotlib as mpl
import matplotlib.cm as cm
from collections import OrderedDict
from sklearn.metrics.cluster import adjusted_rand_score


cmaps=OrderedDict()

def euclidean_dist(p1,p2):
	return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def calculate_correlation(cell_cell_distance,umap):
    dist=np.zeros(len(cell_cell_distance),dtype=np.float)
    for i in range(len(cell_cell_distance)):
        x=int(cell_cell_distance[i][0])
        y=int(cell_cell_distance[i][1])
        #print(x,y)
        dist[i]=euclidean_dist(umap[x], umap[y])



def readfiles(name):
    f1='neighborhoodData/figures_'+name+'_original/umap_coordinate.dat'
    f2='neighborhoodData/figures_'+name+'_avg/umap_coordinate.dat'
    f3='neighborhoodData/figures_'+name+'_sum/umap_coordinate.dat'
    f4='neighborhoodData/figures_'+name+'_normalized/umap_coordinate.dat'


    umap_ori=np.loadtxt(f1,delimiter=',',skiprows=0,usecols=None)
    umap_avg=np.loadtxt(f2,delimiter=',',skiprows=0,usecols=None)
    umap_sum=np.loadtxt(f3,delimiter=',',skiprows=0,usecols=None)
    #umap_normalized=np.loadtxt(f4,delimiter=',',skiprows=0,usecols=None)


    points=np.loadtxt('modified_tissue_positions.dat',delimiter='\t')

    f5='neighobor2_cell-cell_distance.dat'
    cell_cell_distance=np.loadtxt(f5,delimiter='\t',skiprows=0,usecols=None)
    dist_physical=cell_cell_distance[:,2]
    '''
    dist_orig_expression=calculate_correlation(cell_cell_distance,umap_ori)
    dist_avg_expression=calculate_correlation(cell_cell_distance,umap_avg)
    dist_sum_expression=calculate_correlation(cell_cell_distance,umap_sum)

    rho1=np.corrcoef(dist_physical,dist_orig_expression)
    rho2=np.corrcoef(dist_physical,dist_avg_expression)
    rho3=np.corrcoef(dist_physical,dist_sum_expression)

    #plt.hist(dist_physical,bins=30)
    print(len(cell_cell_distance,len(dist_orig_expression)))
    plt.hist(dist_orig_expression,bins=30)

    plt.show()

    print(rho1,rho2,rho3)
    '''

    f5=open('neighborhoodData/figures_'+name+'_original/umap_colors.dat')
    cont=f5.readlines()
    mycolors=[]
    legend=range(1,len(cont)+1)
    for i in range(len(cont)):
        mycolors.append(Color(cont[i][0:-1]))
    combinedcelltype=[]
    for fi in ['original','avg','sum']:
	    #f4=open('neighborhoodData/figures_'+name+'_original/leiden_output.dat')
	    #f4=open('neighborhoodData/figures_'+name+'_avg/leiden_output.dat')
	    f4=open('neighborhoodData/figures_'+name+'_'+fi+'/leiden_output.dat')
	    cont=f4.readlines()
	    celltype=[]
	    for i in range(1,len(cont)):
	        l=cont[i].split(',')
	        id=int(l[1])
	        #celltype.append(colors[id])
	        celltype.append(id)
	    combinedcelltype.append(np.array(celltype))

    print(umap_ori.shape,umap_avg.shape,umap_sum.shape,len(mycolors))

    print('expression and avg', adjusted_rand_score(combinedcelltype[0],combinedcelltype[1]))
    print('expression and sum', adjusted_rand_score(combinedcelltype[0],combinedcelltype[2]))
    print('avg and sum', adjusted_rand_score(combinedcelltype[1],combinedcelltype[2]))


    umap_normalized=0
    celltype=combinedcelltype[0]

    return [umap_ori,umap_avg,umap_sum,umap_normalized,celltype,mycolors,points]


factors='ecm'


[umap_ori,umap_avg,umap_sum,umap_normalized,celltype,mycolors,points]=readfiles(factors)
cmaps['mycolors']=mycolors
#print(cm.gray)
fig,ax=plt.subplots(2,2,figsize=(8,7))
l0=ax[0,0].scatter(umap_ori[:,0],umap_ori[:,1],c=celltype,cmap=mpl.rcParams["image.cmap"])
l1=ax[0,1].scatter(umap_avg[:,0],umap_avg[:,1],c=celltype,cmap=mpl.rcParams["image.cmap"])#,s=8,marker='o')
l2=ax[1,0].scatter(umap_sum[:,0],umap_sum[:,1],c=celltype,cmap=mpl.rcParams["image.cmap"])#,s=8,marker='o')
l3=ax[1,1].scatter(points[:,0],points[:,1],c=celltype,cmap=mpl.rcParams["image.cmap"])#,s=8,marker='o')
#l3=ax[1,1].scatter(umap_normalized[:,0],umap_normalized[:,1],c=celltype,cmap=mpl.rcParams["image.cmap"])#,s=8,marker='o')


#ax[0].legend()
legend1 = ax[0,1].legend(*l1.legend_elements(),
                    ncol=1,  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. , prop={"size":6},title="cell types")#bbox_to_anchor=(1, 0.5)
ax[0,1].add_artist(legend1)


ax[0,0].set_title('original expression (E)')
ax[0,1].set_title(r'(abs coeff. of central cell)')
ax[1,0].set_title(r'(coeff. of central cell)')#abs($\theta$)')
#ax[0,1].set_title(r'(E)*(abs coeff. of central cell)')
#ax[1,0].set_title(r'(E)*(coeff. of central cell)')#abs($\theta$)')
ax[1,1].set_title('cells in tissue')

#plt.axis('off')
for i in range(2):
    for j in range(2):
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        ax[1,j].set_xlabel('UMAP1')
        ax[i,0].set_ylabel('UMAP2')

#text=ax.text(0,0,"umap")
#fig.legend([l1,l2,l3],labels=legend,loc="center right",borderaxespad=0.1,    # Small spacing around legend box
#           title="cell types")
#plt.subplots_adjust(right=0.85)
plt.tight_layout()
fig.savefig('4umap_'+factors+'.png',bbox_inches='tight')
#fig.savefig('4umap_'+factors+'.png', bbox_extra_artists=(legend1), bbox_inches='tight')
