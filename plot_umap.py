
#import scanpy as sc
import numpy as np
import pandas as pd
#import sklearn.datasets
import umap
#import umap.plot
from matplotlib import pyplot as plt

#from sklearn.feature_extraction.text import TfidfVectorizer
#from yellowbrick.datasets import load_hobbies
#from yellowbrick.text import UMAPVisualizer
import seaborn as sns


f=open('./../cellTypeIntegerName_cortex.dat')

d={}
for line in f:
    l=line.split('\t')
    d[int(l[1])]=l[0]


radius=[50,100,150,200,250,300]

for i in range(len(radius)):
    filename='neighbors_logistic_regression_normalized_'+str(radius[i])+'.dat'
    #csvfilename='genebycellmatrix_'+filename+'_full.csv'
    #csvfilename='neighbors_logistic_regression_normalized_200.dat'

    adata = np.genfromtxt(open(filename, "rb"), delimiter='\t', skip_header=0)
    mydata=adata[:,2:]
    target1=adata[:,1].astype(int)
    target=[]
    for j in range(len(target1)):
        target.append(str(target1[j])+':'+d[target1[j]])

    target=np.array(target)


    embedding = umap.UMAP(random_state=42).fit_transform(mydata)
    df = pd.DataFrame(dict(UMAP_1=embedding[:,0], UMAP_2=embedding[:,1], CellType = target))

    #fig,ax=plt.subplots(2,2,figsize=(8,7))
    plt.subplot(2,3,i+1)

    if i==2:
        gs=sns.scatterplot(x='UMAP_1', y='UMAP_2', data=df, s=20,hue='CellType',legend='full', palette='Paired')
        plt.legend(ncol=1,bbox_to_anchor=(1.02, 1), loc=2, prop={"size":6}, borderaxespad=0.)

    else:
        gs=sns.scatterplot(x='UMAP_1', y='UMAP_2', data=df, s=20,hue='CellType',legend=False, palette='Paired')
        #plt.legend(ncol=1,bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)


    #plt.figure(figsize=(13,7))
    #snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classes, yticklabels=classes)
    #plt.subplot(1,1,1)
    #mapper = umap.UMAP().fit(mydata)
    #embedding = umap.UMAP(n_neighbors=15).fit_transform(mydata)

    #if i==0:
    #    p=umap.plot.points(mapper, labels=target, color_key_cmap='Paired', background='white')
    #else:
    #    p=umap.plot.points(mapper,labels=target1,color_key_cmap='Paired', background='white')

    '''
    fig, ax = plt.subplots(1, figsize=(14, 10))
    plt.scatter(*embedding.T, s=0.3, c=target, cmap='Spectral', alpha=1.0)
    plt.setp(ax, xticks=[], yticks=[])
    cbar = plt.colorbar(boundaries=np.arange(11)-0.5)
    cbar.set_ticks(np.arange(10))
    cbar.set_ticklabels(classes)
    '''

    #p.legend(bbox_to_anchor=(1, 1.05), loc=2)
    plt.title('Radius = '+str(radius[i]))

    plt.xticks([])
    plt.yticks([])

    if (i==0)|(i==3):
        pass
    else:
        gs.set(ylabel=None)
    if (i<=2):
        gs.set(xlabel=None)

#plt.tight_layout()
plt.savefig('Neighborhood.png',bbox_inches='tight',dpi=150)
