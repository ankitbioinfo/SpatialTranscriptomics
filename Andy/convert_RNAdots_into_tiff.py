
import sys,os
#from bs4 import BeautifulSoup
from skimage import io,color
from skimage.color import rgb2gray
import skimage as sk
import numpy as np
#import matplotlib.pyplot as plt
#from scipy.interpolate import splprep, splev
import pandas as pd 




name=sys.argv[1]
f=open(name,'r')
cont=f.readlines()
#data = np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
#table.columns
'''
Index(['intensity', 'z', 'y', 'x', 'radius', 'spot_id', 'z_min', 'z_max',
       'y_min', 'y_max', 'x_min', 'x_max', 'features', 'xc', 'yc', 'zc',
       'target', 'distance', 'passes_thresholds'],
      dtype='object')
'''

df=pd.read_excel('Cardiomyocytes_readout_numbering_for_genes.xlsx')
data=df.to_numpy()

genes=data[:,0]


def make_mask(data):
    mat=np.zeros((1648,2056),dtype=np.ubyte)
    mat=np.zeros((1648,935),dtype=np.ubyte)
    mat=np.zeros((1648,856),dtype=np.ubyte)
    mat=np.zeros((1768,1776),dtype=np.ubyte)

    for i in range(len(data)):
        mat[data[i][0]:data[i][1],data[i][2]:data[i][3]]=255
    return mat




na=name.split('_')
na=na[-1].split('.')
maindir='tif/parameter_'+na[0]+'/'
answer=os.path.isdir(maindir)
if answer==True:
	pass
else:
	os.mkdir(maindir)
fw=open(maindir+'temporaryFile.dat','w')

genecount=[]
for j in range(len(genes)):
    fw.write(genes[j]+'\n')
    data1=[]
    for i in range(len(cont)):
        l=cont[i].split(',')
        if l[16]==genes[j]:
            data1.append([int(l[8]),int(l[9]),int(l[10]),int(l[11])])


    genecount.append(len(data1))

    mygene=make_mask(data1)
    io.imsave(maindir+genes[j]+'_'+str(na[0])+'.tif',mygene)


fw.write('all genes'+'\n')
data1=[]
for i in range(len(cont)):
    l=cont[i].split(',')
    if (l[3]!='fake1')|(l[3]!='fake2')|(l[3]!='fake3'):
        data1.append([int(l[8]),int(l[9]),int(l[10]),int(l[11])])
genecount.append(len(data1))
mygene=make_mask(data1)
io.imsave(maindir+'All_genes'+'_'+str(na[0])+'.tif',mygene)



for j in range(len(genecount)):
    fw.write(str(genecount[j])+'\n')
