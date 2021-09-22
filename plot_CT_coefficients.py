

from sklearn.datasets import make_classification
from matplotlib import pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve

from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.pipeline import Pipeline

import pandas as pd
import numpy as np
import seaborn as snn
import os
from random import randrange
from scipy.stats import norm
#import matplotlib.mlab as mlab


def initialize(noct):
        data=[]
        for i in range(noct):
            data.append([])
        return data

R=[50,75,100,150]
C=[1,0.1,10]

filename=['LRF_L2_multi','LRF_L2_ovr','LRF_L1_multi','LRF_L1_ovr','LRF_elasticnet_multi','LRF_elasticnet_ovr']
strategy=['L2_mul','L2_ovr','L1_mul','L1_ovr','EN_mul','EN_ovr']


maindir='coefficients_strategy_wise/'
answer=os.path.isdir(maindir)
if answer==True:
    pass
else:
    os.mkdir(maindir)


f=open('BiologicalNameOfCT.dat')
nameOfCellType={}
for line in f:
    l=line.split('\t')
    nameOfCellType[int(l[0])]=l[1]

noct=len(nameOfCellType)
for j in range(len(R)):
    data=initialize(noct)
    xlabel=[]
    for fi in range(len(filename)):
        for i in range(len(C)):
            name=filename[fi]+'/matrix_coefficients_'+str(C[i])+'_'+str(R[j])+'.dat'
            try:
                d=np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
                flag=True
                #print(d.shape)
            except FileNotFoundError:
                flag=False
            if flag:
                xlabel.append(strategy[fi]+',C='+str(C[i]))
                for k in range(noct):
                    data[k].append(d[k])

    data=np.array(data)
    for i in range(noct):
        plt.figure(figsize=(5,3))
        b=snn.heatmap(data[i],yticklabels=xlabel)
        plt.xlabel('features')
        plt.title('R = '+str(R[j])+', '+ nameOfCellType[i])
        plt.ylabel('strategies')
        _, ylabels= plt.yticks()
        b.set_yticklabels(ylabels, size = 5)
        plt.tight_layout()
        plt.savefig(maindir+'weight_matrix_'+str(nameOfCellType[i])+'_R='+str(R[j])+'.png')
        plt.close()
    print(data.shape,xlabel)
