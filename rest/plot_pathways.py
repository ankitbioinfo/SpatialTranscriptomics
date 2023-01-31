

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, ConvexHull,voronoi_plot_2d, Delaunay
import os
import random
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import seaborn as snn
from matplotlib.lines import Line2D
import scanpy as sc
from matplotlib.cm import ScalarMappable
import matplotlib as mpl
import pickle
from scipy.stats import pearsonr
from gseapy.plot import gseaplot, heatmap
import gseapy




def create_directory(outputFolder):
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)




def find_interest_of_genes(ct1,p1,fname,cutoff):
    xl = pd.ExcelFile(fname)
    sheet=[]
    for i in range(len(xl.sheet_names)):
        if xl.sheet_names[i].find('HVG')!=-1:
            sheet.append(xl.sheet_names[i])

    #print(xl.sheet_names,sheet)
    df=pd.read_excel(fname,sheet_name=sheet[p1-1])
    data1=df.to_numpy()
    #df2=pd.read_excel(fname,sheet_name=sheet[p2-1])
    #data2=df2.to_numpy()

    #print(df[t1])
    header=df.columns
    for i in range(len(header)):
        if header[i]==ct1:
            first=[i,i+p1]
        #if header[i]==ct2:
        #    second=[i,i+p2]

    source=data1[:,first]
    #target=data2[:,second]

    #print(source)

    #cutoff=0.1
    #cutoff=-1000
    #print('first',source)
    #print('seco',target)
    interestofGene1=[]
    v1=[]
    for i in range(len(source)):
        if abs(source[i,1])>cutoff:
            interestofGene1.append(source[i,0])
            v1.append(source[i,1])
            #v1.append('%0.2f'%source[i,1])


    interestofGene2=[]
    v2=[]
    '''
    for i in range(len(target)):
        if (target[i,1])>cutoff:
            interestofGene2.append(target[i,0])
            v2.append(target[i,1])
            #v2.append('%0.2f'%target[i,1])
    '''


    return interestofGene1,interestofGene2,v1,v2



def pathway_analysis(titlename,g1,savename):
    database=['BioPlanet_2019','Reactome_2016']
    pathwayorganism='Mouse'
    for i in range(len(database)):
        titlename1=titlename+'_'+database[i]
        finalsavename=savename+titlename1
        titlename1=titlename1+'_#G='+str(len(g1))
        #print("tt",titlename)
        #print("fin",finalsavename)
        enr_res1 = gseapy.enrichr(gene_list=g1,organism=pathwayorganism,gene_sets=database[i], description='pathway',cutoff = 0.5)
        #enr_res1 = gseapy.enrichr(gene_list=g1,organism='Mouse',gene_sets=background_model,description='pathway',cutoff = 0.5)
        finalsavename.replace(' ','_')


        try:
            gseapy.barplot(enr_res1.res2d,title=titlename1,ofname=finalsavename)#database[i]+titlename
        except Exception as e: #Exception: Error getting the Enrichr libraries
            pass






def main():

    cutoff=0.5
    #mct1=['Paneth','Tuft','Stem_TA','BZE']
    #mp1=[3,1,3,1]    # single cell'B_cells','Basophils',

    celltype=['Central_Vein_EC','Cholangiocytes','Fibroblast','Hep378_MC','Hep145_P',
    'Hep2_MP','Hep6_C','HsPCs','ILC1s','KCs','LSECs','Lymphatic_EC','Macro','Mig_cDCs',
    'Monocytes','NK_cells','Neutrophils','Portain_Vein_EC','Stellate_cells','T_cells','cDC2s','cDC1s','pDCs']

    mct1=['Central_Vein_EC', 'Fibroblast', 'Stellate_cells', #1
    'Cholangiocytes', 'Portain_Vein_EC','Hep145_P','Hep2_MP', #2
    'LSECs','Hep2_MP','Portain_Vein_EC','Hep2_MP', #3
    'LSECs',#'Hep6_C'  #4
     'KCs','Hep378_MC','Hep145_P','Hep145_P',#5
    'LSECs','Hep378_MC', 'KCs', #6
    'KCs','T_cells','Fibroblast', #7
    'LSECs','cDC2s', 'Hep145_P','Hep145_P','Hep2_MP',   #8
    'ILC1s','Monocytes','LSECs',
    'Hep2_MP','Hep6_C','Hep6_C','Hep6_C',
    ]
    mp1=[4,5,4, #1
    5,5,3,5,    #2
    3,4,4,2,    #3
    1,#5,4,     #4
    4,4,4,5,    #5
    4,5,1,      #6
    5,1,4,       #7
    3,4,2,5,4,           #8
    3,4,2,
    2,4,5,3
    ]

    mct1=['Hep145_P','Hep145_P','Hep145_P','Hep145_P','Hep145_P',
    'Hep2_MP','Hep2_MP','Hep2_MP','Hep2_MP','Hep2_MP',
    'Hep378_MC','Hep378_MC','Hep378_MC','Hep378_MC','Hep378_MC']

    #mct1=['Hep124_P']
    mp1=[1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]

    mp1=[]
    mct1=[]
    for i in range(len(celltype)):
        for j in range(5):
            mp1.append(j+1)
            mct1.append(celltype[i])


    print(len(mct1),len(mp1))
    print(mp1)
    print(mct1)

    inputpath='./spatial_ct_ct_interactions_poss4_res0.5_dis15/Regression_on_genePC_loading_0_withoutscaling/'
    fname=inputpath+'gene_correlation5.xlsx'
    #fname='./gene_correlation3_single.xlsx'
    savename=inputpath+'figpathway/'
    create_directory(savename)


    for i in range(59,len(mct1)): #33
        ct1=mct1[i]
        p1=mp1[i]
        # only one is required to plot two is not needed
        # one is dummy plot
        ga1,gb1,va1,vb1=find_interest_of_genes(ct1,p1,fname,cutoff)
        print(i,ct1,p1,len(ga1))
        print(i,ga1[0:5])
        #print('\n\n\n\n')
        pathway_analysis('PC_'+str(p1)+'_'+ct1+'_cut_'+str(int(100*cutoff)),ga1,savename)






main()
