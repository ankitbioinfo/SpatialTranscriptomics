


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


def findNeighbors_in_given_radius(location,radius):
    n=location.shape[0]
    #print('ss',location.shape)
    neighbor={}
    for i in range(n):
        loc1=location[i]
        #print(loc1)
        t1=(loc1[0]-1.1*radius) <= location[:,0]
        t2=location[:,0] <= (loc1[0]+1.1*radius)
        t3=(loc1[1]-1.1*radius) <= location[:,1]
        t4=location[:,1] <= (loc1[1]+1.1*radius)

        index=  np.where ( t1 & t2 & t3 & t4     )
        #print(index[0])
        #neighborIndex= index[0].remove(i)
        #print( neighborIndex)

        #fw=open('ankit'+str(i),'w')
        #for k in range(len(t1)):
        #    fw.write(str(location[k])+'\t'+str(t1[k])+'\t'+str(t2[k])+'\t'+str(t3[k])+'\t'+str(t4[k])+'\n')

        #for j in range(i+1,n):
        count=0

        for k in range(len(index[0])):
            j=index[0][k]
            if j!=i:
                count+=1
                loc2=location[j]
                dist=euclidean_dist(loc1,loc2)
                if dist<radius:
                    if i not in neighbor:
                        neighbor[i]=[j]
                    else:
                        if j not in neighbor[i]:
                            neighbor[i].append(j)

                    if j not in neighbor:
                        neighbor[j]=[i]
                    else:
                        if i not in neighbor[j]:
                            neighbor[j].append(i)

        #print('t',count,len(index[0]))


    newneig=[]
    for i in range(n):
        try:
            l=neighbor[i]
        except KeyError:
            l=[]
        #print(l)
        newneig.append(l)

    return newneig





def findNeighbors_in_given_radius2(location,radius):
    n=location.shape[0]
    neighbor={}
    for i in range(n):
        loc1=location[i]
        for j in range(i+1,n):
            loc2=location[j]
            dist=euclidean_dist(loc1,loc2)
            if dist<radius:
                if i not in neighbor:
                    neighbor[i]=[j]
                else:
                    if j not in neighbor[i]:
                        neighbor[i].append(j)

                if j not in neighbor:
                    neighbor[j]=[i]
                else:
                    if i not in neighbor[j]:
                        neighbor[j].append(i)

    newneig=[]
    for i in range(n):
        newneig.append(neighbor[i])

    return newneig


def euclidean_dist(p1,p2):
	return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)



def reading_data():
    #f=pd.read_excel('aau5324_Moffitt_Table-S7_coordinate.xlsx',engine='openpyxl',sheet_name=0,header=1)
    #data=f.to_numpy()
    name='gene_by_cell_counts.csv'
    data = np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
    print(data.shape)
    #print(data[0,7:])

    '''
    leiden=[]
    for i in range(data.shape[0]):
        name='CT'+str(int(data[i,1]))
        #leiden.append(cname2int[name])
    '''

    f4=open('leiden_output.dat')
    cont=f4.readlines()
    celltype={}
    order_of_cellid={}
    noct=[]
    cellsinCT={}
    total=len(cont)-1
    louvain=[]
    for i in range(1,len(cont)):
        l=cont[i].split(',')
        id=int(l[1])
        louvain.append(id)
        celltype[l[0]]=id
        order_of_cellid[l[0]]=i-1
        if id not in noct:
            noct.append(id)

        if id not in cellsinCT:
            cellsinCT[id]=[l[0]]
        else:
            cellsinCT[id].append(l[0])

    print('no of cell types',len(noct))


    total=len(louvain)



    #CTname=['EryCel4','CD8','EryCel2','EryCel3','CD4','HSCM',
    #'Bcel','EryCel1','AEC','Tcel','Tcel2',
    #'Stel2','EP','Stel1','MK','MacroP','Myeloid','Endo']
    #fw=open('BiologicalNameOfCT.dat','w')
    fraction_CT=[]
    noct=sorted(noct)
    for key in noct:
        fraction_CT.append(len(cellsinCT[key])/float(total))
        #print(CTname[key], key, total,len(cellsinCT[key]),fraction_CT)
        #fw.write(str(key)+'\t'+CTname[key]+'\t'+str('%0.3f'%fraction_CT[key])+'\n')


    f=open('tissue_positions_list.csv')
    cont=f.readlines()
    location_cellname2int={}
    location_int2cellname={}
    #this code got changed due to change in the number of cells
    points=np.zeros((len(order_of_cellid),2),dtype=np.float)
    for i in range(len(cont)):
        l=cont[i].split(',')
        #name=l[0].replace('-','.')
        name=l[0]
        try:
            flag=1
            index=order_of_cellid[name]
        except KeyError:
            flag=0

        if flag==1:
            location_cellname2int[name]=index
            location_int2cellname[index]=name
            points[index]=[float(l[4]), float(l[5])]

    PP=points
    return PP, louvain, noct,fraction_CT



def performa_analysis(radius):
    #radius=100
    tri=Delaunay(PP)
    #neighbors = find_DT_neighbors (tri,PP)
    neighbors=findNeighbors_in_given_radius(PP,radius)
    n=len(neighbors)
    print("total neighbor",n )


    fw=open('neighbors_logistic_regression_normalized_'+str(radius)+'.dat','w')
    expectedNeighbors=[]
    #print(noct)
    for i in range(n):
        cell1=i
        CT1=louvain[i]
        V=neighbors[i]
        CT2=np.zeros(len(noct),dtype=np.float)
        for j in range(len(V)):
            name=louvain[V[j]]
            try:
                CT2[name]+=1.0
            except KeyError:
                pass
        fw.write(str(cell1)+'\t'+str(CT1))
        expected=np.array(fraction_CT)*np.sum(CT2)
        tt=CT1
        #print(np.concatenate(np.array(celltype[key]),CT2))
        expectedNeighbors.append(np.concatenate([np.asarray([tt]),CT2]))
        #print(expectedNeighbors)
        CT2=CT2/expected   #np.sum(CT2) #np.linalg.norm(CT2)
        for j in CT2:
            fw.write('\t'+'%0.5f'%j)
        fw.write('\n')

    expectedNeighbors=np.array(expectedNeighbors)

    M=[]
    for i in range(len(noct)):
        a=np.where(expectedNeighbors[:,0]==i)
        b=np.where(expectedNeighbors[:,0]!=i)
        #print('a',len(a[0]),len(b[0]))
        myCT=np.mean(expectedNeighbors[a[0],1:],axis=0)
        remainCT=np.mean(expectedNeighbors[b[0],1:],axis=0)
        M.append(myCT/remainCT)
        #print(i,M[i])


    M=np.array(M)
    return M




plt.figure(figsize=(13,7))
#snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classes, yticklabels=classes)
PP,louvain,noct,fraction_CT= reading_data()

M=performa_analysis(50)
plt.subplot(2,3,1)
snn.heatmap(M)
plt.title('Radius =50')
plt.ylabel('Avg neighborhood in own CT')


M=performa_analysis(100)
plt.subplot(2,3,2)
snn.heatmap(M)
plt.title('Radius =100')


M=performa_analysis(150)
plt.subplot(2,3,3)
snn.heatmap(M)
plt.title('Radius =150')
#plt.xlabel('Avg neighborhood in remaining CT')

M=performa_analysis(75)
plt.subplot(2,3,6)
snn.heatmap(M)
plt.title('Radius =300')
plt.xlabel('Avg neighborhood in remaining CT')


#plt.savefig('Neighborhood.png')

#make_neighborhood_matrix(df,ecm,neighbors,NeighborDistance,csvFile_cell_name,'ecm',PP,celltype,noct,cellsinCT)
