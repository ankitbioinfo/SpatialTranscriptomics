
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, ConvexHull,voronoi_plot_2d, Delaunay
import os, sys
import random
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from scipy.stats import mode



def euclidean_dist(p1,p2):
	return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

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
        newneig.append(neighbor[i])
    return newneig

'''
def findNeighbors_in_given_radius(location,radius):
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
        try:
            newneig.append(neighbor[i])
        except KeyError:
            newneig.append([])

    return newneig
'''

def reading_clustering():
	f4=open('leiden_output.dat')
	cont=f4.readlines()
	celltype={}
	noct=[]
	cellsinCT={}
	order_of_cellid={}
	total=len(cont)-1
	for i in range(1,len(cont)):
	    l=cont[i].split(',')
	    id=int(l[1])
	    celltype[l[0]]=id
	    order_of_cellid[l[0]]=i-1
	    if id not in noct:
	        noct.append(id)

	    if id not in cellsinCT:
	        cellsinCT[id]=[l[0]]
	    else:
	        cellsinCT[id].append(l[0])

	print('no of cell types',len(noct))
	return celltype,order_of_cellid

def read_tissue_coordinates(radius,order_of_cellid):
	f=open('tissue_positions_list.csv')
	cont=f.readlines()
	location_cellname2int={}
	location_int2cellname={}
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

	#tri=Delaunay(PP)
	#neighbors = find_DT_neighbors (tri,PP)

	neighbors=findNeighbors_in_given_radius(PP,radius)
	print("total neighbor", len(neighbors))



	return neighbors,location_cellname2int,location_int2cellname




def read_expression():
    	df=pd.read_csv('LRgenebycell.csv',sep=',')
    	expression=df.to_numpy()
    	topExpressedGene=[]
    	for i in range(len(expression)):
        	gene=expression[i,0].upper()
        	topExpressedGene.append(i)
    	return expression,df

def make_neighborhood_LRmatrix(radius,df,location_cellname2int,location_int2cellname,neighbors,ligands,celltype):

	fw=open('neighobor.dat','w')
	csvFile_cell_name={}
	pairwise={}
	for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    csvFile_cell_name[key]=i
	    cell1=location_cellname2int[key]
	    #print(cell1)
	    V=neighbors[cell1]
	    fw.write(key+'\t'+str(neighbors[cell1])+'\n')
	    for j in range(len(V)):
	        t=sorted([cell1,V[j]])
	        name=str(t[0])+'-'+str(t[1])
	        pairwise[name]=1


	fL=open('neighbors_weight_'+str(radius)+'.dat','w')
	for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    cell1=location_cellname2int[key]
	    V=neighbors[cell1]
	    #CT2=np.zeros(len(noct),dtype=np.float)
	    lig=np.zeros((len(V),ligands.shape[0]),dtype=np.float)
	    #rec=np.zeros((len(V),receptors.shape[0]),dtype=np.float)
	    good=[]
	    neighbortypes=[]
	    for j in range(len(V)):
	        name=location_int2cellname[V[j]]
	        try:
	        	y=csvFile_cell_name[name]
		        lig[j]=ligands[:,y]
		        #rec[j]=receptors[:,y]
		        #CT2[celltype[name]]+=1.0/prop[celltype[name]]
		        neighbortypes.append(celltype[name])
		        good.append(j)
	        except KeyError:
	            pass
	    unique_nt=np.unique(neighbortypes)
	    #print(neighbortypes,unique_nt)

	    final_lig=np.zeros((len(unique_nt),ligands.shape[0]),dtype=np.float)
	    #final_rec=np.zeros((len(unique_nt),receptors.shape[0]),dtype=np.float)

	    weighted_cell_neighbor=0.0
	    for j in range(len(unique_nt)):
	    	consider=[]
	    	for k in range(len(neighbortypes)):
	    		if neighbortypes[k]==unique_nt[j]:
	    			consider.append(good[k])
	    	final_lig[j]=len(consider)*np.mean(lig[consider],axis=0)
	    	#final_rec[j]=len(consider)*np.mean(rec[consider],axis=0)
	    	weighted_cell_neighbor+=len(consider)

	    if len(V)>0:
		    fL.write(str(cell1)+'\t'+str(celltype[key]))
		    #CT2=CT2/np.sum(CT2)
		    a=(1.0/weighted_cell_neighbor)*np.mean(final_lig,axis=0)
		    #print(a.shape,b.shape)
		    for j in a:
		        fL.write('\t'+str(j))
		    fL.write('\n')


def make_top_20_celltype_expressed_in_cell_types_of_central_cell(df,ligands,celltype):

	ct=[]
	for key in celltype:
		if celltype[key] not in ct:
			ct.append(celltype[key])
	ct=sorted(ct)

	genes=ligands[:,0]

	#print(ct,ligands.shape,genes)

	M=np.zeros((len(ct),ligands.shape[0]),dtype=float)
	counts=np.zeros((len(ct)),dtype=float)

	#print(celltype)
	for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    a=ligands[:,[i]]
	    #print('a',type(a[0]),type(M[0]),celltype[key],a)
	    M[celltype[key]]=M[celltype[key]]+a.T
	    counts[celltype[key]]+=1
	    #print(key,celltype[key],ligands[:,i].shape)

	f=open('BiologicalNameOfCT.dat')
	nameOfCellType={}
	for line in f:
		l=line.split('\t')
		nameOfCellType[int(l[0])]=l[1]


	fw=open('highest_to_lowest_expressed_genes_in_celltypes.dat','w')
	for i in range(len(ct)):
		fw.write(str(ct[i]))
		#fw.write(nameOfCellType[ct[i]] )
		M[ct[i]]=M[ct[i]]/counts[ct[i]]
		index=np.argsort(-M[ct[i]])
		#print(M[i])
		for j in range(len(index)):
			#fw.write(';'+genes[index[j]] + ' (%0.2f'%M[ct[i]][index[j]]+') '      )
			fw.write(';'+genes[index[j]]    )
		fw.write('\n')


#df=pd.read_csv('neighborhoodData1/neighborhood_expression_ecm_sum.csv',sep=',')
#df=pd.read_csv('neighborhoodData/original_ligands.csv',sep=',')
#ligands=df.to_numpy()



#celltype['B1_cell166']=10

def main():
	ligands,df=read_expression()
	print('ligands', ligands.shape)
	#radius=int(sys.argv[1])
	radiusList=[50,75,100,150]
	celltype,order_of_cellid=reading_clustering()

	make_top_20_celltype_expressed_in_cell_types_of_central_cell(df,ligands,celltype)


	for i in radiusList:
		radius=i
		neighbors,location_cellname2int,location_int2cellname=read_tissue_coordinates(radius,order_of_cellid)
		make_neighborhood_LRmatrix(radius,df,location_cellname2int,location_int2cellname,neighbors,ligands,celltype)


main()




'''
noct=sorted(noct)
prop=[]
for i in range(len(noct)):
    print(noct[i],len(cellsinCT[noct[i]]))
    prop.append(1.0*len(cellsinCT[noct[i]])/total)
print(prop)
'''
