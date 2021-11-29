
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
	return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2  + (p1[2]-p2[2])**2)

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
        t5=(loc1[2]-1.1*radius) <= location[:,2]
        t6=location[:,2] <= (loc1[2]+1.1*radius)

        index=  np.where ( t1 & t2 & t3 & t4  & t5 & t6    )
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
    avg_neigh=0.0
    for i in range(n):
        try:
            l=neighbor[i]
        except KeyError:
            l=[]
        #print(l)
        newneig.append(l)
        avg_neigh+=len(l)

    print('average neighbors:',avg_neigh/n)

    return newneig


def reading_data(positionFilename,clusterFilename,celltypeFilename,expressionFilename,radius,saveCommunication):


    f4=open(clusterFilename)
    cont=f4.readlines()
    celltype={}
    noct=[]
    cellsinCT={}
    order_of_cellid={}
    louvain=[]
    total=len(cont)-1
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

    louvain=np.array(louvain)
    #print(louvain[0:10])
    total=len(louvain)

    f=open(positionFilename)
    cont=f.readlines()
    location_cellname2int={}
    location_int2cellname={}
    points=np.zeros((len(cont),3),dtype=float)
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
            location_cellname2int[name]=i
            location_int2cellname[i]=name
            if len(l)==3: # this is 2d
	                points[i]=[float(l[1]), float(l[2]),0]
            else:
	                points[i]=[float(l[1]), float(l[2]),float(l[3])]

    PP=points

    CTname=[]
    with open(celltypeFilename,'r') as f:
        cont = f.read()
        lines=cont.split('\n')
        CTname=[]
        for i in range(len(lines)):
            l=lines[i].split('\t')
            if len(l)>1:
                name=l[1].replace('/','_')
                name=name.replace(' ','_')
                name=name.replace('"','')
                name=name.replace("'",'')
                name=name.replace(')','')
                name=name.replace('(','')
                name=name.replace('+','p')
                name=name.replace('-','n')
                CTname.append(name)

    noct=sorted(noct)
    actual_count=[]
    fraction_CT=[]
    for key in noct:
        actual_count.append(len(cellsinCT[key]))
        fraction_CT.append(len(cellsinCT[key])/float(total))

    temp=np.where(np.array(actual_count)>=5)
    good_index_cell_counts=temp[0]
    less_no_cells_remove=[]
    for i in range(len(good_index_cell_counts)):
        index=np.where(louvain==good_index_cell_counts[i])
        less_no_cells_remove+=list(index[0])

    less_no_cells_remove=sorted(less_no_cells_remove)
    PP=PP[less_no_cells_remove]
    louvain=louvain[less_no_cells_remove]

    new_CT_id={}
    for i in range(len(good_index_cell_counts)):
        new_CT_id[good_index_cell_counts[i]]=i

    for i in range(len(louvain)):
        value=louvain[i]
        louvain[i]=new_CT_id[value]



    fw=open(saveCommunication+'BiologicalNameOfCT.dat','w')
    for i in range(len(new_CT_id)):
            value=fraction_CT[good_index_cell_counts[i]]
            name=CTname[good_index_cell_counts[i]]
            #print(CTname[key], key, total,len(cellsinCT[key]))
            fw.write(str(i)+'\t'+name+'\t'+str('%0.4f'%value)+'\n')
    fw.close()


    df=pd.read_csv(expressionFilename,sep=',')
    expression=df.to_numpy()
    topExpressedGene=[]
    f=open(saveCommunication+'total_LR_genes_exist.dat','w')
    for i in range(len(expression)):
        gene=expression[i,0]
        topExpressedGene.append(i)
        f.write(gene+'\t'+str(i+1)+'\n')
    f.close()


	#return celltype,order_of_cellid


    neighbors=findNeighbors_in_given_radius(PP,radius)
    ct=[]
    for key in celltype:
        if celltype[key] not in ct:
            ct.append(celltype[key])
    ct=sorted(ct)
    genes=expression[:,0]
	#print(ct,ligands.shape,genes)
    M=np.zeros((len(ct),expression.shape[0]),dtype=float)
    counts=np.zeros((len(ct)),dtype=float)
	#print(celltype)
    for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    a=expression[:,[i]]
	    #print('a',type(a[0]),type(M[0]),celltype[key],a)
	    M[celltype[key]]=M[celltype[key]]+a.T
	    counts[celltype[key]]+=1
	    #print(key,celltype[key],ligands[:,i].shape)
    fw=open(saveCommunication+'highest_to_lowest_expressed_genes_in_celltypes.dat','w')
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


    return celltype,neighbors,location_cellname2int,location_int2cellname, expression, df




def make_neighborhood_LRmatrix(radius,df,location_cellname2int,location_int2cellname,neighbors,expression,celltype,saveCommunication):
	#fw=open(saveCommunication+'neighobor'+str(radius)+'.dat','w')
	csvFile_cell_name={}
	pairwise={}
	for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    csvFile_cell_name[key]=i
	    cell1=location_cellname2int[key]
	    #print(cell1)
	    V=neighbors[cell1]
	    #fw.write(key+'\t'+str(neighbors[cell1])+'\n')
	    for j in range(len(V)):
	        t=sorted([cell1,V[j]])
	        name=str(t[0])+'-'+str(t[1])
	        pairwise[name]=1

	fL=open(saveCommunication+'weighted_normalized_comm_neighbors_'+str(radius)+'.dat','w')
	for i in range(1,df.shape[1]):
	    key=df.columns[i]
	    cell1=location_cellname2int[key]
	    V=neighbors[cell1]
	    #CT2=np.zeros(len(noct),dtype=np.float)
	    lig=np.zeros((len(V),expression.shape[0]),dtype=float)
	    #rec=np.zeros((len(V),receptors.shape[0]),dtype=np.float)
	    good=[]
	    neighbortypes=[]
	    for j in range(len(V)):
	        name=location_int2cellname[V[j]]
	        try:
	        	y=csvFile_cell_name[name]
		        lig[j]=expression[:,y]
		        #rec[j]=receptors[:,y]
		        #CT2[celltype[name]]+=1.0/prop[celltype[name]]
		        neighbortypes.append(celltype[name])
		        good.append(j)
	        except KeyError:
	            pass
	    unique_nt=np.unique(neighbortypes)
	    #print(neighbortypes,unique_nt)

	    final_lig=np.zeros((len(unique_nt),expression.shape[0]),dtype=float)
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


def find_ligand_receptor_genes_in_expression(expressionFilename,saveCommunication):
	f=open('sort_3_db_L_R_high_confident.dat')
	totalLR={}
	for line in f:
	    l=line.split()
	    totalLR[l[0]]=1
	    totalLR[l[1]]=1
	print('LR',len(totalLR))

	df=pd.read_csv(expressionFilename,sep=',')
	expression=df.to_numpy()
	print('shape',expression.shape)

	gene=expression[:,0]
	cellname=df.columns[1:]
	print('gene and cell', len(gene),len(cellname))
	index=[]
	LRgenes=[]
	for i in range(len(gene)):
	    x=gene[i].upper()
	    try:
	        totalLR[x]
	        index.append(i)
	        LRgenes.append(gene[i])
	    except KeyError:
	        pass

	print('index',len(index))
	mat=expression[index,1:]
	ndf=pd.DataFrame(mat,index=LRgenes)#print(expression.shape,expression[:,0])
	ndf.to_csv(saveCommunication+'LR_gene_by_cell.csv',index=True, index_label="GENEID",  header=cellname)



def create_directory(outputFolder):
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)

def main():
	#ligands,df=read_expression()
	#print('ligands', ligands.shape)
	#radius=int(sys.argv[1])
	radiusList=[25,50]#100,150,200,250,300]
	inputCluster='louvain_cluster.dat'
	#inputCluster='leiden_cluster.dat'

	inputPositions='tissue_positions_list.csv'
	inputCellTypeName='NameOfCT.dat'
	saveCommunication='communication/'
	expressionFilename='gene_by_cell_counts.csv'
	create_directory(saveCommunication)
	find_ligand_receptor_genes_in_expression(expressionFilename,saveCommunication)


	for i in radiusList:
		radius=i
		expressionFilename=saveCommunication+'LR_gene_by_cell.csv'
		celltype,neighbors,location_cellname2int,location_int2cellname, expression, df=reading_data(inputPositions,inputCluster,inputCellTypeName,expressionFilename,radius,saveCommunication)
		make_neighborhood_LRmatrix(radius,df,location_cellname2int,location_int2cellname,neighbors,expression,celltype,saveCommunication)


main()
