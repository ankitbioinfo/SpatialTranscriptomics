
import pandas as pd
import numpy as np


filename='Liver1Slice1_cell_by_gene.csv'
df=pd.read_csv(filename,sep=',')
expression=df.to_numpy()
print('shape',expression.shape)


cellname=[]
d={}
temp=expression[:,0]
for i in range(len(temp)):
    name='cell'+str(i)
    cellname.append(name)
    d[temp[i]]=name


LRgenes=[]
gene_index=[]
for i in range(1,len(df.columns)):
    if df.columns[i].find('Blank')==-1:
        gene_index.append(i)
        LRgenes.append(df.columns[i])
    #print(i,df.columns[i])

print(len(gene_index))


mat=np.transpose(expression[:,gene_index])

ndf=pd.DataFrame(mat,index=LRgenes)#print(expression.shape,expression[:,0])

ndf.to_csv('Blank_genes_removed.csv',index=True, index_label="GENEID",  header=cellname)


f2=pd.read_csv('Liver1Slice1_cell_metadata.csv',sep=',')
data=f2.to_numpy()

f3=open('tissue_positions_list.csv','w')

for i in range(len(data)):
    #cell15,0,0,0,-3702.522707,-3904.721471

    f3.write(d[data[i][0]]+','+ str(data[i,3])+','+str(data[i,4]) +'\n')
