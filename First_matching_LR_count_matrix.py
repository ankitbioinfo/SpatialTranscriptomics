

import pandas as pd


f=open('sort_3_db_L_R_high_confident.dat')
totalLR={}
for line in f:
    l=line.split()
    totalLR[l[0]]=1
    totalLR[l[1]]=1

print('LR',len(totalLR))



filename='gene_by_cell_counts.csv'
df=pd.read_csv(filename,sep=',')
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

ndf.to_csv('LRgenebycell.csv',index=True, index_label="GENEID",  header=cellname)

f=open('total_LR_genes_exist.dat','w')
for i in range(len(LRgenes)):
    f.write(LRgenes[i].upper()+'\t'+str(i+1)+'\n')
