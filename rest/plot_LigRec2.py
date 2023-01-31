import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np
import seaborn as snn
from matplotlib import gridspec




def triangulation_for_triheatmap(M, N):
    xv, yv = np.meshgrid(np.arange(-0.5, M), np.arange(-0.5, N))  # vertices of the little squares
    xc, yc = np.meshgrid(np.arange(0, M), np.arange(0, N))  # centers of the little squares
    x = np.concatenate([xv.ravel(), xc.ravel()])
    y = np.concatenate([yv.ravel(), yc.ravel()])
    cstart = (M + 1) * (N + 1)  # indices of the centers

    trianglesN = [(i + j * (M + 1), i + 1 + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesE = [(i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesS = [(i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesW = [(i + (j + 1) * (M + 1), i + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    return [Triangulation(x, y, triangles) for triangles in [trianglesN, trianglesE, trianglesS, trianglesW]]


days = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri']
df = pd.DataFrame({'cols': np.random.choice([*'abcdefghij'], 40),
                   'rows': np.random.choice(days, 40),
                   'north': np.random.rand(40),
                   'east': np.random.rand(40),
                   'south': np.random.rand(40),
                   'west': np.random.rand(40)})

#'''
inputpath='./spatial_ct_ct_interactions_poss4_res0.5_dis15/Regression_on_genePC_loading_0_withoutscaling/'
fname=inputpath+'Lig_and_Rec_enrichment_in_interacting_celltypes5.xlsx'
xl = pd.ExcelFile(fname)
sheet=xl.sheet_names
#print(,sheet)
df=pd.read_excel(fname,sheet_name=sheet)
data1=df['LR enrichment'].to_numpy()
data=data1[0:340]
#data=data1[0:16]


AC=data[:,1]
BN=data[:,2]
localized=data[:,3]
AF=data[:,4]
BF=data[:,5]
regCoff=data[:,6]
lig=data[:,7]
rec=data[:,8]
lig_cor=data[:,9]
rec_cor=data[:,10]
#ind1=np.argsort(-abs(lig_cor))
#ind2=np.argsort(-abs(rec_cor))
#print(lig_cor[ind1[0:5]])
#print(rec_cor[ind2[0:5]])

A=[]
B=[]
for i in range(len(data)):
    A.append(AC[i]+'_Fa_'+str(AF[i])+'_'+lig[i])
    B.append(BN[i]+'_Fa_'+str(BF[i])+'_'+rec[i])


A=np.sort(np.unique(A))
B=np.sort(np.unique(B))

ML=np.zeros((len(A), len(B)),dtype=float)
MR=np.zeros((len(A), len(B)),dtype=float)

for i in range(len(data)):
    name1=AC[i]+'_Fa_'+str(AF[i])+'_'+lig[i]
    name2=BN[i]+'_Fa_'+str(BF[i])+'_'+rec[i]
    row=0
    col=0
    for j in range(len(A)):
        if A[j]==name1:
            row=j
    for j in range(len(B)):
        if B[j]==name2:
            col=j
    ML[row,col]=lig_cor[i]
    MR[row,col]=rec_cor[i]



fig=plt.figure(figsize=(25,40))
#gs = gridspec.GridSpec(2, 1, width_ratios=[1, 1])
ax0=plt.subplot(2,1,1)#    plt.subplot(gs[0])
ax1=plt.subplot(2,1,2)#plt.subplot(gs[1])
fig.suptitle('LigRec Plot',fontsize=12)
snn.heatmap(ax=ax0,data=ML,annot=False, fmt='0.2f',   xticklabels=B, annot_kws={"size": 5},yticklabels=A)
snn.heatmap(ax=ax1,data=MR,annot=False, fmt='0.2f',xticklabels=B, annot_kws={"size": 5},yticklabels=A)
#plt.xlabel('Spatial cell clusters')
#plt.ylabel('Single cell clusters')

#plt.title('R = '+str(radius)+', C='+str(lambda_c))
#g.xaxis.set_ticks_position("top")
plt.tight_layout()
fig.savefig('LigRecAnalysis.png',dpi=300)
