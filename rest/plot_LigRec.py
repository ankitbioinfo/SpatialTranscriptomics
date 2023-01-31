import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np


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
ind1=np.argsort(-abs(lig_cor))
ind2=np.argsort(-abs(rec_cor))

print(lig_cor[ind1[0:5]])
print(rec_cor[ind2[0:5]])


how_many_pairs_want2plot=25

n=len(ind1)
limit=min(how_many_pairs_want2plot,n)
flag=1
while flag:
    a=set(ind1[0:limit])
    b=set(ind2[0:limit])
    c=sorted(list(a.intersection(b)))
    print(len(a),len(b),len(c),limit)
    if (len(c)>how_many_pairs_want2plot)|(len(c)==len(a)):
        break
        flag=0
    else:
        limit=limit+10
        if limit>n:
            limit=n

#print(ind1[0:50],ind2[0:50])


#for i in range(20):
#    print(i,ind1[i],ind2[i])

A=[]
B=[]
for i in range(len(c)):
    A.append(AC[c[i]]+'_Fa_'+str(AF[c[i]])+'_'+lig[c[i]]) #str(i)+'_'+
    B.append(str(i)+'_'+BN[c[i]]+'_Fa_'+str(BF[c[i]])+'_'+rec[c[i]])

df = pd.DataFrame({'cols': B,
                   'rows': A,
                   'north': lig_cor[c],
                   'south': rec_cor[c],
                   'east': regCoff[c],
                   'west':localized[c]})
#'''

#print(df)
df['rows'] = pd.Categorical(df['rows'],categories=A)  # fix an ordering,
df_piv = df.pivot_table(index='rows', columns='cols')
#print('\n\n\n\n\n',df_piv)
M = len(df_piv.columns) // 4
N = len(df_piv)
values = [df_piv[dir] for dir in
          ['north', 'east', 'south', 'west']]  # these are the 4 column names in df



print('M',M)
print('N',N)
print('values',values)



triangul = triangulation_for_triheatmap(M, N)
cmaps = ['RdYlBu'] * 4
norms = [plt.Normalize(0, 1) for _ in range(4)]

fig, ax = plt.subplots(figsize=(8, 7))
imgs = [ax.tripcolor(t, np.ravel(val), cmap=cmap, norm=norm, ec='white')
        for t, val, cmap, norm in zip(triangul, values, cmaps, norms)]

#ax.tick_params(length=0)
ax.set_title('localizationCoef='+'%0.3f'%np.unique(localized)+',regressionCoef='+'%0.3f'%np.unique(regCoff))
ax.set_xticks(range(M))
ax.set_xticklabels(df_piv['north'].columns,rotation=90)
ax.set_yticks(range(N))
ax.set_yticklabels(df_piv.index)
ax.invert_yaxis()
ax.margins(x=0, y=0)
#ax.set_aspect('equal', 'box')  # square cells
plt.colorbar(imgs[0], ax=ax)
plt.tight_layout()
fig.savefig('LR1.png',dpi=300)
