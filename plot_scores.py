
from matplotlib import pyplot as plt
import numpy as np



filename=['LRF_L2_multi','LRF_L2_ovr','LRF_L1_multi','LRF_L1_ovr','LRF_elasticnet_multi','LRF_elasticnet_ovr']


C=[1,0.1,10]
R=[50,75,100,150]


fig,axs=plt.subplots(2,2,figsize=(10,7))
pindex=[]
for i in range(2):
    for j in range(2):
        pindex.append([i,j])


color=['r','b','m','k','g','y']
symbol=['s','o','v']

#xlabels=['R=50,C=1','R=75,C=1','R=100,C=1','R=150,C=1','R=50,C=0.1','R=75,C=0.1','R=100,C=0.1','R=150,C=0.1','R=50,C=10','R=75,C=10','R=100,C=10','R=150,C=10',]
strategy=['L2 mul','L2 ovr','L1 mul','L1 ovr','EN mul','EN ovr' ]

xlabels=['accuracy','macro F1','macro precision','macro recall','micro [all]','weighted F1','weighted precision','weighted recall','cohen kappa','mcc','heming loss']
xlabels=['accuracy','macro F1','macro precision','macro recall','micro [all]','weighted F1','weighted precision','weighted recall','cohen kappa','mcc']


index=[0,1,2,3,4,7,8,9,10,12]
for j in range(len(R)):
    data=[]
    for fi in range(len(filename)):
        for i in range(len(C)):
            name=filename[fi]+'/matrix_score_'+str(C[i])+'_'+str(R[j])+'.dat'
            data=np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
            mycolor=color[fi]+symbol[i]+'-'
            axs[pindex[j][0],pindex[j][1]].plot(data[index],mycolor,label=strategy[fi]+','+str(C[i]))
    axs[pindex[j][0],pindex[j][1]].set_title('Radius = '+str(R[j]))

    if j==3:
        legend1= axs[pindex[j][0],pindex[j][1]].legend(loc='upper right',bbox_to_anchor=(1.25, 1),ncol=1, borderaxespad=0., prop={"size":6},fancybox=True, shadow=True)

    axs[pindex[j][0],pindex[j][1]].set_xticks(range(len(index)))
    if j>1:
        axs[pindex[j][0],pindex[j][1]].set_xticklabels(xlabels)
        for tick in axs[pindex[j][0],pindex[j][1]].get_xticklabels():
            tick.set_rotation(90)

    if (j==0)|(j==2):
        axs[pindex[j][0],pindex[j][1]].set_ylabel('score')


fig.tight_layout()

fig.savefig('scores_of_all.png',bbox_inches='tight')
