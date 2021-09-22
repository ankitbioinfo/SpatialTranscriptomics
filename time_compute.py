

import matplotlib.pyplot as plt
import numpy as np

xlabels=['R=50,C=1','R=75,C=1','R=100,C=1','R=150,C=1','R=50,C=0.1','R=75,C=0.1','R=100,C=0.1','R=150,C=0.1','R=50,C=10','R=75,C=10','R=100,C=10','R=150,C=10',]



filename=['LRF_L2_multi','LRF_L2_ovr','LRF_L1_multi','LRF_L1_ovr','LRF_elasticnet_multi','LRF_elasticnet_ovr']
C=[1,0.1,10]
R=[50,75,100,150]

data=[]
for fi in range(len(filename)):
    t=[]
    for i in range(len(C)):
        for j in range(len(R)):
            name=filename[fi]+'/prediction_'+str(C[i])+'_'+str(R[j])+'.dat'
            f=open(name)
            cont=f.readlines()
            l=cont[-1].split('=')
            t.append(float(l[1]))
    data.append(t)




L2multi=data[0]
L2ovr=data[1]
L1multi=data[2]
L1ovr=data[3]
elasticnetmulti=data[4]
elasticnetovr=data[5]


fig, ax = plt.subplots()

ax.plot(np.array(L2multi)/60,'rs-',label='L2 multi lbfgs')
ax.plot(np.array(L1multi)/60,'bs-',label='L1 multi saga')
ax.plot(np.array(elasticnetmulti)/60,'ks-',label='En multi saga')

ax.plot(np.array(L2ovr)/60,'ro:',label='L2 ovr lbfgs')
ax.plot(np.array(L1ovr)/60,'bo:',label='L1 ovr liblinear')
ax.plot(np.array(elasticnetovr)/60,'ko:',label='EN ovr saga')


ax.legend()
ax.set_ylabel('time in minutes')
positions = range(len(xlabels))
plt.xticks(positions, xlabels)
labels = [item.get_text() for item in ax.get_xticklabels()]
#print(labels)

for tick in ax.get_xticklabels():
    tick.set_rotation(90)

fig.tight_layout()

fig.savefig('time_taking.png',bbox_inches='tight')
