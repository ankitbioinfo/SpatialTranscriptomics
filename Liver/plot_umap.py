
import numpy as np
import matplotlib.pyplot as plt



f=open('annot_mouseStStAll.csv','r')
f.readline()

fw=open('sc_cluster.dat','w')
fw.write('cell,cluster\n')
d={}
dn={}
count=0
ctname={}
for line in f:
    l=line.split(',')
    fw.write('cell'+str(count)+','+l[2]+'\n')
    count+=1
    cluster=l[2]

    barcodename=l[5]


    ctname[barcodename]=1
    if cluster in d:
        d[cluster].append([float(l[0]),float(l[1])])
    else:
        d[cluster]=[[float(l[0]),float(l[1])]]
        dn[cluster]=l[3]

print(len(d))
fw=open('NameOfCT.dat','w')
for key in d:
    fw.write(key+'\t'+dn[key]+'\t'+str(len(d[key]))+'\n')
    a=np.array(d[key])
    plt.plot(a[:,0],a[:,1],'.')

plt.savefig('ankit.png')
#plt.show()


f=open('countTable_mouseStSt/barcodes.tsv')
count=0
notpresent=[]
for line in f:
    l=line[0:-1]
    if l in ctname:
        count=count+1
    else:
        print(l)

print('barcode',len(ctname),count)
