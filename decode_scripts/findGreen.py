
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

f=open('rsfish_allspots.csv')
'''
0: R1-G
1: R1-R
2: R1-FR
3: R2-G
4: R2-R
5: R2-FR
6: R3-G
7: R3-R
8: R3-FR
9: R4-G
10: R4-R
11: R4-FR
12: R5-G
13: R5-R
14: R5-FR
15: R6-G
16: R6-R
17: R6-FR

209126 total spots from all 18 channels
1015 is total spots have green channel and at least 3 spots
and from these unique 3 green channels are 711

'''

#Dnm3os=['H1.4',	'H2.4',	'H3.1']
#Igfbp5=['H1.1',	'H2.4','H3.4']
#Lrrc55=['H1.4',	'H2.4',	'H3.4']
#Olfr78=['H1.1',	'H2.1',	'H3.1']
#Osmr=['H1.4','H2.1','H3.4']
#Plin1=['H1.4','H2.1','H3.1']

roundid={}
roundid['H1.1']=0
roundid['H1.4']=3
roundid['H2.1']=6
roundid['H2.4']=9
roundid['H3.1']=12
roundid['H3.4']=15

Dnm3os=['H1.4',	'H2.4',	'H3.1']
Igfbp5=['H1.1',	'H2.4','H3.4']
Lrrc55=['H1.4',	'H2.4',	'H3.4']
Olfr78=['H1.1',	'H2.1',	'H3.1']
Osmr=['H1.4','H2.1','H3.4']
Plin1=['H1.4','H2.1','H3.1']

green=[0,3,6,9,12,15]

ggenes=[Dnm3os,Igfbp5,Lrrc55,Olfr78,Osmr,Plin1]

chid=[]
for i in range(len(ggenes)):
    l=ggenes[i]
    t=[]
    for k in range(len(l)):
        t.append(roundid[l[k]])
    chid.append(t)

print(chid)

d={}
for k in range(len(chid)):
    d[k]=0

samespot={}
fc=0
d2={}
fw=open('rsfish_green.csv','w')
unique_row={}
for line in f:
    l=line[0:-1].split(',')
    unique_row[line]=1
    a=l[0].split('-')
    name=l[0]
    if name not in d2:
        d2[name]=1
    else:
        d2[name]+=1
    b=list(set(a[1:]))

    count=0
    for j in range(len(b)):
        b[j]=int(b[j])
        if b[j] in green:
            count+=1
    #if (count==len(b))&(len(b)==3):
    if (count==len(b))&(len(b)==3)&(len(l)==4):
        nameOfSpots=str(sorted(l[1:]))
        if nameOfSpots not in samespot:
            samespot[nameOfSpots]=1
            for k in range(len(chid)):
                x=np.array(sorted(b))
                y=np.array(sorted(chid[k]))
                value=np.array_equal(x,y)
                if value:
                    d[k]+=1
            #print(x,y,value)
        #print(a,b,len(l))
        fw.write(line)
        fc+=1


print(len(samespot),fc)
print('all keys',len(d2),sum(d2.values()))
print('unique rows',len(unique_row))


ggenes=['Dnm3os','Igfbp5','Lrrc55','Olfr78','Osmr','Plin1']
values=[30481,   55372  ,   564 ,   3773 ,  10091  ,  2014 ]
x=[]
for k in range(len(ggenes)):
    print(ggenes[k],d[k],np.log10(d[k]))
    x.append(np.log10(d[k]))

print(np.log10(values))

y=np.log10(values)
corr,_ = pearsonr(x,y)
plt.plot(x,y,'o')
plt.xlabel('seqFISH')
plt.ylabel('pseudobulk')
plt.title('corr='+str(corr))
plt.savefig('correlation.png')
