import pandas as pd 


#fname1='IntensityTable_R0C0.dat'
#fname2='smFISH localizations_R1_z3_G.csv'

fname1='IntensityTable_R0C2.dat'
fname2='smFISH localizations_R1_z3_FR.csv'
df=pd.read_csv(fname1)
spots=df.to_numpy()


df=pd.read_csv(fname2)#,header=None
andy=df.to_numpy()


print(spots.shape,andy.shape)



seqfish=spots[:,[11,12,13,14]]
rsfish=andy[:,[0,1]]

#y_min,y_max,x_min,x_max,
#print(seqfish[0:5])
#print(rsfish[0:5])


count=0
fw=open('common.dat','w')
for i in range(len(rsfish)):
    x=rsfish[i,0]
    y=rsfish[i,1]
    flag=0
    for j in range(len(seqfish)):
        d=seqfish[j]
        if (d[0]<= y <=d[1]) & (d[2]<= x <=d[3]):
             flag=1
             fw.write(str(i+1)+':'+str(andy[i])+'\t\t')
             fw.write(str(j+1)+':'+str(seqfish[j])+'\n')

    if flag==1:
        count+=1

print("common",count)
