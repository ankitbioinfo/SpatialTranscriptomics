

f=open('cellbarcodeid_and_name.dat','r')
ctname={}
for line in f:
    l=line[0:-1].split('\t')
    ctname[l[0]]=l[1][1:]

fname=open('annot_mouseStStAll_correct_monocytes.csv','w')


f=open('annot_mouseStStAll.csv')
uniquehep={}
fname.write(f.readline())
d={}
for line in f:
        l=line.split(',')
        clusterid=int(l[2])
        barcodeid=l[5]
        celltypename=l[3]

        try:
            myname=ctname[l[5]]
            myname=myname.replace('Capsular Fibroblasts#', '')
            myname=myname.replace('Fibroblasts#', '')
            myname=myname.replace('Mesothelial cells#', '')
            myname=myname.replace('Stellate cells#', '')
            myname=myname.replace('Central Vein Endothelial cells', 'Central Vein EC')
            myname=myname.replace('Lymphatic Endothelial cells', 'Lymphatic EC')
            myname=myname.replace('Portain Vein Endothelial cells', 'Portain Vein EC')
            myname=myname.replace('CD8 Effector Memory T cells', 'CD8 Eff Mem T')
            myname=myname.replace('Hepatocytes NucSeq', 'Hepatocytes')
            myname=myname.replace('Patrolling Monocytes','Monocytes')
            myname=myname.replace('Trans. Monocytes 2','Monocytes')
            myname=myname.replace('Trans. Monocytes','Monocytes')
            myname=myname.replace('Fibroblast 2','Fibroblast')
            myname=myname.replace('Fibroblast 1','Fibroblast')
            myname=myname.replace('Capsule fibroblasts','Fibroblast')
            myname=myname.replace('MoMac1','Macro')
            myname=myname.replace('MoMac2','Macro')
            myname=myname.replace('Peritoneal Macrophages','Macro')
            myname=myname.replace('Mesothelial cells','Fibroblast')


        except KeyError:
            myname=celltypename+''
            myname=myname.replace('Monocytes & Monocyte-derived cells', 'Monocytes')
            myname=myname.replace('Kupffer cells', 'KCs')


        if myname=="Hepatocytes":
            if l[2] not in uniquehep:
                uniquehep[l[2]]=1
            else:
                uniquehep[l[2]]+=1
            #print(l[2])
        if (int(l[2])==3)|(int(l[2])==6)|(int(l[2])==7)|(int(l[2])==8):
                myname='Hep678'
        if (int(l[2])==1)|(int(l[2])==2)|(int(l[2])==4)|(int(l[2])==5):
                myname='Hep12'

        if myname=="Hepatocytes":
            myname='Hep12'        

        d[myname]=1

        fname.write(l[0]+','+l[1]+','+l[2]+','+myname+','+l[4]+','+l[5]+','+l[6]+','+l[7])

print(uniquehep)
a=sorted(list(d.keys()))
for i in range(len(a)):
    print(i,a[i])




'''
cname=sorted(list(d.keys()))

for key in cname:
    a=d[key]
    t={}
    for i in range(len(a)):
        #print(1,i,a[i])
        name=ctname[a[i]]
        if name not in t:
            t[name]=1
        else:
            t[name]+=1

    print('\n',key,t)
'''
