
#features,r,c,intensity,z,y,x,radius,spot_id,z_min,z_max,y_min,y_max,x_min,x_max,IntensityTable
#zmin is 9
#zmax is 10

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

from scipy.stats import pearsonr
import cv2 as cv
import os
import math


def makedir(maindir):
    answer=os.path.isdir(maindir)
    if answer==True:
	    pass
    else:
	    os.mkdir(maindir)

def euclidean_dist(p1,p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2  + (p1[2]-p2[2])**2)

'''
def findNeighbors_in_given_radius(location,radius):
    n=location.shape[0]
    #print('ss',location.shape)
    neighbor={}
    for i in range(n):
        loc1=location[i]
        #print(loc1)
        t1=(loc1[0]-1.1*radius) <= location[:,0]
        t2=location[:,0] <= (loc1[0]+1.1*radius)
        t3=(loc1[1]-1.1*radius) <= location[:,1]
        t4=location[:,1] <= (loc1[1]+1.1*radius)
        t5=(loc1[2]-1.1*radius) <= location[:,2]
        t6=location[:,2] <= (loc1[2]+1.1*radius)

        index=  np.where ( t1 & t2 & t3 & t4 & t5 & t6    )

        count=0

        for k in range(len(index[0])):
            j=index[0][k]
            if j!=i:
                count+=1
                loc2=location[j]
                dist=euclidean_dist(loc1,loc2)
                if dist<radius:
                    if i not in neighbor:
                        neighbor[i]=[j]
                    else:
                        if j not in neighbor[i]:
                            neighbor[i].append(j)

                    if j not in neighbor:
                        neighbor[j]=[i]
                    else:
                        if i not in neighbor[j]:
                            neighbor[j].append(i)

        #print('t',count,len(index[0]))


    newneig=[]
    avg_neigh=0.0
    for i in range(n):
        try:
            l=neighbor[i]
        except KeyError:
            l=[]
        #print(l)
        newneig.append(l)
        avg_neigh+=len(l)

    #print('average neighbors:',avg_neigh/n)

    return newneig
'''


def read_data(fname,d,RCindex):
    df=pd.read_csv(fname)
    #print(df['x_max'])
    #,z_max,y_min,y_max,x_min,x_max,IntensityTable]

    data=df.to_numpy()
    #print(fname,data.shape)
    for i in range(data.shape[0]):
        #print(data[i])
        intensity=data[i,15]
        name=''
        for j in range(9,15):
            name+=str(float(data[i,j]))+'-'
        #print(data[i,9:],name,intensity)
        if (len(name.split('-'))!=7):
            print(fname,name)
        if name not in d:
            #d[name]=np.zeros((18))
            #d[name][RCindex]=intensity
            d[name]=[1]
        else:
            #d[name][RCindex]=intensity
            d[name][0]+=1

    return d
    #print(len(d))


def readdata(rounds,channels,path):
    d={}
    RCindex=0
    #for R in range(6):
    for R in rounds:
        for C in channels:
            fname=path+'/IntensityTable_R'+str(R)+'C'+str(C)+'.dat'
            d=read_data(fname,d,RCindex)
            RCindex+=1
    #print(len(d))
    return d

def read_codebook_file():
    f=open('inputCodebook/retained_intron_count.txt')
    finalgeneinAndyPool=[]
    for line in f:
        l=line.split()
        finalgeneinAndyPool.append(l[0])

    #print(sorted(finalgeneinAndyPool))

    df=pd.read_excel('inputCodebook/Cardiomyocytes_readout_numbering_for_genes.xlsx')
    data=df.to_numpy()

    genename=data[:,0]
    probename=np.unique(data[:,1:])

    #print(probename)

    key_ofProbesRound={}

    #Round and channel set manually for the experiment
    #python numbering start from 0 so [1,2,3] channel becomes [0,1,2]
    key_ofProbesRound['H1.1']=[0]
    key_ofProbesRound['H1.2']=[1]
    key_ofProbesRound['H1.3']=[2]
    key_ofProbesRound['H1.4']=[3]
    key_ofProbesRound['H1.5']=[4]
    key_ofProbesRound['H1.6']=[5]

    key_ofProbesRound['H2.1']=[6]
    key_ofProbesRound['H2.2']=[7]
    key_ofProbesRound['H2.3']=[8]
    key_ofProbesRound['H2.4']=[9]
    key_ofProbesRound['H2.5']=[10]
    key_ofProbesRound['H2.6']=[11]

    key_ofProbesRound['H3.1']=[12]
    key_ofProbesRound['H3.2']=[13]
    key_ofProbesRound['H3.3']=[14]
    key_ofProbesRound['H3.4']=[15]
    key_ofProbesRound['H3.5']=[16]
    key_ofProbesRound['H3.6']=[17]


    #fw=open('inputCodebook/codebook.csv','w')
    #fw.write('gene,I,II,III\n')
    codebook={}
    for i in range(len(genename)):
        gene=genename[i]
        flag=1
        try:
            a1=key_ofProbesRound[data[i,1]]
            a2=key_ofProbesRound[data[i,2]]
            a3=key_ofProbesRound[data[i,3]]

        except KeyError:
            flag=0
        if flag==1:
            #print(a1,a2,a3)
            if gene in finalgeneinAndyPool:
                #fw.write(gene)
                t=[a1[0],a2[0],a3[0]]
                name=''
                for j in range(3):
                    #fw.write(','+str(t[j]))
                    name+='-'+str(t[j])
                #fw.write('\n')
                #print(name)
                codebook[name]=gene
    return codebook,finalgeneinAndyPool


'''
def find_nearest_distance(a1,a2,a3,fname,myradius):
    fw=open(fname,'w')
    a1=list(a1)
    a2=list(a2)
    a3=list(a3)

    location=np.zeros((len(a1)+len(a2)+len(a3),3),dtype=float)


    mid=[]
    all=[]
    for i in range(len(a1)):
        l1=a1[i].split('-')
        all.append(a1[i])
        z=[float(l1[0]),float(l1[1])]
        y=[float(l1[2]),float(l1[3])]
        x=[float(l1[4]),float(l1[5])]
        location[i]=[np.mean(x), np.mean(y), np.mean(z)]
        mid.append(1)

    posid=len(a1)

    for j in range(len(a2)):
        all.append(a2[j])
        l1=a2[j].split('-')
        z=[float(l1[0]),float(l1[1])]
        y=[float(l1[2]),float(l1[3])]
        x=[float(l1[4]),float(l1[5])]
        location[posid+j]=[np.mean(x), np.mean(y), np.mean(z)]
        mid.append(2)

    posid=len(a1)+len(a2)

    for j in range(len(a3)):
        all.append(a3[j])
        l1=a3[j].split('-')
        z=[float(l1[0]),float(l1[1])]
        y=[float(l1[2]),float(l1[3])]
        x=[float(l1[4]),float(l1[5])]
        location[posid+j]=[np.mean(x), np.mean(y), np.mean(z)]
        mid.append(3)


    mid=np.array(mid)

    print('ankit',len(a1)+len(a2)+len(a3),mid.shape)

    neigh=findNeighbors_in_given_radius(location,myradius)

    count=0
    for i in range(len(mid)):
        if mid[i]==1:
            l=neigh[i]
            flag2=0
            flag3=0
            for j in range(len(l)):
                if mid[l[j]]==2:
                    flag2=1
                    neigh1=all[l[j]]
                if mid[l[j]]==3:
                    flag3=1
                    neigh2=all[l[j]]
            if ((flag2==1)&(flag3==1)):
                count+=1
                fw.write(all[i]+'\t'+neigh1+'\t'+neigh2+'\n')
                #print(mid[i],mid[l[j]])

    print(mid.shape,location.shape,count)
'''

def remove_localized_spots(ch,neighborRadius):
    totalsum=len(ch[0])
    for i in range(len(ch)):
        ch[i]=sorted(list(ch[i].keys()))
    com=set(ch[0])
    for j in range(1,len(ch)):
        com=com.intersection(set(ch[j]))
        totalsum+=len(ch[j])
    #left over spots
    leftover=[]
    for i in range(len(ch)):
        leftover.append(sorted(list(set(ch[i])-com)))

    spots=[]
    for i in range(len(ch)):
        for j in range(len(leftover[i])):
            l=leftover[i][j].split('-')
            #print(leftover[i][j],l)
            ymin=0.5*(float(l[2])+float(l[3]))
            ywid=float(l[3])-float(l[2])
            xmin=0.5*(float(l[4])+float(l[5]))
            xwid=float(l[5])-float(l[4])

            zmin=0.5*(float(l[1])+float(l[0]))
            zwid=float(l[1])-float(l[0])
            #print(ymin,ymax,xmin,xmax)

            spots.append([ymin,xmin,ywid/2,xwid/2,zmin,zwid/2,i])
            #print(l,[ymin,xmin,xwid/2,ywid/2])
        #print(i,len(ch[i]),len(leftover[i]))

    spots=np.array(spots)
    data=spots[:,[0,1]]
    #print(data,'\n\n')
    '''
    for i in range(len(data)):
        for j in range(i+1,len(data)):
            d=euclidean_dist(data[i],data[j])
            if d<400:
                print(i,j,d)
    '''

    k_index = cKDTree(data).query_ball_point(x=data,r=neighborRadius, workers=-1)

    neighbors=[]
    for i in range(len(k_index)):
        k_index[i].remove(i)
        neighbors.append(k_index[i])

    registerDots=[]
    for i in range(len(neighbors)):
        if len(neighbors)>0:
            l=neighbors[i]
            neighspots1=spots[i]
            chid1=int(neighspots1[-1])
            for j in range(len(l)):
                neighspots2=spots[l[j]]
                chid2=int(neighspots2[-1])
                if (chid1!=chid2):#&(chid1==0):
                    t=  str(int(neighspots1[4]-neighspots1[5]))+'-'+str(int(neighspots1[4]+neighspots1[5]))+'-'
                    t=t+str(int(neighspots1[0]-neighspots1[2]))+'-'+str(int(neighspots1[0]+neighspots1[2]))+'-'
                    t=t+str(int(neighspots1[1]-neighspots1[3]))+'-'+str(int(neighspots1[1]+neighspots1[3]))+'-'
                    #print(i,j,l[j],neighspots1,neighspots2,t)
                    registerDots.append(t)
                    t=  str(int(neighspots2[4]-neighspots2[5]))+'-'+str(int(neighspots2[4]+neighspots2[5]))+'-'
                    t=t+str(int(neighspots2[0]-neighspots2[2]))+'-'+str(int(neighspots2[0]+neighspots2[2]))+'-'
                    t=t+str(int(neighspots2[1]-neighspots2[3]))+'-'+str(int(neighspots2[1]+neighspots2[3]))+'-'
                    #print(i,j,l[j],neighspots1,neighspots2,t)
                    registerDots.append(t)


    total_registered_spots=com.union(set(registerDots))
    #total=list(com)
    print("common",totalsum,len(com),len(total_registered_spots))

    newch=[]
    unique_spots_in_each_image=[]
    for i in range(len(ch)):
        tt=set(ch[i])-total_registered_spots
        newch=newch+list(tt)
        unique_spots_in_each_image.append(list(tt))

    print(len(newch),len(set(newch)))
    combined_unique_spots_from_all_images=set(newch)

    return [unique_spots_in_each_image,combined_unique_spots_from_all_images]

def merge_spots_from_H1H2H3(un_H1H2H3):
    spots_from_all_R_C={}
    count=0
    for i in range(len(un_H1H2H3)):
        for j in range(len(un_H1H2H3[i])):
            #print(i,j,count)
            for k in range(len(un_H1H2H3[i][j])):
                name=un_H1H2H3[i][j][k]
                if name not in spots_from_all_R_C:
                    spots_from_all_R_C[name]=[count]
                else:
                    spots_from_all_R_C[name].append(count)
            count=count+1
    return  spots_from_all_R_C


def find_3_closest_pairs(neigh_spots,channels):
    r1=[0,1,2,3,4,5]
    r2=[6,7,8,9,10,11]
    r3=[12,13,14,15,16,17]

    t1=[]
    t2=[]
    t3=[]
    sp1=[]
    sp2=[]
    sp3=[]
    for i in range(len(channels)):
        for j in range(len(channels[i])):
            for k in range(len(r1)):
                if channels[i][j]==r1[k]:
                    t1.append(r1[k])
                    sp1.append(neigh_spots[i])
            for k in range(len(r2)):
                if channels[i][j]==r2[k]:
                    t2.append(r2[k])
                    sp2.append(neigh_spots[i])
            for k in range(len(r3)):
                if channels[i][j]==r3[k]:
                    t3.append(r3[k])
                    sp3.append(neigh_spots[i])


    chooseBestNeighbors=[]
    config=[]
    best_spots_return=[]
    for i in range(len(t1)):
        for j in range(len(t2)):
            for k in range(len(t3)):
                d1=sp1[i][0:2]
                d2=sp2[j][0:2]
                d3=sp3[k][0:2]
                best_spots_return.append([sp1[i],sp2[j],sp3[k]])
                data=np.array([d1,d2,d3]).astype(float)
                config.append([t1[i],t2[j],t3[k]])
                #print(data)
                #print(data.shape)
                dist = pdist(data)
                dist=np.mean(dist)
                #print(i,j,k,dist)
                chooseBestNeighbors.append(dist)

    #print(chooseBestNeighbors)
    #print(channels,'\t\t',t1,t2,t3)
    #print(sp1,sp2,sp3)


    #if len(chooseBestNeighbors)>1:
    index=np.argmin(chooseBestNeighbors)
    #print(chooseBestNeighbors,index)
    final_channel=config[index]
    final_spots=best_spots_return[index]


    tfinal=[]
    unique=''
    for i in range(len(final_spots)):
        neighspots1=final_spots[i]
        t=  str('%0.2f'%float(neighspots1[4]-neighspots1[5]))+'-'+str('%0.2f'%float(neighspots1[4]+neighspots1[5]))+'-'
        t=t+str('%0.2f'%float(neighspots1[0]-neighspots1[2]))+'-'+str('%0.2f'%float(neighspots1[0]+neighspots1[2]))+'-'
        t=t+str('%0.2f'%float(neighspots1[1]-neighspots1[3]))+'-'+str('%0.2f'%float(neighspots1[1]+neighspots1[3]))+'-'
        tfinal.append(t)
        unique=unique+'#'+t
        #print(t)
    #print('\n\n\n\n')
    return final_channel,final_spots,tfinal,unique


def find_registered_dots_within_distance(ch,spots_with_ids,neighborRadius):

    for i in range(len(ch)):
        ch[i]=sorted(list(ch[i]))

    com=set(ch[0])
    for j in range(1,len(ch)):
        com=com.intersection(set(ch[j]))

    #left over spots
    leftover=[]
    for i in range(len(ch)):
        leftover.append(sorted(list(set(ch[i])-com)))

    spots=[]
    for i in range(len(ch)):
        for j in range(len(leftover[i])):
            l=leftover[i][j].split('-')
            ymin=0.5*(float(l[2])+float(l[3]))
            ywid=float(l[3])-float(l[2])
            xmin=0.5*(float(l[4])+float(l[5]))
            xwid=float(l[5])-float(l[4])

            zmin=0.5*(float(l[1])+float(l[0]))
            zwid=float(l[1])-float(l[0])

            myid=spots_with_ids[leftover[i][j]]
            spots.append([ymin,xmin,ywid/2,xwid/2,zmin,zwid/2,myid,i])
            #print(l,[ymin,xmin,xwid/2,ywid/2])
        #print(i,len(ch[i]),len(leftover[i]))

    spots=np.array(spots)
    data=spots[:,[0,1]]
    fw=open('final_close.dat','w')
    for i in range(len(data)):
        fw.write(str(spots[i,0])+','+str(spots[i,1])+','+str(spots[i,2])+','+str(spots[i,3])+','+str(spots[i,4])+
        ','+str(spots[i,5])+','+str(spots[i,6][0])+','+str(spots[i,7])+'\n')
    #print(data)
    k_index = cKDTree(data).query_ball_point(x=data,r=neighborRadius, workers=-1)
    #distance_1,k_index1 = cKDTree(data).query(x=data, k=3, p=,workers=n_jobs)

    neighbors=[]
    for i in range(len(k_index)):
        k_index[i].remove(i)
        neighbors.append(k_index[i])

    '''
    fw=open('neighbors.dat','w')
    for i in range(len(neighbors)):
        for j in range(len(neighbors[i])):
            fw.write(str(neighbors[i][j])+':'+ str(spots[neighbors[i][j]]) +  ',')
        fw.write('\n')
    '''

    fw=open('rsfish_allspots.csv','w') #each row is not unique here
    registerDots=[]
    registerDots_ids=[]
    unique_registered_spots={}
    for i in range(len(neighbors)):
    #for i in range(12):
        if len(neighbors)>0:
            l=neighbors[i]
            neighspots1=spots[i]
            chid1=int(neighspots1[-1])
            temp_probe_id=[chid1] #chid1 =0 means belong to H1
            temp_channel_id=[neighspots1[6]]
            all_neighbors_spots=neighspots1
            for j in range(len(l)):
                neighspots2=spots[l[j]]
                #print(i,j,neighspots1.shape,neighspots2.shape,all_neighbors_spots.shape)
                all_neighbors_spots=np.vstack((all_neighbors_spots,neighspots2))
                chid2=int(neighspots2[-1])
                temp_channel_id.append(neighspots2[6])
                #print(i,j,chid1,chid2)
                if chid2 not in temp_probe_id:
                    temp_probe_id.append(chid2) #chid2=1,2 means belong to H2 and H3 so H1 agains should not come


            #print(i,tempid,temp)
            temp_probe_id=np.unique(temp_probe_id)
            if (len(temp_probe_id)==3):
                temp_channel_id,final_spots,final_spots_with_str_name,un_spots=find_3_closest_pairs(all_neighbors_spots,temp_channel_id)
                an=''
                for k in range(len(temp_channel_id)):
                    an+='-'+str(temp_channel_id[k])
                fw.write(an)
                for k in range(len(final_spots)):
                    l=final_spots[k]
                    fw.write(',%0.2f'%l[0]+':'+'%0.2f'%l[1]+':'+str(l[6][0]))
                    #print(,temp)
                fw.write('\n')

                if un_spots not in unique_registered_spots:
                    unique_registered_spots[un_spots]=1
                    registerDots.append(final_spots_with_str_name)
                    registerDots_ids.append(an)


    final_registered_spots={}
    com=list(com)
    for i in range(len(com)):
        myid=spots_with_ids[com[i]]
        name=''
        for j in range(len(myid)):
            name+='-'+str(myid[j])
        if name not in final_registered_spots:
            final_registered_spots[name]=[com[i]]
        else:
            final_registered_spots[name].append(com[i])


        #print(i,com[i],myid,name)
    count=0
    for key in final_registered_spots:
        count+= len(final_registered_spots[key])
    print("length1",len(com),len(final_registered_spots),count)
    print('total spots found',len(registerDots),len(registerDots_ids))
    #print(final_registered_spots['-5-10-17'])
    #print(final_registered_spots)

    for i in range(len(registerDots)):
        '''
        p1=registerDots_ids[i][0]
        p2=registerDots_ids[i][1]
        p3=registerDots_ids[i][2]
        r1=[0,1,2,3,4,5]
        r2=[6,7,8,9,10,11]
        r3=[12,13,14,15,16,17]
        print(p1,p2,p3,registerDots_ids[i])
        flag1=list(set(r1).intersection(set(p1)))
        flag2=list(set(r2).intersection(set(p2)))
        flag3=list(set(r3).intersection(set(p3)))
        if ((len(flag1)==1)&(len(flag2)==1)&(len(flag3)==1)):
            name='-'+str(flag1[0])+'-'+str(flag2[0])+'-'+str(flag3[0])
        '''
        name=registerDots_ids[i]
        if name not in final_registered_spots:
            final_registered_spots[name]=[registerDots[i]]
        else:
            final_registered_spots[name].append(registerDots[i])


    count=0
    for key in final_registered_spots:
        count+= len(final_registered_spots[key])
    print("length2",len(com),len(final_registered_spots),"final registered spots",count)

    print("\nlocalized", len(com)+ len(registerDots))

    return final_registered_spots


def find_index(sp_genename,sc_genename):
    index_sc=[]
    index_sp=[]
    d={}
    for j in range(len(sc_genename)):
        name=sc_genename[j]
        d[name]=j

    for i in range(len(sp_genename)):
        name=sp_genename[i]
        try:
            d[name]
            flag=1
        except KeyError:
            flag=0
        if flag==1:
            index_sc.append(d[name])
            index_sp.append(i)
    return index_sp,index_sc



def convert_spots_into_genes(registered_spots):
    codebook,genes=read_codebook_file()
    #for key in codebook:
    #    print(key,codebook[key]) #genes[i]
    print('list of spots and genes',len(registered_spots),len(genes))


    genecounts={}
    spot_coordinates={}
    for i in range(len(genes)):
        genecounts[genes[i]]=0
        spot_coordinates[genes[i]]=[]


    decode=[]
    genename=[]
    for key in codebook:
        genename.append(codebook[key])
        try:
            value=len(registered_spots[key])

            A=registered_spots[key]
            for i in range(len(A)):
                spot_coordinates[codebook[key]].append(A[i][0]) #0 means only save the first channel
        except KeyError:
            value=0
        decode.append(value)

    print('all spots decoded',sum(decode))



    notdecodable=0
    available=0
    unique_available={}
    for key in registered_spots:
        B=registered_spots[key]
        available+=len(B)
        for j in range(len(B)):
            for k in range(len(B[j])):
                if B[j][k] not in unique_available:
                    unique_available[B[j][k]]=1

        if key not in codebook:
            notdecodable+=len(registered_spots[key])

    total=0
    fw=open('decoded_genes.csv','w')
    for i in range(len(genename)):
        fw.write(genename[i]+'\t'+str(decode[i])+'\n')
        total+=decode[i]
    fw.write('Total'+'\t'+str(total)+'\n')

    print("Total decodable",total, 'not decodable', notdecodable,"available spots", available,'unique',len(unique_available))



    maindir='decoded_tif/'
    #datapath=sys.argv[1]
    #maindir='png_'+datapath
    makedir(maindir)
    dim3=2304
    dim1=2304
    Bl=[255,0,0]
    Gr=[0,255,0]
    Re=[0,0,255]
    purple=[240,32,160]#   rgb(160,32,240)
    yellow=[0,255,255]#    #RGB	rgb(255,255,0)
    magenta=[255,0,255]#RGB	rgb(255,0,255)
    cyan=[255,255,0]    #RGB	rgb(0,255,255)
    mycolor=[Re,Gr,Bl,cyan,magenta,yellow]

    mycolor=[Bl,Gr,Re]


    for mygenename in spot_coordinates:
        img = np.ones((dim1,dim3,3), np.uint8)
        #img.fill(255)
        fname=maindir+mygenename
        data=spot_coordinates[mygenename]
        #print(mygenename,len(data))
        for i in range(len(data)):
            l=data[i].split('-')
            #ya=[int(l[11]),int(l[12])]
            #xa=[int(l[13]),int(l[14])]
            ya=[math.floor(float(l[2])),math.ceil(float(l[3]))]
            xa=[math.floor(float(l[4])),math.ceil(float(l[5]))]
            cv.rectangle(img,(xa[0],ya[0]),(xa[1],ya[1]),mycolor[-1],2) #2 means fill 1 means boundary
        fname=fname+'_'+str(len(data))
        #cv.imwrite(fname+".png",img)
        img = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
        cv.imwrite(fname+".tif",img)

    #print(decode,genename)

    return genename,decode


def main():
    path='./slide4_tile1/IntensityTile1/case9/'
    path='./smfish/'

    mylist1=[0,1]
    mylist2=[2,3]
    mylist3=[4,5]

    channels=[0,1,2]
    d1=[]
    for i in mylist1:
        for j in channels:
            d1.append(readdata([i],[j],path))
    d2=[]
    for i in mylist2:
        for j in channels:
            d2.append(readdata([i],[j],path))
    d3=[]
    for i in mylist3:
        for j in channels:
            d3.append(readdata([i],[j],path))



    #a1=set(list(d1.keys()))
    #a2=set(list(d2.keys()))
    #a3=set(list(d3.keys()))

    neighborRadius=0
    searchRadius=5

    #find_nearest_distance(a1,a2,a3,'dist_H1H2H3_'+str(myradius)+'.dat',myradius)
    #print(len(a4.intersection(a3)))
    un_H1,comb_H1=remove_localized_spots(d1,neighborRadius)
    un_H2,comb_H2=remove_localized_spots(d2,neighborRadius)
    un_H3,comb_H3=remove_localized_spots(d3,neighborRadius)

    print('1',len(un_H1),len(comb_H1),len(np.sum(un_H1)))
    print('2',len(un_H2),len(comb_H2),len(np.sum(un_H2)))
    print('3',len(un_H3),len(comb_H3),len(np.sum(un_H3)))


    spots_from_all_R_C=merge_spots_from_H1H2H3([un_H1,un_H2,un_H3])


    count=0
    for key in spots_from_all_R_C:
        #print(key,spots_from_all_R_C[key])
        count=count+len(spots_from_all_R_C[key])
    print(count)

    #print("all",len(spots_from_all_R_C))



    common=find_registered_dots_within_distance([comb_H1,comb_H2,comb_H3],spots_from_all_R_C,searchRadius)
    genename,gcounts=convert_spots_into_genes(common) #decodable_spots,save_spot_coordinates
    seqFISH_count=np.array(gcounts)


    df=pd.read_csv('./Correlation/data2/expected_sc_counts.dat')
    andydata=df.to_numpy()
    gc=[]
    for i in range(len(genename)):
        ind_andy,ind_star=find_index(andydata[:,0],[genename[i]])
        andy1=andydata[ind_andy,1].astype(int)
        gene=andydata[ind_andy,0]

        #print(i,andy1,gene)
        gc.append(sum(andy1))

    sc_count=np.array(gc)

    green_and_red_genename=['Ackr2', 'Acsm2', 'Ankrd1', 'Ccl11', 'Ccl17', 'Cd28', 'Cd4', 'Cd59a', 'Cd86', 'Ctla4', 'Cx3cr1', 'Cxcl12',
    'Cxcl14', 'Dnm3os', 'Efnb2', 'Esm1', 'Ets1', 'Flt1', 'H19', 'H2-Eb1', 'Hbegf', 'Ifitm6', 'Igfbp5', 'Itgax', 'Kcna1',
    'Kit', 'Lcn2', 'Lrrc55', 'Ly6c1', 'Mrc1', 'Myh11', 'Myl2', 'Nfil3', 'Nkx2-5', 'Nppa', 'Olfr78', 'Osmr', 'Pcna', 'Pcsk1',
    'Plin1', 'Postn', 'Pparg', 'Ptprc', 'Rorc', 'Slc5a3', 'Th', 'Thbs4', 'Timd4', 'Uhrf1', 'Uqcrb', 'Vtn', 'Wnt4', 'Xcl1']

    green_gene=['Dnm3os','Igfbp5','Lrrc55','Olfr78','Osmr','Plin1']
    red_gene=['Cd59a','Cxcl14','Ets1','Flt1','Nppa','Pcna','Pparg','Uhrf1']
    farred_gene=['Myh6','Retnla','Itgam','Ccl5','Lcmt2','Plac8']


    #print(seqFISH_count,'\n',sc_count)
    green_and_red_index=[]
    green_index=[]
    red_index=[]
    farred_index=[]
    for i in range(len(genename)):
        if genename[i] in green_and_red_genename:
            green_and_red_index.append(i)
        if genename[i] in green_gene:
            green_index.append(i)
        if genename[i] in red_gene:
            red_index.append(i)
        if genename[i] in farred_gene:
            farred_index.append(i)



    row=1
    col=1
    fig,ax=plt.subplots(row,col,figsize=(8,5))
    x=np.log10(1+seqFISH_count)
    y=np.log10(1+sc_count)
    #with all the genes
    m, b = np.polyfit(x, y, 1)
    X_plot = np.linspace(min(x),max(x),100)
    corr,_ = pearsonr(x,y)
    ax.plot(X_plot, m*X_plot + b, 'b:',label='all genes r=%0.2f'%corr)
    ax.plot(x,y,'b.',markerfacecolor="None")

    #green_and_red_index
    m, b = np.polyfit(x[green_and_red_index], y[green_and_red_index], 1)
    corr,_ = pearsonr(x[green_and_red_index],y[green_and_red_index])
    ax.plot(X_plot, m*X_plot + b, 'k:',label='GR genes r=%0.2f'%corr,markerfacecolor="None")
    ax.plot(x[green_and_red_index],y[green_and_red_index],'k.',markerfacecolor="k")


    #green_index
    m, b = np.polyfit(x[green_index], y[green_index], 1)
    corr,_ = pearsonr(x[green_index],y[green_index])
    ax.plot(X_plot, m*X_plot + b, 'g-',label='Green r=%0.2f'%corr)
    ax.plot(x[green_index],y[green_index],'g*',markerfacecolor="g",markersize=12)

    #RED_index
    m, b = np.polyfit(x[red_index], y[red_index], 1)
    corr,_ = pearsonr(x[red_index],y[red_index])
    ax.plot(X_plot, m*X_plot + b, 'r-',label='Red r=%0.2f'%corr)
    ax.plot(x[red_index],y[red_index],'r*',markerfacecolor="r",markersize=12)

    #FARRED_index
    m, b = np.polyfit(x[farred_index], y[farred_index], 1)
    corr,_ = pearsonr(x[farred_index],y[farred_index])
    ax.plot(X_plot, m*X_plot + b, 'm-',label='FarRed r=%0.2f'%corr)
    ax.plot(x[farred_index],y[farred_index],'m*',markerfacecolor="m",markersize=12)



    ax.legend( prop={'size': 6},loc='upper left')
    ax.set_xlabel('seqFISH log10')
    ax.set_ylabel('pseudobulk gene exp log10')


    for j in range(len(genename)):
        if genename[j] in green_gene:
            ax.text(x[j],y[j],genename[j],fontsize=8,color='g')
            print('green','%0.2f'%x[j],'%0.2f'%y[j],genename[j], seqFISH_count[j],sc_count[j])
        elif genename[j] in red_gene:
            ax.text(x[j],y[j],genename[j],fontsize=8,color='r')
            #print('red','%0.2f'%x[j],'%0.2f'%y[j],genename[j])
        elif genename[j] in farred_gene:
            ax.text(x[j],y[j],genename[j],fontsize=8,color='m')
        elif genename[j] in green_and_red_genename:
            ax.text(x[j],y[j],genename[j],fontsize=7,color='k')
        else:
            ax.text(x[j],y[j],genename[j],fontsize=4,color='b')



    fig.tight_layout()
    l=path.split('/')
    print(l)
    fig.savefig('new_correlation_'+str(searchRadius)+'.png',bbox_inches='tight',dpi=300)







main()
