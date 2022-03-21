
import pandas as pd 
import numpy as np 

df=pd.read_excel('Cardiomyocytes_readout_numbering_for_genes.xlsx')
data=df.to_numpy()

genename=data[:,0]
probename=np.unique(data[:,1:])

print(probename)


'''
experimental design 
order of channels in the tif file  
1   C640    H1.3
2   C405    DAPI
3   C488    H1.1
4   C561    H1.2


I will rename the tif files after taking the maximum projection
into following order 
python      Experiment 
Round 1     H1 R1   
Round 2     H1 R2 
Round 3     H2 R1 
Round 4     H2 R2 
Round 5     H3 R1 
Round 6     H3 R2 

It is convinient to put DAPI in the 4th channel 
So change the order of channel 
1   C488    H1.1    old3
2   C561    H1.2    old4
3   C640    H1.3    old1
4   C405    DAPI    old2
in following command in format_iss*.py file 
#ch_dict = {0: 'FITC', 1: 'Cy3', 2: 'Cy3 5', 3: 'Cy5',4: }
ch_dict = {0:'C2_FOV1',1:'C3_FOV1',2:'C4_FOV1',3:'C1_FOV1'}

'''





key_ofProbesRound={}

#Round and channel set manually for the experiment 
#python numbering start from 0 so [1,2,3] channel becomes [0,1,2]
key_ofProbesRound['H1.1']=[0,0]
key_ofProbesRound['H1.2']=[0,1]
key_ofProbesRound['H1.3']=[0,2]
key_ofProbesRound['H1.4']=[1,0]
key_ofProbesRound['H1.5']=[1,1]
key_ofProbesRound['H1.6']=[1,2]

key_ofProbesRound['H2.1']=[2,0]
key_ofProbesRound['H2.2']=[2,1]
key_ofProbesRound['H2.3']=[2,2]
key_ofProbesRound['H2.4']=[3,0]
key_ofProbesRound['H2.5']=[3,1]
key_ofProbesRound['H2.6']=[3,2]

key_ofProbesRound['H3.1']=[4,0]
key_ofProbesRound['H3.2']=[4,1]
key_ofProbesRound['H3.3']=[4,2]
key_ofProbesRound['H3.4']=[5,0]
key_ofProbesRound['H3.5']=[5,1]
key_ofProbesRound['H3.6']=[5,2]


fw=open('codebook.json','w')

fw.write('{\n')
fw.write('\t"version": "0.0.0",\n')
fw.write('\t"mappings": [\n')


count=0
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
            #if a1==a2==a3:
            #    pass
            #else:
                #print(a1[0])
                count+=1
                fw.write('\t\t{\n')
                fw.write('\t\t\t"codeword": [\n')
                fw.write('\t\t\t\t{"r": '+str(a1[0])+', "c": '+str(a1[1])+', "v": 1},\n'    )
                fw.write('\t\t\t\t{"r": '+str(a2[0])+', "c": '+str(a2[1])+', "v": 1},\n'    )
                fw.write('\t\t\t\t{"r": '+str(a3[0])+', "c": '+str(a3[1])+', "v": 1}\n'    )
                fw.write('\t\t\t],\n')
                fw.write('\t\t\t"target": "'+gene+'"\n')
                if (i+1)==len(genename):
                    fw.write('\t\t}\n')
                else:
                    fw.write('\t\t},\n')

fw.write('\t]\n')
fw.write('}\n')
print(count)
