


key_ofProbesRound={}

key_ofProbesRound['H1P4']=0
key_ofProbesRound['H2P9']=1
key_ofProbesRound['H3P4']=2
key_ofProbesRound['H3P5']=3
key_ofProbesRound['H3P9']=4
key_ofProbesRound['H3P11']=5


f=open('andy_probe_code.dat','r')
fw=open('codebook.json','w')

fw.write('{\n')
fw.write('\t"version": "0.0.0",\n')
fw.write('\t"mappings": [\n')

cont=f.readlines()
count=0
for i in range(len(cont)):
    l=cont[i].split()
    gene=l[1]
    flag=1
    try:
        a1=str(key_ofProbesRound[l[2]])
        a2=str(key_ofProbesRound[l[3]])
        a3=str(key_ofProbesRound[l[4]])

    except KeyError:
        flag=0

    if flag==1:
            if a1==a2==a3:
                pass
            else:
                count+=1
                fw.write('\t\t{\n')
                fw.write('\t\t\t"codeword": [\n')
                fw.write('\t\t\t\t{"c": 0, "r": '+a1+', "v": 1},\n'    )
                fw.write('\t\t\t\t{"c": 0, "r": '+a2+', "v": 1},\n'    )
                fw.write('\t\t\t\t{"c": 0, "r": '+a3+', "v": 1}\n'    )
                fw.write('\t\t\t],\n')
                fw.write('\t\t\t"target": "'+gene+'"\n')
                fw.write('\t\t},\n')

fw.write('\t]\n')
fw.write('}\n')
print(count)
