


#f=open('gencode.v41.annotation_exon.gff3')
#f=open('gencode.v41.annotation_transcripts.gff3')

#f=open('gencode.v41.annotation_CDS.gff3')
#f=open('Final_merge_transcripts.dat')
f=open('Final_merge_exon2.dat')


fw=open('create_bed_files_for_good_merged_genes_exons.bed','w')

d={}
unique_d={}

total=0

pos=[]
for line in f:
    l=line.split('\t')

    p=line.split('gene_name=')
    temp=[]
    for j in range(1,len(p)):
        q=p[j].split(';')
        gname=q[0]
        if gname not in temp:
            temp.append(gname)


    for j in range(len(temp)):
        d[temp[j]]=1
    if len(temp)==1:
        unique_d[temp[0]]=1

    total+=int(l[2])-int(l[1])

    #merge coordiantes must have unique genes
    if len(temp)==1:
        gname=temp[0]
        strand=list(set(l[4].split(',')))[0]
        mix=gname+'#'+l[3] # genename and their isoform
        data=l[0]+'\t'+l[1]+'\t'+l[2]+'\t'+mix+'\t'+l[3]+'\t'+strand+'\t.'+'\n'
        #data=l[0][3:]+'\t'+l[1]+'\t'+l[2]+'\t'+l[3]+'\t'+gname+'\t'+strand+'\t.'+'\n'
        fw.write(data)

    #if l[0]=='chr1':
        #start=int(l[3])
        #end=int(l[4])
        #if gname=='SAMD11':
        #    pos.append([start,end,gname])

fw.close()
print("unique and all",len(unique_d),len(d),total)
'''
gname=sorted(list(d.keys()))
fw=open('mygene.dat','w')
for i in range(len(gname)):
    fw.write(gname[i]+'\n')


for i in range(len(pos)):
    for j in range(i+1,len(pos)):
         if (pos[j][0]<=pos[i][0]<=pos[j][1])  | (pos[j][0]<=pos[i][1]<=pos[j][1]) :
             print(pos[i],pos[j])
'''
