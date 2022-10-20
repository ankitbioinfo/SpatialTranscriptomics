


f=open('gencode/gencode.v41.annotation.gff3')
for i in range(6):
    f.readline()

category= ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene', 'protein_coding']
#category= ['protein_coding']




fw=open('gencode.v41.annotation_transcripts.gff3','w')

d={}
for line in f:
    l=line.split('\t')
    if len(l)>5:
        #if l[2]=='exon':
        #if (l[2]=='five_prime_UTR')|(l[2]=='three_prime_UTR'):
        if l[2]=='transcript':
        #if l[2]=='CDS':
            flag=0
            for j in range(len(category)):
                if line.find(category[j])!=-1:
                    flag=1


            p=line.split('gene_name=')
            q=p[1].split(';')
            gname=q[0]

            unique_name=l[0]+'#'+l[3]+'#'+l[4]+'#'+l[6]

            #start_pos_in_bed_format=str(int(l[3])-1)
            #data=l[0][3:]+'\t'+start_pos_in_bed_format+'\t'+l[4]+'\t'+gname +'\t'+l[5]+'\t'+l[6]+'\t.'

            if unique_name not in d:
                d[unique_name]=1
                #if flag==1:
                #    fw.write(line)
            else:
                d[unique_name]+=1

            if flag==1:
                fw.write(line)




        name=l[2]
        #d[name]=1

#print(d)
