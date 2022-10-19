
f=open('gencode/gencode.v41.annotation.gff3')
for i in range(6):
    f.readline()

category= ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_V_gene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene', 'protein_coding']


fw=open('gencode.v41.annotation_exon_9cat.dat','w')

d={}
for line in f:
    l=line.split('\t')
    if len(l)>5:
        if l[2]=='exon':
        #if (l[2]=='five_prime_UTR')|(l[2]=='three_prime_UTR'):
        #if l[2]=='transcript':
        #if l[2]=='gene':
            flag=0
            for j in range(len(category)):
                if line.find(category[j])!=-1:
                    flag=1

            if flag==1:
                fw.write(line)

        name=l[2]
        d[name]=1

print(d)
