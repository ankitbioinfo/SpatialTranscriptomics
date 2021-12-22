

f=open('gencode.vM28.annotation.gff3')
fw=open('exon.bed','w')
extract_gene=['Sox9','Sox17']

for line in f:
    l=line.split()
    if len(l)>3:
        if l[2]=='exon':
            for i in range(len(extract_gene)):
                if line.find(extract_gene[i])!=-1:

                    t=l[8].split(';')
                    d={}
                    for j in range(len(t)):
                        a=t[j].split('=')
                        d[a[0]]=a[1]


                    start_pos_in_bed_format=str(int(l[3])-1)
                    name=l[0][3:]+'\t'+start_pos_in_bed_format+'\t'+l[4]+'\t'+d['transcript_id']+'\t'+'0\t'+l[6]+'\t.'

                    #for i in range(8,len(l)):
                    #    name+='\t'+l[i]
                    name+='\n'
                    fw.write(name)
