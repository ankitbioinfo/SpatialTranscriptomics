

f=open('Introns_file.tsv')
fw=open('intron.bed','w')

extract_gene=['Sox9','Sox17']

for line in f:
    l=line.split()
    if len(l)>3:
            for i in range(len(extract_gene)):
                if line.find(extract_gene[i])!=-1:
                    start_pos_in_bed_format=str(int(l[1])-1)
                    name=l[0][3:]+'\t'+start_pos_in_bed_format+'\t'+l[2]+'\t'+l[4]+'\t'+'0\t'+l[3]+'\t.\n'
                    fw.write(name)
