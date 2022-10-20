

sortBed -i gencode.v41.annotation_exon.gff3 > sort.gff 

#bedtools merge -i sort.gff -S + > output.bed 

bedtools merge -i sort.gff -s -c 4,7,9 -o count,collapse,collapse > Final_merge_.dat

rm sort.gff 
 