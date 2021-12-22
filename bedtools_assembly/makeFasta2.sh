

bedtools getfasta -fi Mus_musculus.GRCm39.dna.primary_assembly.fa -bed exon.dat -s -name -fo exon.fa
bedtools getfasta -fi Mus_musculus.GRCm39.dna.primary_assembly.fa -bed intron.dat -s -name -fo intron.fa
