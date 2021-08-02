# SpatialTranscriptomics

To create the input images for transcriptome analysis. 

Run appropriate pipeline either projection one or direct 3d tif files. 

(1) Already projected one 
python37 format_iss_Andy.py fov1 Andy_output 

(2) 3d tif files 
python37 format_seqFISH_Andy.py --input-dir fov1 --output-dir Andy_output --codebook-csv 2017-08-23-10k-gene-barcodes\ round\ I\ to\ V.csv

(3) Run codebook and remove the comma from the output files and then copy into inside folder. 


Interesting read 
https://medium.com/apprentice-journal/evaluating-multi-class-classifiers-12b2946e755b

(4) Scan py 
https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html

https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_04_clustering.html

https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html


(5) Mouse cell atlas 
http://bis.zju.edu.cn/MCA/gallery.html?tissue=Bone-Marrow

(6) https://www.proteinatlas.org/ENSG00000072952-MRVI1/celltype/liver


(7) There is considerable functional overlap and interplay between megakaryocytes and endothelial cells. The ultimate function of platelets is to repair disrupted endothelium and “plug” up minute holes. This occurs via adhesion to exposed subendothelium structures, activation, aggregation, cell flattening, and activation of angiogenesis. Both platelets and endothelial cells utilize prostaglandin signaling pathways, and modulate hemostasis and thrombosis. Both megakaryocytes and endothelial cells synthesize and secrete von Willebrand factor (vWF), which is involved in linking platelets to exposed basement membrane, and P-selectin, which acts as a key adhesion molecule during hemostasis. Activated platelets also secrete a large number of vasoactive and angiogenic modulatory factors.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2741141/

(8) Single cell RNA sequencing of human liver reveals
distinct intrahepatic macrophage populations

(9) check the package version 
packageVersion("RaceID")
