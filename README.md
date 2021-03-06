ad_spatial.write_h5ad('saveall')

adata=sc.read_h5ad("saveall")
sc.pl.umap(adata, color=["leiden","louvain"], wspace=0.4,show=True, save='_spatial_leiden_louvain.png')

https://github.com/alexcwsmith/singleCellTools/blob/master/ACWS_scanPy_MASTER.py

# SpatialTranscriptomics


1)conda create -n starfish "python=3.7"
2)conda activate starfish
3)pip install scikit-image==0.15.0
4)pip install napari 




To create the input images for transcriptome analysis. 

Run appropriate pipeline either projection one or direct 3d tif files. 

(1) Already projected one 
python37 format_iss_Andy.py fov1 Andy_output 

(2) 3d tif files 
python37 format_seqFISH_Andy.py --input-dir fov1 --output-dir Andy_output --codebook-csv 2017-08-23-10k-gene-barcodes\ round\ I\ to\ V.csv

(3) Run codebook and remove the comma from the output files and then copy into inside folder. 

(4) Steps to select the small region in image for starfish.<br/> 
    Suppose the dimension of the image is 2056 X 1648<br/>
    open macro record<br/>
    Select first half region and it gives [0,0,960, 1648]</br>
    Equivalent starfish command is y_slice = slice(0, 1648), x_first = slice(0, 960)</br>
    Select second half regions and it gives [1080 0 976 1648] [xmin ymin xmax ymax]<br/>
    Equivalent starfish command is x_second = slice(1080, 2056)<br/>
    xmin and xmax is the top left corner of the selected region 


(*) Interesting read 
https://medium.com/apprentice-journal/evaluating-multi-class-classifiers-12b2946e755b

(*) Scan py 
https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html

https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_04_clustering.html

https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html


(*) Mouse cell atlas 
http://bis.zju.edu.cn/MCA/gallery.html?tissue=Bone-Marrow

(*) https://www.proteinatlas.org/ENSG00000072952-MRVI1/celltype/liver


(*) There is considerable functional overlap and interplay between megakaryocytes and endothelial cells. The ultimate function of platelets is to repair disrupted endothelium and ???plug??? up minute holes. This occurs via adhesion to exposed subendothelium structures, activation, aggregation, cell flattening, and activation of angiogenesis. Both platelets and endothelial cells utilize prostaglandin signaling pathways, and modulate hemostasis and thrombosis. Both megakaryocytes and endothelial cells synthesize and secrete von Willebrand factor (vWF), which is involved in linking platelets to exposed basement membrane, and P-selectin, which acts as a key adhesion molecule during hemostasis. Activated platelets also secrete a large number of vasoactive and angiogenic modulatory factors.
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2741141/

(*) Single cell RNA sequencing of human liver reveals
distinct intrahepatic macrophage populations

(*) check the package version 
packageVersion("RaceID")

(*) Cell type markers of HSC <br/>
  p1) Immunofluorescence identifies distinct subsets of endothelial cells in the human liver<br/> 
     Type 1 LSEC are CD36hi CD32??? CD14??? LYVE-1??? located in acinar zone 1 of the lobule <br/> 
     Type 2 LSEC are LYVE-1+ CD32hi CD14+ CD54+ CD36mid-lo located in acinar zones 2 and 3 of the lobule <br/> 
     pericyte marker CD146 <br/>
  p2) Prominent Receptors of Liver Sinusoidal Endothelial Cells in Liver Homeostasis and Disease <br/>
     LSECs assist in clearing macromolecular waste (extracellular matrix material and foreign molecules) from the blood and regulate hepatic vascularity. Individual LSEC???s are flat and very small in size, no thicker than 5 ??m at the center and 0.3 ??m at the periphery. <br/>
     LSECs that differentiate them from other endothelial cells is their higher endocytic ability. LSECs only make up about 3% of total liver volume, however, they contribute to about 45% of pinocytic vesicles in the liver
  
  
      A) Parenchymal cells (60-80%) 
          
          
      B) non-parechymal cells (20-40%) 
         LSEC (50%)
         Kupffer cells (20%) 
         stellate cells (1%) 
         lymphocytes (25%) 
         biliary cells (5%)  

