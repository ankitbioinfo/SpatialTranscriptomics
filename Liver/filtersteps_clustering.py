
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


maindir='figures'
answer=os.path.isdir(maindir)
if answer==True:
    pass
else:
    os.mkdir(maindir)


def reading_clustering():
	f4=open('figures1/leiden_output.dat')
	cont=f4.readlines()
	celltype={}
	noct=[]
	cellsinCT={}
	total=len(cont)-1
	authorClustering=[]
	for i in range(1,len(cont)):
	    l=cont[i].split(',')
	    id=int(l[1])
	    authorClustering.append(id)
	    celltype[l[0]]=id
	    if id not in noct:
	        noct.append(id)

	    if id not in cellsinCT:
	        cellsinCT[id]=[l[0]]
	    else:
	        cellsinCT[id].append(l[0])

	#print('no of cell types',len(noct))
	return np.array(authorClustering)

#authorClustering=reading_clustering()
#print(authorClustering)




filename='vizgen'
#csvfilename='Blank_genes_removed.csv'

csvfilename='Blank_genes_removed.csv'

adata=sc.read(csvfilename).transpose()
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#adata.obs_names_make_unique()

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_04_clustering.html
#new_cluster_names = ['AEC', 'SEC','MK', 'Hepatocyte','Macrophage', 'Myeloid','Erythroid progenitor', 'Erythroid cell']
#adata.obs['authorClustering'] = authorClustering.astype(str)
print('original data', adata)

'''
f=open('geneOrdername140.dat','r')
marker_genes=[]
for line in f:
    marker_genes.append(line[0:-1])


remove=['Dll1','Dll4','Dkk3','Dkk1','Ctnnal1','Clec14a','Celsr2','Cd48','Bmp7','B4galt6','Ammecr1','Egfr','Egr1',
'Elk3','Eng','Epcam','Lef1','Lepr','Mecom','Meis1','Meis2','Mertk','Jag1','Jag2',
'Hgf','Hoxb4','Icam1','Igf1','Fzd2','Fzd3','Fzd4','Fzd7','Fzd8','Gata3','Gfap','Fgf1',
'Fgf2','Dkk2','Cxadr','Cspg4','Wnt2','Tox','Vangl2','Vav1','Tmem56','Tgfb2','Tet1','Tet2',
'Tek','Tcf7','Tcf7l1','Tcf7l2','Selp','Sfrp1','Sfrp2','Sgms2','Slamf1','Satb1','Rassf4',
'Procr','Prickle2','Pdpn','Pou2af1','Notch3','Notch4','Olr1','Mpp1','Angpt1','Angptl2',
'Arsb','Axin2','Fbxw7','Flt3','Gca','Maml1','Mrc1','Nes','Nkd2','Rbpj']

marker_genes=sorted(list(set(marker_genes)-set(remove)))

#print("remove non expressed genes",len(remove),len(marker_genes))



#print(marker_genes)
#marker_genes=['Col4a1','Kit','Pdgfra','Epcam','Notch1','Notch2','Notch3','Notch4','Dll1','Dll4','Cd48','Mki67','Pecam1','Kitl','Cd34','E2f2','Cxcl12','Bmp5','Bmp7','Bmp2','Tet1','Il7r','Fzd1']

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
fig.savefig('figures/total_counts_and_n_genes_by_counts.png')

#adata.var is gene
#adata.obs is cells
'''
#sc.pl.highest_expr_genes(adata,n_top=20,show=False,save='.png')
#sc.pp.filter_cells(adata, min_genes=10)
#sc.pp.filter_cells(adata, min_counts=1)
#sc.pp.filter_cells(adata, max_counts=350)

#adata = adata[adata.obs["pct_counts_mt"] < 20]
#print(f"#cells after MT filter: {adata.n_obs}")
#sc.pp.filter_genes(adata, min_cells=3)
#sc.pp.filter_genes(adata,min_counts=500)  # umi counts of genes to filter

#pyscenic
# mito and genes/counts cuts
#mito_genes = adata.var_names.str.startswith('MT-')
#for each cell compute fraction of counts in mito genes vs. all genes
#adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))

print('original data after filter', adata)




#sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#             jitter=0.4, multi_panel=True,show=False,save='.png')
#sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',show=False,save='.png')
#sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',show=False, save='.png')


#adata = adata[adata.obs['n_genes'] < 4000, :]
#adata = adata[adata.obs['percent_mito'] < 0.15, :]
# above two from pyscenic
#adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#adata = adata[adata.obs.pct_counts_mt < 5, :]

#adata = adata[adata.obs.n_genes_by_counts < 4000, :]


#print('\n\nbasic filter cell and gene\n\n', adata.shape)
#df=pd.DataFrame(data=adata.X.transpose(), index=adata.var_names , columns=adata.obs_names)
#df.to_csv("xxoutput_filter_"+filename+".csv")




#normalize counts per cell
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
print('normalize',adata.shape)

sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pl.highly_variable_genes(adata,show=False,save='.png')

#adata.raw = adata
#adata = adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#sc.pp.scale(adata, max_value=10)
#sc.tl.pca(adata,svd_solver='arpack')
#sc.pl.pca(adata,color='n_genes_by_counts',show=False, save='.png')

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata,resolution=1.5)
sc.tl.leiden(adata,resolution=1.5)

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, save='_t-test.png')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,show=False, save='_wilcoxon.png')

#adata.write('results_file')
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,show=False, save='_logreg.png')



plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["leiden","total_counts",  "n_genes_by_counts", "louvain"], wspace=0.4,show=False, save='1.png')

sc.pl.umap(adata, color=["leiden"], wspace=0.4,show=False, save='leiden.png')
sc.pl.umap(adata, color=["louvain"], wspace=0.4,show=False, save='louvain.png')


#sc.pl.dotplot(adata, marker_genes, groupby='louvain',show=False,save='_GM_louvain.png')


cellname=adata.obs_names.to_numpy()
genename=adata.var_names.to_numpy()
marker_genes=['G6pc','Pck1','Hsd17b7','Aldh1b1','Acsl1','Cyp1a2','Cald1','Aldh3a2', 'Eif3a',
'Alas1','Hsd17b6','Epcam','Notch3', 'Itgb2','Slc25a37', 'Fzd7', 'Lyve1','H2afy','Hc',
'Hsd3b3','Csf1r','Alcam','Egr1','Sdc3','Gpr182', 'Runx1', 'Pdpn','Pecam1',
'Kdr','Cd44','Kitl', 'Itgb7', 'Wnt4', 'Kit','Aqp1','Cebpa','Stab2','Cd48','Eng','Cxcl12',
 'Sardh', 'Sdc1','Plvap','Pdgfra','Dnase1l3','Psap','Laptm5','Shisa5','Timp3','Ltbp4', 'Acss2',
 'Gck','Lmna','Mki67','Col1a2','Col6a1', 'Notch2', 'Tcf4','Tcf7','Tpi1','Bcam', 'Fzd5',
 'Il7r','Ep300','Dll4','Slamf1','Eno4','Tet2', 'Tcf7l1', 'Cd34','Tet1', 'Cspg4','Hoxb5',
 'Podxl','Itga2b','Vcam1','Angpt2', 'Gapdhs', 'Kcnj16']

marker_genes=sorted(list(set(marker_genes)))
print(len(marker_genes))
sc.pl.dotplot(adata, marker_genes, groupby='leiden',show=False,save='_dot_plot_GM_leiden.png')


'''
'Hsd3b6','Cdh11','Flt3','Jag2','Mpl','Il6',
array(['Comt', 'Ldha',  'Akr1a1', 'Ugt2b1', 'Acsl5', 'Ugt2a3',
       'Igf1', 'Errfi1', 'Serping1', 'Adh4', 'Hsd17b2',
        'Akr1d1',  'Aldh7a1',  'Hsd17b12', 'Pdhb',
       'Gpd1', 'Cyp7b1', 'Pgam1',  'Dld', 'Cyp2c23', 'Proz',
          'Galm',
         'Pdha1', 'Npc2',
       'Adh7', 'Smpdl3a', 'Egfr', 'Pgm1', 'Fasn', 'Ctsc', 'Abcb4', 'Fyb',
       'Alas2', 'Gpi1', 'Fech', 'Lsr', 'Psmd3', 'Gm2a', 'Pabpc1', 'Cbr4',
       'Tkt', 'Tmem56', 'Eif3f', 'Cxadr', 'Srd5a1', 'Cyp2c55', 'Gnai2',
       'Gimap6', 'Hsd3b2', 'Grn', 'Rpp14', 'Csnk1a1',  'Mpeg1',
       'Acsl4', 'Hmgb1', 'Mpp1', 'Lcp1',   'Oxsm',
       'Dlat', 'Csk', 'Mcat',  'Epas1',  'Nrp1', 'Dek',
        'Bpgm',   'Serpinh1', 'Tinagl1',
       'Aldoc', 'Cyp2c38', 'Dpt', 'Mrc1', 'Minpp1', 'Fgf1',
       'Gimap4', 'Cav2',  'Adgre1',   'Esam',
       'Unc93b1', 'Cnp', 'Clec14a',  'Adpgk', 'Gca', 'Pkm', 'Mkrn1',
       'Acaca',  'Bmp2', 'Tfrc',  'Calcrl',
       'Pfkl', 'Wnt2', 'Cybb', 'Icam1', 'Cdh5', 'Sgms2',  'Stk17b',
       'Tubb6',  'Hgf', 'Ramp1', 'Arsb', 'Pld4', 'Smarca4',
       'Fstl1', 'Pfkm', 'Lhfp',  'Cd300lg',  'Timp2',
        'Acacb', 'Cyp1a1', 'Eno3', 'Cd83',
        'Pgm2', 'Mertk', 'Pth1r',  'Kctd12',
       'Srd5a3', 'Bmp5',  'G6pc3', 'Cyp17a1',  'Cygb',
        'Nid1',  'Ctnnal1', 'Ephb4', 'Elk3', 'Foxq1',
       'Cxcl14', 'Fzd4',   'Srd5a2', 'Aldh3b1', 'Flt4',
       'Selp', 'Rbpj',  'Rhoj', 'Fzd1', 'Tcf7l2', 'Ssh2',
        'Tek', 'Trim47', 'Tent5c', 'Ncf1',
       'Lepr', 'Pck2', 'Lmnb1', 'Selplg', 'Myh10', 'Aldoart1',
       'Tcf3', 'Tspan13', 'Fzd8', 'Lad1', 'Procr', 'Ccr2',
       'Akr1c18', 'Maml1', 'Ms4a1', 'Hk3',  'Dkk3',
       'Bank1', 'Itgal', 'Pgam2', 'Axin2', 'Pfkp', 'Meis2', 'Jag1',
       'Gimap3', 'Rassf4', 'Notch1', 'Cd93',
       'Hvcn1', 'Mal',
       'Tnfrsf13c', 'Hk1',  'Apobec3', 'Slc34a2', 'Vav1',
       'Lamp3', 'Meis1', 'Lck', 'Efnb2', 'Notch4', 'Klrb1c',
       'Vwf', 'E2f2', 'Ccr1', 'Angpt1', 'B4galt6', 'Cyp21a1',
       'Dll1', 'Ammecr1', 'Csf3r', 'Ndn', 'Fgf2',  'Mecom',
       'Itgam', 'Hoxb4', 'Tox', 'Prickle2', 'Acss1', 'Cyp2b9', 'Aldh3a1',
       'Bmp7', 'Gata2',  'Satb1', 'Sfrp1', 'Eno2', 'Mrvi1',
        'Nes', 'Tmod1', 'Ace', 'Gfap', 'Tgfb2', 'Tomt',
       'Sult2b1', 'Hkdc1',   'Hk2', 'Mmrn1',
       'Vangl2', 'Pou2af1',  'Aldh3b2', 'Gypa', 'Lrp2',
       'Lef1', 'Olr1', 'Lox', 'Txlnb', 'Slc12a1', 'Aldh3b3', 'Cxcr2',
       'Nkd2', 'Sult1e1', 'Acsl6', 'Ddx4', 'Ldhc', 'Kcnj1', 'Acsbg1',
       'Fzd3', 'F13a1', 'Hsd11b2', 'Dkk2', 'Hsd17b1', 'Fzd2', 'Cyp2b23',
        'Celsr2', 'Obscn',  'Akap14', 'Gnaz', 'Cd177',
        'Aldoart2', 'Cyp2b19', 'Ryr2', 'Ldhal6b', 'Acsf3',
       'Chodl', 'Ivl', 'Cyp11b1', 'Sfrp2', 'Dkk1', 'Cyp11a1',
       '1700061G19Rik', 'Acsbg2', 'Olah', 'Pdha2', 'Hsd17b3'],
      dtype=object)
'''


plt.rcParams["figure.figsize"] = (4, 4)
new_cluster_names = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24']
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden',legend_loc='on data', title='', frameon=False,show=False, save='new_leiden2.png')
'''
#adata.rename_categories('authorClusterin', new_cluster_names)
#sc.pl.umap(adata, color='authorClustering',legend_loc='on data', title='', frameon=False,show=False, save='2.png')
#sc.pl.dotplot(adata, marker_genes, groupby='authorClustering',show=False,save='_GM_AC.png')
'''

#plot tsne
#sc.tl.tsne(adata, n_pcs = 30)
#sc.pl.tsne(adata, color='louvain',show=False,save='.png')

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_02_dim_reduction.html
#https://scanpy.readthedocs.io/en/stable/api/index.html
#https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html


##Computing the neighborhood graph
#sc.tl.paga(adata)
#sc.pl.paga(adata,plot=False)
#sc.tl.umap(adata,init_pos='paga')

adata.write_h5ad('figures/saveall')

adata.obs.louvain.to_csv('figures/louvain_output.dat',header=True)
adata.obs.leiden.to_csv('figures/leiden_output.dat',header=True)

#adata.obs.authorClustering.to_csv('figures/authorClustering_output.dat',header=True)


print('final', adata.shape)
