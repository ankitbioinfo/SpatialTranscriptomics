
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

authorClustering=reading_clustering()
#print(authorClustering)




filename='vizgen'
#csvfilename='Blank_genes_removed.csv'
csvfilename='gene_by_cell_counts.csv'

adata=sc.read(csvfilename).transpose()
adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#adata.obs_names_make_unique()

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/scanpy/scanpy_04_clustering.html
new_cluster_names = ['AEC', 'SEC','MK', 'Hepatocyte','Macrophage', 'Myeloid','Erythroid progenitor', 'Erythroid cell']
adata.obs['authorClustering'] = authorClustering.astype(str)
print('original data', adata)

f=open('../geneOrdername140.dat','r')
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

#marker_genes=sorted(list(set(marker_genes)-set(remove)))

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

#sc.pl.highest_expr_genes(adata,n_top=20,show=False,save='.png')
#sc.pp.filter_cells(adata, min_genes=10)
#sc.pp.filter_cells(adata, min_counts=1)
#sc.pp.filter_cells(adata, max_counts=350)
adata = adata[adata.obs["pct_counts_mt"] < 20]
print(f"#cells after MT filter: {adata.n_obs}")
sc.pp.filter_genes(adata, min_cells=3)
#sc.pp.filter_genes(adata,min_counts=500)  # umi counts of genes to filter

#pyscenic
# mito and genes/counts cuts
#mito_genes = adata.var_names.str.startswith('MT-')
#for each cell compute fraction of counts in mito genes vs. all genes
#adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))





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


print('\n\nbasic filter cell and gene\n\n', adata.shape)
df=pd.DataFrame(data=adata.X.transpose(), index=adata.var_names , columns=adata.obs_names)
df.to_csv("xxoutput_filter_"+filename+".csv")




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
sc.tl.louvain(adata)
sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, save='_t-test.png')

sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,show=False, save='_wilcoxon.png')

#adata.write('results_file')
sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,show=False, save='_logreg.png')



plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["leiden","total_counts",  "n_genes_by_counts", "louvain"], wspace=0.4,show=False, save='1.png')

sc.pl.dotplot(adata, marker_genes, groupby='louvain',show=False,save='_GM_louvain.png')
sc.pl.dotplot(adata, marker_genes, groupby='leiden',show=False,save='_GM_leiden.png')



plt.rcParams["figure.figsize"] = (4, 4)
adata.rename_categories('authorClustering', new_cluster_names)
sc.pl.umap(adata, color='authorClustering',legend_loc='on data', title='', frameon=False,show=False, save='2.png')

sc.pl.dotplot(adata, marker_genes, groupby='authorClustering',show=False,save='_GM_AC.png')


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


adata.obs.louvain.to_csv('figures/louvain_output.dat')

print('final', adata.shape)
