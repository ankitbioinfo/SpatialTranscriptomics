import Spatial_PathwaysNMF as sppath
#import Spatial_Pathways as sppath
import scanpy as sc
from types import SimpleNamespace
#import gseapy
import xlsxwriter
import numpy as np




inputRadius=[0]

inputDir='./'

spdatapath='./inputSP/'
scdatapath='./inputSC/'  #input_sp_dir  input_sc_dir

#spdatapath='./countdata/'
#scdatapath='./countdata/'

input_spatial_analysis_dir=inputDir+'spatial_ct_ct_interactions_poss4_res0.5_dis15/'
no_of_pc=5
strategy='L2_multi'
BothLinearAndCrossTerms=1

if BothLinearAndCrossTerms==1:
    newstrategy=strategy+'_linear/'
else:
    newstrategy=strategy+'_cross/'

# load single cell data
'''
ad_sc=sc.read_h5ad(scdatapath+'scTransform_singleCell_si1.h5ad')
full_ad_sc=sc.read_h5ad(scdatapath+'HVG_scTransform_singleCell_si_3000.h5ad')
ad_sp=sc.read_h5ad(spdatapath+'scTransform_spatial_si1.h5ad')
'''
ad_sc=sc.read_h5ad(scdatapath+'common_counts_sc.h5ad')
full_ad_sc=sc.read_h5ad(scdatapath+'HVG_counts_10000.h5ad')
#full_ad_sc=sc.read_h5ad(scdatapath+'sc_liver_data.h5ad')
ad_sp=sc.read_h5ad(spdatapath+'common_counts_sp.h5ad')


print('A',ad_sc)
print('B',full_ad_sc)

gn1=ad_sc.var_names.to_numpy()
gn2=full_ad_sc.var_names.to_numpy()

for i in range(len(gn1)):
    if gn1[i]=='Serping1':
        d1=ad_sc[:,i]

for i in range(len(gn2)):
    if gn2[i]=='Serping1':
        d2=full_ad_sc[:,i]


A=d1.X.toarray()
B=d2.X.toarray()

value=np.array_equal(A,B)

print(A.shape,B.shape,value)
print(np.mean(A),np.mean(B),np.std(A),np.std(B))

print('\nCommon genes',len(gn1),len(gn2),len(set(gn1).intersection(set(gn2))))
print('\n\n')



#print("1",full_ad_sc)
#print("2",ad_sp)

#sc.pp.filter_cells(full_ad_sc, min_counts=1000) #min_genes=10,
#sc.pp.filter_cells(ad_sp, min_counts=10) #min_counts=500

#sc.pp.filter_genes(full_ad_sc, min_cells=5) #min_genes=10,
#sc.pp.filter_genes(ad_sp, min_cells=5) #min_counts=500
#sc.pp.filter_genes(full_ad_sc, min_counts=1000) #min_genes=10,
#sc.pp.filter_genes(ad_sp, min_counts=10) #min_counts=500



#print("3",full_ad_sc)
#print("4",ad_sp)

#print(mat)


clusterfname='./inputSC/scdata2_most_recent/cluster_poss4_25.csv'
celltypefname='./inputSC/scdata2_most_recent/nameOfCT_poss4_25.dat'


#clusterfname=scdatapath+'cluster_SI_merge.csv'
#celltypefname=scdatapath+'nameOfCT_SI_merge.dat'


'''
ad_sc=sc.read_h5ad(scdatapath+'scTransform_singleCell.h5ad')
full_ad_sc=sc.read_h5ad(scdatapath+'sc_liver_data_HVG.h5ad')

clusterfname=scdatapath+'scdata2_most_recent/cluster_24.csv'
celltypefname=scdatapath+'scdata2_most_recent/nameOfCT_24.dat'

ad_sp=sc.read_h5ad(spdatapath+'scTransform_spatial.h5ad')
'''

outdir=input_spatial_analysis_dir+'Regression_on_genePC_loading_'
figuresize=[5,5]


#gene_set_names = gseapy.get_library_name(organism='Mouse')
gene_set_names=[]
database=['BioPlanet_2019','Reactome_2016']
f=open(inputDir+'./comb.dat','r')
LRdb=f.readlines()
f.close()


input={}
input['celltypefname']=celltypefname
input['clusterfname']=clusterfname
input['ad_sp']=ad_sp
input['ad_sc']=ad_sc
input['full_ad_sc']=full_ad_sc
input['input_spatial_analysis_dir']=input_spatial_analysis_dir
input['inputRadius']=inputRadius
input['no_of_pc']=no_of_pc
input['strategy']=newstrategy
input['outdir']=outdir

input['LRdb']=LRdb
input['database']=database
input['gene_set_names']=gene_set_names
input['pathwayorganism']='Mouse'
input['pathwayCutoff']=0.5

input['shap_cluster_cutoff']=0.5
input['shap_computation']=0 # 1 means perform 0 means not
input['pathway_plot']=0  #1 means plot 0 means no plot
input['coeff_cutoff_for_rid_reg']=0
input['logistic_coef_cutoff']=0

input['n_jobs']=-1 #no of processors For details see here https://scikit-learn.org/stable/glossary.html#term-n_jobs
input['lambda_c']=list(np.power(2.0, np.arange(-10, 10))) # inverse of lambda regularization
#BothLinearAndCrossTerms=1# If only linearterms then put 1; For both linear and crossterms use 2
input['K_fold']=5 #no of cross folds
input['seed']=3685134
input['n_repeats']=1 #no of times to run the logistic regression after finding the hyperparameters
#coeff_cutoff=25 # No. of coefficient want to print in the figure
#strategy='L2_multi' #Nam

input['gene_correlation_plots_and_sheets']=0    #1 means plot 0 means no plot



inputdata=SimpleNamespace(**input)

sppath.main(inputdata)
