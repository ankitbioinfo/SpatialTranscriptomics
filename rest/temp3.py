import Spatial_Interactions as st
import numpy as np
import time
import os
import pickle
import pandas as pd



spdatapath='./inputSP/'
outputFolder='./'
st.create_directory(outputFolder)
saveSpatial=outputFolder+'spatial_ct_ct_interactions/'
st.create_directory(saveSpatial)


spdatapath2='./caseBasics/find_MNN_poss4_res0.5_dis10_basic/'
spdatapath2='./output_from_linux_cluster/caseB_res_0.6/find_MNN_poss4_res0.6_dis5_/'
spdatapath2='./find_MNN_poss4_res0.5_dis15_/'


neigh=50

clusterFilename=spdatapath2+'3_deg_annotation_spatial_'+str(neigh)+'_50_cluster.csv'
positionFilename=spdatapath+'tissue_positions_list.csv'
celltypeFilename=spdatapath2+'3_deg_annotation_spatial_'+str(neigh)+'_50_ct_name.dat'



inputRadius=[0]


#check your other inputs for optimization
seed=36851234

figuresize1=[4,3] #First and second value is the width and height of the figure
#This figure size use for log reg coeff vs cell type features in processed_glioblastoma1/spatial/L2_multi_linear/TopCoeff_R

n_jobs=4 #no of processors For details see here https://scikit-learn.org/stable/glossary.html#term-n_jobs
lambda_c_ranges=list(np.power(2.0, np.arange(-12, 12))) # inverse of lambda regularization
BothLinearAndCrossTerms=1# If only linearterms then put 1; For both linear and crossterms use 2
K_fold=5 #no of cross folds
n_repeats=1 #no of times to run the logistic regression after finding the hyperparameters
coeff_cutoff=25 # No. of coefficient want to print in the figure
strategy='L2_multi' #Name of the strategy you want to compute the interactions options are [L1_multi, L1_ovr, L2_multi, L2_ovr, elasticnet_multi, elasticnet_ovr]
#strategy='L1_ovr'
#strategy='elasticnet_multi'

removed_CTs_before_finding_CT_CT_interactions=['NM']


for i in range(len(inputRadius)):
    print('\nRadius ',inputRadius[i],'\n')
    fname=saveSpatial+'normalized_spatial_neighbors_'+str(inputRadius[i])+'.dat'
    flag=1
    if os.path.isfile(fname):
        filesize = os.path.getsize(fname)
        if filesize>0: #If file is already exist and have size greater than 0 then no need to run again. It will save some time if you want to run it again with different parameters
            flag=0
            #neighbors=pickle.load( open(saveSpatial+'save_neighbors_'+str(inputRadius[i])+'.p', "rb" ) )
    if flag==1:
        delimiter=','
        PP,cluster,noct,fraction_CT,clusterWithBarcodeId= st.reading_data(positionFilename,clusterFilename,celltypeFilename,saveSpatial,delimiter,removed_CTs_before_finding_CT_CT_interactions)
        M, neighbors=st.create_spatial_CT_feature_matrix(inputRadius[i],PP,cluster,noct,fraction_CT,saveSpatial)
        pickle.dump(neighbors,open(saveSpatial+'save_neighbors_'+str(inputRadius[i])+'.p', 'wb'))
        df=pd.DataFrame(clusterWithBarcodeId)
        df.to_csv(saveSpatial+'save_clusterid_'+str(inputRadius[i])+'.csv',index=False)
        #pickle.dump(cluster,open(saveSpatial+'save_clusterid_'+str(inputRadius[i])+'.p', 'wb'))

    st.run_logistic_regression_on_spatial_features(saveSpatial,n_repeats,K_fold,coeff_cutoff,strategy,seed,n_jobs,lambda_c_ranges,BothLinearAndCrossTerms,inputRadius[i],figuresize1)


#plot evaluation score
figuresize2=[4,3] #This figure size use for plotting of scores
st.plot_evaluation_scores(inputRadius,BothLinearAndCrossTerms,saveSpatial,strategy,figuresize2)

'''
#plot_normalized_coefficients_radius_wise
figuresize3=[10,3] #This figure size use for plotting of radius wise log reg coeff
st.plot_normalized_coefficients_radius_wise(inputRadius,BothLinearAndCrossTerms,saveSpatial,strategy,figuresize3)
'''






#########################################################################
niche_cutoff=0.15 #it is normalized large value has fewer connections and small value has larger connections
figuresize1=[10,7]
for i in range(len(inputRadius)):
    st.plot_niche_interactions(saveSpatial,strategy,BothLinearAndCrossTerms,inputRadius[i],figuresize1,niche_cutoff)
