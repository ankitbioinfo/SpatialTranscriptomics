

from sklearn.datasets import make_classification
from matplotlib import pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import make_scorer,accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import Pipeline
#from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.model_selection import RandomizedSearchCV,GridSearchCV,cross_val_predict, cross_val_score,RepeatedStratifiedKFold,StratifiedShuffleSplit

#from imblearn.over_sampling import SMOTE, SMOTEN,ADASYN, KMeansSMOTE, SVMSMOTE
from sklearn.utils import class_weight
from sklearn.metrics import roc_curve, auc

#Metrics
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import log_loss
from sklearn.metrics import zero_one_loss
from sklearn.metrics import matthews_corrcoef


import pandas as pd
import numpy as np
import seaborn as snn
import os,sys
import random
import warnings
import time
warnings.filterwarnings("ignore")


def plot_multiclass_roc(clf, X_test, y_test, n_classes):
    y_score = clf.decision_function(X_test)
    # structures
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    # calculate dummies once
    y_test_dummies = pd.get_dummies(y_test, drop_first=False).values
    for i in range(n_classes):
        fpr[i], tpr[i], _ = roc_curve(y_test_dummies[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    return fpr, tpr, roc_auc




def create_directory(outputFolder):
    answer=os.path.isdir(outputFolder)
    if answer==True:
        pass
    else:
        os.mkdir(outputFolder)

def plot_results(maindir,nrow,ncol,nameOfCellType,radius,lambda_c,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,BothLinearAndCrossTerms,inputdir,fpr, tpr, roc_auc):

	#print(CTFeatures,inputFeatures,classes,nameOfCellType)
	filename='R'+str(radius)
	fig, ax = plt.subplots(nrow,ncol, figsize=(10, 7))
	plotaxis=[]
	for i in range(nrow):
		for j in range(ncol):
			plotaxis.append([i,j])

	highestROCofcelltype=[]
	for w in sorted(roc_auc, key=roc_auc.get, reverse=True):
		#print(w, roc_auc[w])
		highestROCofcelltype.append(w)


	#for i in range(len(classes)):
	'''
	for i in range(nrow*ncol):
		value=plotaxis[i]
		ax[value[0],value[1]].plot([0, 1], [0, 1], 'k--')
		ax[value[0],value[1]].set_xlim([0.0, 1.0])
		ax[value[0],value[1]].set_ylim([0.0, 1.05])
		if value[0]==(nrow-1):
			ax[value[0],value[1]].set_xlabel('False Positive Rate')
		else:
			ax[value[0],value[1]].set_xticks([])

		if i%ncol==0:
			ax[value[0],value[1]].set_ylabel('True Positive Rate')
		else:
			ax[value[0],value[1]].set_yticks([])

		ax[value[0],value[1]].set_title(str(highestROCofcelltype[i])+' : '+nameOfCellType[highestROCofcelltype[i]])
		ax[value[0],value[1]].plot(fpr[highestROCofcelltype[i]], tpr[highestROCofcelltype[i]], label='ROC(area = %0.2f)' % (roc_auc[highestROCofcelltype[i]]))

		ax[value[0],value[1]].legend(loc="best",fontsize=8)
		#ax[value[0],value[1]].grid(alpha=.4)
	snn.despine()
	#plt.suptitle('Receiver operating characteristic example')
	plt.tight_layout()
	plt.savefig(maindir+'ROC_'+filename+'.png')
	'''


	plt.figure(figsize=(12,10))
	classNames=[]
	for i in range(len(classes)):
		classNames.append(nameOfCellType[classes[i]])
	snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classNames, annot_kws={"size": 5},yticklabels=classNames)
	plt.xlabel('Predicted classes')
	plt.ylabel('Truth classes')
	plt.title('R = '+str(radius)+', C='+str(lambda_c))
	plt.tight_layout()
	plt.savefig(maindir+'Confusing_matrix_'+filename+'.png')

	plt.figure(figsize=(5,8))
	#plt.figure()
	#snn.set(font_scale=0.4)
	b=snn.heatmap(coef.transpose(),xticklabels=classNames)#yticklabels=CTFeatures
	#plt.xticks(rotation=90)
	_, ylabels= plt.yticks()
	b.set_yticklabels(ylabels, size = 5)

	if BothLinearAndCrossTerms==1:
		plt.ylabel('Features linear  terms')
	else:
		plt.ylabel('Features cross terms')

	#plt.xlabel('# of classes (no of cell types)')
	plt.title('R = '+str(radius)+', C='+str(lambda_c))

	plt.tight_layout()
	plt.savefig(maindir+'weight_matrix_'+filename+'.png')

	plt.figure(figsize=(12,6))

	plt.subplot(1,3,1)
	snn.heatmap(np.log(x_train+1))
	plt.xlabel('# of input Features')
	plt.title('training set')
	plt.ylabel('75% of data')

	plt.subplot(1,3,2)
	snn.heatmap(np.log(x_test+1))
	plt.xlabel('# of input Features')
	plt.title('testing set')
	plt.ylabel('25% of data')

	plt.subplot(1,3,3)
	snn.heatmap(predicted_probs,xticklabels=classes)
	plt.title('Predicted probability')
	plt.xlabel('# of classes (no of cell types)')
	plt.tight_layout(.5)
	plt.savefig(maindir+'predicted_probability_'+filename+'.png')
	#print(predicted_probs)
	#prob=sigmoid( np.dot([y_train, y_test,1], log_reg_model.coef_.T) + log_reg_model.intercept_ )
	#print(prob)


def counterPartLigRec(gene,totalLRpairs,lig_or_rec,counterPart):

	#if cc is ligand
	if lig_or_rec=='L':
		a=0
		b=1

	#if cc is receptor
	if lig_or_rec=='R':
		a=1
		b=0

	for key in totalLRpairs:
		l=key.split('--')
		if l[a]==gene:
			if l[b] not in counterPart:
				counterPart.append(l[b])

	return counterPart



def read_communication_processed_data(radius,BothLinearAndCrossTerms,ligand_receptor_check_over_db,inputdir):

	name=inputdir+'weighted_normalized_comm_neighbors_'+str(radius)+'.dat'
	dataLR = np.genfromtxt(open(name, "rb"), delimiter='\t', skip_header=0)
	neighborhoodClassLR= dataLR[:,2:]
	target= dataLR[:,1]


	f=open(inputdir+'total_LR_genes_exist.dat','r')
	cont=f.readlines()
	n=len(cont)

	f=open('sort_3_db_L_R_high_confident.dat','r')
	contdb=f.readlines()
	totalLRpairs={}
	ligand_one={}
	receptor_one={}
	combined_one={}

	for j in range(len(contdb)):
		l=contdb[j].split()
		flag1=0
		flag2=0
		for i in range(n):
			t=cont[i][0:-1].split('\t')
			if (t[0].upper()==l[0].upper()):
				flag1=1
			if (t[0].upper()==l[1].upper()):
				flag2=1

		if flag1+flag2==2:
		#if (flag1+flag2)>0:
			totalLRpairs[l[0]+'--'+l[1]]=1
			ligand_one[l[0]]=1
			receptor_one[l[1]]=1
			combined_one[l[0]]=1
			combined_one[l[1]]=1

			#totalLRpairs[l[1]+'--'+l[0]]=1

	print("3 database LR pairs =",len(totalLRpairs), '; ligands =', len(ligand_one), '; receptors =',len(receptor_one),'; combined =',len(combined_one))
	#print(totalLRpairs,ligand_one,receptor_one)



	#fw=open('dummy.dat','w')
	#for key in totalLRpairs:
	#	fw.write(str(key)+'\n')
	#print(top_expressed_genes,sorted(combined_one),sorted(cont))
	proteins=[]
	ne_=[]
	#check what are the top 5 gene expressed by cell types in csv file
	f=open(inputdir+'highest_to_lowest_expressed_genes_in_celltypes.dat')

	total_gene=[]
	CT_specific_top_exp_LRgenes=[]
	for line in f:
		l=line[0:-1].split(';')
		top_10=len(l)-1
		cc_=[]
		pos_=[]
		for j in range(1,top_10+1):#len(l)):
			gene=l[j].upper()
			total_gene.append(gene)
			# central cell is ligand
			try:
				ligand_one[gene]
				if gene not in cc_:
					cc_.append(gene)
					pos_.append(j)
				ne_=counterPartLigRec(gene,totalLRpairs,'L',ne_)
			except KeyError:
				pass
			# central cell is receptor
			try:
				receptor_one[gene]
				if gene not in cc_:
					cc_.append(gene)
					pos_.append(j)
				ne_=counterPartLigRec(gene,totalLRpairs,'R',ne_)
			except KeyError:
				pass
		CT_specific_top_exp_LRgenes.append([cc_,pos_])

	#print('communication channel cc', len(cc_), ', neighbor ', len(ne_), 'total number of genes',p_ng  )
	#print('x1',cc_ligand, ne_receptor)
	#print('x2',ne_ligand, cc_receptor)

	if ligand_receptor_check_over_db:
		data=ne_
	else:
		data=total_gene

	print('lig rec check over db',len(data))
	gene_index_expressed_in_neighbors=[]
	for i in range(n):
		l1=cont[i].split('\t')
		g1=l1[0]
		#print(i,g1)
		if g1.upper() in data:
			gene_index_expressed_in_neighbors.append(i)
			proteins.append([i,g1])


	for i in range(len(gene_index_expressed_in_neighbors)):
		for j in range(i+1,len(gene_index_expressed_in_neighbors)):
			name=proteins[i][1]+'--'+proteins[j][1]
			if BothLinearAndCrossTerms==2:
				proteins.append([gene_index_expressed_in_neighbors[i],gene_index_expressed_in_neighbors[j],name])

	#fw=open('dummy_coomunicationChannel.dat','w')
	#for i in range(len(proteins)):
	#	fw.write(str(i)+'\t'+str(proteins[i][-1])+'\n')

	#print(all_LR_pairs_found)
	all_LR_pairs_found=[]
	for i in range(len(proteins)):
		all_LR_pairs_found.append(proteins[i][-1])


	neighborhoodClass=np.zeros((neighborhoodClassLR.shape[0],len(proteins)),dtype=np.float)
	neighborhoodClass2=np.zeros((neighborhoodClassLR.shape[0],len(proteins)),dtype=np.float)
	for i in range(len(proteins)):
		if len(proteins[i])==2:
			neighborhoodClass2[:,i]=neighborhoodClassLR[:,proteins[i][0]]
		if BothLinearAndCrossTerms==2:
			if len(proteins[i])==3:
				neighborhoodClass2[:,i]=neighborhoodClassLR[:,proteins[i][0]]*neighborhoodClassLR[:,proteins[i][1]]

	neighborhoodClass=neighborhoodClass2

	#print('data shape',dataLR.shape, 'target', target.shape, "original neighbor shape",neighborhoodClassLR.shape,'features',len(all_LR_pairs_found))
	print('total_len_of_feature',len(proteins),'modified neighbor shape',neighborhoodClass.shape)



	return neighborhoodClass,target,all_LR_pairs_found,CT_specific_top_exp_LRgenes,ne_,totalLRpairs




'''
def calculate_class_weights(vector):
	a=np.unique(vector)
	freq=[]
	for i in range(len(a)):
		freq.append(np.sum(vector==a[i]))

	total=np.sum(freq)
	#print('tot',total)
	cw={}
	for i in range(len(a)):
		cw[a[i]]=freq[i]/float(total)

	return cw
'''






def model_log_regression(K_fold,n_repeats,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms,seed,n_jobs):
	#polynomial = PolynomialFeatures(degree = 2, interaction_only=True, include_bias=False)
	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#log_reg_model = LogisticRegression(max_iter=500,penalty='elasticnet',l1_ratio=0.5,class_weight=classweight,solver='saga')

	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#pipe=Pipeline([('polynomial_features',polynomial), ('logistic_regression',log_reg_model)])
	#X_poly = poly.fit_transform(neighborhoodClass[:,[0,1,2]])
	#print(X_poly.shape,CTFeatures)
	#print(X_poly[0])
	hyperparameter_scoring = {  'f1_weighted': make_scorer(f1_score, average = 'weighted')}

	parameters = {'C':lambda_c }

	if strategy=='L1_multi':
		log_reg_model = LogisticRegression(penalty='l1',multi_class='multinomial',class_weight='balanced',solver='saga',n_jobs=n_jobs)#very slow
	if strategy=='L1_ovr':
		log_reg_model = LogisticRegression(penalty='l1',multi_class='ovr',class_weight='balanced',solver='liblinear',n_jobs=n_jobs)
	if strategy=='L2_multi':
		log_reg_model = LogisticRegression(penalty='l2',multi_class='multinomial',class_weight='balanced',solver='lbfgs',n_jobs=n_jobs)
	if strategy=='L2_ovr':
		log_reg_model = LogisticRegression(penalty='l2',multi_class='ovr',class_weight='balanced',solver='lbfgs',n_jobs=n_jobs)
	if strategy=='elasticnet_multi':
		log_reg_model = LogisticRegression(penalty='elasticnet',multi_class='multinomial',class_weight='balanced',solver='saga',n_jobs=n_jobs)
		parameters = {'C':lambda_c, 'multi_class':['ovr','multinomial'], 'l1_ratio':np.linspace(0,1,10)  }
	if strategy=='elasticnet_ovr':
		log_reg_model = LogisticRegression(penalty='elasticnet',multi_class='ovr',class_weight='balanced',solver='saga',n_jobs=n_jobs)
		parameters = {'C':lambda_c, 'multi_class':['ovr','multinomial'], 'l1_ratio':np.linspace(0,1,10)  }


	'''
	flag=1
	while(flag):
		seed=seed+1
		sss = RepeatedStratifiedKFold(n_splits=K_fold, n_repeats=1 ,random_state=seed)
		gs_grid = GridSearchCV(log_reg_model, parameters, scoring=hyperparameter_scoring, refit='f1_weighted',cv=sss,n_jobs=n_jobs)
		gs_random = RandomizedSearchCV(estimator=log_reg_model, param_distributions=parameters, scoring=hyperparameter_scoring, refit='f1_weighted',cv = sss,n_jobs=n_jobs)
		pipe_grid=Pipeline([    ('StandardScaler',StandardScaler()), ('logistic_regression_grid',gs_grid)])
		pipe_random=Pipeline([   ('StandardScaler',StandardScaler()), ('logistic_regression_random',gs_random)])
		pipe_grid.fit(neighborhoodClass,target)
		pipe_random.fit(neighborhoodClass,target)

		LR_grid= pipe_grid.named_steps['logistic_regression_grid']
		LR_random= pipe_random.named_steps['logistic_regression_random']

		#if LR_grid.best_params_['C']==LR_random.best_params_['C']:
		if True:
			flag=0
			lambda_c=LR_grid.best_params_['C']
			print('Inverse of lambda regularization found', lambda_c)
		else:
			print('Searching hyperparameters ', 'Grid method:', LR_grid.best_params_['C'], ', Randomized method:', LR_random.best_params_['C'])
	#'''

	lambda_c=0.0009765625

	scorecalc=[]
	for i in range(15):
		scorecalc.append([])

	seed=seed+1
	sss = RepeatedStratifiedKFold(n_splits=K_fold, n_repeats=n_repeats ,random_state=seed)
	cmn=[]
	coef=[]
	for train_index, test_index in sss.split(neighborhoodClass,target):
		x_train,x_test=neighborhoodClass[train_index],neighborhoodClass[test_index]
		y_train,y_test=target[train_index],target[test_index]

		if strategy=='L1_multi':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='l1',multi_class='multinomial',class_weight='balanced',solver='saga',n_jobs=n_jobs)#very slow
		if strategy=='L1_ovr':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='l1',multi_class='ovr',class_weight='balanced',solver='liblinear',n_jobs=n_jobs)
		if strategy=='L2_multi':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='l2',multi_class='multinomial',class_weight='balanced',solver='lbfgs',n_jobs=n_jobs)
		if strategy=='L2_ovr':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='l2',multi_class='ovr',class_weight='balanced',solver='lbfgs',n_jobs=n_jobs)
		if strategy=='elasticnet_multi':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='elasticnet',multi_class='multinomial',l1_ratio=0.5,class_weight='balanced',solver='saga',n_jobs=n_jobs)
		if strategy=='elasticnet_ovr':
				log_reg_model = LogisticRegression(C=lambda_c,penalty='elasticnet',multi_class='ovr',l1_ratio=0.5,class_weight='balanced',solver='saga',n_jobs=n_jobs)



		pipe=Pipeline([('StandardScaler',StandardScaler()), ('logistic_regression',log_reg_model)])

		S=pipe.named_steps['StandardScaler']

		pipe.fit(x_train, y_train)
		y_pred=pipe.predict(x_test)
		y_prob = pipe.predict_proba(x_test)


		log_metric=log_loss(y_test,y_prob)
		c_k_s=cohen_kappa_score(y_test,y_pred)
		zero_met=zero_one_loss(y_test,y_pred)
		hl=hamming_loss(y_test,y_pred)
		mcc=matthews_corrcoef(y_test,y_pred)

		scorecalc[0].append(pipe.score(x_test, y_test))
		#precision, recall, fscore, support = score(y_test, predicted)
		scorecalc[1].append(f1_score(y_test, y_pred, average="macro"))
		scorecalc[2].append(precision_score(y_test, y_pred, average="macro"))
		scorecalc[3].append(recall_score(y_test, y_pred, average="macro"))
		scorecalc[4].append(f1_score(y_test, y_pred, average="micro"))
		scorecalc[5].append(precision_score(y_test, y_pred, average="micro"))
		scorecalc[6].append(recall_score(y_test, y_pred, average="micro"))
		scorecalc[7].append(f1_score(y_test, y_pred, average="weighted"))
		scorecalc[8].append(precision_score(y_test, y_pred, average="weighted"))
		scorecalc[9].append(recall_score(y_test, y_pred, average="weighted"))
		scorecalc[10].append(c_k_s)
		scorecalc[11].append(log_metric)
		scorecalc[12].append(mcc)
		scorecalc[13].append(hl)
		scorecalc[14].append(zero_met)

		#poly = pipe.named_steps['polynomial_features']
		#print(poly)
		LR= pipe.named_steps['logistic_regression']
		coef.append(LR.coef_)
		cmn.append(confusion_matrix(y_test,y_pred,normalize='true'))

	cmn_std=np.std(np.array(cmn),axis=0)
	coef_std=np.std(np.array(coef),axis=0)
	comp_score_std=np.std(np.array(scorecalc),axis=1)

	cmn=np.mean(np.array(cmn),axis=0)
	coef=np.mean(np.array(coef),axis=0)
	comp_score=np.mean(np.array(scorecalc),axis=1)

	print('training',x_train.shape,'testing',x_test.shape,'coeff',coef.shape,'Iteration',len(scorecalc[0]))

	classes=LR.classes_.astype(int)

	fpr, tpr, roc_auc=plot_multiclass_roc(pipe, x_test, y_test, n_classes=len(classes))


	#CTFeatures=poly.get_feature_names()
	#print("Features", CTFeatures)
	print("accuracy score\t",np.mean(scorecalc[0]))
	print("\n\nmacro")
	print("f1 score\t",np.mean(scorecalc[1]))
	print("precision score\t",np.mean(scorecalc[2]))
	print("recall score\t",np.mean(scorecalc[3]))

	print("\n\nmicro f1, precision, recall all same")
	print("score\t",np.mean(scorecalc[4]))
	#print("precision score in 10 run\t",np.mean(scorecalc[5]))
	#print("recall score in 10 run\t",np.mean(scorecalc[6]))

	print("\n\nWeighted")
	print("f1 score\t",np.mean(scorecalc[7]))
	print("precision\t",np.mean(scorecalc[8]))
	print("recall score\t",np.mean(scorecalc[9]))


	print('\n\ncohen_kappa_score (best=1): {0:.4f}'.format(np.mean(scorecalc[10])))
	print('log_loss or cross entropy (best=lowest): {0:.4f}'.format(np.mean(scorecalc[11])))
	print('matthews_corrcoef: {0:.4f}'.format( np.mean(scorecalc[12])  ))
	print('hemming_loss (best=lowest): {0:.4f}'.format( np.mean(scorecalc[13] )))
	print('zero_one_loss (best=0): {0:.4f}'.format(np.mean(scorecalc[14])))

	return cmn,coef,comp_score,cmn_std,coef_std,comp_score_std,classes, lambda_c,x_test,x_train,y_prob ,fpr, tpr, roc_auc


	#return cmn,coef,classes,x_test,x_train,y_prob ,fpr, tpr, roc_auc, scorecalc


def plot_lig_rec_with_colors(ax,globalIndex,coeff_of_CT,std_of_coeff,rankL,xcc,mycolor):
    delta=0.1*max(np.abs(coeff_of_CT))

    newxcc=[]
    #xcc=xcc[0:5]
    #mycolor=mycolor[0:5]

    for i in range(len(globalIndex)):
        for j in range(len(xcc)):
            if xcc[j]==globalIndex[i]:
                newxcc.append(i)

    largestCoeff=(np.min(coeff_of_CT)-np.max(std_of_coeff))*np.ones(len(newxcc))

    #print(globalIndex,xcc,newxcc)
    #print(largestCoeff,mycolor)

    ax.plot(largestCoeff,newxcc,mycolor)
    for rank in range(len(newxcc)):
        ax.text(largestCoeff[rank]+delta,newxcc[rank],rankL[rank],fontsize=7)


#def find_interacting_LR(cmn,coef,LRFeatures,nameOfCellType,fw):
def find_interacting_LR(cmn,coef,cmn_std,coef_std,coeff_cutoff,cc,ne,LRpairs,LRFeatures,nameOfCellType,fw,filename,figuresize):

	gene_LRB={} #ligand receptor or both
	for key in LRpairs:
		name=key.split('--')
		if name[0] in gene_LRB:
			if gene_LRB[name[0]]=='R':
				gene_LRB[name[0]]='B'
		else:
				gene_LRB[name[0]]='L'

		if name[1] in gene_LRB:
			if gene_LRB[name[1]]=='L':
				gene_LRB[name[1]]='B'
		else:
				gene_LRB[name[1]]='R'

	#print(LRpairs)
	#print(gene_LRB)
	#for i in range(len(cc)):
		#print(i,cc[i][0][0:4],cc[i][1][0:4])
	a=np.diag(cmn)
	b=np.diag(cmn_std)
	goodPredictedCellType=np.argsort(-a)
	create_directory(filename)
	for k in range(len(a)):
		if a[goodPredictedCellType[k]]>=0:
			top_exp_genes_in_ct=cc[goodPredictedCellType[k]][0]
			expression_rank_in_avg_ct=cc[goodPredictedCellType[k]][1]
			#print(k,goodPredictedCellType[k],top_exp_genes_in_ct,expression_rank_in_avg_ct)
			meanCoefficients=coef[goodPredictedCellType[k]]
			stdCoefficients=coef_std[goodPredictedCellType[k]]
			highestIndex=np.argsort(-abs(meanCoefficients))

			#n=min(coeff_cutoff,len(highestIndex))
			n=len(highestIndex)
			coeff_of_CT=[]
			name_of_the_coeff=[]
			std_of_coeff=[]

			fw.write('\n'+str(k+1)+ ' Largest predicted cell type and their top 5 coefficients : '+
					nameOfCellType[goodPredictedCellType[k]]+' ( id = '+str(goodPredictedCellType[k])+',  confusion score = '+str('%0.2f'%a[goodPredictedCellType[k]])+')\n')
			xccL=[]
			xccR=[]
			xccB=[]
			xneL=[]
			xneR=[]
			xneB=[]
			rankL=[]
			rankR=[]
			rankB=[]
			globalIndex=[]
			xtickcolor=[]
			#lrplotlimitCutoff=50
			for i in range(n):
				temp=LRFeatures[highestIndex[i]]
				symbol=''
				xtickcolor.append('k')
				#print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
				#integerName=LRFeatures[highestIndex[i]].replace('x','')
				fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%meanCoefficients[ highestIndex[i]] ) +'\t'+temp+'\n')
				coeff_of_CT.append(meanCoefficients[ highestIndex[i]])
				#print(temp,cc,ne)
				for ii in range(len(top_exp_genes_in_ct)):
					if temp.upper()==top_exp_genes_in_ct[ii]:
					#if temp.upper() in top_exp_genes_in_ct:
						genetype=gene_LRB[temp.upper()]
						#print(genetype)
						if genetype=='L':
							#if len(xccL)<lrplotlimitCutoff:
								xccL.append(i)
								globalIndex.append(i)
								rankL.append(expression_rank_in_avg_ct[ii])
						if genetype=='R':
							#if len(xccR)<lrplotlimitCutoff:
								xccR.append(i)
								globalIndex.append(i)
								rankR.append(expression_rank_in_avg_ct[ii])
						if genetype=='B':
							#if len(xccB)<lrplotlimitCutoff:
								xccB.append(i)
								rankB.append(expression_rank_in_avg_ct[ii])
								globalIndex.append(i)
						#print(temp)
				if temp.upper() in ne:
					genetype=gene_LRB[temp.upper()]
					if genetype=='L':
						#if len(xneL)<lrplotlimitCutoff:
							xneL.append(i)
							xtickcolor[i]='b'
							symbol=r'$\diamondsuit$'
							globalIndex.append(i)
					if genetype=='R':
						#if len(xneR)<lrplotlimitCutoff:
							xneR.append(i)
							xtickcolor[i]='r'
							symbol=r'$\bigstar$'
							globalIndex.append(i)
					if genetype=='B':
						#if len(xneB)<lrplotlimitCutoff:
							xneB.append(i)
							xtickcolor[i]='g'
							symbol=r'$\bigcirc$'
							globalIndex.append(i)
					#print(temp)
				name_of_the_coeff.append(symbol+' '+temp)
				std_of_coeff.append(stdCoefficients[ highestIndex[i]])



			globalIndex=list(set(globalIndex))
			flag=0
			for i in range(min([coeff_cutoff,n])):
				if i not in globalIndex:
					if len(globalIndex)<coeff_cutoff:
						globalIndex.append(i)
					else:
						flag=1

			if flag==1:
				globalIndex=range(min([coeff_cutoff,n]))
				#print('flag1',len(globalIndex) , n , min([coeff_cutoff,n]))
			else:
				globalIndex=sorted(globalIndex)
				#print('flag0',len(globalIndex))

			#print(len(globalIndex),globalIndex)

			fig,ax=plt.subplots( figsize=(figuresize[0],figuresize[1]))
			xx=np.arange(len(globalIndex))#np.arange(len(coeff_of_CT))
			yy=np.zeros(len(globalIndex))

			plot_lig_rec_with_colors(ax,globalIndex,coeff_of_CT,std_of_coeff,rankL,xccL,'bD')
			plot_lig_rec_with_colors(ax,globalIndex,coeff_of_CT,std_of_coeff,rankR,xccR,'r*')
			plot_lig_rec_with_colors(ax,globalIndex,coeff_of_CT,std_of_coeff,rankB,xccB,'go')

			coeff_of_CT=np.array(coeff_of_CT)
			std_of_coeff=np.array(std_of_coeff)
			name_of_the_coeff=np.array(name_of_the_coeff)
			xtickcolor=np.array(xtickcolor)
			#print('a',len(coeff_of_CT),len(xx),len(std_of_coeff),n,len(globalIndex))

			ax.errorbar(coeff_of_CT[globalIndex],xx, xerr=std_of_coeff[globalIndex],fmt='o',capsize=5)
			ax.plot(yy,xx,'k-',linewidth=0.2)


			ax.xaxis.tick_top()
			titlename=nameOfCellType[goodPredictedCellType[k]]+', id = '+str(goodPredictedCellType[k])+', conf. = {0:.3f}'.format(a[goodPredictedCellType[k]]) +'$\pm$'+str('%0.3f'%b[goodPredictedCellType[k]])
			ax.set_title(titlename,fontsize=7)
			ax.set_yticks(xx)
			ax.set_yticklabels(name_of_the_coeff[globalIndex],fontsize=5)

			for ticklabel, tickcolor in zip(ax.get_yticklabels(),xtickcolor[globalIndex]):
				ticklabel.set_color(tickcolor)

			ax.set_ylim(ax.get_ylim()[::-1])

			fig.tight_layout()
			fig.savefig(filename+'/Rank'+str(k+1)+'_'+nameOfCellType[goodPredictedCellType[k]],dpi=300)
			fig.clf()



def run_logistic_regression_on_communication_features(inputdir,n_repeats,K_fold,coeff_cutoff,strategy,seed,n_jobs,lambda_c_ranges,BothLinearAndCrossTerms,radius,figuresize):
		#mypath='/home/ext/gruenlab3/neighbor_analysis/ScienceJeffrey2018/'

		f=open(inputdir+'BiologicalNameOfCT.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[0])]=l[1]

		if BothLinearAndCrossTerms==1:
			maindir=inputdir+strategy+'_linear/'
		else:
			maindir=inputdir+strategy+'_cross/'

		create_directory(maindir)


		#top_expressed_genes=50 # Put 0 if you want to check in full transcriptome; this is irrelevant if you deselect the db
		ligand_receptor_check_over_db=0  # 1 means check in the lig rec db, 0 means use the full transcriptome

		if True:
				start_time = time.time()
				fw=open(maindir+'prediction_R'+ str(radius)+'.dat','w')
				fw.write('\nRadius = '+ str(radius)+   '\n')
				neighborhoodClass,target,LRFeatures,cc,ne,LRpairs=read_communication_processed_data(radius,BothLinearAndCrossTerms,ligand_receptor_check_over_db,inputdir)
				#print(target[0:5],neighborhoodClass[0,0:5])
				#cmn,coef,classes,x_test,x_train,predicted_probs,fpr, tpr, roc_auc,scorecalc=model_log_regression(no_of_times_to_run_logistic_regression,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms)
				cmn,coef,comp_score,cmn_std,coef_std,comp_score_std,classes,lambda_c,x_test,x_train,predicted_probs,fpr, tpr, roc_auc=model_log_regression(K_fold, n_repeats,neighborhoodClass,target,lambda_c_ranges,strategy,BothLinearAndCrossTerms,seed,n_jobs)

				np.savetxt(maindir+'matrix_avg_coefficients_R'+str(radius)+'.dat', coef,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_avg_confusion_R'+str(radius)+'.dat', cmn,fmt='%0.6f',delimiter=',')
				score=np.array([comp_score, comp_score_std]).T
				np.savetxt(maindir+'matrix_std_coefficients_R'+str(radius)+'.dat', coef_std,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_std_confusion_R'+str(radius)+'.dat', cmn_std,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_score_R'+str(radius)+'.dat',score ,fmt='%0.4f',delimiter=',')
				#np.savetxt(maindir+'matrix_Features'+str(radius)+'.dat',CTFeatures ,fmt='%s',delimiter=',')
				np.savez(maindir+'save_numpy_array_'+str(radius)+'.npz',cmn=cmn,coef=coef,cmn_std=cmn_std,coef_std=coef_std)

				find_interacting_LR(cmn,coef,cmn_std,coef_std,coeff_cutoff,cc,ne,LRpairs,LRFeatures,nameOfCellType,fw,maindir+'/TopCoeff_R'+str(radius),figuresize)
				plot_results(maindir,2,3,nameOfCellType,radius,lambda_c,cmn,coef,classes,0,x_test,x_train,predicted_probs,BothLinearAndCrossTerms,inputdir,fpr, tpr, roc_auc)
				finish_time=time.time()
				fw.write('\n\nTotal time to compute = '+ str(finish_time-start_time)+'\n')


def plot_evaluation_scores(inputRadius,BothLinearAndCrossTerms,inputdir,strategy,figuresize):
    if BothLinearAndCrossTerms==1:
        maindir=inputdir+strategy+'_linear/'
    else:
        maindir=inputdir+strategy+'_cross/'

    ## The order of all 15 scores are following
    #1-4 'accuracy','macro F1','macro precision','macro recall',
    #5-7 'micro [all]',
    #8-11 'weighted F1','weighted precision','weighted recall','cohen kappa',
    #12=15 'cross entropy', 'matthew correlation coefficient','heming loss', 'zero one loss'

    xlabels=['accuracy','macro F1','macro precision','macro recall','micro [all]','weighted F1','weighted precision','weighted recall','cohen kappa','mcc']
    index=[0,1,2,3,4,7,8,9,10,12]

    for radius in inputRadius:
        fig,axs=plt.subplots(1,1,figsize=(figuresize[0],figuresize[1]))
        name=maindir+'matrix_score_R'+str(radius)+'.dat'
        data=np.genfromtxt(open(name, "rb"), delimiter=',', skip_header=0)
        yt=data[index,0]
        xt=range(len(yt))
        axs.plot(xt,yt,'b.-',label=strategy)
        lowRg=yt-data[index,1]
        highRg=yt+data[index,1]
        axs.fill_between(xt, lowRg, highRg,facecolor='b',alpha=0.2)
        legend1= axs.legend(loc='lower left',bbox_to_anchor=(0.0, 0.96),ncol=1, borderaxespad=0., prop={"size":6},fancybox=True, shadow=True)

        axs.set_xticks(range(len(index)))
        ytstd=max(data[index,1])
        axs.set_yticks(np.linspace(min(yt)-ytstd,max(yt)+ytstd,4))

        axs.set_xticklabels(xlabels)
        for tick in axs.get_xticklabels():
            tick.set_rotation(90)

        axs.set_ylabel('score')
        fig.tight_layout()
        fig.savefig(maindir+'scores_'+str(radius)+'.png',bbox_inches='tight',dpi=300)
        fig.clf()


def main():
    #check your other inputs for optimization
    seed=36851234
    inputRadius=[25]#,150,200,250,300]
    inputdir='communication/'#'communicationLR120/'
    figuresize1=[3,6] #First and second value is the width and height of the figure
    #This figure size use for log reg coeff vs cell type features in processed_glioblastoma1/spatial/L2_multi_linear/TopCoeff_R

    n_jobs=1 #no of processors For details see here https://scikit-learn.org/stable/glossary.html#term-n_jobs
    lambda_c_ranges=list(np.power(2.0, np.arange(-12, 12))) # inverse of lambda regularization
    BothLinearAndCrossTerms=1# If only linearterms then put 1; For both linear and crossterms use 2
    K_fold=5 #no of cross folds
    n_repeats=1 #no of times to run the logistic regression after finding the hyperparameters
    coeff_cutoff=50 # No. of coefficient want to print in the figure
    strategy='L2_multi' #Name of the strategy you want to compute the interactions options are [L1_multi, L1_ovr, L2_multi, L2_ovr, elasticnet_multi, elasticnet_ovr]
    #strategy='L1_ovr'
    #strategy='elasticnet_multi'
    outputFolder=inputdir
    for i in range(len(inputRadius)):
        print('\nRadius ',inputRadius[i],'\n')
        run_logistic_regression_on_communication_features(inputdir,n_repeats,K_fold,coeff_cutoff,strategy,seed,n_jobs,lambda_c_ranges,BothLinearAndCrossTerms,inputRadius[i],figuresize1)

    figuresize2=[4,3] #This figure size use for plotting of scores
    plot_evaluation_scores(inputRadius,BothLinearAndCrossTerms,outputFolder,strategy,figuresize2)


main()
