

from sklearn.datasets import make_classification
from matplotlib import pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import Pipeline
#from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.model_selection import StratifiedShuffleSplit
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



def plot_results(maindir,nrow,ncol,nameOfCellType,radius,lambda_c,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,inputdir,fpr, tpr, roc_auc):

	#print(CTFeatures,inputFeatures,classes,nameOfCellType)
	filename=str(radius)+'_'+str(lambda_c)
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


	plt.figure(figsize=(8,6))
	snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classes, annot_kws={"size": 5},yticklabels=classes)
	plt.xlabel('Predicted classes')
	plt.ylabel('Truth classes')
	plt.title('R = '+str(radius)+', C='+str(lambda_c))
	plt.tight_layout()
	plt.savefig(maindir+'Confusing_matrix_'+filename+'.png')

	plt.figure(figsize=(5,8))
	#plt.figure()
	#snn.set(font_scale=0.4)
	b=snn.heatmap(coef.transpose(),xticklabels=classes)#yticklabels=CTFeatures
	#plt.xticks(rotation=90)
	_, ylabels= plt.yticks()
	b.set_yticklabels(ylabels, size = 5)

	plt.ylabel('Features cross terms')
	plt.xlabel('# of classes (no of cell types)')
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



def readdata(radius,top_10_expressed_genes,BothLinearAndCrossTerms,ligand_receptor_check_over_db,inputdir,filename):

	name=inputdir+'neighbors_'+filename+'_'+str(radius)+'.dat'
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

	print("3 database LR pairs =",len(totalLRpairs), '; ligands =', len(ligand_one),'; receptors =',len(receptor_one),'; combined =',len(combined_one))



	#fw=open('dummy.dat','w')
	#for key in totalLRpairs:
	#	fw.write(str(key)+'\n')

	#print(sorted(combined_one),sorted(cont))



	proteins=[]
	cc_=[]
	ne_=[]


	#check what are the top 5 gene expressed by cell types in csv file
	f=open('highest_to_lowest_expressed_genes_in_celltypes.dat')


	total_gene=[]
	for line in f:
		l=line[0:-1].split(';')

		p_ng=len(l)-1
		lig_rec_present_in_ct=[]
		if top_10_expressed_genes==0:
			top_10=p_ng
		else:
			top_10=	top_10_expressed_genes
		for j in range(1,top_10+1):#len(l)):
			gene=l[j].upper()
			total_gene.append(gene)
			# central cell is ligand
			try:
				ligand_one[gene]
				if gene not in cc_:
					cc_.append(gene)
				ne_=counterPartLigRec(gene,totalLRpairs,'L',ne_)
			except KeyError:
				pass

			# central cell is receptor
			try:
				receptor_one[gene]
				if gene not in cc_:
					cc_.append(gene)
				ne_=counterPartLigRec(gene,totalLRpairs,'R',ne_)
			except KeyError:
				pass

	print('communication channel cc', len(cc_), ', neighbor ', len(ne_), 'total number of genes',p_ng  )
	#print('x1',cc_ligand, ne_receptor)
	#print('x2',ne_ligand, cc_receptor)

	if ligand_receptor_check_over_db:
		data=ne_
	else:
		data=total_gene

	gene_index_expressed_in_neighbors=[]
	for i in range(n):
		l1=cont[i].split('\t')
		g1=l1[0]
		if g1 in data:
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



	return neighborhoodClass,target,all_LR_pairs_found




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






def model_log_regression(no_of_times_to_run_logistic_regression,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms):
	#polynomial = PolynomialFeatures(degree = 2, interaction_only=True, include_bias=False)
	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#log_reg_model = LogisticRegression(max_iter=500,penalty='elasticnet',l1_ratio=0.5,class_weight=classweight,solver='saga')

	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#pipe=Pipeline([('polynomial_features',polynomial), ('logistic_regression',log_reg_model)])
	#X_poly = poly.fit_transform(neighborhoodClass[:,[0,1,2]])
	#print(X_poly.shape,CTFeatures)
	#print(X_poly[0])

	scorecalc=[]
	for i in range(15):
		scorecalc.append([])

	sss = StratifiedShuffleSplit( n_splits=no_of_times_to_run_logistic_regression, test_size=0.25, random_state=0)
	cmn=[]
	coef=[]
	for train_index, test_index in sss.split(neighborhoodClass,target):
		x_train,x_test=neighborhoodClass[train_index],neighborhoodClass[test_index]
		y_train,y_test=target[train_index],target[test_index]
		#cw=calculate_class_weights(y_train)

		#class_weights= list(class_weight.compute_class_weight('balanced',np.unique(y_train),y_train))
		#print(cw,class_weights)

		if strategy=='L1_multi':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='l1',multi_class='multinomial',class_weight='balanced',solver='saga')#very slow
		if strategy=='L1_ovr':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='l1',multi_class='ovr',class_weight='balanced',solver='liblinear')
		if strategy=='L2_multi':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='l2',multi_class='multinomial',class_weight='balanced',solver='lbfgs')
		if strategy=='L2_ovr':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='l2',multi_class='ovr',class_weight='balanced',solver='lbfgs')
		if strategy=='elasticnet_multi':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='elasticnet',multi_class='multinomial',l1_ratio=0.5,class_weight='balanced',solver='saga')
		if strategy=='elasticnet_ovr':
				log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='elasticnet',multi_class='ovr',l1_ratio=0.5,class_weight='balanced',solver='saga')



		pipe=Pipeline([('StandardScaler',StandardScaler()), ('logistic_regression',log_reg_model)])

		#print('a',len(y_train),cw)
		#x_train, x_test, y_train, y_test=train_test_split(neighborhoodClass,target,test_size=0.25)
		#x_train = poly.fit_transform(x_train)
		#x_test= poly.transform(x_test)
		#log_reg_model.fit(x_train,y_train)
		#score.append( log_reg_model.score(x_test,y_test))
		#coef=log_reg_model.coef_

		#oversample = SMOTE(k_neighbors=3)
		#oversample = ADASYN(n_neighbors=3)
		#oversample = KMeansSMOTE(k_neighbors=2)
		#oversample = SMOTEN(k_neighbors=3)
		#oversample = SVMSMOTE(m_neighbors=5)
		#print('a',y_train.shape)
		#x_train, y_train  = oversample.fit_resample(x_train, y_train)
		#print('b',y_train.shape)

		S=pipe.named_steps['StandardScaler']
	    #print('mean', S.mean_)
	    #print('scale',S.scale_)

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



	cmn=np.mean(np.array(cmn),axis=0)
	coef=np.mean(np.array(coef),axis=0)
	print('training',x_train.shape,'testing',x_test.shape,'coeff',coef.shape)


	#cmn=confusion_matrix(y_test,y_pred,normalize='true')
	#cmn = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
	classes=LR.classes_.astype(int)


	#modifiedFeatures=range(1,len(CTFeatures)+1)
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





	return cmn,coef,classes,x_test,x_train,y_prob ,fpr, tpr, roc_auc, scorecalc



def find_interacting_LR(cmn,coef,LRFeatures,nameOfCellType,fw):

	a=np.diag(cmn)
	goodPredictedCellType=np.argsort(-a)
	#print('coeff',coef.shape)
	#print(len(CTFeatures))

	# top 3 cell type in confusion matrix
	for k in range(3):
		goodCoefficients=coef[goodPredictedCellType[k]]
		highestIndex=np.argsort(-abs(goodCoefficients))
		fw.write('\n'+str(k+1)+ ' Largest predicted cell type and their top 5 coefficients : '+
				nameOfCellType[goodPredictedCellType[k]]+' ( id = '+str(goodPredictedCellType[k])+',  confusion score = '+str('%0.2f'%a[goodPredictedCellType[k]])+')\n')
		for i in range(5):
		#for i in range(len(highestIndex)):
			temp=LRFeatures[highestIndex[i]]
			#print(temp)
			#print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
			#integerName=LRFeatures[highestIndex[i]].replace('x','')
			fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%goodCoefficients[ highestIndex[i]] ) +'\t'+temp+' ('+ 'X'  +')\n')



def main():
		#mypath='/home/ext/gruenlab3/neighbor_analysis/ScienceJeffrey2018/'
		filename='weight'
		inputdir='./'
		#25 is pretty standard
		no_of_times_to_run_logistic_regression=1
		'''
		f=open('./../cellTypeIntegerName_ob.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[1])]=l[0]
		'''
		f=open(inputdir+'BiologicalNameOfCT.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[0])]=l[1]

		no_of_times_to_run_logistic_regression=1
		#strategy='L1_ovr'
		strategy='L2_multi'
		#strategy='elasticnet_multi'

		maindir=inputdir+'LRF_'+strategy+'/'
		answer=os.path.isdir(maindir)
		if answer==True:
			pass
		else:
			os.mkdir(maindir)

		lambda123=[1,0.1,10]
		lambda123=[1]
		top_10_expressed_genes=10 # Put 0 if you want to check in full transcriptome; this is irrelevant if you deselect the db
		BothLinearAndCrossTerms=1# if both linear and crossterms then 2, if only linearterm then 1

		ligand_receptor_check_over_db=0  # 1 means check in the lig rec db, 0 means use the full transcriptome

		radiusList=[50,75,100,150]
		radiusList=[50]
		for i in radiusList:
			print('\n\n\n')
			for j in lambda123:
				start_time = time.time()
				radius=i
				lambda_c=j
				fw=open(maindir+'prediction_'+str(lambda_c)+'_'+ str(radius)+'.dat','w')
				print('\n\nRadius = ',radius,', c = ',lambda_c,'\n')
				fw.write('\nRadius = '+ str(radius)+  '\t lambda C = '+ str(lambda_c) + '\n')
				neighborhoodClass,target,LRFeatures=readdata(radius,top_10_expressed_genes,BothLinearAndCrossTerms,ligand_receptor_check_over_db,inputdir,filename)
				cmn,coef,classes,x_test,x_train,predicted_probs,fpr, tpr, roc_auc,scorecalc=model_log_regression(no_of_times_to_run_logistic_regression,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms)

				np.savetxt(maindir+'matrix_coefficients_'+str(lambda_c)+'_'+ str(radius)+'.dat', coef,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_confusion_'+str(lambda_c)+'_'+ str(radius)+'.dat', cmn,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_score_'+str(lambda_c)+'_'+ str(radius)+'.dat', scorecalc,fmt='%0.4f',delimiter=',')

				find_interacting_LR(cmn,coef,LRFeatures,nameOfCellType,fw)
				plot_results(maindir,2,3,nameOfCellType,radius,lambda_c,cmn,coef,classes,0,x_test,x_train,predicted_probs,inputdir,fpr, tpr, roc_auc)
				finish_time=time.time()
				fw.write('\n\nTotal time to compute = '+ str(finish_time-start_time)+'\n')

main()
