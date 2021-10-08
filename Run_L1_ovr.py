

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
from sklearn.model_selection import RepeatedStratifiedKFold

#Metrics
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import hamming_loss
from sklearn.metrics import log_loss
from sklearn.metrics import zero_one_loss
from sklearn.metrics import matthews_corrcoef


import pandas as pd
import numpy as np
import seaborn as snn
import os
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



def plot_results(maindir,nrow,ncol,nameOfCellType,radius,lambda_c,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,inputFeatures,inputdir,fpr, tpr, roc_auc):
	filename='R'+str(radius)+'_C'+str(lambda_c)
	#print(nameOfCellType)
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


	plt.figure(figsize=(8,6))
	snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classes, annot_kws={"size": 5},yticklabels=classes)
	plt.xlabel('Predicted classes')
	plt.ylabel('Truth classes')
	plt.title('R = '+str(radius)+', C='+str(lambda_c))
	plt.tight_layout()
	plt.savefig(maindir+'Confusing_matrix_'+filename+'.png',dpi=300)

	plt.figure(figsize=(5,8))
	#plt.figure()
	#snn.set(font_scale=0.4)
	b=snn.heatmap(coef.transpose(),yticklabels=CTFeatures,xticklabels=classes)
	#plt.xticks(rotation=90)
	_, ylabels= plt.yticks()
	b.set_yticklabels(ylabels, size = 5)

	plt.ylabel('Features cross terms')
	plt.xlabel('# of classes (no of cell types)')
	plt.title('R = '+str(radius)+', C='+str(lambda_c))
	plt.tight_layout()
	plt.savefig(maindir+'weight_matrix_'+filename+'.png',dpi=300)

	plt.figure(figsize=(12,6))

	plt.subplot(1,3,1)
	snn.heatmap(x_train,xticklabels=inputFeatures)
	plt.xlabel('# of input Features')
	plt.title('training set')
	plt.ylabel('75% of data')

	plt.subplot(1,3,2)
	snn.heatmap(x_test,xticklabels=inputFeatures)
	plt.xlabel('# of input Features')
	plt.title('testing set')
	plt.ylabel('25% of data')

	plt.subplot(1,3,3)
	snn.heatmap(predicted_probs,xticklabels=classes)
	plt.title('Predicted probability')
	plt.xlabel('# of classes (no of cell types)')
	plt.tight_layout(.5)
	plt.savefig(maindir+'predicted_probability_'+filename+'.png')

	plt.close('all')
	#print(predicted_probs)
	#prob=sigmoid( np.dot([y_train, y_test,1], log_reg_model.coef_.T) + log_reg_model.intercept_ )
	#print(prob)


def readdata(radius,inputdir):

	name=inputdir+'neighbors_logistic_regression_normalized_'+str(radius)+'.dat'
	#data = np.loadtxt(name, delimiter=',', usecols=(), unpack=True)
	#data = np.loadtxt(open(path_to_data, "rb"), delimiter=",", skiprows=1, usecols=np.arange(1,n))
	data1 = np.genfromtxt(open(name, "rb"), delimiter='\t', skip_header=0)
	ind=~np.isnan(data1).any(axis=1)
	data=data1[ind,:]

	prop={}
	for i in range(len(data)):
	    mytype=data[i,1]
	    if mytype in prop:
	        prop[mytype]+=1
	    else:
	        prop[mytype]=1

	#print('cell type proportion')
	total=sum(prop.values())
	keys=sorted( list( prop.keys())  )
	#classweight={}
	#for key in keys:
	    #print(key,'\t', 100.0*prop[key]/total)
	#    classweight[key]=(prop[key]/total)

	nct=len(prop)
	featureVector=range(2,2+nct) # #just neighborhood
	#filename='justLR'; featureVector=range(2+nct, 2+2*nct) #     2+2*nct) # just LR pairs
	#filename='both'; featureVector=range(2, 2+2*nct) #        2+2*nct) # all
	neighborhoodClass= data[:,featureVector]
	#mu=np.sum(neighborhoodClass,axis=1)

	#for i in range(len(mu)):
	#	neighborhoodClass[i]=neighborhoodClass[i]/mu[i]
	#mu=np.sum(neighborhoodClass,axis=1)
	#print('mu',mu)

	target= data[:,1]
	print('data shape',data.shape, target.shape, "neighbor shape",neighborhoodClass.shape)
	inputFeatures=range(nct)


	return neighborhoodClass,target,inputFeatures





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
for i in range(len(target)):
    r=random.randint(0,4)
    if target[i]<4:
        target[i]=r
    elif (target[i]>=4)&(target[i]<=10):
        target[i]=r
    else:
        target[i]=r
'''


def model_log_regression(K_fold,n_repeats,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms):
	polynomial = PolynomialFeatures(degree = BothLinearAndCrossTerms, interaction_only=True, include_bias=False)
	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#log_reg_model = LogisticRegression(max_iter=500,penalty='elasticnet',l1_ratio=0.5,class_weight=classweight,solver='saga')

	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#pipe=Pipeline([('polynomial_features',polynomial), ('logistic_regression',log_reg_model)])
	#X_poly = poly.fit_transform(neighborhoodClass[:,[0,1,2]])
	#print(X_poly.shape,CTFeatures)
	#print(X_poly[0])
	#targetFreq=pd.Series(target).value_counts(normalize=True)


	scorecalc=[]
	for i in range(15):
		scorecalc.append([])

	sss = RepeatedStratifiedKFold(n_splits=K_fold, n_repeats=n_repeats ,random_state=36851234)
	#sss = StratifiedShuffleSplit( n_splits=no_of_times_to_run_logistic_regression, test_size=0.25, random_state=0)
	cmn=[]
	coef=[]
	for train_index, test_index in sss.split(neighborhoodClass,target):
		x_train,x_test=neighborhoodClass[train_index],neighborhoodClass[test_index]
		y_train,y_test=target[train_index],target[test_index]
		#cw=calculate_class_weights(y_train)

		#class_weights= list(class_weight.compute_class_weight('balanced',np.unique(y_train),y_train))
		#print(cw,class_weights)

		#trainFreq=pd.Series(y_train).value_counts(normalize=True)
		#testFreq=pd.Series(y_test).value_counts(normalize=True)
		#print('match')
		#print('index',sorted(test_index))
		#for i in range(len(targetFreq)):
		#	print(targetFreq[i], trainFreq[i], testFreq[i])

		#weights={}
		#for index,weight in enumerate(class_weights):
		#	weights[index]=weight

		if strategy=='L1_multi':
			log_reg_model = LogisticRegression(max_iter=50000,C=lambda_c,penalty='l1',multi_class='multinomial',class_weight='balanced',solver='saga')#very slow
		if strategy=='L1_ovr':
			log_reg_model = LogisticRegression(max_iter=50000,C=lambda_c,penalty='l1',multi_class='ovr',class_weight='balanced',solver='liblinear')
		if strategy=='L2_multi':
			log_reg_model = LogisticRegression(max_iter=50000,C=lambda_c,penalty='l2',multi_class='multinomial',class_weight='balanced',solver='lbfgs')
		if strategy=='L2_ovr':
			log_reg_model = LogisticRegression(max_iter=50000,C=lambda_c,penalty='l2',multi_class='ovr',class_weight='balanced',solver='lbfgs')
		if strategy=='elasticnet_multi':
			log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='elasticnet',multi_class='multinomial',l1_ratio=0.5,class_weight='balanced',solver='saga')
		if strategy=='elasticnet_ovr':
			log_reg_model = LogisticRegression(max_iter=5000,C=lambda_c,penalty='elasticnet',multi_class='ovr',l1_ratio=0.5,class_weight='balanced',solver='saga')


		pipe=Pipeline([  ('polynomial_features',polynomial),   ('StandardScaler',StandardScaler()), ('logistic_regression',log_reg_model)])

		#pipe=Pipeline([    ('StandardScaler',StandardScaler()), ('logistic_regression',log_reg_model)])


		#pipe2=StandardScaler()



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

		'''
		pipe2.fit(x_train)

		print('\nxtrain',pipe2.mean_[0:5])
		print('S',S.mean_[0:5])

		print('\npipe2 train\t',np.std(pipe2.transform(x_train),axis=0)[0:10]   )
		print('pipe 2 test\t',np.std(pipe2.transform(x_test),axis=0)[0:10]   )
		print('pipe train\t',np.std(S.transform(x_train),axis=0)[0:10]   )
		print('pipe test\t',np.std(S.transform(x_test),axis=0)[0:10]   )
        	'''

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

		poly = pipe.named_steps['polynomial_features']
		LR= pipe.named_steps['logistic_regression']
		coef.append(LR.coef_)
		cmn.append(confusion_matrix(y_test,y_pred,normalize='true'))

	cmn_std=np.std(np.array(cmn),axis=0)
	coef_std=np.std(np.array(coef),axis=0)
	comp_score_std=np.std(np.array(scorecalc),axis=1)


	cmn=np.mean(np.array(cmn),axis=0)
	coef=np.mean(np.array(coef),axis=0)
	comp_score=np.mean(np.array(scorecalc),axis=1)

	#print('a',comp_score)
	#print('b',comp_score_std)


	print('training',x_train.shape,'testing',x_test.shape,'coeff',coef.shape,'Iteration',len(scorecalc[0]))


	#cmn=confusion_matrix(y_test,y_pred,normalize='true')
	#cmn = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
	classes=LR.classes_.astype(int)


	#modifiedFeatures=range(1,len(CTFeatures)+1)
	fpr, tpr, roc_auc=plot_multiclass_roc(pipe, x_test, y_test, n_classes=len(classes))


	CTFeatures=poly.get_feature_names()
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




	return cmn,coef,comp_score,cmn_std,coef_std,comp_score_std,classes, CTFeatures,x_test,x_train,y_prob ,fpr, tpr, roc_auc



def find_interacting_cell_types(cmn,coef,cmn_std,coef_std,confusion_cutoff,coeff_cutoff,CTFeatures,nameOfCellType,fw,filename):


	answer=os.path.isdir(filename)
	if answer==True:
		pass
	else:
		os.mkdir(filename)

	a=np.diag(cmn)
	b=np.diag(cmn_std)
	goodPredictedCellType=np.argsort(-a)


	#for i in range(len(goodPredictedCellType)):
	#	print(i,a[goodPredictedCellType[i]])

	# top 3 cell type in confusion matrix
	for k in range(len(a)):
		if a[goodPredictedCellType[k]]>=confusion_cutoff:
			meanCoefficients=coef[goodPredictedCellType[k]]
			stdCoefficients=coef_std[goodPredictedCellType[k]]
			highestIndex=np.argsort(-abs(meanCoefficients))

			n=min(coeff_cutoff,len(highestIndex))
			coeff_of_CT=[]
			name_of_the_coeff=[]
			std_of_coeff=[]

			fw.write('\n'+str(k+1)+ ' Largest predicted cell type and their top 5 coefficients : '+
					nameOfCellType[goodPredictedCellType[k]]+' ( id = '+str(goodPredictedCellType[k])+',  confusion score = '+str('%0.2f'%a[goodPredictedCellType[k]])+')\n')
			for i in range(n):
			#for i in range(len(highestIndex)):
				l=CTFeatures[highestIndex[i]].split()
				temp=''
				for j in range(len(l)):
					temp+=nameOfCellType[int(l[j][1:])]
					if j!=(len(l)-1):
						temp+='--'
				#print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
				integerName=CTFeatures[highestIndex[i]].replace('x','')
				fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%meanCoefficients[ highestIndex[i]] ) +'\t'+temp+' ('+ integerName  +')\n')
				coeff_of_CT.append(meanCoefficients[ highestIndex[i]])
				name_of_the_coeff.append(temp)
				std_of_coeff.append(stdCoefficients[ highestIndex[i]])


			fig,ax=plt.subplots( figsize=(5,3))
			xx=range(len(coeff_of_CT))
			yy=np.zeros((len(coeff_of_CT)))
			ax.errorbar(xx, coeff_of_CT, yerr=std_of_coeff,fmt='o',capsize=5)
			ax.plot(xx,yy,'k-',linewidth=0.2)
			ax.set_ylabel('value of coeff.')
			ax.set_xlabel('name of the coeff.')
			titlename=nameOfCellType[goodPredictedCellType[k]]+', id = '+str(goodPredictedCellType[k])+', confusion score = {0:.3f}'.format(a[goodPredictedCellType[k]]) +'$\pm$'+str('%0.3f'%b[goodPredictedCellType[k]])
			ax.set_title(titlename,fontsize=7)


			ax.set_xticks(xx)
			ax.set_xticklabels(name_of_the_coeff)
			for tick in ax.get_xticklabels():
				tick.set_rotation(90)


			fig.tight_layout()
			fig.savefig(filename+'/Rank'+str(k)+'_'+nameOfCellType[goodPredictedCellType[k]],dpi=300)
			fig.clf()



def main():
		#inputpdir='/home/ext/gruenlab3/neighbor_analysis/ScienceJeffrey2018/'
		inputpdir='./../'
		inputdir='./'  #inputpdir+'4expected/'

		'''
		f=open('./../cellTypeIntegerName_cortex.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[1])]=l[0]
		'''

		f=open(inputdir+'BiologicalNameOfCT.dat')
		nameOfCellType={}
		for line in f:
		    l=line[0:-1].split('\t')
		    nameOfCellType[int(l[0])]=l[1]

		#no_of_times_to_run_logistic_regression=1
		n_repeats=20
		K_fold=5
		confusion_cutoff=0.2
		coeff_cutoff=25
		strategy='L1_ovr'
		#strategy='L2_multi'
		#strategy='elasticnet_multi'

		lambda123=[1]
		#radiusList=[50,100,150,200,250,300]
		radiusList=[25,50,75,100]

		BothLinearAndCrossTerms=2# if both linear and crossterms then 2, if only linearterm then 1









		if BothLinearAndCrossTerms==1:
			maindir=inputdir+'LRF_'+strategy+'_linear/'
		else:
			maindir=inputdir+'LRF_'+strategy+'_cross/'
		answer=os.path.isdir(maindir)
		if answer==True:
			pass
		else:
			os.mkdir(maindir)	
		for i in radiusList:
			print('\n\n\n')
			for j in lambda123:
				start_time = time.time()
				radius=i
				lambda_c=j
				fw=open(maindir+'prediction_R'+str(radius)+'_C'+ str(lambda_c)+'.dat','w')
				print('\n\nRadius and c',radius,lambda_c,'\n')
				fw.write('\nRadius = '+ str(radius)+  '\t lambda C = '+ str(lambda_c) + '\n')
				neighborhoodClass,target,inputFeatures=readdata(radius,inputdir)
				cmn,coef,comp_score,cmn_std,coef_std,comp_score_std,classes,CTFeatures,x_test,x_train,predicted_probs,fpr, tpr, roc_auc=model_log_regression(K_fold, n_repeats,neighborhoodClass,target,lambda_c,strategy,BothLinearAndCrossTerms)

				np.savetxt(maindir+'matrix_avg_coefficients_R'+str(radius)+'_C'+ str(lambda_c)+'.dat', coef,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_avg_confusion_R'+str(radius)+'_C'+ str(lambda_c)+'.dat', cmn,fmt='%0.6f',delimiter=',')

				score=np.array([comp_score, comp_score_std]).T

				np.savetxt(maindir+'matrix_std_coefficients_R'+str(radius)+'_C'+ str(lambda_c)+'.dat', coef_std,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_std_confusion_R'+str(radius)+'_C'+ str(lambda_c)+'.dat', cmn_std,fmt='%0.6f',delimiter=',')
				np.savetxt(maindir+'matrix_score_R'+str(radius)+'_C'+ str(lambda_c)+'.dat',score ,fmt='%0.4f',delimiter=',')

				find_interacting_cell_types(cmn,coef,cmn_std,coef_std,confusion_cutoff,coeff_cutoff,CTFeatures,nameOfCellType,fw,maindir+'/TopCoeff_R'+str(radius)+'_C'+str(lambda_c))

				plot_results(maindir,2,3,nameOfCellType,radius,lambda_c,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,inputFeatures,inputdir,fpr, tpr, roc_auc)
				finish_time=time.time()

				fw.write('\n\nTotal time to compute = '+ str(finish_time-start_time)+'\n')

main()
