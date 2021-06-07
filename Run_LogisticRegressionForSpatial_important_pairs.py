

from sklearn.datasets import make_classification
from matplotlib import pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, classification_report,confusion_matrix,roc_curve, roc_auc_score, precision_score, recall_score, precision_recall_curve
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import Pipeline
#from sklearn.metrics import precision_recall_fscore_support as score
from sklearn.model_selection import StratifiedShuffleSplit
#from imblearn.over_sampling import SMOTE, SMOTEN,ADASYN, KMeansSMOTE, SVMSMOTE
from sklearn.utils import class_weight


import pandas as pd
import numpy as np
import seaborn as snn
import os
import random


def plot_results(radius,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,inputFeatures,inputdir):

	filename=str(radius)

	maindir=inputdir+'LogisticRegressionFigures/'
	answer=os.path.isdir(maindir)
	if answer==True:
		pass
	else:
		os.mkdir(maindir)


	plt.figure(figsize=(8,6))
	snn.heatmap(cmn,annot=True, fmt='.2f',xticklabels=classes, annot_kws={"size": 5},yticklabels=classes)
	plt.xlabel('Predicted classes')
	plt.ylabel('Truth classes')
	plt.tight_layout()
	plt.savefig(maindir+'Confusing_matrix_'+filename+'.png')

	plt.figure(figsize=(5,8))
	#plt.figure()
	#snn.set(font_scale=0.4)
	b=snn.heatmap(coef.transpose(),yticklabels=CTFeatures,xticklabels=classes)
	#plt.xticks(rotation=90)
	_, ylabels= plt.yticks()
	b.set_yticklabels(ylabels, size = 5)

	plt.ylabel('Features cross terms')
	plt.xlabel('# of classes (no of cell types)')
	plt.tight_layout()
	plt.savefig(maindir+'weight_matrix_'+filename+'.png')

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
	#print(predicted_probs)
	#prob=sigmoid( np.dot([y_train, y_test,1], log_reg_model.coef_.T) + log_reg_model.intercept_ )
	#print(prob)


def readdata(radius,inputdir):

	name=inputdir+'neighbors_logistic_regression_normalized_'+str(radius)+'.dat'
	#data = np.loadtxt(name, delimiter=',', usecols=(), unpack=True)
	#data = np.loadtxt(open(path_to_data, "rb"), delimiter=",", skiprows=1, usecols=np.arange(1,n))
	data = np.genfromtxt(open(name, "rb"), delimiter='\t', skip_header=0)

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
	classweight={}
	for key in keys:
	    #print(key,'\t', 100.0*prop[key]/total)
	    classweight[key]=(prop[key]/total)

	#print(classweight)


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


def model_log_regression(no_of_times_to_run_logistic_regression,neighborhoodClass,target):
	polynomial = PolynomialFeatures(degree = 2, interaction_only=True, include_bias=False)
	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#log_reg_model = LogisticRegression(max_iter=500,penalty='elasticnet',l1_ratio=0.5,class_weight=classweight,solver='saga')

	#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',class_weight=classweight,solver='newton-cg')
	#pipe=Pipeline([('polynomial_features',polynomial), ('logistic_regression',log_reg_model)])
	#X_poly = poly.fit_transform(neighborhoodClass[:,[0,1,2]])
	#print(X_poly.shape,CTFeatures)
	#print(X_poly[0])

	scorecalc=[]
	for i in range(10):
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

		#weights={}
		#for index,weight in enumerate(class_weights):
		#	weights[index]=weight

		log_reg_model = LogisticRegression(max_iter=500,penalty='l2',multi_class='multinomial',class_weight='balanced',solver='newton-cg')
		#log_reg_model = LogisticRegression(max_iter=500,penalty='elasticnet',l1_ratio=0.5,class_weight=weights,solver='saga')
		#log_reg_model = LogisticRegression(max_iter=500,penalty='l2',solver='newton-cg')

		pipe=Pipeline([('polynomial_features',polynomial), ('logistic_regression',log_reg_model)])

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

		pipe.fit(x_train, y_train)
		y_pred=pipe.predict(x_test)
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
		#print(precision,recall,fscore,support)
		poly = pipe.named_steps['polynomial_features']
		LR= pipe.named_steps['logistic_regression']
		coef.append(LR.coef_)
		cmn.append(confusion_matrix(y_test,y_pred,normalize='true'))

	y_pred = pipe.predict(x_test)
	#print(y_pred.shape,y_test.shape)

	cmn=np.mean(np.array(cmn),axis=0)
	coef=np.mean(np.array(coef),axis=0)
	print('training',x_train.shape,'testing',x_test.shape,'coeff',coef.shape)


	#cmn=confusion_matrix(y_test,y_pred,normalize='true')
	#cmn = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
	classes=LR.classes_.astype(int)
	predicted_probs = pipe.predict_proba(x_test)

	#modifiedFeatures=range(1,len(CTFeatures)+1)

	CTFeatures=poly.get_feature_names()
	#print("Features", CTFeatures)
	print("accuracy score in 10 run\t",np.mean(scorecalc[0]))
	print("\n\nmacro")
	print("f1 score in 10 run\t",np.mean(scorecalc[1]))
	print("precision score in 10 run\t",np.mean(scorecalc[2]))
	print("recall score in 10 run\t",np.mean(scorecalc[3]))

	print("\n\nmicro")
	print("f1 score in 10 run\t",np.mean(scorecalc[4]))
	print("precision score in 10 run\t",np.mean(scorecalc[5]))
	print("recall score in 10 run\t",np.mean(scorecalc[6]))

	print("\n\nWeighted")
	print("f1 score in 10 run\t",np.mean(scorecalc[7]))
	print("precision score in 10 run\t",np.mean(scorecalc[8]))
	print("recall score in 10 run\t",np.mean(scorecalc[9]))

	return cmn,coef,classes, CTFeatures,x_test,x_train,predicted_probs



def find_interacting_cell_types(cmn,coef,CTFeatures,nameOfCellType,fw):

	a=np.diag(cmn)
	goodPredictedCellType=np.argmax(a)
	goodCoefficients=coef[goodPredictedCellType]
	#print(a,goodPredictedCellType,goodCoefficients.shape)
	highestIndex=np.argsort(-abs(goodCoefficients))


	fw.write('Top 5 coefficient used for predicted cell type : '+
				nameOfCellType[goodPredictedCellType]+' ( int id = '+str(goodPredictedCellType)+',  confusion score = '+str('%0.2f'%a[goodPredictedCellType])+')\n')
	for i in range(5):
	#for i in range(len(highestIndex)):
		l=CTFeatures[highestIndex[i]].split()

		temp=''
		for j in range(len(l)):
			temp+=nameOfCellType[int(l[j][1:])]
			if j!=(len(l)-1):
				temp+='--'

		#print(temp,highestIndex[i],CTFeatures[highestIndex[i]],goodCoefficients[ highestIndex[i]   ])
		integerName=CTFeatures[highestIndex[i]].replace('x','')
		fw.write(str(highestIndex[i])+'\t'+str('%0.2f'%goodCoefficients[ highestIndex[i]] ) +'\t'+temp+' ('+ integerName  +')\n')




def main():
		inputpdir='/home/ext/gruenlab3/neighbor_analysis/ScienceJeffrey2018/'
		#inputpdir='./../'
		inputdir=inputpdir+'4expected/'

		#25 is pretty standard
		no_of_times_to_run_logistic_regression=25
		'''
		f=open('./../cellTypeIntegerName_cortex.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[1])]=l[0]
		'''

		f=open(inputpdir+'BiologicalNameOfCT.dat')
		nameOfCellType={}
		for line in f:
		    l=line.split('\t')
		    nameOfCellType[int(l[0])]=l[1]

		fw=open(inputdir+'prediction_of_cell_types'+str(no_of_times_to_run_logistic_regression)+'.dat','w')

		radiusList=[50,100,150,200,250,300]
		#radiusList=[50,100]
		for i in radiusList:
			radius=i
			print('\n\nRadius',radius,'\n')
			fw.write('\nRadius = '+ str(radius)+'\n')
			neighborhoodClass,target,inputFeatures=readdata(radius,inputdir)
			cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs=model_log_regression(no_of_times_to_run_logistic_regression,neighborhoodClass,target)
			find_interacting_cell_types(cmn,coef,CTFeatures,nameOfCellType,fw)
			plot_results(radius,cmn,coef,classes,CTFeatures,x_test,x_train,predicted_probs,inputFeatures,inputdir)
main()
