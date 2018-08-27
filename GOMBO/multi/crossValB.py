# These are all the modules we'll be using later. Make sure you can import them
# before proceeding further.
from __future__ import print_function
import numpy as np
import os
import sys
import tarfile
import math
import random
import sys
from funcCNN import *
# used for normalization
from sklearn.preprocessing import StandardScaler

def get_Info(indexVar,labelSize,vectorSize):

	data=[]
	data = np.genfromtxt('./data/'+'read_count.txt', delimiter=' ')
	data=np.array(data)
	print('data set', data.shape)
	
	
	labels=openVector('./data/labels.txt')
	testIndex=openVector('./data/index/'+str(indexVar)+'test_index.txt')
	valIndex=openVector('./data/index/'+str(indexVar)+'val_index.txt')
	trainIndex=openVector('./data/index/'+str(indexVar)+'train_index.txt')
	
	testIndex=testIndex.astype(int)
	valIndex=valIndex.astype(int)
	trainIndex=trainIndex.astype(int)
		
	train=[]
	test=[]
	valid=[]
	
	
	trainLabels=[]
	testLabels=[]
	validLabels=[]
	
	scaler_sample = StandardScaler()
	scaler_sample2 = StandardScaler()
	data = scaler_sample.fit_transform(data.T).T
	#test***************************************************************************
	for i in range (0,len(testIndex)):
		testLabels.append(labels[testIndex[i]])
		temp=[]
		for j in range (0,len(data[0])):
			if(data[testIndex[i]][j]==-1):
				temp.append(0)
			else:
				temp.append(data[testIndex[i]][j])
		test.append(temp)
	#valid***************************************************************************
	for i in range (0,len(valIndex)):
		validLabels.append(labels[valIndex[i]])
		temp=[]
		for j in range (0,len(data[0])):
			if(data[valIndex[i]][j]==-1):
				temp.append(0)
			else:
				temp.append(data[valIndex[i]][j])
		valid.append(temp)
	#train***************************************************************************
	for i in range (0,len(trainIndex)):
		trainLabels.append(labels[trainIndex[i]])
		temp=[]
		for j in range (0,len(data[0])):
			if(data[trainIndex[i]][j]==-1):
				temp.append(0)
			else:
				temp.append(data[trainIndex[i]][j])
		train.append(temp)
	
	test=np.array(test)
	testLabels=np.array(testLabels)

	valid=np.array(valid)
	validLabels=np.array(validLabels)
	
	train=np.array(train)
	trainLabels=np.array(trainLabels)
	
	print(train.shape)
	print(trainLabels.shape)
	print(valid.shape)
	print(validLabels.shape)
	print(test.shape)
	print(testLabels.shape)
	

	oneHot_train_labels=oneHot(trainLabels,labelSize)
	print(oneHot_train_labels.shape)

	oneHot_valid_labels=oneHot(validLabels,labelSize)
	print(oneHot_valid_labels.shape)

	oneHot_test_labels=oneHot(testLabels,labelSize)
	print(oneHot_test_labels.shape)


	return(test,oneHot_test_labels,valid,oneHot_valid_labels,train,oneHot_train_labels)





