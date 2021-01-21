#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 22:50:49 2020

@author: lily
"""
import pandas as pd
import os
import import_data
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
## import data 
x_train, y_train, x_test, y_test = import_data.load_data()

# divide the original training data into taining set and validation set,80%-training, 20%-validation
index = np.arange(0,x_train.shape[0])
train_length = int(0.8 *index.shape[0])
train_indices = index[:train_length]
vali_indices = index[train_length:]
SVM_x_train = x_train[train_indices]
SVM_y_train = y_train[train_indices]
SVM_x_vali = x_train[vali_indices]
SVM_y_vali = y_train[vali_indices]

#1(a)training for different C
#training for different C
# from sklearn import svm
from sklearn.svm import LinearSVC

C_values = np.array([10**-4,10**-3,10**-2,10**-1, 10**0,10**1,10**2,10**3,10**4])
train_acc = np.zeros((len(C_values),1))
vali_acc = np.zeros((len(C_values),1))
test_acc = np.zeros((len(C_values),1))
n_support_vec = np.zeros((len(C_values),1))
for i in range(0,len(C_values)):
    linear_svc = LinearSVC(C= C_values[i], max_iter= 200)
#     linear_svc = svm.SVC(kernel = 'linear', C= C_values[i], max_iter= 200)
    linear_svc.fit(SVM_x_train, SVM_y_train)
    train_acc[i] = linear_svc.score(SVM_x_train, SVM_y_train)
#     n_support_vec[i] = sum(linear_svc.n_support_) 
    
    vali_acc[i] = linear_svc.score(SVM_x_vali, SVM_y_vali)
    
    test_acc[i] = linear_svc.score(x_test, y_test)
    
   # print timestamp
    timestamp = datetime.now()    
    print("progress:", int(100*(i+1)/len(C_values)), "%","Time Stamp:",timestamp)
# plot
plt.figure(num="Accuracy Curves for Different C Parameter")
lb1 = "Training Accuracy"
lb2 = "Validation Accuracy"
lb3 = "Testing Accuracy"

plt.plot(np.arange(-4, 5, 1), train_acc, 'x-k', label = lb1)
plt.plot(np.arange(-4, 5, 1), vali_acc, 'o-r', label = lb2)
plt.plot(np.arange(-4, 5, 1), test_acc, '*-g', label = lb3)
plt.xlabel("C's exponent")
plt.ylabel("Accuracy")
plt.legend() 
plt.savefig('1_a.png')


#1(b)    
linear_b= LinearSVC(C= 10**-2, max_iter= 200)
linear_b.fit(x_train, y_train)
train_acc_b = linear_b.score(x_train, y_train)
test_acc_b = linear_b.score(x_test, y_test)
# plot non_normalized confusion matrix
from sklearn import metrics
y_predict = linear_b.predict(x_test)
cm_test = metrics.confusion_matrix(y_test, y_predict)
print(cm_test)   

 #1（c）
from sklearn import SVC
deg_values = [2,3,4]
train_acc = np.zeros((len(deg_values),1))
vali_acc = np.zeros((len(deg_values),1))
test_acc = np.zeros((len(deg_values),1))
n_support_vec = np.zeros((len(deg_values),1))
for i in range(0,len(deg_values)):
    poly_svc = SVC(kernel = 'poly', C= 10**-2, degree=deg_values[i] ,max_iter= 200)
#     linear_svc = svm.SVC(kernel = 'linear', C= C_values[i], max_iter= 200)
    poly_svc.fit(SVM_x_train, SVM_y_train)
    train_acc[i] = linear_svc.score(SVM_x_train, SVM_y_train)
#     n_support_vec[i] = sum(linear_svc.n_support_) 
    
    vali_acc[i] = linear_svc.score(SVM_x_vali, SVM_y_vali)
    
    test_acc[i] = linear_svc.score(x_test, y_test)
 # 1(c) plot_accuracy
plt.figure(num="Accuracy Curves for Different Degrees of Polynomial Kernel")
lb1 = "Training Accuracy"
lb2 = "Validation Accuracy"
lb3 = "Testing Accuracy"

plt.plot(np.arange(2, 5, 1), train_acc, 'x-k', label = lb1)
plt.plot(np.arange(2, 5, 1), vali_acc, 'o-r', label = lb2)
plt.plot(np.arange(2, 5, 1), test_acc, '*-g', label = lb3)
plt.xlabel("Degrees")
plt.ylabel("Accuracy")
plt.legend() 
plt.savefig('1_c.png')   
 # 1(c) plot_number of vectors
plt.figure(num="Support Vectors Number for Different Degrees of Polynomial Kernel")
plt.plot(np.arange(2, 5, 1), n_support_vec, 'x-k', label = lb1)
plt.xlabel("Degrees")
plt.ylabel("Number of Support Vectors")
plt.legend() 
plt.savefig('1_c_2.png')


##3. 
 
def load_data():
    base_dir = os.getcwd()
    path = os.path.join(base_dir, "breast_cancer/wdbc.data")
    df_data = pd.read_csv(path, header = None)  
    ar_data = np.zeros((df_data.shape[0], df_data.shape[1]-2)) # transform dataframe to array
    ar_class = np.zeros((df_data.shape[0],1))
    for i in range(len(df_data)):
        ar_data[i] = df_data.iloc[i,2:]
        if df_data.loc[i,1]=='M':
            ar_class[i]=1        
        else:
            ar_class[i]=0            
    return  ar_data, ar_class


def data_partition(data,cls):
    train_length = int(0.7*len(data)) # 70%data for training 
    vali_length = int(0.1*len(data))  # 10% for validation
    x_train = data[0:train_length]
    x_vali = data[train_length:(train_length+vali_length)]
    x_test = data[(train_length+vali_length): ]
    y_train = cls[0:train_length]
    y_vali = cls[train_length:(train_length+vali_length): ]
    y_test = cls[(train_length+vali_length): ]
    return x_train, x_vali, x_test,  y_train, y_vali, y_test

ar_data,ar_class = load_data()
x_train, x_vali, x_test,  y_train, y_vali, y_test = data_partition(ar_data,ar_class)


# entropy of a feature
# get the unique values of a feature
def entropy_f(feature, labels):  
    # uniq values of a feature
    uniq = []   
    for i in range(len(feature)): 
        if feature[i] not in uniq:
            uniq.append(feature[i])    
    uniq.sort()   
        
    ## conditional entropy--for a feature--imput uniq of the feature
    thresholds = []
    Nz0_x1=0 
    Nz0_x0=0 
    Nz1_x0=0 
    Nz1_x1=0 
    min_entro_z = 9999
    min_threshold =0
    for val in range(len(uniq)-1):
        thresholds.append(uniq[val]+(uniq[val+1]-uniq[val])/2)
    for z in range(len(thresholds)):
        th_temp = thresholds[z]
        for idx in range(len(feature)):            
            if feature[idx] >=th_temp:
                if labels[idx] ==1:
                    Nz0_x1 +=1
                else:
                    Nz0_x0 +=1                    
            else:
                if labels[idx] ==1:
                    Nz1_x1 +=1
                else:
                    Nz1_x0 +=1         
        Pz0 = (Nz0_x1+Nz0_x0)/len(labels)
        Pz1 = (Nz1_x1+Nz1_x0)/len(labels)
        if (Nz0_x0+Nz0_x1)==0:
            Pz0_x0=0
            Pz0_x1=0
        elif (Nz1_x0+Nz1_x1)==0:
            Pz1_x0=0
            Pz1_x1=0
        else:  
            Pz0_x0 = Nz0_x0/(Nz0_x0+Nz0_x1)
            Pz0_x1 = Nz0_x1/(Nz0_x0+Nz0_x1)
            Pz1_x0 = Nz1_x0/(Nz1_x0+Nz1_x1)
            Pz1_x1 = Nz1_x1/(Nz1_x0+Nz0_x1)
        if Pz0_x0==0 or Pz0_x1==0:
            Hx_z0 = 0
            Hx_z1 = -Pz1_x0*np.log(Pz1_x0)- Pz1_x1*np.log(Pz1_x1)
       
        elif Pz1_x0==0 or Pz1_x1==0:
            Hx_z1 = 0
            Hx_z0 = -Pz0_x0*np.log(Pz0_x0)- Pz0_x1*np.log(Pz0_x1)
       
        else:
            Hx_z0 = -Pz0_x0*np.log(Pz0_x0)- Pz0_x1*np.log(Pz0_x1)
            Hx_z1 = -Pz1_x0*np.log(Pz1_x0)- Pz1_x1*np.log(Pz1_x1)        
        entro_z = Pz0*Hx_z0 + Pz1* Hx_z1
        if entro_z < min_entro_z:
            min_entro_z = entro_z
            min_threshold = th_temp 
    return  min_entro_z, min_threshold, uniq

##test entropy function

min_entro_z, min_threshold,uniq = entropy_f(x_train[:,21], y_train)
print(min_entro_z, min_threshold,uniq)
        
  