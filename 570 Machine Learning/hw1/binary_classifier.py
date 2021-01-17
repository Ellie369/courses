## multiclass classifier


import import_data
import numpy as np
import matplotlib.pyplot as plt

## import data 
x_train, y_train, x_test, y_test = import_data.load_data()
y_label_train, y_label_test = import_data.y_label(y_train, y_test)

##
def perceptron(images,labels, itr):
    learning_curve = [0]*itr 
    accuracy_curve = [0]*itr 
    weight_list = [0]*itr
    weight = np.zeros(len(images[0])) # initialize weight as zeros
    for i in range(itr):
        errors = 0 #number of mistakes
        for j in range(len(images)):
            y_cap = np.sign(np.dot(weight,images[j])) # prediction of yt
            if labels[j] != y_cap:
                weight = np.add(weight,np.dot(labels[j],images[j]))
                errors += 1
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(images))
        weight_list[i] = weight
    return weight,weight_list, learning_curve, accuracy_curve

def aver_perceptron(images,labels, itr):
    learning_curve = [0]*itr 
    accuracy_curve = [0]*itr 
    aver_weight_list = [0]*itr
    weight = np.zeros(len(images[0])) # initialize weight as zeros
    weight_aver = np.zeros(len(images[0]))
    for i in range(itr):
        errors = 0 #number of mistakes
        for j in range(len(images)):
            y_cap = np.sign(np.dot(weight,images[j])) # prediction of yt
            if labels[j] != y_cap:
                weight = np.add(weight,np.dot(labels[j],images[j]))
                errors += 1
            weight_aver += weight 
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(images))
        aver_weight_list[i] = weight_aver
    return weight, aver_weight_list, learning_curve, accuracy_curve

def PA(images,labels, itr):
    learning_curve = [0]*itr 
    accuracy_curve = [0]*itr 
    weight_list = [0]*itr
    weight = np.zeros(len(images[0])) # initialize weight as zeros   
    for i in range(itr):
        errors = 0 #number of mistakes
        for j in range(len(images)):
            y_cap = np.sign(np.dot(weight,images[j])) # prediction of yt
            if labels[j] != y_cap:
                t = (1-labels[j]*np.dot(weight,images[j]))/(np.linalg.norm(images[j]))**2 #learning rate
                weight = np.add(weight,t*np.dot(labels[j],images[j]))
                errors +=1
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(images))
        weight_list[i] = weight
    return weight, weight_list, learning_curve, accuracy_curve

def test(weight, image_test, label_test):
    right = 0   
    for t in range(len(image_test)):
        y_cap = np.sign(np.dot(weight,image_test[t]))
        if label_test[t]==y_cap:
            right +=1         
    accuracy = (right/len(image_test))
    return accuracy


def test_accu_curve(weight_list, x_test, y_test):
    test_accuracy_curve =[]
    for wt in range(weight_list):
        accu = test(wt,x_test, y_test)
        test_accuracy_curve.append(accu)
    return test_accuracy_curve
    



## 5.1a 
[w_train_pt,wlist_train_pt, learning_curve_pt,accuracy_curve_pt] = perceptron(x_train,y_label_train, 50)
[w_train_PA,wlist_train_PA, learning_curve_PA,accuracy_curve_PA] = PA(x_train,y_label_train, 50)

plt.figure(num="Binary Classifier Online Learning Curves")
lb1 = "Perceptron Algorithm"
lb2 = "Passive Aggressive Algorithm"
plt.plot(np.arange(1, 51, 1), learning_curve_pt, 'x-k', label = lb1)
plt.plot(np.arange(1, 51, 1), learning_curve_PA, 'o-r', label = lb2)
plt.xlabel("Number of Training Iterations")
plt.ylabel("Number of Mistakes")
plt.legend() 

## 5.1b
[w_train_pt,wlist_train_pt, learning_curve_pt,accuracy_curve_pt] = perceptron(x_train,y_label_train, 20)
[w_train_PA,wlist_train_PA, learning_curve_PA,accuracy_curve_PA] = PA(x_train,y_label_train, 20)
test_accu_pt = test_accu_curve(wlist_train_pt,x_test,y_label_test)
test_accu_PA = test_accu_curve(wlist_train_PA,x_test,y_label_test)

plt.figure(num="Online Learning Accuracy Curves")
lb1 = "Perceptron-Training"
lb2 = "PA-Training"
lb3 = "Perceptron-Testing"
lb4 = "PA-Testing"
plt.plot(np.arange(1, 21, 1), accuracy_curve_pt, 'x-b', label = lb1)
plt.plot(np.arange(1, 21, 1), accuracy_curve_PA, 'o-r', label = lb2)
plt.plot(np.arange(1, 21, 1), test_accu_pt, 'd-g', label = lb3)
plt.plot(np.arange(1, 21, 1), test_accu_PA, '*-y', label = lb4)
plt.xlabel("Number of Training Iterations")
plt.ylabel("Accuracy")
plt.legend() 

## 5.1c

[w_train_pt,wlist_train_pt, learning_curve_pt,accuracy_curve_pt] = perceptron(x_train,y_label_train, 20)
test_accu_pt = test_accu_curve(wlist_train_pt,x_test,y_label_test)
[aver_w_train_pt,aver_wlist_train_pt, learning_curve_pt,accuracy_curve_pt] = aver_perceptron(x_train,y_label_train, 20)
aver_test_accu_pt = test_accu_curve(aver_wlist_train_pt,x_test,y_label_test)

plt.figure(num="Accuracy Curves of Plain Perceptron and Averaged Perceptron-Binary Classifier")
lb1 = "Plain Perceptron"
lb2 = "Averaged Perceptron"
plt.plot(np.arange(1, 21, 1), test_accu_pt, 'x-b', label = lb1)
plt.plot(np.arange(1, 21, 1),aver_test_accu_pt, 'o-r', label = lb2)
plt.xlabel("Number of Iterations")
plt.ylabel("Accuracy")
plt.legend() 

# 5.1d

weight_num=[]
for num in range(20):    
    train_length = int(100+100*num)    
    x_train_num = x_train[0:train_length]
    y_train_num = y_label_train[0:train_length]
    [w_train_pt,wlist_train_pt, mistake_train_pt,accuracy_train_pt] = perceptron(x_train_num,y_train_num, 20)
    weight_num.append(w_train_pt)

test_accu_pt_num = test_accu_curve(weight_num, x_test, y_label_test)

plt.figure(num="General Learning Curve for Perceptron-Binary")
lb1 = "Perceptron"
plt.plot(np.arange(100, 2100,100), test_accu_pt_num, 'x-b', label = lb1)
plt.xlabel("Number of Training Examples")
plt.ylabel("Accuracy")
plt.legend() 


