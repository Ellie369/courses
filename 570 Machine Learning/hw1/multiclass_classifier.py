## multiclass classifier

import import_data
import numpy as np
import matplotlib.pyplot as plt

## import data 
x_train, y_train, x_test, y_test = import_data.load_data()


# augmented feature function

def Fun_aug(xt,y):
    num_k =10  #In the dataset of 'Fashion', there are 10 classes
    aug_x = np.zeros(num_k*len(xt))
    aug_x[(y*len(xt)):(len(xt)+y*len(xt))] = xt
    return aug_x

   
def mc_perceptron(x,y, itr): # x:examples; y:labels; itr:iteration
    learning_curve = [0]*itr #for 5.2a, count num of mistakes at each itr
    accuracy_curve = [0]*itr # for 5.2b, accuracy of each iteration  
    weight_list = [0]*itr
    num_k =10  #num of classes
    weight = np.zeros(num_k*len(x[0])) #initialize weight as zeros       
    for i in range(itr):
        errors = 0 #number of mistakes
        for t in range(len(x)):
            y_temp = np.zeros(num_k)
            for k in range(num_k):                              
                y_temp[k] = np.dot(weight,Fun_aug(x[t],k))
            y_cap = np.argmax(y_temp) # prediction of yt
            if y_cap != y[t]:
                dif = Fun_aug(x[t],y[t])-Fun_aug(x[t],y_cap)                
                weight = np.add(weight,dif)
                errors +=1
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(x))
        weight_list[i] = weight
    return weight,weight_list, learning_curve, accuracy_curve

def mc_aver_perceptron(x, y, itr): # x:examples; y:labels; itr:iteration
    learning_curve = [0]*itr #for 5.2a, count num of mistakes at each itr
    accuracy_curve = [0]*itr # for 5.2b, accuracy of each iteration  
    aver_weight_list = [0]*itr
    num_k =10  #num of classes
    weight = np.zeros(num_k*len(x[0])) #initialize weight as zeros     
    weight_aver = np.zeros(num_k*len(x[0]))
    for i in range(itr):
        errors = 0 #number of mistakes
        for t in range(len(x)):
            y_temp = np.zeros(num_k)
            for k in range(num_k):                              
                y_temp[k] = np.dot(weight,Fun_aug(x[t],k))
            y_cap = np.argmax(y_temp) # prediction of yt
            if y_cap != y[t]:
                dif = Fun_aug(x[t],y[t])-Fun_aug(x[t],y_cap)                
                weight = np.add(weight,dif)
                errors +=1
            weight_aver += weight            
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(x))
        aver_weight_list[i] = weight_aver
    return weight, aver_weight_list, learning_curve, accuracy_curve   

def mc_pa(x,y, itr): # x:examples; y:labels; itr:iteration
    learning_curve = [0]*itr #for 5.2a, count num of mistakes at each itr
    accuracy_curve = [0]*itr # for 5.2b, accuracy of each iteration  
    weight_list = [0]*itr
    num_k =10  #num of classes
    weight = np.zeros(len(x[0])*num_k) #initialize weight as zeros       
    for i in range(itr):
        errors = 0 #number of mistakes   
        for t in range(len(x)):
            y_temp = np.zeros(num_k)
            for k in range(num_k):               
                y_temp[k] = np.dot(weight,Fun_aug(x[t],k))
            y_cap = np.argmax(y_temp) # prediction of yt
            if y_cap != y[t]:
                dif = Fun_aug(x[t],y[t])-Fun_aug(x[t],y_cap)  
                # learning rate
                rate_a = 1-np.subtract(np.dot(weight,Fun_aug(x[t],y[t])),np.dot(weight,Fun_aug(x[t],y_cap)))
                rate_b = np.linalg.norm(dif)**2
                rate = rate_a/rate_b   
                # update weight
                weight = np.add(weight,np.dot(rate,dif))
                errors +=1
        learning_curve[i] = errors         
        accuracy_curve[i]=1-(errors/len(x))
        weight_list[i] = weight
    return weight, weight_list, learning_curve, accuracy_curve

def test(weight, x_test, y_test):    
    num_k =10  #num of classes    
    right = 0
    for t in range(len(x_test)):
        y_temp = np.zeros(num_k)
        for k in range(num_k):
            y_temp[k] = np.dot(weight,Fun_aug(x_test[t],k))
            y_cap = np.argmax(y_temp) # prediction of yt  
            if y_test[t]==y_cap:
                right +=1 
        accuracy = (right/len(x_test))          
    return accuracy

def test_accu_curve(weight_list, x_test, y_test):
    test_accuracy_curve =[]
    for wt in range(weight_list):
        accu = test(wt,x_test, y_test)
        test_accuracy_curve.append(accu)
    return test_accuracy_curve
    

## 5.2a
[w_train_mpt,wlist_train_mpt, learning_curve_mpt,accuracy_curve_mpt] = mc_perceptron(x_train,y_train, 50)
[w_train_mpa,wlist_train_mpa, learning_curve_mpa,accuracy_curve_mpa] = mc_pa(x_train,y_train, 50)
plt.figure(num="Multi-Class Classifier Online Learning Curve")
lb1 = "Perceptron Algorithm"
lb2 = "Passive Aggressive Algorithm"
plt.plot(np.arange(1, 51, 1), learning_curve_mpt, 'x-k', label = lb1)
plt.plot(np.arange(1, 51, 1), learning_curve_mpa, 'o-r', label = lb2)
plt.xlabel("Number of Training Iterations")
plt.ylabel("Number of Mistakes")
plt.legend()  

#5.2b
[w_train_mpt,wlist_train_mpt, learning_curve_mpt,accuracy_curve_mpt] = mc_perceptron(x_train,y_train, 20)
test_accu_mpt = test_accu_curve(wlist_train_mpt,x_test,y_test)
[w_train_mpa,wlist_train_mpa, learning_curve_mpa,accuracy_curve_mpa] = mc_pa(x_train,y_train, 20)
test_accu_mpa = test_accu_curve(wlist_train_mpa,x_test,y_test)
plt.figure(num="Multi-Class Classifier Online Accuracy Curves")
lb1 = "Perceptron-Training"
lb2 = "PA-Training"
lb3 = "Perceptron-Testing"
lb4 = "PA-Testing"
plt.plot(np.arange(1, 21, 1), accuracy_curve_mpt, 'x-b', label = lb1)
plt.plot(np.arange(1, 21, 1),accuracy_curve_mpa, 'o-r', label = lb2)
plt.plot(np.arange(1, 21, 1), test_accu_mpt, 'd-g', label = lb3)
plt.plot(np.arange(1, 21, 1), test_accu_mpa, '*-y', label = lb4)
plt.xlabel("Number of Iterations")
plt.ylabel("Accuracy")
plt.legend() 
    
#5.2c
[w_train_mpt,wlist_train_mpt, learning_curve_mpt,accuracy_curve_mpt] = mc_perceptron(x_train,y_train, 20)
test_accu_mpt = test_accu_curve(wlist_train_mpt,x_test,y_test)
[aver_w_train_mpt,aver_wlist_train_mpt, learning_curve_mpt,accuracy_curve_mpt] = mc_aver_perceptron(x_train,y_train, 20)
test_accu_mpt = test_accu_curve(aver_wlist_train_mpt,x_test,y_test)

plt.figure(num="Accuracy Curves of Plain Perceptron and Averaged Perceptron")
lb1 = "Plain Perceptron"
lb2 = "Averaged Perceptron"
plt.plot(np.arange(1, 21, 1), test_accu_mpt, 'x-b', label = lb1)
plt.plot(np.arange(1, 21, 1),test_accu_mpa, 'o-r', label = lb2)
plt.xlabel("Number of Iterations")
plt.ylabel("Accuracy")
plt.legend() 

#5.2d

weight_num=[]
index_list = np.arange(0, x_train.shape[0])
for num in range(20):    
    train_length = int(100+100*num) 
    train_indices = index_list[:train_length]
    x_train_num = x_train[train_indices]
    y_train_num = y_train[train_indices]
    [w_train_pt,wlist_train_pt, mistake_train_pt,accuracy_train_pt] = mc_perceptron(x_train_num,y_train_num, 20)
    weight_num.append(w_train_pt)
    
test_accu_mpt_num = test_accu_curve(weight_num, x_test, y_test)
plt.figure(num="General Learning Curve for Perceptron")
lb1 = "Perceptron"
plt.plot(np.arange(100, 2100,100), test_accu_mpt_num, 'x-b', label = lb1)
plt.xlabel("Number of Training Examples")
plt.ylabel("Accuracy")
plt.legend() 



