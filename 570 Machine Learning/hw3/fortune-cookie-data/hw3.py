
"""
hw3-programming 1
"""

import numpy as np
# input data
def load_data(filename):
    dataset = list()
    with open(filename, "r") as file:
        for row in file:   
            dataset.append(row.split()) 
        
    return dataset
#split data
def load_label(filename):
    dataset = list()
    with open(filename, "r") as file:
        for row in file:   
            dataset.append(row.strip())
    return dataset
#Form the vocabulary
def vocabulary_build(dataset):
    vocabulary = list()
    for row in range(len(dataset)):
        vector = dataset[row]     
        for i in range(len(vector)):
            if (vector[i] not in vocabulary) and (vector[i] not in stoplist):
                vocabulary.append(vector[i])
    return vocabulary        
#convert training data into features
def convert_feature(data, vocabulary):
    features = list()
    for row in range(len(data)):  
        temp = np.zeros(len(vocabulary))
        for i in range(len(vocabulary_sorted)):
            if vocabulary_sorted[i] in data[row]:
                temp[i]=1
        features.append(temp)                     
    return features
## sort data into classes,x=features, y=labels
def sort_data(x, y):
    xy0 = list()
    xy1 = list()
    for i in range(len(y)):
        if y[i]=='0':
            xy0.append(x[i])
        else:
            xy1.append(x[i])
    return xy0, xy1
#probability distribution of features
def pdf(data):
    F = len(data[0]) #num of features
    px1_y =np.zeros(F) # probability of xi=1|y=0    
    for i in range(len(data)):
        for j in range(F):
            if data[i][j] == 1:
                px1_y[j]+=1
            
    px1_y = np.add(px1_y,1)/(len(data)+F)
    return px1_y

def probability(prior, sample,cdp_x1): # cdp-conditional probability
    prob_xi = np.zeros(len(sample))
    for i in range(len(sample)):
        if sample[i]==1:
            prob_xi[i] = cdp_x1[i]
        else:
            prob_xi[i] =1- cdp_x1[i]
    
    prob = prior*(np.prod(prob_xi))
    return prob

#accuracy of prediction
def accuracy(data, labels):
    y_hat = np.array(['0','1']) #classes of labels
    accu =0
    for i in range(len(data)):
        predict0 = probability(py0,data[i], px1_y0)
        predict1 = probability(py1,data[i], px1_y1)
        idx = np.argmax([predict0,predict1])
        if y_hat[idx] == labels[i]:
            accu +=1
    accu = accu/len(data)
    return accu
    
# input data
traindata = load_data('traindata.txt')
trainlabels = load_label('trainlabels.txt')
testdata = load_data('testdata.txt')
testlabels = load_label('testlabels.txt')
stoplist = load_label('stoplist.txt')
## training
vocabulary = vocabulary_build(traindata)
vocabulary_sorted = sorted(vocabulary)
train_feature = convert_feature(traindata, vocabulary_sorted)
xy0, xy1 = sort_data(train_feature, trainlabels)

#calculate priors
py0 = len(xy0)/len(trainlabels)
py1 = len(xy1)/len(trainlabels)
# conditional probability
px1_y0 = pdf(xy0)  
px1_y1 = pdf(xy1) 
# accuracy of training data
accu_train = accuracy(train_feature, trainlabels)
# accuracy of testing data
test_feature = convert_feature(testdata, vocabulary_sorted)
accu_test = accuracy(test_feature, testlabels)

print('predicting accuracies for training and testing data are %.3f and %.3f' %(accu_train, accu_test))

