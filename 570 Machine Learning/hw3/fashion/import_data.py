# import data
import os
import mnist_reader


def load_data():
    base_dir = os.getcwd()
    path = os.path.join(base_dir)
    print(path)
    x_train,y_train = mnist_reader.load_mnist(path,kind='train') # x,images; y,labels
    x_test, y_test = mnist_reader.load_mnist(path, kind='t10k')
    return x_train, y_train, x_test, y_test


# define labels, even= 1; odd= -1
def y_even_odd(y_data):
    y_label = [0]*len(y_data)
    for i in range(len(y_data)):
        if (y_data[i]%2)==0:
            y_label[i] = 1
        else:
            y_label[i] = -1
    return y_label

def y_label(y_train, y_test):
    y_label_train = y_even_odd(y_train)
    y_label_test = y_even_odd(y_test)
    return y_label_train, y_label_test
