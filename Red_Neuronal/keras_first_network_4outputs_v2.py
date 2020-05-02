# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:00:10 2020

@author: MRS
"""

# first neural network with keras tutorial
#import numpy as np
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense
from random import sample
import matplotlib.pyplot as plt
import numpy as np

# load the dataset
dataset = loadtxt('data_neural_network_csv.csv', delimiter=',')
# split into input (X) and output (y) variables
X = dataset[:,0:6]
y = dataset[:,6:10]

#Select 70% to fit the network and 30% to test it

numSimulations = X.shape[0]
numDataFit = round(numSimulations*0.7)
#Generate a vector of round(numSimulations*0.7) random numbers 
L=sample(range(0,numSimulations), numDataFit)

X_fit = X[L]
y_fit = y[L]


# define the keras model
model = Sequential()
model.add(Dense(12, input_dim=6, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(4, activation='relu'))

# compile the keras model
model.compile(loss='mean_squared_error', optimizer='adam', 
              metrics=['accuracy'])

# fit the keras model on the dataset
history = model.fit(X, y, validation_split=0.33, epochs=280, batch_size=10, verbose=0)


# evaluate the keras model
_, accuracy = model.evaluate(X, y)
print('Accuracy: %.2f' % (accuracy*100))

# summarize history for accuracy
f = plt.figure(1)
plt.plot(history.history['accuracy'])
plt.plot(history.history['val_accuracy'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
f.show()
#plt.show()
# summarize history for loss
g = plt.figure(2)
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
g.show()
#plt.show()

#input()

#test the network with the rest of the data set
index_fit_set = {l for l in L }
all_index = list(range(0, numSimulations))
all_index_set = {l for l in all_index }
index_test_set = all_index_set - index_fit_set

X_test = X[list(index_test_set)]
y_test = y[list(index_test_set)]

# make class predictions with the model
predictions = model.predict(X_test,batch_size=1)


#plot results
# f = open("results70_30.txt", "a")

# for i in range(len(index_test_set)):
#     f.write('%s => %s (expected %s)\n' % (X_test[i].tolist(), 
#                                      predictions[i].tolist(), 
#                                     y_test[i].tolist()))
#     if(i%10 == 0):
#         f.write('-------------------------------------------\n')
#     # print('%s => %s (expected %s)' % (X_test[i].tolist(), 
#     #                                 predictions[i].tolist(), 
#     #                                 y_test[i].tolist()))
    
# f.close()
print('-----------------------------------------------------')


# # load the dataset
# dataset_prueba = loadtxt('data_neural_network_csv_prueba.csv', delimiter=',')
# X_prueba = dataset_prueba[:,0:6]
# y_prueba = dataset_prueba[:,6:10]

# # make class predictions with the model
# predictions = model.predict(X_prueba,batch_size=1)

# #summarize the first 15 cases
# for i in range(15):
#     print('%s => %s (expected %s)' % (X_prueba[i].tolist(), 
#                                     predictions[i].tolist(), 
#                                     y_prueba[i].tolist()))
# print('-----------------------------------------------------')

# for i in range(15):
#     print('%s => %s (expected %s)' % (X_prueba[200-i].tolist(), 
#                                     predictions[200-i].tolist(), 
#                                     y_prueba[200-i].tolist()))
# print('-----------------------------------------------------')


# for i in range(15):
#     print('%s => %s (expected %s)' % (X_prueba[100-i].tolist(), 
#                                     predictions[100-i].tolist(), 
#                                     y_prueba[100-i].tolist()))
# print('-----------------------------------------------------')

# for i in range(15):
#     print('%s => %s (expected %s)' % (X_prueba[400-i].tolist(), 
#                                     predictions[400-i].tolist(), 
#                                     y_prueba[400-i].tolist()))
    
# print('-----------------------------------------------------')    
# for i in range(15):
#     print('%s => %s (expected %s)' % (X_prueba[900-i].tolist(), 
#                                     predictions[900-i].tolist(), 
#                                     y_prueba[900-i].tolist()))

    
# print('------------------------------------------------------')
myValue = np.array([[116.68,89,4.5, 6.6,6.9,8.7],
                    [74.4,92,3.15,4.8,3.75,6.45]])

myPrediction = model.predict(myValue, batch_size=1)
print('Prediccion:')
print(myPrediction.tolist())