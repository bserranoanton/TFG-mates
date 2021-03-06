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
import numpy as np


# load the dataset
dataset = loadtxt('data_neural_network_csv.csv', delimiter=',')
# split into input (X) and output (y) variables
X = dataset[:,0:6]
y = dataset[:,6:10]



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
#model.fit(X, y, epochs=150, batch_size=10, verbose=0)
model.fit(X, y, epochs=250, batch_size=10, verbose=0)


# evaluate the keras model
_, accuracy = model.evaluate(X, y)
print('Accuracy: %.2f' % (accuracy*100))

# make class predictions with the model
#predictions = model.predict(X,batch_size=1)

# summarize the first 5 cases
# for i in range(5):
#     print('%s => %s (expected %s)' % (X[i].tolist(), 
#                                     predictions[i].tolist(), 
#                                     y[i].tolist()))
    
# for i in range(5):
#     print('%s => %s (expected %s)' % (X[100-i].tolist(), 
#                                 predictions[100-i].tolist(), 
#                                 y[100-i].tolist()))
    
# for i in range(5):
#     print('%s => %s (expected %s)' % (X[300-i].tolist(), 
#                                 predictions[300-i].tolist(), 
#                                 y[300-i].tolist()))


# load the dataset
# dataset_prueba = loadtxt('data_neural_network_csv_prueba.csv', delimiter=',')
# X_prueba = dataset_prueba[:,0:6]
# y_prueba = dataset_prueba[:,6:10]

# make class predictions with the model
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

    
myValue = np.array([[116.68,89,4.5, 6.6,6.9,8.7],
                    [450.086,342,3.3,8.4,4.8,11.7]])

myPrediction = model.predict(myValue, batch_size=1)
print('Prediccion:')
print(myPrediction.tolist())
