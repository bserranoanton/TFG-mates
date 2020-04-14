# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 12:00:10 2020

@author: MRS
"""

# first neural network with keras tutorial
from numpy import loadtxt
from keras.models import Sequential
from keras.layers import Dense

# load the dataset
dataset = loadtxt('pima-indians-diabetes.csv', delimiter=',')
# split into input (X) and output (y) variables
X = dataset[:,0:5]
y = dataset[:,5:9]



# define the keras model
model = Sequential()
model.add(Dense(12, input_dim=5, activation='relu'))
model.add(Dense(8, activation='relu'))
#model.add(Dense(8, activation='relu'))
model.add(Dense(4, activation='linear'))

# compile the keras model
model.compile(loss='mean_squared_error', optimizer='adam', 
              metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(X, y, epochs=150, batch_size=10, verbose=0)


# evaluate the keras model
_, accuracy = model.evaluate(X, y)
print('Accuracy: %.2f' % (accuracy*100))

# make class predictions with the model
predictions = model.predict(X)
# summarize the first 5 cases
for i in range(5):
    print('%s => %s (expected %s)' % (X[i].tolist(), 
                                   predictions[i].tolist(), 
                                   y[i].tolist()))
    