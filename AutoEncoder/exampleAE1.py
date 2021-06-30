# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 12:06:03 2021

@author: Bruin
"""

import pandas as pd
import tensorflow as tf
import keras as ks
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV 
import keras
from matplotlib import pyplot as p

df = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.hd5') 
df2 = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.msk90.hd5') 
colnames=df.columns
colnames2=df2.columns
df_array=np.array(df)
df2_array=np.array(df2)
df_log=np.array([np.log10(df_array[:,i]+1) for i in np.arange(df_array.shape[1])]).transpose()
#df_log.columns=colnames
df2_log=np.array([np.log10(df2_array[:,i]+1) for i in np.arange(df2_array.shape[1])]).transpose()
#df2_log.columns=colnames2

pca=PCA(n_components=10)
PCA_matrix=pca.fit_transform(df_log)

fig = p.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
ax.scatter(PCA_matrix[:,1], PCA_matrix[:,2])

train_set_idx=[]
props=[1]*df2.shape[0]
#training_data=

for i in np.arange(df2_log.shape[0]):
    props[i]=sum(df2_log[i,:]>0)/len(df2_log[i,:])
    if props[i]>0.1:
        train_set_idx.append(i)
        
p.hist(props)

training_data=df2_log[train_set_idx,:]
    

from sklearn.model_selection import train_test_split

validate,train=train_test_split(training_data, 
                                test_size=0.7, 
                                random_state=2)

validate, test=train_test_split(validate,
                                test_size=0.5, 
                                random_state=2)

# validate2,train2=train_test_split(df2_log, 
#                                 test_size=0.7, 
#                                 random_state=2)

# validate2, test2=train_test_split(validate2,
#                                 test_size=0.5, 
#                                 random_state=2)


#=============================================================================

def create_encoder(hidden_units=[400], act='relu', opt='adam', ID=949, OD=200):
    
    '''creates the encoder'''
    
    myAE_encoder = keras.models.Sequential()
    
    myAE_encoder.add(keras.Input(shape=(ID,)))
        
    #2.2 build the output layer and use the softmax activation  
    myAE_encoder.add(keras.layers.Dense(
        units=hidden_units[0],
        input_dim=ID,
        kernel_initializer='glorot_uniform',
        bias_initializer='zeros',
        activation=act))
    
    myAE_encoder.add(keras.layers.Dense(
        units=OD,
        input_dim=hidden_units[0],
        kernel_initializer='glorot_uniform',
        bias_initializer='zeros',
        activation=act))
        
        
    #2.3 choose the optimizer, compile the network and return it. Use 'accuracy' as the metrics
    if(opt=="rmsprop"):
        optimizer = keras.optimizers.RMSprop(lr=0.003, decay=1e-7, momentum=.9)
        
    elif(opt=="adam"):
        optimizer = keras.optimizers.Adam(lr=0.003)
        
    else:
        optimizer = keras.optimizers.SGD(lr=0.003, decay=1e-7, momentum=.9)
        
    myAE_encoder.compile(optimizer=optimizer, loss=None,
                 metrics=None)
    
    
    return myAE_encoder


# myAE_encoder=create_encoder()
# myAE_encoder.summary()
# myAE_encoder.fit(train,
#                  train,
#                  batch_size=256, 
#                  epochs=100)

#=============================================================================


def create_AE(hidden_layers = [400,200,400], act = 'relu', opt = 'adam', 
              ID=949, OD=949): 
    '''create a deep feedforwad neural network using keras
    
    Parameters
    -----------
    hidden_layers: a list that defines the numbers of hidden nodes for all hidden layers, e.g., [1000] indicates
    the nn has only one hidden layer with 1000 nodes, while [1000, 500] defines two hidden layers and the first
    layer has 1000 nodes and the second has 500 nodes.
    act: activation function for all hidden layers
    opt: optimizer
    
    Returns
    -------
    myAE: the neural network autoencoder model
    
    '''
    in_dim = ID
    out_dim = OD
    
    ## add your code here
    myAE = keras.models.Sequential()
    
    #2.1 build all hidden layers
    for i in np.arange(len(hidden_layers)):
        myAE.add(keras.layers.Dense(
            units=hidden_layers[i],
            input_dim=in_dim,
            use_bias=True,
            kernel_initializer='glorot_uniform',
            bias_initializer='Zeros',
            activation=act))
        
        in_dim=hidden_layers[i]
        print(in_dim)
    
        
    #2.2 build the output layer and use the softmax activation  
    myAE.add(keras.layers.Dense(
        units=out_dim,
        input_dim=hidden_layers[len(hidden_layers)-1],
        use_bias=True,
        kernel_initializer='glorot_uniform',
        bias_initializer='Zeros',
        activation=act))
    
    #2.3 choose the optimizer, compile the network and return it. Use 'accuracy' as the metrics
    if(opt=="rmsprop"):
        optimizer = keras.optimizers.RMSprop(lr=0.003, decay=1e-7, momentum=.9)
        
    elif(opt=="adam"):
        optimizer = keras.optimizers.Adam(lr=0.003)
        
    else:
        optimizer = keras.optimizers.SGD(lr=0.003, decay=1e-7, momentum=.9)
        
    myAE.compile(optimizer=optimizer, loss='mse',
                 metrics='mse')
    
    
    return myAE

myAE = create_AE()
myAE.summary()


myAE.fit(train,
         train,
         batch_size=256, 
         epochs=1000)












# #=============================================================================
# from keras.wrappers.scikit_learn import KerasClassifier

# def AE_params_search(ae, X, y,param_grid): # 30 points
#     '''Search best paramaters
    
#     Parameters
#     ----------
#     X_train: features
#     y_train: target of the input
#     param_grid: a dict that defines the parameters
    
#     Returns
#     -------
#     best_params_
        
#     '''
#     ## add your code here. set cv = 3, scoring = 'accuracy', and verbose = 2

#     aeps = GridSearchCV(estimator = ae, cv = 2, param_grid = param_grid , 
#                        scoring = 'mse', verbose = 2)
    
#     aeps.fit(X,y)
    
#     return aeps

# ae = KerasClassifier(build_fn = create_AE, batch_size = 128, epochs = 10) # using the keras wapper

# #param_grid = {'batch_size': [128, 256], 
#    #           'epochs':[10,15,20]}

# #aetuned = AE_params_search(ae, train, param_grid = param_grid)







