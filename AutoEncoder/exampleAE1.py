# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 12:06:03 2021

@author: Bruin
"""
#import necessary modules
import pandas as pd
import tensorflow as tf
import keras as ks
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV 
import keras
from matplotlib import pyplot as p

#read in .HD5 files: the PBMC_G949_10k and PBMC_G949_10k_90msk datasets
df = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.hd5') 
df2 = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.msk90.hd5') 
#retain the Gene Names (column names)
colnames=df.columns
colnames2=df2.columns
#convert to array
df_array=np.array(df)
df2_array=np.array(df2)
#perform log10(x+1) transformation of the data to reduce variance in measure
df_log=np.array([np.log10(df_array[:,i]+1) for i in np.arange(df_array.shape[1])]).transpose()
df2_log=np.array([np.log10(df2_array[:,i]+1) for i in np.arange(df2_array.shape[1])]).transpose()

#conduct PCA on unmasked data and view first two components
pca=PCA(n_components=10)
PCA_matrix=pca.fit_transform(df_log)

fig = p.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
ax.scatter(PCA_matrix[:,1], PCA_matrix[:,2])

#create histogram for sparsity of matrix (proportions of nonzero elements in each
#sample)
#train_set_idx=[]
props=[1]*df2.shape[0]


for i in np.arange(df2_log.shape[0]):
    props[i]=sum(df2_log[i,:]>0)/len(df2_log[i,:])
    # if props[i]>0.1:
    #     train_set_idx.append(i)

figH = p.figure(figsize = (8,8))
ax = figH.add_subplot(1,1,1) 
ax.set_xlabel('proportion nonzero', fontsize = 15)
ax.set_ylabel('density', fontsize = 15)
ax.set_title('hist of nonzero density', fontsize = 20)
ax.hist(props)

#training_data=df2_log[train_set_idx,:]
    
#partition the masked data into training, testing, and validation sets 70:15:15
from sklearn.model_selection import train_test_split

#first random split into 70:30
validate,train=train_test_split(df2_log, 
                                test_size=0.7, 
                                random_state=2)

#split remaining 30 into halves
validate, test=train_test_split(validate,
                                test_size=0.5, 
                                random_state=2)


#=============================================================================
#----------------------Specialized-Non-Zero-Loss-Function---------------------
#=============================================================================

def MSE_nz(input_mat, output_mat):

    '''
    
    Loss Function described in Badsha et. al., 2020 : "Imputation of single-cell 
    gene expression with an autoencoder neural network"  which describes the 
    loss in terms of the mse between all non-zero values in the imputed and 
    training matrices
    
    Parameters
    ----------
    input_mat : The input data: an n x m dataset
    output_mat : The predicted data: an imputed n x m data set

    Returns
    -------
    the m x 1 loss array.
    
    Dependencies
    ------------
    This function relies on tensorflow functions because keras specifically
    uses the tensor object. Thus the input matrix with be an n x m tensor

    '''
    # input_mat=input_mat.numpy()
    # output_mat=output_mat.numpy()
    omega = tf.sign(input_mat)
    diff = tf.subtract(input_mat, output_mat)
    square_err_ = tf.pow(diff, 2)
    square_err_nz = tf.reduce_sum(tf.multiply(square_err_, omega), axis=1)
    
    return square_err_nz
        

#=============================================================================
#-----------------------------Build-AutoEncoder-Model-------------------------
#=============================================================================


def create_AE(hidden_layers = [400,200,400], act = 'relu', opt = 'adam', 
              ID=949, OD=949, BI="Zeros", KI='glorot_uniform', LR=0.0003):
    
    
    '''
    creates keras neural network model
    
    Parameters
    -----------
    hidden_layers: a list that defines the numbers of hidden nodes for all hidden 
        layers, e.g., [1000] indicates the nn with only one hidden layer and 1000 
        nodes has hidden_layers=[1000], while [1000, 500] defines two hidden 
        layers and the first layer has 1000 nodes and the second has 500 nodes.
    act: activation function for all hidden layers
    opt: optimization function to be passed to keras.compile()
    ID: the input dimension of the data (# of columns of input)
    OD: the output dimension of the model (# of nodes in output layer)
    BI: the bias initialization function: see keras API documentation
    KI: the kernal initialization function: see keras API documentation
    LR: the learning rate of the optimization function
    
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
            kernel_initializer=KI,
            bias_initializer=BI,
            activation=act))
        
        in_dim=hidden_layers[i]
        print(in_dim)
    
        
    #2.2 build the output layer and use the softmax activation  
    myAE.add(keras.layers.Dense(
        units=out_dim,
        input_dim=hidden_layers[len(hidden_layers)-1],
        use_bias=True,
        kernel_initializer=KI,
        bias_initializer=BI,
        activation=act))
    
    #2.3 choose the optimizer, compile the network and return it. Use 'accuracy' as the metrics
    if(opt=="rmsprop"):
        optimizer = keras.optimizers.RMSprop(lr=LR, decay=1e-7, momentum=.9)
        
    elif(opt=="adam"):
        optimizer = keras.optimizers.Adam(lr=LR)
        
    else:
        optimizer = keras.optimizers.SGD(lr=LR, decay=1e-7, momentum=.9)
        
    myAE.compile(optimizer=optimizer, loss=MSE_nz,
                 metrics='mse')
    
    
    return myAE


#create and fit model to training data
myAE = create_AE(hidden_layers=[800,400,800], LR=0.003)
myAE.summary()


myAE.fit(train,
         train,
         batch_size=256, 
         validation_data=(validate,validate),
         epochs=20)

w1_e=myAE.layers[0].get_weights()[0]
w2_e=myAE.layers[1].get_weights()[0]
w1_d=myAE.layers[2].get_weights()[0]
w2_d=myAE.layers[3].get_weights()[0]

#plot and view the model weight matrices:

#plot the first hidden layer weight matrix (encoder)
# fig = p.figure(figsize=(9, 9))
# p.imshow(w1_e, cmap="rainbow", 
#            vmin=np.min(w1_e), 
#            vmax=np.max(w1_e), 
#            aspect='auto')
# p.colorbar()

# #plot the second hidden layer weight matrix
# fig = p.figure(figsize=(9, 9))
# p.imshow(w2_e, cmap="rainbow", 
#            vmin=np.min(w2_e), 
#            vmax=np.max(w2_e), 
#            aspect='auto')
# p.colorbar()
    
#load in the imputed matrix from running example 1

imputed_data=pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\imputation.step2.hd5')
imputed_data=np.array(imputed_data)

#use trained keras model to impute full 90% masked dataset
my_imputed=myAE.predict(df2_log, batch_size=64)

#head map of their imputed, my imputed, and original data
fig = p.figure(figsize=(9, 9))
p.imshow(df_log, cmap="rainbow", 
           vmin=np.min(df_log), 
           vmax=np.max(df_log), 
           aspect='auto')
p.title("Original Log10 data (unmasked)", fontsize=20)
p.colorbar()
#save
p.savefig(r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\DiagnosticPlots\ordatamatrix.png')

fig = p.figure(figsize=(9, 9))
p.imshow(my_imputed, cmap="rainbow", 
           vmin=np.min(my_imputed), 
           vmax=np.max(my_imputed), 
           aspect='auto')
p.title("Keras Imputed Result", fontsize=20)
p.colorbar()
p.savefig(r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\DiagnosticPlots\KerasImpute.png')

fig = p.figure(figsize=(9, 9))
p.imshow(imputed_data, cmap="rainbow", 
           vmin=np.min(imputed_data), 
           vmax=np.max(imputed_data), 
           aspect='auto')
p.title("Example 1 Imputed Result", fontsize=20)
p.colorbar()
p.savefig(r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\DiagnosticPlots\Ex1Imputed.png')


#plot few of the genes in pairwise mannor to see differences between methods

for i in (np.arange(10)+1):

    fig = p.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel(colnames[i-1], fontsize = 15)
    ax.set_ylabel(colnames[i], fontsize = 15)
    ax.set_title("Gene biplot comparing Keras solution to Badsha-Fu", fontsize = 20)
    ax.scatter(my_imputed[:,i-1],
               my_imputed[:,i],                
               c='b', 
               marker="s",
               label="Keras_W1")
    ax.scatter(imputed_data[:,i-1],
               imputed_data[:,i],
               c='r', 
               marker="o",
               label="Their_W1")
    p.legend(loc='upper left')

#save model weight matrices to .csv file using pandas










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







