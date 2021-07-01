# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:30:34 2021

@author: Bruin
"""



import numpy as np
import pandas as pd

#read in weights, bias, and imputed data (encoder)
data_w2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_w2.step2.npy')
data_b2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_b2.step2.npy')
data_w1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_w1.step2.npy')
data_b1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_b1.step2.npy')
imputed_data=pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\imputation.step2.hd5')


#read in weights and bias's (decoder)
data_dw2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_w2.step2.npy')
data_db2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_b2.step2.npy')
data_dw1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_w1.step2.npy')
data_db1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_b1.step2.npy')

from matplotlib import pyplot as plt

#plot the first hidden layer weight matrix (encoder)
fig = plt.figure(figsize=(9, 9))
plt.imshow(data_w2, cmap="rainbow", 
           vmin=np.min(data_w1), 
           vmax=np.max(data_w1), 
           aspect='auto')
plt.colorbar()

#plot the second hidden layer weight matrix
fig = plt.figure(figsize=(9, 9))
plt.imshow(data_w2, cmap="rainbow", 
           vmin=np.min(data_w2), 
           vmax=np.max(data_w2), 
           aspect='auto')
plt.colorbar()


#look at first 12 pairwise plots of the columns of W1 and W2
for i in (np.arange(10)+1):
    fig2=plt.figure(figsize=(8,8))
    ax=fig2.add_subplot(1,1,1)
    ax.scatter(data_w1[:,i-1], 
               data_w1[:,i], 
               c='b', 
               marker="s",
               label="W1")
    ax.scatter(data_w2[:,i-1], 
               data_w2[:,i], 
               c='r', 
               marker="o",
               label="W2")
    
    plt.legend(loc='upper left')
