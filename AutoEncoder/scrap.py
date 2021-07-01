# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 15:13:16 2021

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

imputed_data=pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\imputation.step2.hd5')
df = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.hd5') 
df2 = pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\pbmc.g949_c10k.msk90.hd5') 

imputed_data=np.array(imputed_data)
df=np.array(df)
df2=np.array(df2)

            

mse=MSE_nz(df2, imputed_data)