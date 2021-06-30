# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:30:34 2021

@author: Bruin
"""
import numpy as np
import pandas as pd


data_w2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_w2.step2.npy')
data_b2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_b2.step2.npy')
data_w1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_w1.step2.npy')
data_b1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\e_b1.step2.npy')
imputed_data=pd.read_hdf (r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\imputation.step2.hd5')



data_dw2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_w2.step2.npy')
data_db2=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_b2.step2.npy')
data_dw1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_w1.step2.npy')
data_db1=np.load(file=r'C:\Users\Bruin\Documents\GitHub\MRPC_support\AutoEncoder\step2\d_b1.step2.npy')




