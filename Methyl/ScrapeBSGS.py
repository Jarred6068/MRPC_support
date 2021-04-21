# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:41:51 2021

@author: Bruin
"""

import requests
import numpy as np
from bs4 import BeautifulSoup
import re
import pandas as pd

dftxt = pd.read_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Methyl/GSE53195_GSE56105_ID_equivalence.txt',
                    delimiter = " ")
a=np.arange(1286226, 1286578).tolist()
b=[1286580, 1286582, 1286584, 1286586, 1286588]
c=np.arange(1286590, 1286731).tolist()
numrange=a+b+c
patients=['GSM'+str(numrange[i]) for i in np.arange(0,len(numrange))]
processed_extracts=np.array([1]*len(patients)*6).reshape(len(patients),6)
pe_df = pd.DataFrame(processed_extracts, 
                     columns = ['ID','Gender','Twin','Age','GSE53195_ID', 'GSE56105_ID'])

for i in np.arange(0, len(patients)):

    URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+ patients[i]
    page = requests.get(URL)

    soup = BeautifulSoup(page.content, 'html.parser')

    results = soup.find('table')
    #print(results.text)

    extracted_text=results.text.replace('\n', ' ').strip()
    
    print(numrange[i])
    m = re.search('Characteristics(.+?)Treatment', extracted_text)
    if m:
        found = m.group(1)
        
    gender=re.search('gender:(.+?)twin', found).group(1)
    twin=re.search('twin:(.+?)age', found).group(1)
    age=re.search('age:(.+?)yrs', found).group(1)
    
    pe_df['Gender'][i]=gender
    pe_df['Twin'][i]=twin
    pe_df['Age'][i]=age
    pe_df['ID'][i]=patients[i]
    pe_df['GSE53195_ID'][i]=dftxt['GSE53195_ID'][i]
    pe_df['GSE56105_ID'][i]=dftxt['GSE56105_ID'][i]
    
pe_df.to_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Methyl/GEO_Accession_data.csv')
    
    