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

dfgeo = pd.read_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE116754_METH_stem/GSE116754/sample.accessions.116754.csv') #place "r" before the path string to address special character, such as '\'. Don't forget to put the file name at the end of the path + '.xlsx'

patients=dfgeo['Accession'].tolist()
processed_extracts=np.array([1]*len(patients)*5).reshape(len(patients),5)
pe_df = pd.DataFrame(processed_extracts, 
                     columns = ['Tissue', 'Sample_Type' ,'Developement_Stage',
                                'Gender','GSE116754_Sample_ID'])

for i in np.arange(0, len(patients)):

    URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+ patients[i]
    page = requests.get(URL)

    soup = BeautifulSoup(page.content, 'html.parser')

    results = soup.find('table')
    #print(results.text)

    extracted_text=results.text.replace('\n', ' ').strip()
    
    m = re.search('Characteristics(.+?)protocol', extracted_text)
    print(i)
    if m:
        found = m.group(1)
        
    tiss=re.search('tissue:(.+?)sample', found).group(1)
    check=re.search('type:(.+?)stage', found)
    if check is None:
        samptype=None
    else:
        samptype=re.search('type:(.+?)stage', found).group(1)
    check2=re.search('stage:(.+?)gender', found)
    if check2 is None:
        dev_stage=None
    else:
        dev_stage=re.search('stage:(.+?)gender', found).group(1)
    check3=re.search('gender:(.+?)Growth', found)
    if check3 is None:
        gender=None
    else:
        gender=re.search('gender:(.+?)Growth', found).group(1)

    
    
    pe_df['Tissue'][i]=tiss
    pe_df['Sample_Type'][i]=samptype
    pe_df['Developement_Stage'][i]=dev_stage
    pe_df['Gender'][i]=gender
    pe_df['GSE116754_Sample_ID'][i]=patients[i]
    
    
pe_df.to_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE116754_METH_stem/GSE116754/GSE116754_meta_info.csv')
    
    