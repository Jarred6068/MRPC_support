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

dfgeo = pd.read_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE101961_normal_breast/sample.accessions.csv') #place "r" before the path string to address special character, such as '\'. Don't forget to put the file name at the end of the path + '.xlsx'

patients=dfgeo['Accession'].tolist()
processed_extracts=np.array([1]*len(patients)*4).reshape(len(patients),4)
pe_df = pd.DataFrame(processed_extracts, 
                     columns = ["Age", "Race", "BMI" , 'GSE101961_Sample_ID'])

for i in np.arange(0, len(patients)):

    URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+ patients[i]
    page = requests.get(URL)

    soup = BeautifulSoup(page.content, 'html.parser')

    results = soup.find('table')
    #print(results.text)

    extracted_text=results.text.replace('\n', ' ').strip()
    
    m = re.search('Characteristics(.+?)molecule', extracted_text)
    print(i)
    if m:
        found = m.group(1)
        
    age=re.search('age:(.+?)race', found).group(1)
    race=re.search('race:(.+?)bmi', found).group(1)
    
    bmi=re.search('bmi:(.+?)Extracted', found).group(1)
    
    pe_df['Age'][i]=age
    pe_df['Race'][i]=race
    pe_df['BMI'][i]=bmi
    pe_df['GSE101961_Sample_ID'][i]=patients[i]
    
    
pe_df.to_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE101961_normal_breast/GSE101961_meta_info.csv')
    
    