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

dfgeo = pd.read_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE116754_METH_stem/sample.accessions.59091.csv') #place "r" before the path string to address special character, such as '\'. Don't forget to put the file name at the end of the path + '.xlsx'

patients=dfgeo['Accession'].tolist()
processed_extracts=np.array([1]*len(patients)*8).reshape(len(patients),8)
pe_df = pd.DataFrame(processed_extracts, 
                     columns = ['Donor_ID', 'Cell_Line_Num' ,'Donor_Cell_Type',
                                'Reprogram_Method','Sex', 'Replicate_Num',
                                'Diff_Capacity','GSE59091_Sample_ID'])

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
        
    donor_id=re.search('id:(.+?)cell', found).group(1)
    cell_line=re.search('line #(.+?)donor cell', found).group(1)
    cell_type=re.search('cell type:(.+?)reprogramming', found).group(1)
    reprog_meth=re.search('method:(.+?)Sex', found).group(1)
    sex=re.search('Sex:(.+?)replicate', found).group(1)
    rep_num=re.search('replicate #:(.+?)differentiation', found).group(1)
    diff_cap=re.search('capacity:(.+?)Treatment', found).group(1)
    
    
    pe_df['Donor_ID'][i]=donor_id
    pe_df['Cell_Line_Num'][i]=cell_line
    pe_df['Donor_Cell_Type'][i]=cell_type
    pe_df['Reprogram_Method'][i]=reprog_meth
    pe_df['Sex'][i]=sex
    pe_df['Replicate_Num'][i]=rep_num
    pe_df['Diff_Capacity'][i]=diff_cap
    pe_df['GSE59091_Sample_ID'][i]=patients[i]
    
    
pe_df.to_csv('C:/Users/Bruin/Documents/GitHub/MRPC_support/Cancer_Selection/GSE116754_METH_stem/GSE59091_meta_info.csv')
    
    