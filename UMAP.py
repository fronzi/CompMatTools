#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 10:40:19 2021

@author: Marco
"""

import pandas as pd
import umap.umap_ as umap
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import sys 

umap_test = pd.read_csv("ALMAP_test.csv")
#umap_test

umap_test.info()



umap_test = umap_test.drop(columns=['Unnamed: 0','Bilayer'])
print(umap_test)

umap_test.isna().sum()   #counts NaN 



#umap_test = umap_test.fillna(0)
umap_test = umap_test.dropna(0)

#sys.exit(-1)

umap_test.isna().sum()

reducer = umap.UMAP(random_state=42)


umap_data = umap_test[umap_test.columns[0:50]].values
scaled_umap_data = StandardScaler().fit_transform(umap_data)

umap_com = umap_test[umap_test.columns[0:49]].values
scaled_umap_com = StandardScaler().fit_transform(umap_com)

embedding = reducer.fit_transform(scaled_umap_com)

embedding.shape

plt.figure(figsize=(16,8))
plt.title("UMAP plot of Band Gap")

### for continuos values 

#plt.scatter(x=embedding[:,0],y=embedding[:,1],s=0.3,c=umap_test['HSE'].values)       

####   for discrete values
 
umap_test['HSE'] = umap_test['HSE'].replace({'L':'0', 'AL1':'1','AL2':'2','AL3':'3'})
plt.scatter(x=embedding[:,0],y=embedding[:,1],c=umap_test['HSE'].astype(int),cmap='Spectral', s=5)
plt.show()

