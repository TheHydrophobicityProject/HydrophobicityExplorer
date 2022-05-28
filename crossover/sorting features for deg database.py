#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from scipy.stats import pearsonr


# In[2]:


df = pd.read_csv(r'C:\Users\rtm11\Box\rtm11_files\Documents\python\vol_datav4.csv')
df.head()


# In[5]:


p = df.corr(method = 'pearson')

import seaborn as sb
sb.heatmap(p, 
            xticklabels=p.columns,
            yticklabels=p.columns,
            cmap='RdBu_r',
            annot=True,
            linewidth=0.5)


# In[4]:


p_abs = abs(p["krel"])
rel_feat = p_abs[p_abs >0.4].drop('krel')
rel_feat.sort_values(ascending=False)

