import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import seaborn as sb

df = pd.read_csv(r'C:\Users\rtm11\Box\rtm11_files\Documents\python\vol_datav4.csv')
df.head()

p = df.corr(method = 'pearson')

sb.heatmap(p, 
            xticklabels=p.columns,
            yticklabels=p.columns,
            cmap='RdBu_r',
            annot=True,
            linewidth=0.5)

p_abs = abs(p["krel"])
rel_feat = p_abs[p_abs >0.4].drop('krel')
rel_feat.sort_values(ascending=False)