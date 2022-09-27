import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


dict = {'MHP': [0.012376251236427096, 0.013573507027943717, 0.014126042433896363, 0.01443698774172117, 0.014628104033066558, 
0.014986073983959734, 0.015104605019622333, 0.015056990576087965, 0.0149831465679757, 0.015223976222613582], 
'N': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 
'smi': ['CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 
'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 
'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 
'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 
'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)', 
'CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)CC(c1ccccc1)']}

def crossover(calc,*,data_dict=None,data_csv=None):
    if data_csv is not None:
        df = pd.read_csv(data_csv)
    elif data_dict is not None:
        df = pd.DataFrame.from_dict(data_dict)
    else:
        raise Exception("No data provided")

    X=df[["N"]]
    Y=df[[calc]]

    #calculate lin reg for each 
    
#crossover("MHP", data_csv="data.csv")
crossover("MHP", data_dict=dict)