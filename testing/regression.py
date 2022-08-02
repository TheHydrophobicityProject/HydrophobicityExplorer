import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def func_exp(x, a, b, c):
    #c = 0
    return a * (x**b) + c

def exponential_regression (x_data, y_data):
    popt, _ = curve_fit(func_exp, x_data, y_data)
    print(popt)
    points = plt.plot(x_data, y_data, 'x', color='xkcd:maroon', label = "data")
    regresion_curve = plt.plot(x_data, func_exp(x_data, *popt), color='xkcd:teal', label = "fit: {:.3f}*x^{:.3f}+{:.3f}".format(*popt))
    plt.legend()
    plt.show()
    return func_exp(x_data, *popt)

data = pd.read_csv('RG.csv', usecols=["RG", "N"])

rg=data["RG"]
n=data["N"]

exponential_regression(n, rg)

# print(data)
# print(data["RG"])
