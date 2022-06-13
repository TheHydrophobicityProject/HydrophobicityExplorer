#data organizer loop
#this simply plots data
import matplotlib.pyplot as plt
import numpy as n

def graphdata(datasetx, datasety, xlabel=None,ylabel=None):
    plt.plot(datasetx,datasety,'red') #ro is scatterplot
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

def linreg(xdata,ydata,p):
    #calculating mean
    m=n.mean(ydata)

    #calculating linear regression coeficcients
    xy=0
    for i in range(0,len(xdata)):
        x=xdata[i]
        y=ydata[i]
        xy+=x*y
    
    #must have the same amount of entries in x and y
    b = float(((len(xdata) * xy) - (sum(xdata) * sum(ydata))) / ((len(xdata) * (sum(n.square(xdata)))) - (sum(xdata)**2)))
    bo = float((sum(ydata) - b * sum(xdata)) / len(xdata))

    #finding y values for linear trendline
    estvalues=[]
    lindata={}
    for x in xdata:
        y=(b*x)+bo
        estvalues.append(y)
    lindata["estvalues"]=[estvalues]
    lindata["m"]=b
    lindata["b"]=bo
    if p==1:
        plt.plot(xdata, estvalues)
    return lindata

#logrithmic regression
def logreg(xdata,ydata,p):
    xplotted=xdata
    xdata=n.log(xdata)
    #calculating logrithmic regression coefficients
    xy=0
    for i in range(0,len(xdata)):
        x=xdata[i]
        y=ydata[i]
        xy+=x*y
    
    #must have the same amount of entries in x and y
    b=float(((len(xdata)*xy)-(sum(xdata)*sum(ydata)))/((len(xdata)*(sum(n.square(xdata))))-(sum(xdata)**2)))
    bo= float((sum(ydata)-(b*sum(xdata)))/len(xdata))

    #finding y values for log trendline
    estvalues=[]
    for x in xdata:
        y=(b*x)+bo
        estvalues.append(y)
    if p==1:
        plt.plot(xplotted, estvalues, 'red')

    return estvalues

#calculating r2
def R2(xdata,ydata,yest):
    ybar=n.mean(ydata)

    #estimated values and normalvalues minus mean of y
    estd=[(y-ybar) for y in yest]
    normd=[(y-ybar) for y in ydata]
    
    #summing and squaring differences
    estdsum=sum(n.square(estd))
    normdsum=sum(n.square(normd))

    #r2 equation
    r2=estdsum/normdsum

    return r2

def main():
    x=[2,4,6,8,10,12,14]
    y=[0.0137,0.0145,0.0151,0.0151,0.0152,0.0154,0.0153]

    logreg(x,y,0)
    #saving data as png and compiling graphs 

if __name__ == "__main__":
    main()
