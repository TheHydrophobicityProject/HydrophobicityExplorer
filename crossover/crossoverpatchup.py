def crossover(xdata,ydata,plotlines=False,r2front=0.9,r2back=0.9,returnPlot=False,ax=None):
    import functionality as f
    import matplotlib.pyplot as plt
    
    frontxaxisbin=[]  #coordinates for linear reg line near orgin
    frontyaxisbin=[]
    
    backxaxisbin=[]   #coordinates for second linear reg line
    backyaxisbin=[]
    
    #calculating r2 values of linear from front and finding breakoff point 
    for i in range(0, len(xdata)):
        #building up the data
        #getting the ith datapoint
        xi=xdata[i]
        yi=ydata[i]
        
        #dropping it in the bin
        frontxaxisbin.append(xi)
        frontyaxisbin.append(yi)
        
        
        #ith r2 value
        #i>1 for defined r2 value
        if i>=0:
            # get the R^2 value for the line up to this point
            R2i=f.R2(frontxaxisbin,frontyaxisbin,f.linreg(frontxaxisbin,frontyaxisbin,0)["estvalues"][0])
            #checking for drop in r2 - if the new R^2 value is less than 0.9, stop here
            if R2i<r2front:
                frontxaxisbin.remove(xi)
                frontyaxisbin.remove(yi)
                break
            
    #TODO make more efficinent
    for i in range(0, len(xdata)):
        #building up the data (but backwards this time)
        #getting the ith datapoint
        # using j as an index starting from the END of the data
        j=(len(xdata)-1)-i
        xi=xdata[j]
        yi=ydata[j]
        
        #dropping it in the bin
        backxaxisbin.append(xi)
        backyaxisbin.append(yi)
        
        
        #ith r2 value
        #i>1 for defined r2 value

        if i>=1:
            R2i=f.R2(backxaxisbin,backyaxisbin,f.linreg(backxaxisbin,backyaxisbin,0)["estvalues"][0])
            #checking for drop in r2
            if R2i<r2back:
                break
    
    print(frontxaxisbin,backxaxisbin)
    #getting the linear regression outputs to graph
    frontregline=f.linreg(frontxaxisbin, frontyaxisbin,0)["estvalues"][0]
    backregline=f.linreg(backxaxisbin, backyaxisbin,0)["estvalues"][0]
    
    #linear regression equation coeficients - m for both front and back lines
    frontm=f.linreg(frontxaxisbin, frontyaxisbin,0)["m"]
    backm=f.linreg(backxaxisbin, backyaxisbin,0)["m"]
    
    # b for both front and back lines
    frontb=f.linreg(frontxaxisbin, frontyaxisbin,0)["b"]
    backb=f.linreg(backxaxisbin, backyaxisbin,0)["b"]
    
    frontline=[]
    backline=[]
    for i in xdata:
        yf=(frontm*i)+frontb
        frontline.append(yf)
    for i in xdata:
        yb=(backm*i)+backb
        backline.append(yb)
        
    #plotting lines (if plotting lines on the graph is enabled)
    if plotlines is True:
        plt.plot(xdata,frontline)
        plt.plot(xdata,backline)
        plt.plot(xdata,ydata,'o')
        plt.show()
    
    # modify the plot itself if the user asked for that
    if returnPlot is True:
        ax = ax or plt.gca()
        ax.plot(xdata,frontline)
        ax.plot(xdata,backline)
        ax.plot(xdata,ydata,'o')


#x=[2,4,6,8,10,12,14]
x=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #number of monomer units

#y=[0.0137,0.0145,0.0151,0.0151,0.0152,0.0154,0.0153]
y=[0.0124, 0.0136, 0.014, 0.0143, 0.0145, 0.0148, 0.0149, 0.0148, 0.015, 0.0149] #LogP/SA values for polystyrene

crossover(x,y,1)  #1 enables plotting of lines
        
        