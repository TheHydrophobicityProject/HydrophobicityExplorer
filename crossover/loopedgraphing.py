import numpy as np
import datacatch as data
import matplotlib.pyplot as plt
import functionality as f
import crossoverpatchup as c

def loop(x,yddict,file=0):    
    fig, axs = plt.subplots(len(yddict))
    fig.tight_layout()
    
    fig.text(0.5, 0.04, 'common xlabel', ha='center', va='center',color='red')
    fig.text(0.06, 0.5, 'common ylabel', ha='center', va='center', rotation='vertical',color='red')

    #looping through 2d list graphing the data crossovers
    #variable keeping track of index
    indextrack=0
    for key in yddict:
        #crossover unit to split trendline
        crossunit=c.crossover(x,yddict[key],0)

        #dataslpit trendlines
        linedata1=f.linreg(x[:crossunit],yddict[key][:crossunit],0)
        linedata2=f.linreg(x[crossunit-1:],yddict[key][crossunit-1:],0)

        #plotting on a subplot
        axs[indextrack].title.set_text(key)
        axs[indextrack].plot(x[:crossunit],linedata1)
        axs[indextrack].plot(x[crossunit-1:],linedata2)
        # axs[indextrack].set_ylabel('yaxis')
        # # axs[indextrack].set_xlabel('xaxis')
        axs[indextrack].plot(x,yddict[key],'o')

        #iterating to a newsubplot
        indextrack+=1
    
    if file==1:
        filename=input('Enter file title: ')
        plt.savefig(filename+'.pdf', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format='pdf',
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)
