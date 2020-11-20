#/usr/bin/env python
#############################################################################
# Program : plot_timeseries.py
# Author  : Sarah Sparrow
# Date    : 31/01/2020
# Purpose : Plot ERA5 monthly mean time series
#############################################################################

import sys, os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from return_time_plot import *
from calc_PR import *

def read_data(fname,restore):
    # Read or restore the data
    if restore:
        vals=np.load(fname)
        # adjust tempture thresholds into C from K
        vals=vals-273.14
    else:
        # Read in the dataset
        infile = open(fname, 'r')
        data = [float(line) for line in infile.readlines()]
        infile.close()
        vals=np.array(data,float)

    return vals.flatten()

def plot_return_time(variable,dataAct,dataNat,threshold,dirn,fig):
    # Setup the plot parameters and axis limits
    ax = plt.subplot2grid((1,3),(0,0))
    plt.title("(a) "+variable)
    ax.set_ylabel(variable,fontsize=16)
    ax.set_xlabel("Chance of event occurring in a given year",fontsize=16)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)
    ax.set_xlim(1,1e2)
    ax.set_ylim(28.0, 38.0)

    # Plot the return time for the historical and historicalNat simulations
    plot_rt(dataAct,["royalblue","cornflowerblue","mediumblue"],"Historical",ax,"both",dirn,threshold,"",'--')
    plot_rt(dataNat,["orange","gold","darkorange"],"HistoricalNat",ax,"both",dirn,threshold,"Threshold",'--')
    labels=['','1/1','1/10','1/100']
    ax.set_xticklabels(labels)

    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)

    ll=ax.legend(loc="upper left",prop={"size": 14},fancybox=True,numpoints=1)


def plot_rt(data,cols,plabel,ax,errb,dirn,threshold,tlabel,tstyle):
    # Plot the return times with bootstrap 9-95% confience intervals.
    # Calculate the return times
    y_data_all, x_data_all = calc_return_times(data,direction=dirn,period=1)
    
    # Calculate the bootstrap confidences in both the x and y directions
    conf_all = calc_return_time_confidences(data,direction=dirn,bsn=1000)
    conf_all_x = calc_return_time_confidences(x_data_all,direction="descending",bsn=1000)

    # Plot  the return time curve
    l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=4,
                       linestyle='None',mec=cols[0],mfc=cols[0],
                       color=cols[0],fillstyle='full',
                       label=plabel,zorder=2)
    
    conf_all_5=conf_all[0,:].squeeze()
    conf_all_95=conf_all[1,:].squeeze()

    conf_all_x_5=conf_all_x[0,:].squeeze()
    conf_all_x_95=conf_all_x[1,:].squeeze()

    ax.grid(b=True,which='major')
    ax.grid(b=True, which='minor',linestyle='--')

    # Plot the error bars onn the return times
    if errb=="both":
    	cl0=ax.fill_between(x_data_all,conf_all_5,conf_all_95,color=cols[1],alpha=0.2,linewidth=1.,zorder=0)
    if errb=="magnitude" or errb=="both":
    	cl1=ax.semilogx([x_data_all,x_data_all],[conf_all_5,conf_all_95],color=cols[1],linewidth=1.,zorder=1)
    if errb=="return_time" or errb=="both":
	    cl2=ax.semilogx([conf_all_x_5,conf_all_x_95],[y_data_all,y_data_all],color=cols[1],linewidth=1.,zorder=1)

    # Calculate GEV fit to data 
    shape,loc,scale=stats.genextreme.fit(data)
    T=np.r_[1:10000]*0.01
    # Perform K-S test and print goodness of fit parameters
    D, p = stats.kstest(data.flatten(), 'genextreme', args=(shape, loc,scale));
    print(plabel+' GEV fit, K-S test parameters p: '+str(p)+" D: "+str(D))
    # Plot fit line
    if dirn=="ascending":
        rt=stats.genextreme.ppf(1./T,shape,loc=loc,scale=scale)
    else:
        rt=stats.genextreme.isf(1./T,shape,loc=loc,scale=scale)
    l1=ax.semilogx(T,rt,color=cols[2],label="GEV fit")
    
    # Highlight where the threshold is and where the return time curve bisects this threshold
    xmin,xmax=ax.get_xlim()
    ymin,ymax=ax.get_ylim()
    ax.semilogx([xmin,xmax],[threshold,threshold],color="Grey",linestyle=tstyle,linewidth=2.5,label=tlabel, zorder=2)
    nidx=find_nearest(y_data_all,threshold)
    ax.axvspan(conf_all_x_5[nidx],conf_all_x_95[nidx],ymin=0,ymax=(threshold-ymin)/(ymax-ymin),facecolor='silver',edgecolor=cols[2],linewidth=2.,alpha=0.3,zorder=0)

    # Reeturn the fit parameters
    return [shape,loc,scale]

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx


def plot_PR_rt(PR_median, PR_max,PR_min, ev_periods,act_params,nat_params,fig):
    ax = plt.subplot2grid((1,3),(0,1))
    plt.title("(b) PR return times")

    cols=["Indigo","MediumSlateBlue"]

    ax.set_ylabel("Probability Ratio",fontsize=16)
    ax.set_xlabel("Chance of event occurring in a given year",fontsize=16)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)
    ax.set_xlim(1,1e2)
    ax.set_ylim(1,32)

    ax.grid(b=True,which='major')
    ax.grid(b=True, which='minor',linestyle='--')

    l1=ax.loglog(ev_periods,PR_median,basey=2, marker='o',markersize=3,
                       linestyle='None',mec=cols[0],mfc=cols[0],
                       color=cols[0],fillstyle='full',zorder=2)

    cl1=ax.semilogx([ev_periods,ev_periods],[PR_min,PR_max],color=cols[1],linewidth=1.,zorder=1,alpha=0.3)

    # Add fit line
    T=np.r_[1:1001]
    rt=stats.genextreme.isf(1./T,act_params[0],loc=act_params[1],scale=act_params[2])
    cdf=stats.genextreme.cdf(rt,act_params[0],loc=act_params[1],scale=act_params[2])
    probAct=1-cdf
    cdf=stats.genextreme.cdf(rt,nat_params[0],loc=nat_params[1],scale=nat_params[2])
    probNat=1-cdf
    PR=probAct/probNat
    l2=ax.loglog(T,PR,basey=2,color=cols[1],linewidth=2,zorder=1, label="GEV fit")

    ll=ax.legend(loc="upper left",prop={"size": 14},fancybox=True,numpoints=1)

    labels=['','1/1','1/10','1/100']
    ax.set_xticklabels(labels)

    plt.tight_layout()

def plot_PR_threshold(PR_median, PR_max,PR_min, thresholds,variable,act_params, nat_params,flevel,fig):
    ax = plt.subplot2grid((1,3),(0,2))
    plt.title("(c) PR thresholds")

    cols=["#c80000","LightPink"]

    ax.set_ylabel("Probability Ratio",fontsize=16)
    ax.set_xlabel(variable,fontsize=16)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)
    ax.set_xlim(28.0,38.0)
    ax.set_ylim(1,32)

    ax.grid(b=True,which='major')
    ax.grid(b=True, which='minor',linestyle='--')

    l1=ax.semilogy(thresholds,PR_median, basey=2,marker='o',markersize=3,
                       linestyle='None',mec=cols[0],mfc=cols[0],
                       color=cols[0],fillstyle='full',zorder=2)

    cl1=ax.semilogy([thresholds,thresholds],[PR_min,PR_max],color=cols[1],linewidth=1.,zorder=1,alpha=0.3)

    # Add fit line
    T=np.r_[1:1001]
    rt=stats.genextreme.isf(1./T,act_params[0],loc=act_params[1],scale=act_params[2])
    cdf=stats.genextreme.cdf(rt,act_params[0],loc=act_params[1],scale=act_params[2])
    probAct=1-cdf
    cdf=stats.genextreme.cdf(rt,nat_params[0],loc=nat_params[1],scale=nat_params[2])
    probNat=1-cdf
    PR=probAct/probNat
    l2=ax.semilogy(rt,PR,basey=2,color=cols[1],linewidth=2,zorder=1, label='GEV fit')

    xmin,xmax=ax.get_xlim() 
    ymin,ymax=ax.get_ylim()
    ax.set_ylim(ymin,ymax)
    ax.semilogy([flevel,flevel],[ymin,ymax],basey=2,color="Grey",linestyle="--",linewidth=2.5,label="Threshold", zorder=2)
   
    nidx=find_nearest(thresholds,flevel)
    ax.axhspan(ymin=PR_min[nidx],ymax=PR_max[nidx],xmin=0, xmax=(flevel-xmin)/(xmax-xmin),facecolor='silver',edgecolor=cols[1],linewidth=2.,alpha=0.3,zorder=0)
    print("PR median (5% 95%):",PR_median[nidx],PR_min[nidx],PR_max[nidx])
    ll=ax.legend(loc="upper left",prop={"size": 14},fancybox=True,numpoints=1)

    plt.tight_layout()


#Main controling function
def main():
    #Read in the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--fname_act", help="The filename including path to the historical (actual) data.")
    parser.add_argument("--fname_nat", help="The filename including path to the historicalNat (natural) data.")
    parser.add_argument("--variable",  help="The variable name")
    parser.add_argument("--threshold", type=float, help="The year of interest that the extreme event ocurred in")
    parser.add_argument("--dirn", default="descending", help="The direction of the variable.  'Descending' for a variable where the threshold is exceeded (default) and 'ascending' where the threshold is less than a given value.")
    parser.add_argument("--restore", action="store_true", help="Restore the data")
    args = parser.parse_args()

    # Set up the plot
    font = {'family' : 'sans-serif',
            'size'   : 20}

    matplotlib.rc('font', **font)
    fig = plt.figure()
    fig.set_size_inches(21,6)

    dataAct=read_data(args.fname_act,args.restore)
    dataNat=read_data(args.fname_nat,args.restore)

    plot_return_time(args.variable,dataAct,dataNat,args.threshold,args.dirn,fig)

    PR_median, PR_max,PR_min, ev_periods, thresholds=get_PR_data(dataAct,dataNat,'example_data',args.restore)
    # adjust tempture thresholds into C from K
    thresholds=thresholds-273.14

    act_params=stats.genextreme.fit(dataAct)
    nat_params=stats.genextreme.fit(dataNat)

    plot_PR_rt(PR_median, PR_max,PR_min, ev_periods,act_params,nat_params,fig)
    plot_PR_threshold(PR_median, PR_max,PR_min, thresholds,args.variable,act_params, nat_params,args.threshold,fig)


    plt.tight_layout()
    fig.savefig("Return_time_PR_"+args.variable+".png",dpi=28.75*2)
	
    print('Finished!')

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
