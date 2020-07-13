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
import matplotlib.dates as mdates
import seaborn as sns
from scipy import stats
import netCDF4
import calendar
from return_time_plot import *


def read_data(fname,variable):
    # Read in the dataset
    data=netCDF4.Dataset(fname)
    times=data.variables['time']
    time=netCDF4.num2date(times[:],times.units,times.calendar)
    vals=data.variables[variable.lower()][:]
    vals=np.array(vals.flatten())
    # Below included for monthly datetime time variable and was included because initial read of time variable did not work correclty. 
    # Change date ranges below as required to fit input data
    time=np.arange('1979-01', '2020-01', dtype='datetime64[M]')

    return vals,time

def plot_tseries(variable,data,time,fig):
    ax = plt.subplot2grid((3,4),(0,1),colspan=3)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)
    plt.title("(b) "+variable)
    decades=mdates.YearLocator(5,month=1,day=1)
    decades_fmt=mdates.DateFormatter('%Y')
    years=mdates.YearLocator()
    x=np.arange(0,len(time))
    ## Plot the data
    ax.plot(time,data,linewidth=2.5,label=variable)
    ax.xaxis.set_major_locator(decades)
    ax.xaxis.set_major_formatter(decades_fmt)
    ax.xaxis.set_minor_locator(years)
    x=np.arange(0,len(time))
    z = np.polyfit(x,data,1)
    y=z[0]*x+z[1]
    ax.plot(time,y,color="black")


def plot_distribution(variable,data,sm,em,fig):
    ax = plt.subplot2grid((3,4),(0,0)) 
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)

    # Convert month abbreviation to index
    s=list(calendar.month_abbr).index(sm)
    e=list(calendar.month_abbr).index(em)

    # Reshape the data and mean the chosen months.  
    # Note this reshaping will be different if daily data is used and will need to be amended accordingly.
    vals=np.reshape(data,(-1,12))
    
    if s!=e:
        plt.title("(a) "+variable+" "+sm+"-"+em+" mean")
        sel_vals=vals[:,s:e]
        # Change np.mean to np.max or np.min as required.
        sel_plot=np.mean(sel_vals,1)
    else:
        plt.title("(a) "+variable+" "+sm+" mean")
        sel_plot=vals[:,s]


    # Plot the data
    sns.distplot(np.array(sel_plot), kde=False,fit=stats.genextreme, color="#00b478",label=variable,fit_kws={"linewidth":2.5,"color":"#00b478"})

def plot_seasonal_cycle(variable,data,time,year,fig):
    ax = plt.subplot2grid((3,4),(1,0),colspan=2,rowspan=2)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)

    # Note this reshaping will be different if daily data is used and will need to be amended accordingly.
    vals=np.reshape(data,(-1,12))
    times=np.reshape(time,(-1,12))

    mm_val=np.mean(vals[:-1,:],0)
    print(mm_val.shape)
    mm_std=np.std(vals[:-1,:],0)
    
    ## Plot the data
    months=[calendar.month_abbr[i] for i in range(1,13)]
    plt.title("(c) "+variable)
    ax.plot(months,mm_val,linewidth=3,color="orange",zorder=1)
    ax.fill_between(months,mm_val+mm_std,mm_val-mm_std,color="gold",alpha=0.3)

    for i in range(0,40):
	# Adjust 1979 to reflect the start year of the data
        if i+1979 in [year]:
                col="black"
                width=1
        else:
                col="silver"
                width=0.5
        ax.plot(months,vals[i,:].flatten(),linewidth=width,color=col)

def plot_return_time(variable,data,time,year,dirn,sm,em,fig):
    # Plot the return time for the selected period
    ax = plt.subplot2grid((3,4),(1,2),colspan=2,rowspan=2)
    plt.title("(d) "+variable+" "+sm+"-"+em+" mean")
    ax.set_ylabel(variable,fontsize=16)
    ax.set_xlabel("Chance of event occurring in a given year",fontsize=16)
    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)
    ax.set_xlim(1,1e2)

    # Note this reshaping will be different if daily data is used and will need to be amended accordingly.
    vals=np.reshape(data,(-1,12))
    times=np.reshape(time,(-1,12))

    # Convert month abbreviation to index
    s=list(calendar.month_abbr).index(sm)
    e=list(calendar.month_abbr).index(em)

    if s!=e:
        plt.title("(d) "+variable+" "+sm+"-"+em+" mean")
        sel_vals=vals[:,s:e]
        # Change np.mean to np.max or np.min as required.
        sel_plot=np.mean(sel_vals,1)
    else:
        plt.title("(d) "+variable+" "+sm+" mean")
        sel_plot=vals[:,s]

    # Calculate the threshold for the given year.  Adjust the start year of 1979 in the data as required
    threshold=sel_plot[year-1979]
    print(str(year)+" Threshold: "+str(threshold))

    plot_rt(sel_plot.flatten(),["Orange","Gold"],"ERA5",ax,"both",dirn,[threshold],[str(year)+" event"])
    labels=['','1/1','1/10','1/100']
    ax.set_xticklabels(labels)

    plt.setp(ax.get_xticklabels(),fontsize=16)
    plt.setp(ax.get_yticklabels(),fontsize=16)

def plot_rt(data,cols,plabel,ax,errb,dirn,thresholds,tlabels):
    # Plot the return times with bootstrap 9-95% confience intervals.
    y_data_all, x_data_all = calc_return_times(data,direction=dirn,period=1)
    
    conf_all = calc_return_time_confidences(data,direction=dirn,bsn=1000)
    conf_all_x = calc_return_time_confidences(x_data_all,direction="descending",bsn=1000)

    l1=ax.semilogx(x_data_all,y_data_all, marker='o',markersize=4,
                       linestyle='None',mec=cols[0],mfc=cols[0],
                       color=cols[0],fillstyle='full',
                       label=plabel,zorder=2)
    
    conf_all_5=conf_all[0,:].squeeze()
    conf_all_95=conf_all[1,:].squeeze()

    conf_all_x_5=conf_all_x[0,:].squeeze()
    conf_all_x_95=conf_all_x[1,:].squeeze()

    if errb=="both":
    	cl0=ax.fill_between(x_data_all,conf_all_5,conf_all_95,color=cols[1],alpha=0.2,linewidth=1.,zorder=0)
    if errb=="magnitude" or errb=="both":
    	cl1=ax.semilogx([x_data_all,x_data_all],[conf_all_5,conf_all_95],color=cols[1],linewidth=1.,zorder=1)
    if errb=="return_time" or errb=="both":
	    cl2=ax.semilogx([conf_all_x_5,conf_all_x_95],[y_data_all,y_data_all],color=cols[1],linewidth=1.,zorder=1)

    # Plot GEV fit to data
    shape,loc,scale=stats.genextreme.fit(data)
    T=np.r_[1:10000]*0.01
    if dirn=="ascending":
        rt=stats.genextreme.ppf(1./T,shape,loc=loc,scale=scale)
    else:
        rt=stats.genextreme.isf(1./T,shape,loc=loc,scale=scale)
    l1=ax.semilogx(T,rt,color="black",label="GEV fit")
    
    xmin,xmax=ax.get_xlim()
    ymin,ymax=ax.get_ylim()

    styles=['--',':']
    for i,threshold in enumerate(thresholds):
        ax.semilogx([xmin,xmax],[threshold,threshold],color="Silver",linestyle=styles[i],linewidth=2.5,label=tlabels[i], zorder=2)
        nidx=find_nearest(y_data_all,threshold)
        ax.axvspan(conf_all_x_5[nidx],conf_all_x_95[nidx],ymin=0,ymax=(threshold-ymin)/(ymax-ymin),facecolor='silver',alpha=0.5,zorder=0)

    ll=ax.legend(loc="upper right",prop={"size": 14},fancybox=True,numpoints=1)

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx

#Main controling function
def main():
    #Read in the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--fname", help="The filename including path of the 1D data.")
    parser.add_argument("--variable", help="The variable name")
    parser.add_argument("--year", type=int, help="The year of interest that the extreme event ocurred in")
    parser.add_argument("--dirn", default="descending", help="The direction of the variable.  'Descending' for a variable where the threshold is exceeded (default) and 'ascending' where the threshold is less than a given value.")
    parser.add_argument("--start_month", default="Jan", help="The start month for the distribution and return time plot, default 'Jan'")
    parser.add_argument("--end_month", default="Dec", help="The end month for the distribution and return time plot, default 'Dec'")
    args = parser.parse_args()

    # Set up the plot
    font = {'family' : 'sans-serif',
            'size'   : 20}

    matplotlib.rc('font', **font)
    fig = plt.figure()
    fig.set_size_inches(21,15)

    data,time=read_data(args.fname,args.variable)

    plot_distribution(args.variable,data,args.start_month, args.end_month,fig)
    plot_tseries(args.variable,data,time,fig)
    plot_seasonal_cycle(args.variable,data,time,args.year,fig)
    plot_return_time(args.variable,data,time,args.year,args.dirn,args.start_month,args.end_month,fig)
    plt.tight_layout()
    fig.savefig("Timeseries_analysis_"+args.variable+".png",dpi=28.75*2)
	
    print('Finished!')

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()
