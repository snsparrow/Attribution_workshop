#!/usr/bin/env python2.7
#############################################################################
# Program : plot_q_q.py
# Author  : Sarah Sparrow
# Date    : 26/10/2020
# Purpose : Quantile-quantile plot adapted from Neil Massey q_q_plot.py (in https://github.com/CPDN-git/cpdn_analysis)
#############################################################################

import sys, os
import argparse
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import netCDF4
import calendar

def read_data(var,fname):
    print("Reading "+fname)
    # Read in the dataset
    data=netCDF4.Dataset(fname)
    vals=data.variables[var][:]
    vals=np.array(vals.flatten())
    return vals

def q_q_plot(obs, ens, p_tiles,ptitle):
    sp = plt.subplot2grid((1,1),(0,0))
    plt.title(ptitle)
    # plot the (non-conditional) quantiles
    obs_quantiles = []
    ens_quantiles = []
    for i in range(1, 99):
        obs_quantiles.append(stats.scoreatpercentile(obs.flatten(), i))
        ens_quantiles.append(stats.scoreatpercentile(ens.flatten(), i))
    sp.plot(obs_quantiles, ens_quantiles, 'k--', lw=1.5, zorder=1)
    ens_range = ens_quantiles[-1] - ens_quantiles[0]
    obs_range = obs_quantiles[-1] - obs_quantiles[0]
    # plot the percentiles in the p_tiles list)
    for p in p_tiles:
        obs_p = stats.scoreatpercentile(obs.flatten(), p)
        ens_p = stats.scoreatpercentile(ens.flatten(), p)
        sp.plot(obs_p, ens_p, 'r+', ms=12, zorder=1)
        sp.text(obs_p+0.03*obs_range, ens_p-0.03*ens_range, str(p),c='r', ha='center', va='bottom',
                zorder=1, fontsize=18)

    smallest_x = np.min([obs_quantiles[0], ens_quantiles[0]])
    largest_x = np.max([obs_quantiles[-1], ens_quantiles[-1]])

    # 1:1 line - lowest z order
    sp.plot([smallest_x, largest_x], [smallest_x, largest_x], 'k', lw=2.0,
            zorder=0)

    sp.set_xlabel("Observations",fontsize=22)
    sp.set_ylabel("Model",fontsize=22)
    # limits
    sp.set_xlim([smallest_x, largest_x])
    sp.set_ylim([smallest_x, largest_x])
    sp.set_aspect(1.0)
    return sp

#Main controling function
def main():
    #Read in the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--obs_var", help="The observations variable name")
    parser.add_argument("--model_var", help="The model variable name")
    parser.add_argument("--obs_file", help="The netCDF file containing the observations")
    parser.add_argument("--model_file", help="The netCDF file containing the model output") 
    parser.add_argument("--plot_title", default="", help="Title string for plot")
    args = parser.parse_args()

    # Set up the plot
    font = {'family' : 'sans-serif',
            'size'   : 20}

    matplotlib.rc('font', **font)
    fig = plt.figure()
    fig.set_size_inches(8,8)

    dataObs=read_data(args.obs_var,args.obs_file)
    dataModel=read_data(args.model_var, args.model_file)
    
    q_q_plot(dataObs, dataModel, [5,10,25,50,75,90,95],args.plot_title)
    plt.tight_layout()
    fig.savefig(args.obs_var+"_q_q_plot.png",dpi=28.75*2)

    print('Finished!')

#Washerboard function that allows main() to run on running this file
if __name__=="__main__":
  main()

