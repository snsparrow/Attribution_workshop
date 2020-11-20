#!/usr/bin/env python2.7
#############################################################################
# Program : calc_PR.py
# Author  : Sarah Sparrow
# Date    : 05/11/2020
# Purpose : Functions for calculating probability ratios
#############################################################################

import sys
import numpy as np
from return_time_plot import *


def probability_ratio(dataHist,dataNat,threshold):
    HistCount=sum(np.array(dataHist)>threshold)
    NatCount=sum(np.array(dataNat)>threshold)
    P_hist=float(HistCount)/float(len(dataHist))
    P_nat=float(NatCount)/float(len(dataNat))
    try:
        PR=P_hist/P_nat
    except:
        PR=None

    return PR

def calc_PR_conf(actBoot,natBoot,percentile,threshold,bsn=10000):
    PR=[]
    for ib in range(0,bsn):
        PR.append(probability_ratio(actBoot[ib,:].flatten(),natBoot[ib,:].flatten(),threshold))
    PR_5=get_val_percentile(PR,percentile[0])
    PR_50=get_val_percentile(PR,percentile[1])
    PR_95=get_val_percentile(PR,percentile[2])
    return [PR_5,PR_50,PR_95]

def get_val_percentile(data, percentile):
    if len(data) < 1:
        value = None
    elif (percentile >= 100):
        sys.stderr.write('ERROR: percentile must be < 100.  you supplied: %s\n'%
 percentile)
        value = None
    else:
        element_idx = int(len(data) * (percentile / 100.0))
        data.sort()
        value = data[element_idx]
    return value

def get_thresholds(data):
    y_data, x_data = calc_return_times(data,direction="descending",period=1)
    return y_data,x_data

def get_PR_data(dataAct,dataNat,expt_dir,restore):
    
    #  Calculate the range of thresholds and return times
    thresholds,ev_periods=get_thresholds(dataAct)
    
    # Calculate the bootstrap ensembles
    actBoot = calc_bootstrap_ensemble(dataAct,direction="descending",bsn=10000, slen=525)
    natBoot = calc_bootstrap_ensemble(dataNat,direction="descending",bsn=10000, slen=525)

    # Initialise some blank  arrays
    PR_median=[]
    PR_min=[]
    PR_max=[]

    # Calculating confidence intervals by bootstrap on ensembles can be intensive so only recalculate if needed
    if restore==False:
        for threshold in thresholds:
                try:
                        [PR_5,PR_50,PR_95]=calc_PR_conf(actBoot,natBoot,[5,50,95],threshold)
                        PR_median.append(PR_50)
                        PR_min.append(PR_5)
                        PR_max.append(PR_95)
                except:
                        PR_median.append(0)
                        PR_min.append(0)
                        PR_max.append(0)

        if not os.path.exists(expt_dir):
            os.makedirs(expt_dir)
        
        np.save(expt_dir+'/PR_median.npy',np.array(PR_median))
        np.save(expt_dir+'/PR_min.npy',np.array(PR_min))
        np.save(expt_dir+'/PR_max.npy',np.array(PR_max))
        np.save(expt_dir+'/PR_ev_periods.npy',np.array(ev_periods))
        np.save(expt_dir+'/PR_thresholds.npy',np.array(thresholds))
    else:
         PR_median=np.load(expt_dir+'/PR_median.npy')
         PR_max=np.load(expt_dir+'/PR_max.npy')
         PR_min=np.load(expt_dir+'/PR_min.npy')
         ev_periods=np.load(expt_dir+'/PR_ev_periods.npy')
         thresholds=np.load(expt_dir+'/PR_thresholds.npy')

    return PR_median, PR_max, PR_min, ev_periods, thresholds

