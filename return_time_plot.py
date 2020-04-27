#! /usr/bin/env python
###############################################################################
# Program : return_time_plot.py
# Author  : Neil Massey
# Date    : 23/01/12
# Purpose : plot return time periods
###############################################################################

import sys, os
EU = os.path.expanduser
sys.path.append(EU("~massey/DPhil/Coding/python_lib"))

import matplotlib.pyplot as plt
from scipy import stats
import numpy
import random

###############################################################################

def calc_return_times(em, direction="ascending", period=1):
	ey_data = em.flatten()
	ey_data.sort()
	# reverse if necessary
	if direction == "descending":	# being beneath a threshold value
		ey_data = ey_data[::-1]
	# create the n_ens / rank_data
	val = float(len(ey_data)) * 1.0/period
	end = float(len(ey_data)) * 1.0/period
	start = 1.0
	step = (end - start) / (len(ey_data)-1)
	ranks = [x*step+start for x in range(0, len(ey_data))]
	ex_data = val / numpy.array(ranks, dtype=numpy.float32)
	return ey_data, ex_data

###############################################################################

def calc_return_time_confidences(em, direction="ascending", c=[0.05, 0.95], bsn=1e5):
	# c = confidence intervals (percentiles) to calculate
	# bsn = boot strap number, number of times to resample the distribution
	ey_data = em.flatten()
	# create the store
	sample_store = numpy.zeros((int(bsn), ey_data.shape[0]), 'f')
	# do the resampling
	for s in range(0, int(bsn)):
		t_data = numpy.zeros((ey_data.shape[0]), 'f')
		for y in range(0, ey_data.shape[0]):
			x = random.uniform(0, ey_data.shape[0])
			t_data[y] = ey_data[int(x)]
		t_data.sort()
		# reverse if necessary
		if direction == "descending":
			t_data = t_data[::-1]
		sample_store[s] = t_data
	# now for each confidence interval find the  value at the percentile
	conf_inter = numpy.zeros((len(c), ey_data.shape[0]), 'f')
	for c0 in range(0, len(c)):
		for y in range(0, ey_data.shape[0]):
			data_slice = sample_store[:,y]
			conf_inter[c0,y] = stats.scoreatpercentile(data_slice, c[c0]*100)
	return conf_inter

###############################################################################

def return_time_plot(plt_obj, ens, direction="ascending", period=1):
	# first calculate the return time data
	y_data, x_data = calc_return_times(ens, direction, period)
	l1 = plt_obj.semilogx(x_data, y_data, 'ko', marker='o', mec='k', mfc='w')
	return l1

###############################################################################

def multi_return_time_plot(plt_obj, ens, direction="ascending", highlight=-2e20, period=1):
	# first calculate the return time data
	y_data = []
	x_data = []
	for em in ens:
		ey_data, ex_data = calc_return_times(em, direction, period)
		y_data.append(ey_data)
		x_data.append(ex_data)
	# draw each return time plot
	lines = []
	colors = ['b', 'r', 'k']
	for e in range(0, len(y_data)):
		l = plt_obj.semilogx(x_data[e], y_data[e], 'ko', marker='o', mec=colors[e], mfc='w',
							 ms = 4.0, zorder=1)
		# add an extra marker at the highlight
		if highlight != -2e20:
			# search for the nearest location
			if direction == "descending":
				fi = numpy.interp([highlight], y_data[e][::-1], x_data[e][::-1])
			else:
				fi = numpy.interp([highlight], y_data[e], x_data[e])
			print(fi)
			# plot the highlight
			plt_obj.plot(fi, highlight, 'ko', marker='o', 
						 mec=colors[e], mfc=colors[e], ms = 8.0, zorder=2,
						 fillstyle='full')

		lines.append(l)

	return lines

###############################################################################

def multi_return_time_plot_confidence(ens, direction="ascending", c=0.05, period=1):
	# multi return time plot with confidence interval between c and 1.0-c
	# first calculate the actual real return data
	y_data = []
	x_data = []
	c_data = []
	for em in ens:
		ey_data, ex_data = calc_return_times(em, direction, period)
		confs = calc_return_time_confidences(ey_data, direction, [c,1-c], 1e4)
		y_data.append(ey_data)
		x_data.append(ex_data)
		c_data.append(confs)
	# draw each return time plot
	sp = plt.subplot(111)
	lines = []
	colors = ['b', 'r', 'k']
	for e in range(0, len(y_data)):
		l = sp.semilogx(x_data[e], y_data[e], color=colors[e], lw=2.0)
		sp.fill_between(x_data[e], c_data[e][0], c_data[e][1], facecolor=colors[e],
						edgecolor=colors[e], alpha=0.5)
		lines.append(l)

	return sp, lines

###############################################################################

def mega_multi_return_time_plot(rt_pairs, direction="ascending", period=1):
	# produce comparative return time plots between two sets of data,
	# e.g. bias-corrected and non-bias-corrected data
	# first calculate the return time data
	y_data_pairs = []
	x_data_pairs = []
	for rt_p in rt_pairs:
		y_data = []
		x_data = []
		for em in rt_p:
			ey_data, ex_data = calc_return_times(em, direction, period)
			y_data.append(ey_data)
			x_data.append(ex_data)
		y_data_pairs.append(y_data)
		x_data_pairs.append(x_data)
	# draw each return time plot
	sp = plt.subplot(111)
	lines = []
	colors = ['k', 'r', 'b']
	p = 0
	for y_data in y_data_pairs:
		sp.fill_between(x_data[0], y_data[0], y_data[1], facecolor=colors[p], 
						edgecolor=colors[p], alpha=0.5)
		l = sp.plot(x_data[0], y_data[0], color=colors[p], alpha=1.0)
		sp.plot(x_data[0], y_data[1], color=colors[p], alpha=1.0)
		lines.append(l)
		p += 1
	sp.set_xscale("log")

	return sp, lines

###############################################################################

def find_event(ens, ev_period, direction="ascending", period=1):
	# find the event value in the first ensemble and the equivalent return period in
	# the other ensembles
	# first calculate the return time data
	y_data = []
	x_data = []
	for em in ens:
		ey_data, ex_data = calc_return_times(em, direction, period)
		y_data.append(ey_data)
		x_data.append(ex_data)

	# find the ev_period event in the x_data[0] and y_data[0]
	ev_pos = numpy.intersect1d(numpy.where(x_data[0] <= ev_period)[0],
							   numpy.where(x_data[0] >= ev_period)[0])
	ev_val = y_data[0][ev_pos]
	ev_yr_equivs = [ev_period]
	# now search in the other ensembles for this value
	for e in range(1, len(ens)):
		ev_pos = -1
		for i in range(0, y_data[e].shape[0]-1):
			if direction == "descending":
				if y_data[e][i] >= ev_val and y_data[e][i+1] < ev_val:
					ev_pos = i
					break
			elif direction == "ascending":
				if y_data[e][i] <= ev_val and y_data[e][i+1] > ev_val:
					ev_pos = i
					break
		if ev_pos == -1:
			yr_equiv = x_data[e][0]
		else:
			yr_equiv = x_data[e][ev_pos]
		ev_yr_equivs.append(yr_equiv)

	return ev_val, ev_yr_equivs

###############################################################################

def get_return_time(ens, value, direction="ascending", period=1):
	# get the return time for a given value
	y_data, x_data = calc_return_times(ens, direction, period)
	ev_pos = -1
	for e in range(0, y_data.shape[0]-1):
		if direction == "descending":
			if y_data[e] >= value and y_data[e+1] < value:
				ev_pos = e
				break
		elif direction == "ascending":
			if y_data[e] <= value and y_data[e+1] > value:
				ev_pos = e
				break
	ev_yr = x_data[ev_pos]
	return ev_yr
