import os
import sys
import h5py as h5
import lsqfit
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt
import string
import fitters



def plot_data(y_data,Pcm,parameter):
	fig = plt.figure(figsize=(7,4.326237))
	ax = plt.axes([0.15,0.15,0.8,0.8])
	x = np.arange(0,48,1)
	if parameter == "M_eff":
		meff = np.log(y_data/np.roll(y_data,-1))
		y = [i.mean for i in meff]
		err = [i.sdev for i in meff]
		ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',fillstyle='none',markersize='5',elinewidth=1,capsize=2)
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlim([0, 24])
		ax.set_ylim([0.35, 0.55])
		ax.set_xlabel('t', fontsize=20)
		ax.set_ylabel("$M_{eff}$", fontsize=20)
		plt.draw()
		#plt.show()
		plt.savefig(Pcm+"_"+parameter+".eps")
		plt.close()
	elif parameter == "A_eff":
		meff = np.log(y_data/np.roll(y_data,-1))
		E0 = [i.mean for i in meff]
		y_const = y_data*np.exp(E0*x)
		y = [i.mean for i in y_const]
		err = [i.sdev for i in y_const]
		ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',fillstyle='none',markersize='5',elinewidth=1,capsize=2)
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlabel('t', fontsize=20)
		ax.set_ylabel("$A_{eff}$", fontsize=20)
		ax.set_xlim([0, 24])
		ax.set_ylim([0.0, 0.0000025])
		plt.draw()
		#plt.show()
		plt.savefig(Pcm+"_"+parameter+".eps")
		plt.close()
	elif parameter == "E_1":
		meff = np.log(y_data/np.roll(y_data,-1))
		E0 = [i.mean for i in meff]
		A_0 = y_data*np.exp(E0*x)
		C_1 = y_data - A_0*(np.exp(-meff*x)+np.exp(-meff*(48-x)))
		E_1 = np.log(C_1/np.roll(C_1,-1))
		y = [i.mean for i in E_1]
		err = [i.sdev for i in E_1]
		ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',fillstyle='none',markersize='5',elinewidth=1,capsize=2)
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlabel('t', fontsize=20)
		ax.set_ylabel("$E_{1}$", fontsize=20)
		ax.set_xlim([0, 24])
		ax.set_ylim([-1.0, 1.0])
		plt.draw()
		#plt.show()
		plt.savefig(Pcm+"_"+parameter+".eps")
		plt.close()
	elif parameter == "A_1":
		meff = np.log(y_data/np.roll(y_data,-1))
		E0 = [i.mean for i in meff]
		A_0 = y_data*np.exp(E0*x)
		C_1 = y_data - A_0*(np.exp(-meff*x)+np.exp(-meff*(48-x)))
		E_1 = np.log(C_1/np.roll(C_1,-1))
		y_const = C_1*np.exp(E_1*x)
		y = [i.mean for i in y_const]
		err = [i.sdev for i in y_const]
		ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',fillstyle='none',markersize='5',elinewidth=1,capsize=2)
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlabel('t', fontsize=20)
		ax.set_ylabel("$A_{eff}$", fontsize=20)
		ax.set_xlim([0, 24])
		#ax.set_ylim([0.0008, 0.001])
		ax.set_ylim([0.00, 0.000005])
		plt.draw()
		#plt.show()
		plt.savefig(Pcm+"_"+parameter+".eps")
		plt.close()
	elif parameter == "Corr":
		y = [i.mean for i in y_data]
		err = [i.sdev for i in y_data]
		ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',fillstyle='none',markersize='5',elinewidth=1,capsize=2)
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlabel('t', fontsize=20)
		ax.set_ylabel("$C_{\pi\pi}$", fontsize=20)
		ax.set_xlim([15, 34])
		ax.set_ylim([0.000, 0.00000000025])
		plt.draw()
		#plt.show()
		plt.savefig(Pcm+"_"+parameter+".eps")
		plt.close()


def plot_fit(data,fit_range,params,Pcm,shell,bs,nstates,particles):
	#plots the data alongside the fit with fit error filled in
	print "plotting fit"
	y = [i.mean for i in data]
	x = fit_range
	err = [i.sdev for i in data]
	fig = plt.figure(figsize=(7,4.326237))
	ax = plt.axes([0.15,0.15,0.8,0.8])
	#plot the data:
	plot_label = Pcm+" bs "+bs+" fit"
	ax.errorbar(x=x,y=y,yerr=err,ls='None',color = 'black',marker='o',label=Pcm)
	#start = fit_range[0]
	start = 0
	#end = fit_range[-1]
	end = 48
	xc = dict()
	xc["x"] = np.arange(start,end,0.1)
	xc["sink"] = "s"
	errc = np.zeros_like(xc["x"])
	ft = []
	if particles == "pion":
		#fitc = fitters.pion_fit_function(nstates=nstates,T=len(data),sinks = ["s"])
		fitc = fitters.pion_fit_function(nstates=nstates,T=48,sinks = ["s"])
		ft = fitc.func(xc, params)
		ft_avg = [i.mean for i in ft]
		ft_std = [i.sdev for i in ft]
		ft_min = np.subtract(ft_avg,ft_std)
		ft_max = np.add(ft_avg, ft_std)
		ax.errorbar(x=xc["x"],y=ft_avg,yerr = errc, color='blue',label=plot_label)
		ax.fill_between(xc["x"], ft_min, ft_max, color = 'blue', alpha = '0.25')
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlim([17, 24])
		ax.set_ylim([0.000, 0.00002])
	elif particles == "pipi":
		#fitc = fitters.pipi_fit_function(nstates=nstates,T=len(data),sinks = ["s"])
		fitc = fitters.pipi_fit_function(nstates=nstates,T=48,sinks = ["s"])
		ft = fitc.func(xc, params)
		ft_avg = [i.mean for i in ft]
		ft_std = [i.sdev for i in ft]
		ft_min = np.subtract(ft_avg,ft_std)
		ft_max = np.add(ft_avg, ft_std)
		ax.errorbar(x=xc["x"],y=ft_avg,yerr = errc, color='blue',label=plot_label)
		ax.fill_between(xc["x"], ft_min, ft_max, color = 'blue', alpha = '0.25')
		ax.xaxis.set_tick_params(labelsize=16)
		ax.yaxis.set_tick_params(labelsize=16)
		ax.set_xlim([17, 24])
		ax.set_ylim([0.000, 0.00000000025])
	ax.set_xlabel('$t$', fontsize=20)
	ax.set_ylabel('$Corr$', fontsize=20)
	ax.legend(loc='upper left')
	plt.draw()
	#plt.show()
	plt.savefig("Pcmsq_"+Pcm+"_prelsq_"+str(shell)+"_"+str(particles)+"_"+str(bs)+"_fit_n_states_"+str(nstates)+".eps")
	plt.close()


def plot_single_time_slice(data, time, Pcm):
	data = data[:,time]
	#hist = np.histogram(data,bins='sqrt')
	start = np.amin(data)
	stop = np.amax(data)
	plt.hist(data, bins=100, range=[start, stop], histtype='step')
	plt.draw()
	plt.savefig(Pcm+"_t_"+str(time)+".eps")
	plt.close()





