import os
import sys
import h5py as h5
import lsqfit
import gvar as gv
import numpy as np
import datetime
import string



def fit_data(data,fit_range,particles,Pcm,nstates,nsinks,cov=None, print_fit = True):
	if nsinks ==1:
		if print_fit:
			print "fitting with "+str(nstates)+" states and smeared data"
		else:
			pass
	else:
		if print_fit:
			print "fitting with "+str(nstates)+" states and point and smeared data"
		else:
			pass
	x = fit_range
	y = data
	if cov is None:
		cov_matrix = gv.evalcov(y)
		#cov_matrix = gv.evalcov(y)/len(y)
		#cov_matrix = np.identity(len(x))
		if print_fit:
			print np.shape(cov_matrix)
		else:
			pass
	else:
		cov_matrix = cov
	y = gv.mean(y)
	#you'll need to update the pipi function to deal with number of sinks as well as nstates
	if particles == "pipi":
		p = pipi_priors(nstates,nsinks,Pcm)
		#fitc = pipi_fit_function(nstates=nstates,T=len(data))
		if nsinks == 1:
			fitc = pipi_fit_function(nstates=nstates,T=48,sinks = ["s"])
		else:
			fitc = pipi_fit_function(nstates=nstates,T=48,sinks = ["s","p"])
	elif particles == "pion":
		p = pion_priors(nstates,nsinks,Pcm)
		p0 = pion_priors(nstates,nsinks,Pcm)
		if nsinks == 1:
			fitc = pion_fit_function(nstates=nstates,T=48,sinks = ["s"])
		else:
			fitc = pion_fit_function(nstates=nstates,T=48,sinks = ["s","p"])
	#fit = lsqfit.nonlinear_fit(data=(x,y,cov_matrix),prior=p0,fcn=fitc.full_func)
	fit = lsqfit.nonlinear_fit(data=(x,y,cov_matrix),prior=p,fcn=fitc.full_func)
	if print_fit:
		print fit
	else:
		pass
	return fit


def pion_priors(nstates,nsinks,Pcm): 
	p = dict()
	#p['Cs0'] = [0.0,1.0]
	p['Es0'] = [0.237,0.002]
	p['Bs1'] = [-1000.,1005.0]
	p['Bs2'] = [-1.0E6,5.0]
	p['As0'] = [0.001,0.0001]
	p['As1'] = [0.,8E-4]
	p['As2'] = [0.,8E-6]
	#p['Cp0'] = [0.0,1.0]
	p['Ep0'] = [0.25,0.05]
	p['Bp1'] = [-100,50]
	p['Bp2'] = [-100,50]
	p['Ap0'] = [9E-4,0.001]
	p['Ap1'] = [0.,8E-4]
	p['Ap2'] = [0.,8E-6]
	prior = dict()
	if nsinks == 1:
		#sinks = ['s']
		sinks = ['s']
	else:
		sinks = ['s','p']
	for n in range(nstates):
		for snk in sinks:
			for k in p.keys():
				if int(k[-1]) == n and k[-2:-1] == snk:
					prior[k] = gv.gvar(p[k][0], p[k][1])
				else: pass
	return prior

class pion_fit_function():
	def __init__(self,nstates,T,sinks):
		self.sinks = sinks
		self.T = T
		self.nstates = nstates
		return None
	def A(self,n,snk,p):
		return p['A%s%s' %(snk,n)] 
	def E(self,n,snk,p):
		E = p['E%s0' %snk] 
		for ns in range(1,n+1):
			E += np.exp(p['B%s%s' %(snk,ns)])
		return E
	def C(self,n,snk,p):
		return p['C%s%s' %(snk,n)] 
	def func(self,t_dict,p):
		r = 0
		t = t_dict["x"]
		snk = t_dict["sink"]
		#r += self.C(0,snk,p)*np.ones_like(t)
		for n in range(0,self.nstates):
			An = self.A(n,snk,p)
			En = self.E(n,snk,p)
			r += An * np.exp(-En*t) + An*np.exp(-En*(self.T-t))
		return r
	def full_func(self,t,p):
		r = []
		for snk in self.sinks:
			t_dict = {"x":t,"sink":snk}
			r = np.concatenate((np.array(r),self.func(t_dict,p)))
		#print r
		return r

def pipi_priors(nstates,nsinks,Pcm): 
	p = dict()
	p['Cs0'] = [1E-10,5E-11]
	p['Es0'] = [0.45,0.05]
	p['Bs1'] = [-5.0,3.0]
	p['Bs2'] = [-5,5.0]
	p['As0'] = [0.0000008,0.0000008]
	p['As1'] = [0.00000,0.0000008]
	p['As2'] = [0.,0.0000008]
	#p['Cp0'] = [0.0,1.0]
	p['Ep0'] = [0.5,0.2]
	p['Bp1'] = [-100,50]
	p['Bp2'] = [-100,50]
	p['Ap0'] = [9E-4,0.001]
	p['Ap1'] = [0.,8E-4]
	p['Ap2'] = [0.,8E-6]
	prior = dict()
	if nsinks == 1:
		#sinks = ['s']
		sinks = ['s']
	else:
		sinks = ['s','p']
	for n in range(nstates):
		for snk in sinks:
			for k in p.keys():
				if int(k[-1]) == n and k[-2:-1] == snk:
					prior[k] = gv.gvar(p[k][0], p[k][1])
				else: pass
	return prior

class pipi_fit_function():
	def __init__(self,nstates,T,sinks):
		self.sinks = sinks
		self.T = T
		self.nstates = nstates
		return None
	def A(self,n,snk,p):
		return p['A%s%s' %(snk,n)] 
	def E(self,n,snk,p):
		E = p['E%s0' %snk] 
		for ns in range(1,n+1):
			E += np.exp(p['B%s%s' %(snk,ns)])
		return E
	def C(self,n,snk,p):
		return p['C%s%s' %(snk,n)] 
	def func(self,t_dict,p):
		t = t_dict["x"]
		snk = t_dict["sink"]
		r = self.C(0,snk,p)*np.ones_like(t)
		for n in range(0,self.nstates):
			An = self.A(n,snk,p)
			En = self.E(n,snk,p)
			r += An * np.exp(-En*t) + An*np.exp(-En*(self.T-t))
		return r
	def full_func(self,t,p):
		r = []
		for snk in self.sinks:
			t_dict = {"x":t,"sink":snk}
			r = np.concatenate((np.array(r),self.func(t_dict,p)))
		#print r
		return r