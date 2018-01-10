import os
import data_management_tools as dmt
import data_analysis_tools as analyze
import plotting_tools as pt
import fitters
import gvar as gv
import numpy as np



nsinks = 1
nstates = 2
t_min = 14
pion_filename = "pion_CM_averaged.h5"
pion_CM_avg, N_cfgs = dmt.read_in_CM_averaged(pion_filename)
pion_formatted, cfg_list = dmt.format_data(pion_CM_avg)
'''
for c in pion_formatted.keys():
	print c
	for s in pion_formatted[c].keys():
		print s
		for p in pion_formatted[c][s].keys():
			print p
			print pion_formatted[c][s][p].keys()
		print np.shape(pion_formatted[c][s]["Re"]["p"])
'''
#now data is in format:
#dict[P_cm][shell][Re/Im][s/p] = np.array(N_cfgs,t)
fit_range = np.arange(t_min,24)
for Pcm in pion_formatted:
	for shell in pion_formatted[Pcm]:
		sources_dict = dict()
		sources_dict["s"] = pion_formatted[Pcm][shell]["Re"]["s"]
		sources_dict["p"] = pion_formatted[Pcm][shell]["Re"]["p"]
		s_data = sources_dict["s"]
		#pt.plot_single_time_slice(s_data, 2, Pcm)
		s_data = dmt.fold_data(s_data)
		raw_data = gv.dataset.avg_data(s_data)
		#print np.shape(raw_data)
		#pt.plot_data(raw_data,Pcm,"Corr")
		#pt.plot_data(raw_data,Pcm,"M_eff")
		#pt.plot_data(raw_data,Pcm,"A_eff")
		#pt.plot_data(raw_data,Pcm,"E_1")
		#pt.plot_data(raw_data,Pcm,"A_1")
		s_data = s_data[:,fit_range]
		p_data = sources_dict["p"]
		p_data = dmt.fold_data(p_data)
		p_data = p_data[:,fit_range]
		if nsinks == 1:
			b0_data = s_data
		else:
			b0_data = np.concatenate((s_data,p_data),axis=1)
		gv_b0_data = gv.dataset.avg_data(b0_data)
		cov_matrix = gv.evalcov(gv_b0_data)/(N_cfgs-1) #if you want to freeze the covariance matrix
		fit0 = fitters.fit_data(gv_b0_data,fit_range,"pion",Pcm,nstates,nsinks,cov=None)
		params = fit0.p
		E0 = fit0.p["Es0"]
		print "E_{0}(t_min = %s, nstates = %s) = %s"%(str(t_min),str(nstates),E0)
		if nsinks == 1:
			fit_plot_data = gv_b0_data
		else:
			fit_plot_data = gv_b0_data[0:len(gv_b0_data)/2]
		pt.plot_fit(fit_plot_data,fit_range,params,Pcm,shell,"b0",nstates,"pion")
		#bootstrap
		fit_dict = dict()
		chi2_dict = dict()
		bs_data_dict = dict()
		fit_dict["b0"] = fit0.p
		chi2_dict["b0"] = fit0.chi2/fit0.dof
		all_data = []
		#bs_list = analyze.bs_corr_list(b0_data,N_cfgs,seed=1,return_b0 = True)
		bs_list = analyze.bs_corr_list(b0_data,5,seed=1,return_b0 = True)
		for bs in range(bs_list.shape[0]):
			bs_data = b0_data[bs_list[bs]]
			#print np.shape(bs_data)
			gv_bs_data = gv.dataset.avg_data(bs_data)
			#print np.shape(gv_bs_data)
			gv_bs_means_only = [i.mean for i in gv_bs_data]
			#print np.shape(gv_bs_means_only)
			#print gv_bs_means_only[0]
			bs_data_dict[bs] = gv_bs_data
			#print np.shape(gv_bs_data)
			fit_bs = fitters.fit_data(gv_bs_data,fit_range,"pion",Pcm,nstates, nsinks, cov=None, print_fit = False)#cov=None if you don't want to freeze the covariance matrix
			fit_dict[bs] = fit_bs.p
			chi2_dict[bs] = fit_bs.chi2/fit_bs.dof
			all_data.append(gv_bs_means_only)
			#plot the bs fit
			if nsinks == 1:
				fit_plot_data = gv_bs_data
			else:
				fit_plot_data = gv_bs_data[0:len(gv_bs_data)/2]
			#pt.plot_fit(fit_plot_data,fit_range,fit_bs.p,Pcm,shell,"bs"+str(bs),nstates,"pion")
		#print np.shape(all_data)
		#print all_data[0][0]
		#print all_data[0][1]
		gv_all_bs = gv.dataset.avg_data(all_data, bstrap = True)
		fit_all = fitters.fit_data(gv_all_bs,fit_range,"pion",Pcm,nstates,nsinks,cov=None)
		bs_params = fit_all.p
		if nsinks == 1:
			all_bs_data = gv_all_bs
		else:
			all_bs_data = gv_all_bs[0:len(gv_all_bs)/2]
		pt.plot_fit(all_bs_data,fit_range,fit_bs.p,Pcm,shell,"bs",nstates,"pion")
		#save all the energy fits from the bootstraps for later data analysis