import os
import sys
import h5py as h5
import lsqfit
import gvar as gv
import numpy as np
import datetime
import string
import tables as pyt
import weights as wgt




def make_cfg_list(data_dir,system):
	cfg_count = 0
	cfgs = []
	if system == "two pion":
		print "acquiring two pion data from .h5 files"
		for file in os.listdir(data_dir):
			if (file.endswith('.h5') and file.startswith('mm')):
				cfg_count += 1
				cfgs.append(file)
		print str(cfg_count)+" files to compile"
	elif system == "one pion":
		print "acquiring one pion data from .h5 files"
		for file in os.listdir(data_dir):
			if (file.endswith('.h5') and file.startswith('l')):
				cfg_count += 1
				cfgs.append(file)
		print str(cfg_count)+" files to compile"
	return cfgs, cfg_count

def get_h5_group(path,file,filename):
	group_list = []
	if path not in file:
		print path+" does not exist in file "+filename
		pass
	else:
		for g in file.get_node(path,''):
			temp = str(g)
			g_list = temp.split('/')
			group_temp = g_list[-1].split('(')
			group_name = group_temp[0].strip(' ')
			group_list.append(group_name)
	return group_list

def get_h5_array_names(path,file,filename):
	array_list = []
	for g in file.get_node(path,''):
		temp = str(g)
		#print temp
		g_list = temp.split('/')
		x = g_list[-1]
		x_list = x.split('(')
		array_list.append(x_list[0].strip(' '))
	return array_list

def get_real_part(data):
	smearing0 = np.squeeze(data[:,0,:])
	smearing1 = np.squeeze(data[:,1,:])
	Re0 = []
	Re1 = []
	for i,d0 in enumerate(smearing0):
		Re0.append(d0.real)
	for j,d1 in enumerate(smearing1):
		Re1.append(d1.real)
	real_data = np.asarray([Re0,Re1])
	#print np.shape(real_data)
	return real_data

def get_imag_part(data):
	smearing0 = np.squeeze(data[:,0,:])
	smearing1 = np.squeeze(data[:,1,:])
	Im0 = []
	Im1 = []
	for i,d0 in enumerate(smearing0):
		Im0.append(d0.imag)
	for j,d1 in enumerate(smearing1):
		Im1.append(d1.imag)
	imag_data = np.asarray([Im0,Im1])
	#print np.shape(imag_data)
	return imag_data

def get_data(data_dir,system):
	#extracts data from h5 file and returns it in a dictionary
	#pipi data is returned with structure:
	#	dict[P_cm][p_rel][cfg][source]["Re"/"Im"] = np.array(2,48) (two smearings)
	if system == "two pion":
		data_dict = dict()
		cfgs, cfg_count = make_cfg_list(data_dir,system)
		cfg_list = []
		for filename in cfgs:
			file = data_dir+filename
			if not os.path.isfile(file):
				print 'file does not exist: '+filename
				pass
			else:
				f = h5.File(file,'r')
				temp_list = filename.split('_')
				config = temp_list[-1][:-3]
				cfg_list.append(config)	
				fread = pyt.open_file(file)
				path = '/wf1p0_m51p3_l512_a52p0_smrw4p2_n60/mm/ml0p0158_ms0p0902/pipi/'
				COM_momenta = get_h5_group(path,fread,filename)
				for i,c in enumerate(COM_momenta):
					if c not in data_dict:
						data_dict[c] = dict()
					c_path = path+str(c)
					rel_momenta = get_h5_group(c_path,fread,filename)
					for j,r in enumerate(rel_momenta):
						if r not in data_dict[c]:
							data_dict[c][r] = dict()
						r_path = c_path+'/'+str(r)
						sources = get_h5_array_names(r_path,fread,filename)
						if config not in data_dict[c][r]:
							data_dict[c][r][config] = dict()
						for x in sources:
							data = f["wf1p0_m51p3_l512_a52p0_smrw4p2_n60"]["mm"]["ml0p0158_ms0p0902"]["pipi"][c][r][x].value
							#print np.shape(data)
							#raw_input("proceed")
							real_data = get_real_part(data)
							imag_data = get_imag_part(data)
							if x in data_dict[c][r][config]:
								print "duplicate data in "+str(c)+'/'+str(r)+'/'+str(x)+" for config "+config
								pass
							else:
								data_dict[c][r][config][x] = dict()
								data_dict[c][r][config][x]["Re"] = real_data
								data_dict[c][r][config][x]["Im"] = imag_data
				f.close()
				fread.close()
		return data_dict, cfg_count, cfg_list
	elif system == "one pion":
		data_dict = dict()
		cfgs, cfg_count = make_cfg_list(data_dir,system)
		cfg_list = []
		for filename in cfgs:
			file = data_dir+filename
			if not os.path.isfile(file):
				print 'file does not exist: '+filename
				pass
			else:
				f = h5.File(file,'r')
				temp_list = filename.split('_')
				config = temp_list[-1][:-3]
				cfg_list.append(config)	
				fread = pyt.open_file(file)
				for pi_type in ["spectrum","dispersion"]:
					path = '/wf1p0_m51p3_l512_a52p0_smrw4p2_n60/'+pi_type+'/ml0p0158_ms0p0902/pion'
					if pi_type == "spectrum":
						c = "Px0Py0Pz0"
						if c not in data_dict:
							data_dict[c] = dict()
						sh = "0.0"
						if sh not in data_dict[c]:
							data_dict[c][sh] = dict()
						if config not in data_dict[c][sh]:
							data_dict[c][sh][config] = dict()
						sources = get_h5_array_names(path,fread,filename)
						for x in sources:
							data = f["wf1p0_m51p3_l512_a52p0_smrw4p2_n60"][pi_type]["ml0p0158_ms0p0902"]["pion"][x].value
							real_data = get_real_part(data)
							imag_data = get_imag_part(data)
							if x in data_dict[c][sh][config]:
								print "duplicate data in "+str(c)+'/'+str(sh)+'/'+str(x)+" for config "+config
								pass
							else:
								data_dict[c][sh][config][x] = dict()
								data_dict[c][sh][config][x]["Re"] = real_data
								data_dict[c][sh][config][x]["Im"] = imag_data
					elif pi_type == "dispersion":
						#data_dict["boosted"] = dict()
						path = '/wf1p0_m51p3_l512_a52p0_smrw4p2_n60/'+pi_type+'/ml0p0158_ms0p0902/pion/'
						COM_momenta = get_h5_group(path,fread,filename)
						temp = [i.replace("_", "") for i in COM_momenta]
						CM_momenta = [i.replace("p", "P") for i in temp]
						#print CM_momenta
						for i,c in enumerate(COM_momenta):
							c2 = CM_momenta[i]
							if c2 not in data_dict:
								data_dict[c2] = dict()
								r = 0
								data_dict[c2][r] = dict()
							#print c
							c_path = path+str(c)
							sources = get_h5_array_names(c_path,fread,filename)
							if config not in data_dict[c2][r]:
								data_dict[c2][r][config] = dict()
								for x in sources:
									data = f["wf1p0_m51p3_l512_a52p0_smrw4p2_n60"][pi_type]["ml0p0158_ms0p0902"]["pion"][c][x].value
									real_data = get_real_part(data)
									imag_data = get_imag_part(data)
									if x in data_dict[c2][r][config]:
										print "duplicate data in "+str(c)+'/'+str(x)+" for config "+config
										pass
									else:
										data_dict[c2][r][config][x] = dict()
										data_dict[c2][r][config][x]["Re"] = real_data
										data_dict[c2][r][config][x]["Im"] = imag_data
					else:
						pass
				fread.close()
				f.close()
		return data_dict, cfg_count, cfg_list

def output_data(data_dict, filename,data_type):
	print "saving "+data_type+" data to file "+filename
	fout = h5.File(filename,'w')
	for c in data_dict:
		CM_group = fout.create_group(c)
		for r in data_dict[c]:
			r_group = CM_group.create_group(str(r))
			for cfg in data_dict[c][r]:
				cfg_group = r_group.create_group(str(cfg))
				for x in data_dict[c][r][cfg]:
					x_group = cfg_group.create_group(str(x))
					for part in data_dict[c][r][cfg][x]:
						smearing0 = data_dict[c][r][cfg][x][part][0] #[i for i in data_dict[c][r][cfg][x][part][0]] 
						smearing1 = data_dict[c][r][cfg][x][part][1] #[i for i in data_dict[c][r][cfg][x][part][1]]
						temp = [smearing0,smearing1]
						data = np.asarray(temp)
						dset = x_group.create_dataset(str(part),data = data)
						#dset = p_group.create_dataset(str(part),data = data,dtype = float)
	fout.close()

def get_compiled_data(filename):
	#extracts data from h5 file and returns it in a dictionary
	#pipi data is returned with structure:
	#dict[P_cm][p_rel][cfg][source]["Re"/"Im"] = np.array(2,48) (two smearings)def read_in_summed_shells(filename):
	data_dict = dict()
	work_dir  = os.getcwd()
	all_cfgs = []
	file = work_dir+'/'+filename
	if not os.path.isfile(file):
		print 'file does not exist: '+filename
		pass
	else:
		print "reading in data from file "+filename
		h5f = h5.File(file,'r')
		fread = pyt.open_file(file)
		path = '/'
		COM_momenta = get_h5_group(path,fread,filename)
		for c in COM_momenta:
			if c not in data_dict:
				data_dict[c] = dict()
			r_path = path+str(c)
			prel = get_h5_group(r_path,fread,filename)
			for r in prel:
				if r not in data_dict[c]:
					data_dict[c][r] = dict()
				cfg_path = r_path+'/'+str(r)
				config_list = get_h5_group(cfg_path,fread,filename)
				for config in config_list:
					if config not in all_cfgs:
						all_cfgs.append(config)
					if config not in data_dict[c][r]:
						data_dict[c][r][config] = dict()
					x_path = cfg_path+'/'+str(config)
					sources = get_h5_group(x_path, fread, filename)
					for x in sources:
						if x not in data_dict[c][r][config]:
							data_dict[c][r][config][x] = dict()
						p_path = x_path+'/'+str(x)
						parts = get_h5_array_names(p_path,fread,filename)
						for part in parts:
							if part not in data_dict[c][r][config][x]:
								data_dict[c][r][config][x][part] = h5f[c][r][config][x][part].value
		h5f.close()
		fread.close()
	N_cfgs = len(all_cfgs)
	return data_dict, N_cfgs, all_cfgs

def average_sources(data_dict, pi_type):
	#pipi data is returned with structure:
	#	dict[P_cm][p_rel][cfg]["Re"/"Im"] = np.array(2,48) (two smearings)
	#the data is gvar type
	print "averaging sources for "+pi_type+" data"
	avg_data = dict()
	for c in data_dict.keys():
		avg_data[c] = dict()
		for r in data_dict[c].keys():
			avg_data[c][r] = dict()
			for f in data_dict[c][r].keys():
				avg_data[c][r][f] = dict()
				real_sources = []
				imag_sources = []
				for s in data_dict[c][r][f].keys():
					real_sources.append(data_dict[c][r][f][s]["Re"])
					imag_sources.append(data_dict[c][r][f][s]["Im"])
				avg_data[c][r][f]["Re"] = gv.gvar(np.mean(real_sources,axis=0),np.std(real_sources,axis=0))
				avg_data[c][r][f]["Im"] = gv.gvar(np.mean(imag_sources,axis=0),np.std(imag_sources,axis=0))
	return avg_data

def output_averaged_sources(avg_data, filename, data_type, cfg_array = False):
	print "saving "+data_type+" averaged data to .h5 file"
	fout = h5.File(filename,'w')
	for c in avg_data:
		CM_group = fout.create_group(c)
		for r in avg_data[c]:
			r_group = CM_group.create_group(str(r))
			if cfg_array:
				ReTemp = []
				ImTemp = []
				for f in avg_data[c][r]:
					smearing0R = [i.mean for i in avg_data[c][r][f]["Re"][0]] 
					smearing1R = [i.mean for i in avg_data[c][r][f]["Re"][1]]
					ReTemp.append(np.asarray([float(f)*np.ones_like(smearing0R),smearing0R,smearing1R]))
					smearing0I = [i.mean for i in avg_data[c][r][f]["Im"][0]] 
					smearing1I = [i.mean for i in avg_data[c][r][f]["Im"][1]]
					ImTemp.append(np.asarray([float(f)*np.ones_like(smearing0I),smearing0I,smearing1I]))
				ReData = np.asarray(ReTemp)
				ImData = np.asarray(ImTemp)
				dset = r_group.create_dataset("Re",data = ReData)
				dset = r_group.create_dataset("Im",data = ImData)
			else:
				for f in avg_data[c][r]:
					cfg_group = r_group.create_group(str(f))
					#print avg_data[c][r][f].keys()
					for part in avg_data[c][r][f].keys():
						#print part
						temp = []
						smearing0 = [i.mean for i in avg_data[c][r][f][part][0]] 
						smearing1 = [i.mean for i in avg_data[c][r][f][part][1]]
						temp.append([smearing0,smearing1])
						data = np.asarray(temp)
						dset = cfg_group.create_dataset(part,data = data)
	fout.close()

def get_averaged_data(filename, cfg_array = False):
	#extracts data from h5 file and returns it in a dictionary
	#data is returned with structure:
	#dict[P_cm][shell][cfg]["Re"/"Im"] = np.array(2,48) (two smearings)
	#this is true even if cfg_array is flagged as True!
	data_dict = dict()
	work_dir  = os.getcwd()
	all_cfgs = []
	file = work_dir+'/'+filename
	if not os.path.isfile(file):
		print 'file does not exist: '+filename
		pass
	else:
		print "reading in data from file "+filename
		h5f = h5.File(file,'r')
		fread = pyt.open_file(file)
		path = '/'
		COM_momenta = get_h5_group(path,fread,filename)
		for c in COM_momenta:
			if c not in data_dict:
				data_dict[c] = dict()
			r_path = path+str(c)
			prel = get_h5_group(r_path,fread,filename)
			for r in prel:
				if r not in data_dict[c]:
					data_dict[c][r] = dict()
				if cfg_array:
					config_list = []
					p_path = r_path+'/'+str(r)
					parts = get_h5_array_names(p_path,fread,filename)
					#print parts
					for p in parts:
						temp = h5f[c][r][p].value
						#print np.shape(temp)
						config_list = [i[0][0] for i in temp]
						#print config_list
						for i,f in enumerate(config_list):
							if f not in data_dict[c][r]:
								data_dict[c][r][f] = dict()
							data_dict[c][r][f][p] = [temp[i][1],temp[i][2]]
				else:
					cfg_path = r_path+'/'+str(r)
					config_list = get_h5_group(cfg_path,fread,filename)
					for config in config_list:
						if config not in all_cfgs:
							all_cfgs.append(config)
						if config not in data_dict[c][r]:
							data_dict[c][r][config] = dict()
						p_path = cfg_path+'/'+str(config)
						#print p_path
						parts = get_h5_array_names(p_path,fread,filename)
						for p in parts:
							if p not in data_dict[c][r][config]:
								data_dict[c][r][config][p] = h5f[c][r][config][p].value
		h5f.close()
		fread.close()
	N_cfgs = len(all_cfgs)
	return data_dict, N_cfgs, all_cfgs

def add_rel_momenta(data, system):
	#pipi data is returned with structure:
	#	dict[P_cm][shell][cfg]["Re"/"Im"] = np.array(2,48) (two smearings)
	#print "PxPyPz,P_{CM}^{2}"
	#print "PxPyPz,Pcm^2,pxpypx_qxqyqz,p_{rel}^{2}"
	if system == "two pion":
		print "performing weighted relative momentum sum for "+system+" data"
		rm_summed_data = dict()
		momentum_shells = dict()
		for c in data:
			Pcm,k = extract_CM(c)
			#print str(c)+','+str(Pcm)
			if c not in rm_summed_data:
				rm_summed_data[c] = dict()
			for r in data[c]:
				s = extract_shell(r)
				#print str(c)+','+str(Pcm)+','+str(r)+','+str(s)
				if s not in rm_summed_data[c]:
					rm_summed_data[c][s] = dict()
					for f in data[c][r]:
						#w = wgt.compute_weight_rel(Pcm,s,k)
						w = wgt.compute_weight_rel(Pcm,s)*np.ones_like(data[c][r][f]["Re"])
						#print data[c][r][f].keys()
						if f not in rm_summed_data[c][s]:
							rm_summed_data[c][s][f] = dict()
						for part in data[c][r][f]:
							if part not in rm_summed_data[c][s][f]:
								rm_summed_data[c][s][f][part] = np.zeros_like(w)
							rm_summed_data[c][s][f][part] = np.add(rm_summed_data[c][s][f][part],w*data[c][r][f][part])
		return rm_summed_data
	elif system == "one pion":
		print "consolidating dictionary for "+system+" data"
		rm_summed_data = dict()
		for c in data:
			if c not in rm_summed_data:
				rm_summed_data[c] = dict()
				s = 'na'
				rm_summed_data[c][s] = dict()
			for r in data[c]:
				for f in data[c][r]:
					if f not in rm_summed_data[c][s]:
						rm_summed_data[c][s][f] = dict()
					for part in data[c][r][f]:
						if part not in rm_summed_data[c][s][f]:
							rm_summed_data[c][s][f][part] = dict()
						temp = data[c][r][f][part]
						#print np.shape(temp)
						rm_summed_data[c][s][f][part] = temp
		return rm_summed_data

def extract_CM(s):
	if not os.path.isfile("Momenta.txt"):
		f = open("Momenta.txt",'w')
		f.write("Momentum shells\n")
		f.close()
	f = open("Momenta.txt",'a')
	x = s.find("Px")+2
	x_sgn = 1.0
	if s[x] == 'm':
		x_sgn = -1.0
		x += 1 
	y = s.find("Py")+2
	y_sgn = 1.0
	if s[y] == 'm':
		y_sgn = -1.0
		y += 1 
	z = s.find("Pz")+2
	z_sgn = 1.0
	if s[z]== 'm':
		z_sgn = -1.0
		z += 1 
	x_val = x_sgn*float(s[x])
	y_val = y_sgn*float(s[y])
	z_val = z_sgn*float(s[z])
	k = 3
	if x_val == 0:
		k = k-1
	if y_val == 0:
		k = k-1
	if z_val == 0:
		k = k-1
	Pcm = x_val*x_val + y_val*y_val + z_val*z_val
	f.write(s+'\n')
	f.write(
		"Pcm^2 = {0} ({1}, {2}, {3},)\n".format(
			str(Pcm).ljust(2),
			str(x_val).rjust(2),
			str(y_val).rjust(2),
			str(z_val).rjust(2)
			)
		)
	f.close()
	return Pcm,k

def extract_shell(s):
	if not os.path.isfile("Momenta.txt"):
		f = open("Momenta.txt",'w')
		f.write("Momentum shells\n")
		f.close()
	f = open("Momenta.txt",'a')
	x = s.find("px")+2
	x_sgn = 1.0
	if s[x] == 'm':
		x_sgn = -1.0
		x += 1 
	y = s.find("py")+2
	y_sgn = 1.0
	if s[y] == 'm':
		y_sgn = -1.0
		y += 1 
	z = s.find("pz")+2
	z_sgn = 1.0
	if s[z]== 'm':
		z_sgn = -1.0
		z += 1 
	x_val = x_sgn*float(s[x])
	y_val = y_sgn*float(s[y])
	z_val = z_sgn*float(s[z])
	Ptot1 = x_val*x_val + y_val*y_val + z_val*z_val
	f.write(s+'\n')
	f.write(
		"p^2 = {0} ({1}, {2}, {3},)\n".format(
			str(Ptot1).ljust(2),
			str(x_val).rjust(2),
			str(y_val).rjust(2),
			str(z_val).rjust(2)
			)
		)
	x = s.find("qx")+2
	x_sgn = 1.0
	if s[x] == 'm':
		x_sgn = -1.0
		x += 1 
	y = s.find("qy")+2
	y_sgn = 1.0
	if s[y] == 'm':
		y_sgn = -1.0
		y += 1 
	z = s.find("qz")+2
	z_sgn = 1.0
	if s[z]== 'm':
		z_sgn = -1.0
		z += 1 
	x_val = x_sgn*float(s[x])
	y_val = y_sgn*float(s[y])
	z_val = z_sgn*float(s[z])
	Ptot2 = x_val*x_val + y_val*y_val + z_val*z_val
	f.write(
		"q^2 = {0} ({1}, {2}, {3},)\n".format(
			str(Ptot2).ljust(2),
			str(x_val).rjust(2),
			str(y_val).rjust(2),
			str(z_val).rjust(2)
			)
		)
	shell = Ptot1+Ptot2 
	f.write("shell = "+str(shell)+'\n')
	f.close()
	#return "Shell "+str(shell)
	return shell

def output_summed_shells(rm_summed_data, filename,data_type):
	print "saving "+data_type+" summed shell data to .h5 file"
	fout = h5.File(filename,'w')
	for c in rm_summed_data:
		CM_group = fout.create_group(c)
		for s in rm_summed_data[c]:
			s_group = CM_group.create_group(str(s))
			for f in rm_summed_data[c][s]:
				temp = []
				for part in rm_summed_data[c][s][f]:
					#print np.shape(rm_summed_data[c][s][f][part])
					smearings = np.squeeze(rm_summed_data[c][s][f][part])
					temp.append([smearings[0],smearings[1]])
				#print np.shape(temp)
				data = np.asarray(temp)
				dset = s_group.create_dataset(str(f),data = data)
	fout.close()

def read_in_summed_shells(filename):
	#this produces a dictionary with structure:
	#data[pcm][shell][cfg][Re/Im][p/s]
	print "reading in summed shells from file "+filename
	rm_summed_data = dict()
	work_dir  = os.getcwd()
	file = work_dir+'/'+filename
	if not os.path.isfile(file):
		print 'file does not exist: '+filename
		pass
	else:
		h5f = h5.File(file,'r')
		fread = pyt.open_file(file)
		path = '/'
		COM_momenta = get_h5_group(path,fread,filename)
		N_cfgs = 0
		for c in COM_momenta:
			if c not in rm_summed_data:
				rm_summed_data[c] = dict()
			s_path = path+str(c)
			shells = get_h5_group(s_path,fread,filename)
			for s in shells:
				if s not in rm_summed_data[c]:
					rm_summed_data[c][s] = dict()
				N_cfgs = len(rm_summed_data[c][s].keys())
				cfg_path = s_path+'/'+str(s)
				config_list = get_h5_array_names(cfg_path,fread,filename)
				for config in config_list:
					if config not in rm_summed_data[c][s]:
						rm_summed_data[c][s][config] = dict()
						rm_summed_data[c][s][config]['Re'] = dict()
						rm_summed_data[c][s][config]['Im'] = dict()
					temp = h5f[c][s][config].value
					rm_summed_data[c][s][config]['Re']['p'] = temp[0][0]
					rm_summed_data[c][s][config]['Re']['s'] = temp[0][1]
					rm_summed_data[c][s][config]['Im']['p'] = temp[1][0]
					rm_summed_data[c][s][config]['Im']['s'] = temp[1][1]
		h5f.close()
		fread.close()
	return rm_summed_data, N_cfgs

def average_CM(data, system):
	#this returns data[pcmsq][prelsq][cfg][Re/Im][p/s]
	print "averaging COM contributions for "+system+" data"
	CM_avg_data = dict()
	for c in data:
		Pcm,k = extract_CM(c)
		if Pcm not in CM_avg_data:
			CM_avg_data[Pcm] = dict()
		for s in data[c]:
			if s not in CM_avg_data[Pcm]:
				CM_avg_data[Pcm][s] = dict()
			for f in data[c][s]:
				if f not in CM_avg_data[Pcm][s]:
					CM_avg_data[Pcm][s][f] = dict()
				for part in data[c][s][f]:
					if part not in CM_avg_data[Pcm][s][f]:
						CM_avg_data[Pcm][s][f][part] = dict()
						CM_avg_data[Pcm][s][f][part]['p'] = np.zeros_like(data[c][s][f][part]['p'])
						CM_avg_data[Pcm][s][f][part]['s'] = np.zeros_like(data[c][s][f][part]['s'])
					CM_avg_data[Pcm][s][f][part]['p'] = np.add(CM_avg_data[Pcm][s][f][part]['p'],data[c][s][f][part]['p'])
					CM_avg_data[Pcm][s][f][part]['s'] = np.add(CM_avg_data[Pcm][s][f][part]['s'],data[c][s][f][part]['s'])
	return CM_avg_data

def output_CM_averaged(CM_avg_data, filename,data_type):
	print "saving "+data_type+" averaged CM data to .h5 file"
	fout = h5.File(filename,'w')
	for c in CM_avg_data:
		CM_group = fout.create_group(str(c))
		for s in CM_avg_data[c]:
			s_group = CM_group.create_group(str(s))
			for f in CM_avg_data[c][s]:
				f_group = s_group.create_group(str(f))
				for part in CM_avg_data[c][s][f]:
					p_group = f_group.create_group(str(part))
					for smearing in ['p','s']:
						#smearing_group = p_group.create_group(str(smearing))
						data = np.asarray(CM_avg_data[c][s][f][part][smearing])
						#print data
						#print np.shape(data)
						#dset = smearing_group.create_dataset(smearing,data = data)
						dset = p_group.create_dataset(smearing,data = data)
	fout.close()

def read_in_CM_averaged(filename):
	#this produces a dictionary with structure:
	#data[pcm][shell][cfg][Re/Im][p/s]
	print "reading in CM averaged data from file "+filename
	CM_avg_data = dict()
	work_dir  = os.getcwd()
	file = work_dir+'/'+filename
	if not os.path.isfile(file):
		print 'file does not exist: '+filename
		pass
	else:
		h5f = h5.File(file,'r')
		fread = pyt.open_file(file)
		path = '/'
		COM_momenta = get_h5_group(path,fread,filename)
		cfgs = []
		for c in COM_momenta:
			if c not in CM_avg_data:
				CM_avg_data[c] = dict()
			s_path = path+str(c)
			shells = get_h5_group(s_path,fread,filename)
			for s in shells:
				if s not in CM_avg_data[c]:
					CM_avg_data[c][s] = dict()
				N_cfgs = len(CM_avg_data[c][s].keys())
				cfg_path = s_path+'/'+str(s)
				config_list = get_h5_group(cfg_path,fread,filename)
				for config in config_list:
					if config not in cfgs:
						cfgs.append(config)
					if config not in CM_avg_data[c][s]:
						CM_avg_data[c][s][config] = dict()
					p_path = cfg_path+'/'+str(config)
					parts = get_h5_group(p_path, fread, filename)
					for part in parts:
						if part not in CM_avg_data[c][s][config]:
							CM_avg_data[c][s][config][part] = dict()
						smearing_path = p_path+'/'+str(part)
						smearings = get_h5_array_names(smearing_path, fread, filename)
						for smearing in smearings:
							if smearing not in CM_avg_data[c][s][config][part]:
								CM_avg_data[c][s][config][part][smearing] = h5f[c][s][config][part][smearing].value
		h5f.close()
		fread.close()
		N_cfgs = len(cfgs)
	return CM_avg_data, N_cfgs

def format_data(data_dict):
	#returns data as a dictionary in format:
	#	dict[Pcmsq][shell][Re/Im][s/p] = np.array(N_cfgs,t)
	formatted_data_dict = dict()
	cfg_list = []
	for c in data_dict:
		if c not in formatted_data_dict:
			formatted_data_dict[c] = dict()
		for s in data_dict[c]:
			if s not in formatted_data_dict[c]:
				formatted_data_dict[c][s] = dict()
			Re_s_data = []
			Re_p_data = []
			Im_s_data = []
			Im_p_data = []
			for f in data_dict[c][s]:
				cfg_list.append(f)
				Re_s_data.append(data_dict[c][s][f]["Re"]["s"])
				Re_p_data.append(data_dict[c][s][f]["Re"]["p"])
				Im_s_data.append(data_dict[c][s][f]["Im"]["s"])
				Im_p_data.append(data_dict[c][s][f]["Im"]["p"])
			formatted_data_dict[c][s]["Re"] = dict()
			formatted_data_dict[c][s]["Im"] = dict()
			formatted_data_dict[c][s]["Re"]["s"] = np.array(Re_s_data)
			formatted_data_dict[c][s]["Re"]["p"] = np.array(Re_p_data)
			formatted_data_dict[c][s]["Im"]["s"] = np.array(Im_s_data)
			formatted_data_dict[c][s]["Im"]["p"] = np.array(Im_p_data)
	cfg_set = set(cfg_list)
	return formatted_data_dict, cfg_set

def fold_data(data):
	#takes array of shape (N_cfgs, N_t) and folds it
	flipped = np.flip(data,1)
	for i in range(len(flipped[1])):
		if (flipped[1][-(i+1)] == data[1][i]):
			pass
		else:
			print "entry %s does not match"%i
	folded = 0.5*(np.roll(data,1,axis=1) + flipped)
	return folded