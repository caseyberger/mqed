import numpy as np



def compute_weight_CM(Pcmsq):
	r3n = [1.,1./6.,1./12.,1./8.,1./6.]
	weight = r3n[int(Pcmsq)]
	return weight

def compute_weight_rel(Pcmsq,nsq):
	weight = 0.
	if Pcmsq == 0 and nsq == 0:
		weight = 1.
	elif Pcmsq == 0 and nsq == 2:
		weight = 1./np.sqrt(6.)
	elif Pcmsq == 1 and nsq == 1:
		weight = 1./np.sqrt(12.)
	elif Pcmsq == 1 and nsq == 3:
		weight = 1./np.sqrt(48.)
	elif Pcmsq == 2 and nsq == 2:
		weight = 1./np.sqrt(48.)
	elif Pcmsq == 2 and nsq == 4:
		weight = 1./np.sqrt(96.)
	elif Pcmsq == 3 and nsq == 3:
		weight = 1./np.sqrt(64.)
	elif Pcmsq == 4 and nsq == 4:
		weight = 1./np.sqrt(6.)
	elif Pcmsq == 4 and nsq == 4:
		weight = 1./np.sqrt(36.)	
	else:
		weight = 0.
	return weight
