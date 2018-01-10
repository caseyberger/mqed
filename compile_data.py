import data_management_tools as dmt
import os
import numpy as np




work_dir  = os.getcwd()
data_dir = work_dir+'/Data/' #set this for LLNL later

'''
compile all cfgs and make one .h5 file 
'''

pipi_dict, pipi_cfg_count, pipi_cfgs = dmt.get_data(data_dir,"two pion")
dmt.output_data(pipi_dict,"pipi_all_cfgs.h5","two pion")
pion_dict, pion_cfg_count, pion_cfgs = dmt.get_data(data_dir,"one pion")
dmt.output_data(pion_dict,"pion_all_cfgs.h5","one pion")

'''
read in new .h5 files and average sources, then save to a new .h5 file
'''

pipi_dict, pipi_cfg_count, pipi_cfgs = dmt.get_compiled_data("pipi_all_cfgs.h5")
pipi_avg = dmt.average_sources(pipi_dict,"two pion")
dmt.output_averaged_sources(pipi_avg,"pipi_averaged_sources.h5","two pion", cfg_array = True)
pion_dict, pion_cfg_count, pion_cfgs = dmt.get_compiled_data("pion_all_cfgs.h5")
pion_avg = dmt.average_sources(pion_dict,"one pion")
dmt.output_averaged_sources(pion_avg,"pion_averaged_sources.h5","one pion", cfg_array = True)



'''
read in data with averaged sources, sum the relative momenta into shells, 
	and output to a new file
'''

pipi_avg, pipi_cfg_count, pipi_cfgs = dmt.get_averaged_data("pipi_averaged_sources.h5", cfg_array = True)
pipi_summed_shells = dmt.add_rel_momenta(pipi_avg,"two pion")
dmt.output_summed_shells(pipi_summed_shells,"pipi_summed_shells.h5","two pion")
pion_avg, pion_cfg_count, pion_cfgs = dmt.get_averaged_data("pion_averaged_sources.h5", cfg_array = True)
pion_summed_shells = dmt.add_rel_momenta(pion_avg,"one pion")
dmt.output_summed_shells(pion_summed_shells,"pion_summed_shells.h5","one pion")

'''
read in data with shells, average P_CM contributions, and output to a new file
'''
pipi_summed_shells, pipi_cfg_count = dmt.read_in_summed_shells("pipi_summed_shells.h5")
pipi_CM_averaged = dmt.average_CM(pipi_summed_shells, "two pion")
dmt.output_CM_averaged(pipi_CM_averaged,"pipi_CM_averaged.h5", "two pion")
pion_summed_shells, pion_cfg_count = dmt.read_in_summed_shells("pion_summed_shells.h5")
pion_CM_averaged = dmt.average_CM(pion_summed_shells, "one pion")
dmt.output_CM_averaged(pion_CM_averaged,"pion_CM_averaged.h5", "one pion")

