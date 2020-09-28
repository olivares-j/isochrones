import sys
import os
import numpy as np
import pandas as pd
import h5py

from Globals import *

#------------ Files -------------------------------------------------
file_samples = dir_base + "samples_all.h5"
file_stats   = dir_base + "statistics_all.csv"

h5 = h5py.File(file_samples,'w')

#++++++++++++++++++ Loop over chunks +++++++++++++++++++++++++++
dfs = []
for i in range(1,size+1):
	#------------- Files ---------------------------------------------
	file_samp  = dir_base + "samples_{0}_of_{1}.h5".format(i,size)
	file_stat  = dir_base + "statistics_{0}_of_{1}.csv".format(,size)
	#-----------------------------------------------------------------

	#------------- Load data frames into list --------------------------------
	dfs.append(pd.read_csv(file_stat))
	#---------------------------------------------------------

	#---------------- Read and save samples --------------------------
	with h5py.File(file_samp,'r') as hf:
		for ID in hf.keys():
			h5.create_dataset(ID,data=np.array(hf.get(ID)))
	#-----------------------------------------------------------------
h5.close()

#------- Join data frames -------------------------
df = reduce(lambda left,right: pd.merge(left,right,on=identifier), dfs)

#------------- Save statistics ---------------------
df.set_index(identifier,inplace=True)
df.to_csv(file_stats)
