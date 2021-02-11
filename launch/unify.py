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
	file_stat  = dir_base + "statistics_{0}_of_{1}.csv".format(i,size)
	#-----------------------------------------------------------------

	#------------- Load data frames into list --------------------------------
	dfi = pd.read_csv(file_stat)
	dfi.set_index(identifier,inplace=True)
	dfs.append(dfi)
	#---------------------------------------------------------

	#---------------- Read and save samples --------------------------
	store = pd.HDFStore(file_samp)
	for ID in store.keys():
		store.get(ID).to_hdf(file_samples,key=ID,mode="a")
	store.close()
	#-----------------------------------------------------------------

h5.close()

#------- Concatenate data frames -------------------------
df = pd.concat(dfs)

#------------- Save statistics ---------------------
df.to_csv(file_stats)
