
import sys
import numpy as np
import pandas as pd
from isochrones import get_ichrone,SingleStarModel
from isochrones.priors import GaussianPrior
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mist = get_ichrone('mist', bands=['G','BP','RP','J','H','K'])
parameters = ["age","mass","distance","AV"]



#--------------- Files ---------------------
identifier = "ID"
gaia       = ["parallax","G","BP","RP"]
gaia_obs   = ["parallax","g","bp","rp"]
gaia_unc   = ["parallax_error","g_error","bp_error","rp_error",
					]
dosmass      = ["J","H","K"]
dosmass_obs  = ["Jmag","Hmag","Kmag"]
dosmass_unc  = ["e_Jmag","e_Hmag","e_Kmag"]




dir_base   = "/raid/jromero/OCs/Taurus/"
dir_out    = dir_base + "Isochrones/"
file_data  = dir_base + "members.csv"
file_samp  = dir_out  + "samples.h5"
file_stat  = dir_out  + "statistics.csv"

columns_data = sum([[identifier],gaia_obs,dosmass_obs,gaia_unc,dosmass_unc],[])
stats_names = sum([[identifier],["mean_"+p for p in parameters],["std_"+p for p in parameters]],[])


#------------- Load data --------------------------------
df = pd.read_csv(file_data,usecols=columns_data)
df.set_index(identifier,inplace=True)

stats = pd.DataFrame(data=None,columns=stats_names)

for ID,datum in df.iterrows():
	#-------- Name ----------------
	params = {'name':ID}

	#=============== Observations =======================
	for true,obs,unc in zip(gaia,gaia_obs,gaia_unc):
		params[true] = (datum[obs],datum[unc])

	for true,obs,unc in zip(dosmass,dosmass_obs,dosmass_unc):
		if np.isfinite(datum[obs]):
			if not np.isfinite(datum[unc]):
				datum[unc] = 0.1*datum[obs]

		#-------- Add only if observed -------------
			params[true] = (datum[obs],datum[unc])

	#------- Use BP only in bright sources --------------
	if datum["bp"] > 15.:
		del params["BP"]
	#======================================================

	#------ Start the model ---------------------
	mod = SingleStarModel(mist, **params)

	#--------- Prior -------------------
	mod.set_prior(age=GaussianPrior(7, 1, bounds=(6,9)))

	#------ Fit -------
	mod.fit()

	#------ Statistics -------
	mean  = mod.derived_samples.mean()
	stds  = mod.derived_samples.std()

	row = {identifier:str(ID)}
	for par in parameters:
		row["mean_"+par] = mean[par]
		row["std_"+par]  = stds[par]

	stats = stats.append(row,ignore_index=True)
	#-----------------------------------------------

	#------ Save samples -------
	mod.derived_samples.to_hdf(file_samp,key=str(ID),mode="a")

	#------------- Plots -----------------------------------
	plt.figure()
	mod.corner_params()
	plt.savefig(dir_out+"Plots/{}_par.png".format(ID))
	plt.close()

	plt.figure()
	mod.corner_observed()
	plt.savefig(dir_out+"Plots/{}_obs.png".format(ID))
	plt.close()
	#---------------------------------------------------------

#------------- Save statistics ---------------------
stats.set_index(identifier,inplace=True)
stats.to_csv(file_stat)
