
import sys
import os
import numpy as np
import pandas as pd
from isochrones import get_ichrone,SingleStarModel
from isochrones.priors import GaussianPrior,AgePrior
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from Globals import *


#-------------- Get the isochrones ----------------------------
mist = get_ichrone('mist', bands=bands)

#------- This prints the available photometry ------------------
# mass, age, feh = (1.03, 9.72, -0.11)
# print(mist.generate(mass, age, feh, return_dict=True))
# sys.exit()
#---------------------------------------------------------------

#------------ Files -------------------------------------------------
file_chunk = dir_data + "data_{0}_of_{1}.csv".format(XXX,size)
file_samp  = dir_base + "samples_{0}_of_{1}.h5".format(XXX,size)
file_stat  = dir_base + "statistics_{0}_of_{1}.csv".format(XXX,size)
#-------------------------------------------------------------------

#------------- Load data --------------------------------
df = pd.read_csv(file_chunk,usecols=columns_data)
df.set_index(identifier,inplace=True)
#---------------------------------------------------------

#---------------- Initialize statistics --------------------------
stats = pd.DataFrame(data=None,columns=stats_names)
#-----------------------------------------------------------------
#======================================================================

#==================== Loop over stars ==========================
for ID,datum in df.iterrows():
	#-------- Name ----------------
	params = {'name':ID}

	#=============== Observations =======================
	for true,obs,unc in zip(list_mod,list_obs,list_unc):
		if np.isfinite(datum[obs]):
			params[true] = (datum[obs],datum[unc])
	#======================================================
	
	#------ Start the model ---------------------
	model = SingleStarModel(mist, **params)

	model.mnest_basename = dir_chain+str(ID)+"-"

	#--------- Prior -------------------
	model.set_prior(
			age=AgePrior(),
			AV=GaussianPrior(prior_Av["loc"],prior_Av["scale"], 
			 	bounds=(prior_Av["lower"],prior_Av["upper"])),
			distance=GaussianPrior(prior_distance["loc"],prior_distance["scale"],
				bounds=(prior_distance["lower"],prior_distance["upper"]))
			)

	#------ Fit ---------------------------------
	model.fit(n_live_points=1000,verbose=False)
	#-------------------------------------------

	#------ Statistics -------
	mean  = model.derived_samples.mean()
	stds  = model.derived_samples.std()

	row = {identifier:str(ID)}
	for par in parameters:
		row["mean_"+par] = mean[par]
		row["std_"+par]  = stds[par]

	stats = stats.append(row,ignore_index=True)
	#-----------------------------------------------

	#------ Save samples -------
	model.derived_samples.to_hdf(file_samp,key=str(ID),mode="a")

	#------------- Plots -----------------------------------
	plt.figure()
	model.corner_params()
	plt.savefig(dir_plots+"{}_par.png".format(ID))
	plt.close()

	plt.figure()
	model.corner_observed()
	plt.savefig(dir_plots+"{}_obs.png".format(ID))
	plt.close()

	# ------- Photometric fit -----------
	true_phot       = model.derived_samples[bands_mag]
	true_phot_mean  = true_phot.mean()
	delta_true      = true_phot - true_phot_mean 
	
	obs_phot   = pd.DataFrame([datum[phot_obs].to_numpy()],
								columns=bands_mag)
	delta_obs  = (obs_phot - true_phot_mean).to_numpy().flatten()
	#---------------------------------------------------------------

	fig,ax = plt.subplots()
	ax.set_title(str(ID))
	ax = sns.violinplot(data=delta_true,color="gray",
		inner=None,zorder=0)
	ax.errorbar(bands_mag,delta_obs,yerr=datum.loc[phot_unc],
					fmt="o",color="blue",zorder=1,label="Observed")
	ax.set_xlabel('Band')
	ax.set_ylabel('Difference from best model [mag]')
	ax.set_xticklabels(bands)
	#--------- Legend --------------------------------
	handles, labels = ax.get_legend_handles_labels()
	patch = mpatches.Patch(color='grey', label="Samples from model")
	handles.append(patch) 
	plt.legend(handles=handles, loc='upper left')
	#---------------------------------------------------------------

	plt.savefig(dir_plots+"{}_phot.png".format(ID))
	plt.close()
	#---------------------------------------------------------

#------------- Save statistics ---------------------
stats.set_index(identifier,inplace=True)
stats.to_csv(file_stat)
