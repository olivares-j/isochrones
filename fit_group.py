
import sys
import os
import numpy as np
import pandas as pd
from isochrones import get_ichrone,SingleStarModel
from isochrones.priors import GaussianPrior,AVPrior,DistancePrior
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


#================== Load everything only in master ===================
#--------------- Observables ----------------------------------
identifier   = "ID"
gaia         = ["parallax","BP","G","RP"]
gaia_obs     = ["parallax","bp","g","rp"]
gaia_unc     = ["parallax_error","bp_error","g_error","rp_error"]
dosmass      = ["J","H","K"]
dosmass_obs  = ["Jmag","Hmag","Kmag"]
dosmass_unc  = ["e_Jmag","e_Hmag","e_Kmag"]
allwise      = ["W1","W2","W3"]
allwise_obs  = ["W1mag","W2mag","W3mag"]
allwise_unc  = ["e_W1mag","e_W2mag","e_W3mag"]
bands        = ['BP','G','RP','J','H','K','W1','W2','W3']
parameters   = ["age","mass","distance","AV"]
#---------------------------------------------------------------

#-------------- Get the isochrones ----------------------------
mist = get_ichrone('mist', bands=bands)
#---------------------------------------------------------------

#-------------- Transform lists -------------------------------
bands_mag = [ band+"_mag" for band in bands]
phot_obs   = sum([gaia_obs,dosmass_obs,allwise_obs],[])
phot_unc   = sum([gaia_unc,dosmass_unc,allwise_unc],[])
columns_data = sum([[identifier],gaia_obs,dosmass_obs,allwise_obs,
							     gaia_unc,dosmass_unc,allwise_unc],[])
stats_names = sum([[identifier],["mean_"+p for p in parameters],
					["std_"+p for p in parameters]],[])
phot_obs.pop(0)
phot_unc.pop(0)
#--------------------------------------------------------------

#------------ Files -------------------------------
dir_base   = "/raid/jromero/OCs/Taurus/"
dir_base   = "/home/javier/Cumulos/Taurus/run_1/"
file_data  = dir_base + "members.csv"
dir_out    = dir_base + "Isochrones/"
dir_chain  = dir_out  + "chains/"
dir_plots  = dir_out  + "plots/"
file_samp  = dir_out  + "samples.h5"
file_stat  = dir_out  + "statistics.csv"
os.makedirs(dir_plots,exist_ok=True)
#---------------------------------------------------

#------------- Load data --------------------------------
df = pd.read_csv(file_data,usecols=columns_data)
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
	#---------------- Gaia ---------------------------
	for true,obs,unc in zip(gaia,gaia_obs,gaia_unc):
		params[true] = (datum[obs],datum[unc])

	#------- Use BP only in bright sources --------------
	if datum["bp"] > 15.:
		del params["BP"]
	#--------------------------------------------------

	#---------------- 2MASS ----------------------------------
	for true,obs,unc in zip(dosmass,dosmass_obs,dosmass_unc):
		if np.isfinite(datum[obs]):
			if not np.isfinite(datum[unc]):
				datum[unc] = 0.1*datum[obs]

		#-------- Add only if observed -------------
			params[true] = (datum[obs],datum[unc])
	#---------------------------------------------------------

	#---------------- AllWISE ----------------------------------
	for true,obs,unc in zip(allwise,allwise_obs,allwise_unc):
		if np.isfinite(datum[obs]):
			if not np.isfinite(datum[unc]):
				datum[unc] = 0.2*datum[obs]

		#-------- Add only if observed -------------
			params[true] = (datum[obs],datum[unc])
	#---------------------------------------------------------

	
	#======================================================
	

	#------ Start the model ---------------------
	model = SingleStarModel(mist, **params)

	model.mnest_basename = dir_chain+str(ID)+"-"

	#--------- Prior -------------------
	model.set_prior(age=GaussianPrior(7, 1, bounds=(6,10)),
					AV=AVPrior(bounds=(0,10)),
					distance=DistancePrior(500))

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
