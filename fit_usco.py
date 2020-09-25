
import sys
import os
import numpy as np
import pandas as pd
from isochrones import get_ichrone,SingleStarModel
from isochrones.priors import GaussianPrior,FlatPrior
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


#================== Load everything only in master ===================
#--------------- Observables ----------------------------------
identifier   = "GDR2_ID"
gaia         = ["parallax","BP","G","RP"]
gaia_obs     = ["parallax","BP","G","RP"]
gaia_unc     = ["parallax_error","e_BP","e_G","e_RP"]
dosmass      = ["J","H","K"]
dosmass_obs  = ["J","H","Ks"]
dosmass_unc  = ["e_J","e_H","e_Ks"]
panstar      = ["PS_g", "PS_r", "PS_i", "PS_z","PS_y"]
panstar_obs  = ["g_sdss","r_sdss","i_sdss","z_sdss","Y"]
panstar_unc  = ["e_g_sdss","e_r_sdss","e_i_sdss","e_z_sdss","e_Y"]
bands        = ['BP','G','RP','J','H','K',"PS_g", "PS_r", "PS_i", "PS_z","PS_y"]
parameters   = ["age","mass","distance","AV"]
#---------------------------------------------------------------

#-------------- Get the isochrones ----------------------------
mist = get_ichrone('mist', bands=bands)

#------- This prints the available photometry ------------------
# mass, age, feh = (1.03, 9.72, -0.11)
# print(mist.generate(mass, age, feh, return_dict=True))
# sys.exit()
#---------------------------------------------------------------

#-------------- Transform lists -------------------------------
bands_mag = [ band+"_mag" for band in bands]
phot_obs  = sum([gaia_obs,dosmass_obs,panstar_obs],[])
phot_unc  = sum([gaia_unc,dosmass_unc,panstar_unc],[])
columns_data = sum([[identifier],gaia_obs,dosmass_obs,panstar_obs,
						gaia_unc,dosmass_unc,panstar_unc],[])

stats_names = sum([[identifier],["mean_"+p for p in parameters],
					["std_"+p for p in parameters]],[])
phot_obs.pop(0)
phot_unc.pop(0)
#--------------------------------------------------------------

#------------ Files -------------------------------
dir_base   = "/raid/jromero/OCs/USco/"
file_data  = dir_base + "members_GDR2_p05_remaining.csv"
dir_out    = dir_base + "Isochrones/"
dir_chain  = dir_out  + "chains/"
dir_plots  = dir_out  + "plots/"
file_samp  = dir_out  + "samples_rem.h5"
file_stat  = dir_out  + "statistics_rem.csv"
os.makedirs(dir_plots,exist_ok=True)
#---------------------------------------------------

#------------- Load data --------------------------------
df = pd.read_csv(file_data,usecols=columns_data)
df.replace(to_replace=99.0,value=np.nan,inplace=True)
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
	if datum["BP"] > 15.:
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
	for true,obs,unc in zip(panstar,panstar_obs,panstar_unc):
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
	model.set_prior(age=FlatPrior(bounds=(0,15)),
			 AV=GaussianPrior(1,4, bounds=(0,10)),
			distance=GaussianPrior(145,40,bounds=(0,200)))

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
